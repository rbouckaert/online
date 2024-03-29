package online;


import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beastfx.app.tools.LogCombiner;
import beastfx.app.tools.Application;
import beastfx.app.util.XMLFile;
import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.Logger;
import beast.base.inference.Logger.LogFileMode;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.core.Log;
import beast.base.util.Randomizer;
import beast.base.parser.XMLParserException;
import online.math.DistributionComparator;
import online.math.DistributionComparator.ConvergenceCriterion;

@Description("Create tree and trace files extending an input multiple-state file with different set of taxa")
public class TraceExpander extends BaseStateExpander {
	final public Input<XMLFile> xml2Input = new Input<>("xml2", "BEAST XML file with expanded state. If not specified, xml2 is set equal to xml1, and we are resuming", new XMLFile("[[none]]"));
	final public Input<File> multiStateFileInput = new Input<>("multiStateFile", "state file containing multiple states associated with initial XML file (xml1)."
			+ "If not specified, use xml1+\".state.multi\"", new File("[[none]]"));
	
	final public Input<File> stateFileInput = new Input<>("stateFile", "state file containing operarator optimisation information associated with initial XML file (xml1). "
			+ "If not specified, use xml1+\".state\". If that does not exist, ignore.", new File("[[none]]"));

    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads","maximum number of threads to use, if not specified the number of available cores is used (default)");
    final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of states in multiStateFile to be used as burn-in (these will be ignored)", 10);
    final public Input<Boolean> overwriteInput = new Input<>("overwrite", "overwrite existing tree and trace files", true);
    final public Input<Boolean> verboseInput = new Input<>("verbose", "show more information if true", false);
	final public Input<Integer> sampleCountInput = new Input<>("sampleCount", "determines number of samples taken from multiStateFile evenly "
			+ "selected from multiStateFile. If number of states in multiStateFile > sampleCount some entries may be taken multiple times."
			+ "If <= 0, use every sample in the multiStateFile (minus burn-in)", 200);

	/**
	 * https://www.stata.com/new-in-stata/gelman-rubin-convergence-diagnostic/
	 * " Gelman and Rubin (1992) and Brooks and Gelman (1998) suggest that diagnostic Rc values greater than 1.2 for 
	 *   any of the model parameters should indicate nonconvergence. In practice, a more stringent rule of Rc < 1.1 is 
	 *   often used to declare convergence."
	 *   But even 1.1 is often too low (https://arxiv.org/pdf/1812.09384.pdf) and 1.01 makes more sense. 
	 */
    final public Input<Double> thresholdInput = new Input<>("threshold", "threshold appropriate for convergence criterion, "
    		+ "e.g. maximum acceptable value of Gelman Rubin statistic, or minimum p-value for KS test. "
			+ "Set 'criterion' to 'none' to stop after first cycle. "
			+ "This is a comma separated list matching with one threshold for each convergence criterion", 0.4);
    
	final public Input<String> tempDirInput = new Input<>("tempDir","directory where temporary files are written.", "/tmp/");
    final public Input<String> criterionInput = new Input<>("criterion",
			"Comma separated list of convergence criteria. "
			+ "If not set to 'none', the chain keeps afterburning (with chainLength steps) till all items in trace log converge according to all criteria. " +
    		DistributionComparator.convergenceCriterionDescription, 
    		ConvergenceCriterion.never+"");
    final public Input<Integer> maxCycleInput = new Input<>("maxCycle", "maximum number of cycles before stopping. Ignored if negative (4 by default)", 4);
    
    
    private int nrOfThreads;
    private static ExecutorService exec;
    private static CountDownLatch countDown;
	private List<Logger> loggers;
	private Model model1, model2;
	private long sampleNr;
	private BufferedReader fin;
	private String multiStateInputFile;
	private PrintStream multiStateOut;
	private boolean autoConverge;
	private int cycle;
	private int currentSample = 0, sampleCount, availableSamples, currentSampleNr = 0;
	private ConvergenceCriterion [] criteria;

	@Override
	public void initAndValidate() {
	}
	
	@Override
	public void run() throws Exception {
		inputCheck();
		
		
		Long start = System.currentTimeMillis();

//		Log.setLevel(Log.Level.debug);
		initialise();
		
		
		cycle = 0;
		boolean isResuming = importModels();
		sampleCount = sampleCountInput.get();
		int n = getStateCount();
		if (n == 0) {
			throw new IllegalArgumentException("Cannot find any states in the multi state file");
		}
		int burnIn = skipBurnin(isResuming, n);
		if (sampleCount <= 0) {
			sampleCount = n - burnIn;
		}
		availableSamples = n - burnIn;
		
		initLoggers(cycle);
		
		sampleNr = 0;
		if (nrOfThreads == 1) {
			processUnThreaded(xml2Input.get() == null || xml2Input.get().getName().equals("[[none]]"));
		} else {
			// split work over k threads, with work for (n - burnIn)/k states each
			processThreaded(xml2Input.get() == null || xml2Input.get().getName().equals("[[none]]"));
		}
		
		String xml2Path = xml2Input.get().getAbsolutePath();
		close(cycle, xml2Path);

		Log.info("Cycle " + cycle + " done in " + (System.currentTimeMillis()-start)/1000 + " seconds with " + nrOfThreads + " threads");

		if (autoConverge) {
			availableSamples = sampleCount;
			autoConverge(isResuming, burnIn, xml2Path);

		}
		
		Long end = System.currentTimeMillis();
		Log.info("Done in " + (end-start)/1000 + " seconds with " + nrOfThreads + " threads");
		System.exit(0);
	}

	
	private void inputCheck() {
		if (!xml1Input.get().exists()) {
			throw new IllegalArgumentException("File " + xml1Input.get().getName() + " does not exist.");
		}
		if (xml2Input.get() != null && !xml2Input.get().getName().equals("[[none]]") &&
				!xml2Input.get().exists()) {
			throw new IllegalArgumentException("File " + xml2Input.get().getName() + " does not exist.");
		}
		if (stateFileInput.get() != null && !stateFileInput.get().getName().equals("[[none]]") &&
				!stateFileInput.get().exists()) {
			throw new IllegalArgumentException("File " + stateFileInput.get().getName() + " does not exist.");
		}
		if (multiStateFileInput.get() != null && !multiStateFileInput.get().getName().equals("[[none]]") &&
				!multiStateFileInput.get().exists()) {
			throw new IllegalArgumentException("File " + multiStateFileInput.get().getName() + " does not exist.");
		}
		if (!new File(tempDirInput.get()).exists()) {
			throw new IllegalArgumentException("Directory " + tempDirInput.get() + " does not exist.");
		}
	}

	private void autoConverge(boolean isResuming, int burnIn, String xml2Path) throws IOException, XMLParserException, SAXException, ParserConfigurationException, InterruptedException {
		boolean converged;
		do {
			// prep for next cycle
			long cycleStart = System.currentTimeMillis();
			cycle++;
			isResuming = true;
			multiStateOut = new PrintStream(xml2Path + ".state.multi.tmp");
			multiStateInputFile = xml2Path + ".state.multi";
	        fin = new BufferedReader(new FileReader(multiStateInputFile));

			initLoggers(cycle);
			
			sampleNr = 0; currentSample = 0; currentSampleNr = 0;
			if (nrOfThreads == 1) {
				processUnThreaded(true);
			} else {
				// split work over k threads, with work for (n - burnIn)/k states each
				processThreaded(true);
			}
			
			close(cycle, xml2Path);
			converged = true;
			for (ConvergenceCriterion criterion : criteria) {
				converged = converged && converged(cycle, criterion);
			}
			Log.info("Cycle " + cycle + " done in " + (System.currentTimeMillis()-cycleStart)/1000 + " seconds with " + nrOfThreads + " threads "
					+ "has " + (converged ? "indeed" : "not") + " converged");
			
		} while (cycle != maxCycleInput.get() && !converged);
		
		combine(cycle, xml2Path);		
	}

	private boolean converged(int cycle, ConvergenceCriterion criterion) throws IOException {
		if (!autoConverge) {
			return true;
		}
		if (cycle == 0) {
			// need at least two log files
			return false;
		}
		
		double maxStat = Double.MIN_VALUE, minStat = Double.MAX_VALUE;
		DistributionComparator comparator = new DistributionComparator();
		comparator.setVerbose(!(criterion.equals(ConvergenceCriterion.always) || criterion.equals(ConvergenceCriterion.never)));
		for (Logger logger : loggers) {
			if (!logger.isLoggingToStdout() && logger.mode == Logger.LOGMODE.compound) {
				String fileName1 = getFilename(logger.fileNameInput.get(), cycle);
				String fileName2 = getFilename(logger.fileNameInput.get(), cycle/2);
				double stat = comparator.calcStats(fileName1, fileName2, criterion);
				maxStat = Math.max(maxStat, stat);
				minStat = Math.min(minStat, stat);
			}
 		}
		
		switch (criterion) {
			case always:
				return true;
			case never:
				return false;
			case GR:
			case SplitR:
			case KDE:
			case mean:
			case interval:
				return maxStat < thresholdInput.get();
			case corr:
				return maxStat < thresholdInput.get() && Math.abs(minStat) < thresholdInput.get();
			case KS:
				return minStat > thresholdInput.get();
		}
		return true;
	}

	
	private String getFilename(String fileName, int cycle) {
		if (fileName != null && autoConverge) {
			if (fileName.lastIndexOf(File.separator) > 0) {
				fileName = fileName.substring(fileName.lastIndexOf(File.separator) + 1);
				if (fileName.startsWith("cycle") && fileName.indexOf("_") > 0) {
					fileName = fileName.substring(fileName.indexOf("_") + 1);
				}
			}
			fileName = tempDirInput.get() + File.separator + "cycle" + cycle + "_" + fileName;
		}
		return fileName;
	}

	@SuppressWarnings("unchecked")
	private void initLoggers(int cycle) throws IOException {
		if (loggers == null) {
			Object o = model2.mcmc.getInput("logger").get();
			if (o instanceof List<?>) {
				loggers = new ArrayList<>();
				for (Logger logger : (List<Logger>) o) {
					if (!logger.isLoggingToStdout()) {
						loggers.add(logger);
					}
				}
			} else {
				throw new IllegalArgumentException("Expected list of loggers in XML2");
			}
		}

		for (Logger logger : loggers) {
			logger.everyInput.setValue(1, logger);
			if (autoConverge && !logger.isLoggingToStdout()) {
				String fileName = getFilename(logger.fileNameInput.get(), cycle);
				logger.fileNameInput.setValue(fileName, logger);
				logger.initAndValidate();
				logger.init();
			}
		}
	}

	private void initialise() {
		// set up file mode
		if (overwriteInput.get()) {
			Logger.FILE_MODE = LogFileMode.overwrite;
		} else {
			Logger.FILE_MODE = LogFileMode.only_new_or_exit;
		}
	
		// set up threads
		nrOfThreads = Runtime.getRuntime().availableProcessors();
		if (maxNrOfThreadsInput.get() != null && maxNrOfThreadsInput.get() > 0) {
			nrOfThreads = Math.min(maxNrOfThreadsInput.get(), nrOfThreads);
		}
		if (nrOfThreads > 1) {
		     exec = Executors.newFixedThreadPool(nrOfThreads);
		}

		if (seedInput.get() != null) {
			Randomizer.setSeed(seedInput.get());
		}
		
		// set up criteria and autoConverge
		String [] criterion = criterionInput.get().split(",");
		criteria = new ConvergenceCriterion[criterion.length];
		int i = 0;
		autoConverge = true;
		for (String c : criterion) {
			criteria[i] = ConvergenceCriterion.valueOf(c);
			if (criteria[i] == ConvergenceCriterion.always && criterion.length > 1) {
				throw new IllegalArgumentException("criterion 'always' cannot be combined with other criteria");
			}
			if (criteria[i] == ConvergenceCriterion.never && criterion.length > 1) {
				autoConverge = false;
				throw new IllegalArgumentException("criterion 'never' cannot be combined with other criteria");
			}
		}		
	}

	private boolean importModels() throws IOException, SAXException, ParserConfigurationException, XMLParserException {
		multiStateInputFile = multiStateFileInput.get().getPath();
		if (multiStateInputFile == null || multiStateInputFile.equals("[[none]]")) {
			multiStateInputFile = xml1Input.get().getPath() + ".state.multi"; 
		}
		if (!new File(multiStateInputFile).exists()) {
			throw new IllegalArgumentException("Could not find state file " + multiStateInputFile);
		}
		// import models
		boolean isResuming = false;
		if (xml2Input.get() == null || xml2Input.get().getName().equals("[[none]]")) {
			// if not specified, assume we are resuming
			model2 = getModelFromFile(xml1Input.get());
			isResuming = true;
		} else {
			model2 = getModelFromFile(xml2Input.get());
		}
		multiStateOut = new PrintStream(xml2Input.get().getAbsolutePath() + ".state.multi.tmp");
		return isResuming;
	}

	private void close(int cycle, String xml2Path) throws IOException {
		fin.close();

		for (Logger logger : loggers) {
			logger.close();
		}
					
		Files.move(
				new File(xml2Path + ".state.multi.tmp").toPath(), 
				new File(xml2Path + ".state.multi").toPath(),
				java.nio.file.StandardCopyOption.REPLACE_EXISTING);
	}
	
	private void combine(int cycle, String xml2Path) throws IOException {
		PrintStream stdout = System.out;
		PrintStream stderr = System.err;
        if (!verboseInput.get()) {
			Log.setLevel(Log.Level.error);
			System.setOut(nullStream);
			System.setErr(nullStream);
        }
		if (autoConverge && cycle > 0) {
			for (Logger logger : loggers) {
				if (!logger.isLoggingToStdout()) {
					String from1 = getFilename(logger.fileNameInput.get(), cycle/2);
					String from2 = getFilename(logger.fileNameInput.get(), cycle);
					String to = logger.fileNameInput.get();
					to = to.substring(to.lastIndexOf(File.separator) + 1);
					if (to.startsWith("cycle") && to.indexOf("_") > 0) {
						to = to.substring(to.indexOf("_") + 1);
					}
					Log.err.println("Combining " + from1 + " + " + from2 + " => " + to);
					LogCombiner.main(new String[]{
							"-b", "0", "-log", from1, "-log", from2, "-o", to
					   });
				}
			}
		}
        if (!verboseInput.get()) {
			Log.setLevel(Log.Level.info);
        	System.setOut(stdout);
        	System.setErr(stderr);        	
        }
		
	}

	// skip burn-in states
    private int skipBurnin(boolean isResuming, int n) throws IOException {
        fin = new BufferedReader(new FileReader(multiStateInputFile));
		int burnIn = isResuming ? 0 : burnInPercentageInput.get() * n / 100;
		if (burnIn > 0) {
			Log.info("Skipping " + burnIn + " states as burn-in");
		}
		for (int i = 0; i < burnIn; i++) {
			while (fin.ready()) {
				String str = fin.readLine();
				if (str.startsWith("</itsabeastystatewerein>")) {
					break;
				}
			}
		}
		return burnIn;
	}


	class CoreRunnable implements Runnable {
    	int from; 
    	int to;
    	boolean afterBurnOnly;
    	Model model1, model2;
    	BaseStateExpander expander;

        CoreRunnable(int from, int to, boolean afterBurnOnly) {
    		this.from = from;
    		this.to = to;
    		this.afterBurnOnly = afterBurnOnly;
    		
        	expander = new BaseStateExpander(chainLengthInput.get());
        	

    		try {
	    		this.model1 = getModelFromFile(xml1Input.get());
	    		this.model2 = getModelFromFile(
	    				xml2Input.get() == null || xml2Input.get().getName().equals("[[none]]") ? 
	    				xml1Input.get():
	    				xml2Input.get());
	    		
	    		// restore operator settings, if possible
	    		String stateFile = stateFileInput.get().getPath();
	    		if (stateFile == null || stateFile.equals("[[none]]")) {
	    			stateFile = xml1Input.get().getPath() + ".state";
	    		}		
	    		if (new File(stateFile).exists()) {
	    	        model1.operatorSchedule.setStateFileName(stateFile);
	    	        model1.operatorSchedule.restoreFromFile();	
	    	        model2.operatorSchedule.setStateFileName(stateFile);
	    	        model2.operatorSchedule.restoreFromFile();	
	    			stateFile = xml2Input.get().getPath() + ".state";
	    	        model2.operatorSchedule.setStateFileName(stateFile);
	    		}
    		} catch (IOException | SAXException | ParserConfigurationException | XMLParserException e) {
    			
    		}
        }

        @Override
		public void run() {
            try {
            	for (int i = from; i < to; i++) {
        			String xml = nextState();
        			
        			if (!afterBurnOnly) {
        				model1.state.fromXML(xml);
        				List<String> additions = expander.step1UpdateState(model1, model2);
        				expander.step2OptimiseState(model2, additions);
            			expander.step3RunMCMC(model2);
        			} else {
        				model2.state.fromXML(xml);
            			expander.step3RunMCMC(model2);
        			}

        			logState(model2);
            	}
            	
    	        // release memory
    			expander.optimiser = null;
            } catch (Exception e) {
                Log.err.println("Something went wrong in a calculation of " + from + " " + to + ": " + e.getMessage());
                e.printStackTrace();
                System.exit(1);
            }
            countDown.countDown();
        }
    } // CoreRunnable
    
	private CoreRunnable [] coreRunnable;
	private static   PrintStream nullStream = new PrintStream(new ByteArrayOutputStream() {
		@Override
	    public void write(int b) {
	      // do nothing
	    }

	    @Override
	    public void write(byte[] b, int off, int len) {
	      // do nothing
	    }

	    @Override
	    public void writeTo(OutputStream out) throws IOException {
	      // do nothing
	    }
    });
	
	private void processThreaded(boolean afterBurnOnly)  throws IOException, InterruptedException {
        countDown = new CountDownLatch(nrOfThreads);
        // kick off the threads
        int from = 0;
        int delta = Math.max(1, sampleCount / nrOfThreads); 
        int to = delta;
        if (coreRunnable == null) {
        	coreRunnable = new CoreRunnable[nrOfThreads];
        }
		PrintStream stdout = System.out;
		PrintStream stderr = System.err;
        if (!verboseInput.get()) {
			Log.setLevel(Log.Level.error);
			System.setOut(nullStream);
			System.setErr(nullStream);
        }
        for (int i = 0; i < nrOfThreads; i++) {
        	if (coreRunnable[i] == null) {
        		coreRunnable[i] = new CoreRunnable(from, to, afterBurnOnly);
        	} else {
        		coreRunnable[i].afterBurnOnly = afterBurnOnly;
        	}
            exec.execute(coreRunnable[i]);
            from = to;
            to += delta;
            if (to > sampleCount) {
            	to = sampleCount;
            }
        }
        countDown.await();
        if (!verboseInput.get()) {
    		Log.setLevel(Log.Level.info);
        	System.setOut(stdout);
        	System.setErr(stderr);
        }
        
	}

	private void processUnThreaded(boolean afterBurnOnly)  throws IOException, SAXException, ParserConfigurationException, XMLParserException {
		if (model1 == null) {
			String stateFile = stateFileInput.get().getPath();
			if (stateFile == null || stateFile.equals("[[none]]")) {
				stateFile = xml1Input.get().getPath() + ".state";
			}		

			if (xml2Input.get() == null || xml2Input.get().getName().equals("[[none]]")) {
				model1 = model2;
				if (new File(stateFile).exists()) {
			        model1.operatorSchedule.setStateFileName(stateFile);
			        model1.operatorSchedule.restoreFromFile();	
	    			stateFile = xml2Input.get().getPath() + ".state";
	    	        model2.operatorSchedule.setStateFileName(stateFile);
				}
			} else {
				model1 = getModelFromFile(xml1Input.get());
				if (new File(stateFile).exists()) {
			        model1.operatorSchedule.setStateFileName(stateFile);
			        model1.operatorSchedule.restoreFromFile();	
			        model2.operatorSchedule.setStateFileName(stateFile);
			        model2.operatorSchedule.restoreFromFile();	
	    			stateFile = xml2Input.get().getPath() + ".state";
	    	        model2.operatorSchedule.setStateFileName(stateFile);
				}
			}

		}

        for (int i = 0; i < sampleCount; i++) {
			// get state from file
			String xml = nextState();
			
			if (!afterBurnOnly) {
				model1.state.fromXML(xml);
				List<String> additions = step1UpdateState(model1, model2);
				step2OptimiseState(model2, additions);
				step3RunMCMC(model2);
			} else {
				model2.state.fromXML(xml);
				step3RunMCMC(model2);
			}

			logState(model2);
		}
        
        // release memory
        optimiser = null;
	}

	
	private String prevStartState = null;
	synchronized private String nextState() throws IOException {
		int target = (currentSampleNr + 1) * availableSamples / sampleCount;
		if (target > availableSamples) {
			target = availableSamples;
		}
		if (currentSample == target) {
			if (prevStartState == null) {
				StringBuilder b = new StringBuilder();
				while (fin.ready()) {
					String str = fin.readLine();
					b.append(str);
					if (str.startsWith("</itsabeastystatewerein>")) {
						prevStartState = b.toString();
						currentSampleNr++;
						return b.toString();
					}
				}				
			}
			return prevStartState;
		}
		
		StringBuilder b = new StringBuilder();
		while (fin.ready()) {
			String str = fin.readLine();
			b.append(str);
			if (str.startsWith("</itsabeastystatewerein>")) {
				currentSample++;
				if (currentSample == target) {
					prevStartState = b.toString();
					currentSampleNr++;
					return b.toString();
				} else {
					b.delete(0, b.length());
				}
			}
		}
		if (prevStartState != null) {
			Log.warning("Ran out of states in state file ");
			Log.warning("target = " + target);
			Log.warning("currentSampleNr = " + currentSampleNr);
			Log.warning("availableSamples = " + availableSamples);
			Log.warning("sampleCount = " + sampleCount);
			return prevStartState;
		}
		throw new RuntimeException("Ran out of states in state file");
	}

	/** counts number of states in state file **/
	private int getStateCount() throws IOException {
        BufferedReader fin = new BufferedReader(new FileReader(multiStateInputFile));
        String str = null;
        int n = 0;
        while (fin.ready()) {
            str = fin.readLine();
            if (str.startsWith("<itsabeastystatewerein")) {
            	n++;
            }
        }
        fin.close();
		return n;
	}

	synchronized private void logState(Model model) throws IOException {
		
		
		State other = model.mcmc2.startStateInput.get();
		PrintStream sf = new PrintStream(new File(xml2Input.get().getPath()+".state"));
		sf.println("<itsabeastystatewerein>");
		sf.println("</itsabeastystatewerein>");
		sf.close();
		// TODO: set operator schedule
		model.operatorSchedule.storeToFile();

		State state = model2.state;
		if (state != other) {
			for (int i = 0; i < state.getNrOfStateNodes(); i++) {
				StateNode s1 = state.getStateNode(i);
				StateNode s2 = other.getStateNode(i);
				s1.assignFrom(s2);
			}
		}
		
		Distribution p = model2.mcmc.posteriorInput.get();
//		double logP2 = 
		model2.state.robustlyCalcPosterior(p);
		
//		double logP = p.getCurrentLogP();
//		System.err.println(logP + " - " + logP2 + " = " + (logP - logP2));

		for (Logger logger : loggers) {
			if (logger.getFileName() != null) {
				logger.log(sampleNr);
			}
		}
		
		if (multiStateOut != null) {
			multiStateOut.println(other.toXML(sampleNr));
		}
		sampleNr++;
		
		// print progress bar:
		if (sampleNr == 1) {
			Log.err.print("Cycle " + cycle + ": ");
		}
		if (sampleNr % 1 == 0) {
			if (sampleNr % 10 == 0) {
				Log.err.print("|");
			} else {
				Log.err.print(".");
			}
		}
	} // logState
		
	public static void main(String[] args) throws Exception {
		new Application(new TraceExpander(), "Trace Expander", args);
	}

} // class TraceExpander
