package online;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.app.tools.LogCombiner;
import beast.app.util.Application;
import beast.app.util.XMLFile;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Logger;
import beast.core.Logger.LogFileMode;
import beast.core.State;
import beast.core.StateNode;
import beast.core.util.Log;
import beast.util.Randomizer;
import beast.util.XMLParserException;
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

	/**
	 * https://www.stata.com/new-in-stata/gelman-rubin-convergence-diagnostic/
	 * " Gelman and Rubin (1992) and Brooks and Gelman (1998) suggest that diagnostic Rc values greater than 1.2 for 
	 *   any of the model parameters should indicate nonconvergence. In practice, a more stringent rule of Rc < 1.1 is 
	 *   often used to declare convergence."
	 *   But even 1.1 is often too low (https://arxiv.org/pdf/1812.09384.pdf) and 1.01 makes more sense. 
	 */
    final public Input<Double> thresholdInput = new Input<>("threshold", "threshold appropriate for convergence criterion, "
    		+ "e.g. maximum acceptable value of Gelman Rubin statistic, or minimum p-value for KS test. "
			+ "Set 'criterion' to 'none' to stop after first cycle.", 1.01);
    
	final public Input<String> tempDirInput = new Input<>("tempDir","directory where temporary files are written."
			+ "(Ignored if maxRInput < 1).", "/tmp/");
    final public Input<ConvergenceCriterion> criterionInput = new Input<>("criterion",
			"If not set to 'none', the chain keeps afterburning (with chainLength steps) till all items in trace log converge. " +
    		DistributionComparator.convergenceCriterionDescription, 
    		ConvergenceCriterion.KDE, ConvergenceCriterion.values());
    final public Input<Integer> maxCycleInput = new Input<>("maxCycle", "maximum number of cycles before stopping. Ignored if negative (which is the default)", -1);
    
    
    private int nrOfThreads;
    private static ExecutorService exec;
    private static CountDownLatch countDown;
	private List<Logger> loggers;
	private Model model2;
	private long sampleNr;
	private BufferedReader fin;
	private String multiStateInputFile;
	private PrintStream multiStateOut;
	private boolean autoConverge;
	private ConvergenceCriterion criterion;
	private int cycle;

	@Override
	public void initAndValidate() {
	}
	
	@Override
	public void run() throws Exception {
		Long start = System.currentTimeMillis();
		criterion = criterionInput.get();
		
//		Log.setLevel(Log.Level.debug);
		initialise();
		
		
		cycle = 0;
		boolean isResuming = importModels();
		int n = getStateCount();
		int burnIn = skipBurnin(isResuming, n);
		
		initLoggers(cycle);
		
		sampleNr = 0;
		if (nrOfThreads == 1) {
			processUnThreaded(n - burnIn, false);
		} else {
			// split work over k threads, with work for (n - burnIn)/k states each
			processThreaded(n - burnIn, false);
		}
		
		String xml2Path = xml2Input.get().getAbsolutePath();
		close(cycle, xml2Path);

		Log.info("Cycle " + cycle + " done in " + (System.currentTimeMillis()-start)/1000 + " seconds with " + nrOfThreads + " threads");

		
		do {
			// prep for next cycle
			long cycleStart = System.currentTimeMillis();
			cycle++;
			isResuming = true;
			multiStateOut = new PrintStream(xml2Path + ".state.multi.tmp");
			multiStateInputFile = xml2Path + ".state.multi";
	        fin = new BufferedReader(new FileReader(multiStateInputFile));

			initLoggers(cycle);
			
			sampleNr = 0;
			if (nrOfThreads == 1) {
				processUnThreaded(n - burnIn, true);
			} else {
				// split work over k threads, with work for (n - burnIn)/k states each
				processThreaded(n - burnIn, true);
			}
			
			close(cycle, xml2Path);
			Log.info("Cycle " + cycle + " done in " + (System.currentTimeMillis()-cycleStart)/1000 + " seconds with " + nrOfThreads + " threads");
		} while (cycle != maxCycleInput.get() && !converged(cycle));
		
		Long end = System.currentTimeMillis();
		Log.info("Done in " + (end-start)/1000 + " seconds with " + nrOfThreads + " threads");
		System.exit(0);
	}

	
	private boolean converged(int cycle) throws IOException {
		if (!autoConverge) {
			return true;
		}
		if (cycle == 0) {
			// need at least two log files
			return false;
		}
		
		double maxStat = Double.MIN_VALUE, minStat = Double.MAX_VALUE;
		DistributionComparator comparator = new DistributionComparator();
		comparator.setVerbose(true);
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
			case none:
				return true;
			case GR:
			case SplitR:
			case KDE:
			case mean:
			case corr:
				return maxStat < thresholdInput.get();
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
		Object o = model2.mcmc.getInput("logger").get();
		if (o instanceof List<?>) {
			loggers = (List<Logger>) o;
		} else {
			throw new IllegalArgumentException("Expected list of loggers in XML2");
		}

		for (Logger logger : loggers) {
			logger.everyInput.setValue(1, logger);
			if (autoConverge && !logger.isLoggingToStdout()) {
				String fileName = getFilename(logger.fileNameInput.get(), cycle);
				logger.fileNameInput.setValue(fileName, logger);
			}
			logger.initAndValidate();
			logger.init();
		}
	}

	private void initialise() {
		if (overwriteInput.get()) {
			Logger.FILE_MODE = LogFileMode.overwrite;
		} else {
			Logger.FILE_MODE = LogFileMode.only_new_or_exit;
		}
	
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
		
		autoConverge = !criterionInput.get().equals(ConvergenceCriterion.none);
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
					LogCombiner.main(new String[]{
							"-b", "0", "-log", from1, "-log", from2, "-o", to
					   });
				}
			}
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
			nextState();
		}
		return burnIn;
	}


	class CoreRunnable implements Runnable {
    	int from; 
    	int to;
    	boolean afterBurnOnly;

        CoreRunnable(int from, int to, boolean afterBurnOnly) {
    		this.from = from;
    		this.to = to;
    		this.afterBurnOnly = afterBurnOnly;
        }

        @Override
		public void run() {
            try {
            	BaseStateExpander expander = new BaseStateExpander(chainLengthInput.get());
        		Model model1 = getModelFromFile(xml1Input.get());
        		Model model2 = getModelFromFile(
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
        		}

            	for (int i = from; i < to; i++) {
        			String xml = nextState();
        			if (!afterBurnOnly) {
        				model1.state.fromXML(xml);
        				expander.updateState(model1, model2);
        			} else {
        				model2.state.fromXML(xml);
        				expander.afterBurner(model2, new ArrayList<>(), 0.0);
        			}
        			logState(model2.state);
            	}
            } catch (Exception e) {
                Log.err.println("Something went wrong in a calculation of " + from + " " + to);
                e.printStackTrace();
                System.exit(1);
            }
            countDown.countDown();
        }
    } // CoreRunnable
    
	private void processThreaded(int n, boolean afterBurnOnly)  throws IOException, InterruptedException {
        countDown = new CountDownLatch(nrOfThreads);
        // kick off the threads
        int from = 0;
        int delta = Math.max(1, n / nrOfThreads); 
        int to = delta;
        for (int i = 0; i < nrOfThreads; i++) {
            CoreRunnable coreRunnable = new CoreRunnable(from, to, afterBurnOnly);
            exec.execute(coreRunnable);
            from = to;
            to += delta;
            if (to > n) {
            	to = n;
            }
        }
        countDown.await();
	}

	private void processUnThreaded(int n, boolean afterBurnOnly)  throws IOException, SAXException, ParserConfigurationException, XMLParserException {
		Model model1 = getModelFromFile(xml1Input.get());
		
		String stateFile = stateFileInput.get().getPath();
		if (stateFile == null || stateFile.equals("[[none]]")) {
			stateFile = xml1Input.get().getPath() + ".state";
		}		
		if (new File(stateFile).exists()) {
	        model1.operatorSchedule.setStateFileName(stateFile);
	        model1.operatorSchedule.restoreFromFile();	
	        model2.operatorSchedule.setStateFileName(stateFile);
	        model2.operatorSchedule.restoreFromFile();	
		}
		

        for (int i = 0; i < n; i++) {
			// get state from file
			String xml = nextState();
			
			if (!afterBurnOnly) {
				model1.state.fromXML(xml);
				updateState(model1, model2);
			} else {
				model2.state.fromXML(xml);
				afterBurner(model2, new ArrayList<>(), 0);
			}

//			Distribution p = model2.mcmc.posteriorInput.get();
//			double logP2 = model2.state.robustlyCalcPosterior(p);
			logState(model2.state);

//			double logP = p.getCurrentLogP();
//			System.err.println(logP + " - " + logP2 + " = " + (logP - logP2));
		}
	}

	synchronized private String nextState() throws IOException {
		StringBuilder b = new StringBuilder();
		while (fin.ready()) {
			String str = fin.readLine();
			b.append(str);
			if (str.startsWith("</itsabeastystatewerein>")) {
				return b.toString();
			}
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

	synchronized private void logState(State other) throws IOException {
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
		if (sampleNr == 0) {
			System.err.print("Cycle " + cycle + ": ");
		}
		if (sampleNr % 1 == 0) {
			if (sampleNr % 10 == 0) {
				System.err.print("|");
			} else {
				System.err.print(".");
			}
		}
	} // logState
		
	public static void main(String[] args) throws Exception {
		new Application(new TraceExpander(), "Trace Expander", args);
	}

} // class TraceExpander
