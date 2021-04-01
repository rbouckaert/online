package online;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
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
import beast.util.LogAnalyser;
import beast.util.Randomizer;
import beast.util.XMLParserException;

@Description("Create tree and trace files extending an input multiple-state file with different set of taxa")
public class TraceExpander extends BaseStateExpander {
	final public Input<XMLFile> xml2Input = new Input<>("xml2", "BEAST XML file with expanded state. If not specified, xml2 is set equal to xml1, and we are resuming", new XMLFile("[[none]]"));
	final public Input<File> multiStateFileInput = new Input<>("multiStateFile", "state file containing multiple states associated with initial XML file (xml1)."
			+ "If not specified, use xml1+\".state.multi\"", new File("[[none]]"));
	
	final public Input<File> stateFileInput = new Input<>("stateFile", "state file containing operarator optimisation information associated with initial XML file (xml1). "
			+ "If not specified, use xml1+\".state\". If that does not exist, ignore.", new File("[[none]]"));

    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads","maximum number of threads to use, if not specified the number of available cores is used (default)");
    final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);
    final public Input<Boolean> overwriteInput = new Input<>("overwrite", "overwrite existing tree and trace files", true);

	public final Input<String> tempDirInput = new Input<>("tempDir","directory where temporary files are written", "/tmp/");
	/**
	 * https://www.stata.com/new-in-stata/gelman-rubin-convergence-diagnostic/
	 * " Gelman and Rubin (1992) and Brooks and Gelman (1998) suggest that diagnostic Rc values greater than 1.2 for 
	 *   any of the model parameters should indicate nonconvergence. In practice, a more stringent rule of Rc < 1.1 is 
	 *   often used to declare convergence."
	 *   But even 1.1 is often too low (https://arxiv.org/pdf/1812.09384.pdf) and 1.01 makes more sense. 
	 */
	public final Input<Double> maxRInput = new Input<>("maxR", "maximum acceptable value of Gelman Rubin statistic."
			+ "The chain keeps afterburning (with chainLength steps) till all items in trace log converge. "
			+ "Set to less than 1 to stop after first cycle.", 1.01);
    final public Input<Boolean> useSplitRInput = new Input<>("useSplitR", "use split-R estimate instead of original Gelman-Rubin statistic", true);
    
    
    private int nrOfThreads;
    private static ExecutorService exec;
    private static CountDownLatch countDown;
	private List<Logger> loggers;
	private Model model2;
	private long sampleNr;
	private BufferedReader fin;
	private String multiStateFile;
	private PrintStream multiStateOut;
	private boolean autoConverge;

	@Override
	public void initAndValidate() {
	}
	
	@Override
	public void run() throws Exception {
		Long start = System.currentTimeMillis();

//		Log.setLevel(Log.Level.debug);
		initialise();
		
		
		int cycle = 0;
		do {
			boolean isResuming = importModels();
	
			int n = getStateCount();
			int burnIn = skipBurnin(isResuming, n);
			
			initLoggers(cycle);
			
			sampleNr = 0;
			if (nrOfThreads == 1) {
				processUnThreaded(n - burnIn);
			} else {
				// split work over k threads, with work for (n - burnIn)/k states each
				processThreaded(n - burnIn);
			}
			
			close(isResuming, cycle);
			
			// prep for next cycle
			cycle++;
			xml2Input.setValue(null, this);
		} while (!converged(cycle));
		
		Long end = System.currentTimeMillis();
		Log.info("Done in " + (end-start)/1000 + " seconds with " + nrOfThreads + " threads");
		System.exit(0);
	}

	
	private boolean converged(int cycle) throws IOException {
		if (!autoConverge) {
			return true;
		}
		if (cycle == 1) {
			// need at least two log files
			return false;
		}
		
		double maxR = Double.MIN_VALUE;
		for (Logger logger : loggers) {
			if (!logger.isLoggingToStdout() && logger.mode == Logger.LOGMODE.compound) {
				String fileName1 = getFilename(logger.fileNameInput.get(), cycle-1);
				String fileName2 = getFilename(logger.fileNameInput.get(), cycle-2);
				maxR = Math.max(maxR, calcGRStats(fileName1, fileName2));
			}
 		}
		
		return maxR < maxRInput.get();
	}

	private double calcGRStats(String fileName1, String fileName2) throws IOException {
		Double [][] trace1 = new LogAnalyser(fileName1, 0, true, false).getTraces();
		Double [][] trace2 = new LogAnalyser(fileName2, 0, true, false).getTraces();
		double maxR = Double.MIN_VALUE;
		for (int i = 1; i < trace1.length; i++) {
			double R = useSplitRInput.get()? 
					calcSplitGRStat(trace1[i], trace2[i]):
					calcGRStat(trace1[i], trace2[i]);
			maxR = Math.max(maxR, R);
		}
		return maxR;
	}

	/** original Gelman Rubin statistic for 2 chains **/
	private double calcGRStat(Double[] trace1, Double[] trace2) {
		double sampleCount = trace1.length;
		if (trace2.length != sampleCount) {
			throw new IllegalArgumentException("Expected traces of the same length");
		}
		
		// calc means and squared means
		double mean1 = 0, mean2 = 0, sumsq1 = 0, sumsq2 = 0;
		for (Double d : trace1) {
			mean1 += d;
			sumsq1 += d * d;
		}
		mean1 /= sampleCount;
		for (Double d : trace2) {
			mean2 += d;
			sumsq2 += d * d;
		}
		mean2 /= sampleCount;

		// calculate variances for both chains
		double var1 = (sumsq1 - mean1 * mean1 * sampleCount)/(sampleCount - 1);
		double var2 = (sumsq2 - mean2 * mean2 * sampleCount)/(sampleCount - 1);
		
		// average variance for this item
		double fW = (var1 + var2) / 2;

		// sum to get totals
		double totalMean = (mean1 + mean2) / 2;
		double totalSq = mean1*mean1 + mean2*mean2;
		
		// variance for joint
		double fB = (totalSq - totalMean * totalMean * 2);
		
		
		double varR = ((sampleCount - 1.0)/sampleCount) + (fB/fW)*(1.0/sampleCount);
		double R = Math.sqrt(varR);
		return R;
	}

	/** Split-R, following 
	 * Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A. and Rubin, D.B.. 
	 * Bayesian data analysis. CRC press. 2013.
	 */
	private double calcSplitGRStat(Double[] trace1, Double[] trace2) {
		if (trace2.length != trace1.length) {
			throw new IllegalArgumentException("Expected traces of the same length");
		}
		int sampleCount = trace1.length/2;
		int sampleCountb = trace1.length - sampleCount;
		
		// calc means and squared means
		double mean1a = 0, mean2a = 0, sumsq1a = 0, sumsq2a = 0;
		double mean1b = 0, mean2b = 0, sumsq1b = 0, sumsq2b = 0;
		for (int i = 0; i < sampleCount; i++) {
			final double d = trace1[i];
			mean1a += d;
			sumsq1a += d * d;
			final double d2 = trace2[i];
			mean2a += d2;
			sumsq2a += d2 * d2;
		}
		mean1a /= sampleCount;
		mean2a /= sampleCount;
		for (int i = sampleCount; i < trace1.length; i++) {
			final double d = trace1[i];
			mean1b += d;
			sumsq1b += d * d;
			final double d2 = trace2[i];
			mean2b += d2;
			sumsq2b += d2 * d2;
		}
		mean1b /= sampleCountb;
		mean2b /= sampleCountb;
		

		// calculate variances for both chains
		double var1a = (sumsq1a - mean1a * mean1a * sampleCount)/(sampleCount - 1);
		double var2a = (sumsq2a - mean2a * mean2a * sampleCount)/(sampleCount - 1);
		double var1b = (sumsq1b - mean1b * mean1b * sampleCountb)/(sampleCountb - 1);
		double var2b = (sumsq2b - mean2b * mean2b * sampleCountb)/(sampleCountb - 1);
		
		// average variance for this item
		double fW = (var1a + var2a + var1b + var2b) / 4;

		// sum to get totals
		double totalMean = (mean1a + mean2a + mean1b + mean2b) / 4;
		double totalSq = mean1a*mean1a + mean2a*mean2a + mean1b*mean1b + mean2b*mean2b;
		
		// variance for joint
		double fB = (totalSq - totalMean * totalMean * 4)/3.0;
		
		
		double varR = ((sampleCount - 1.0)/sampleCount) + (fB/fW)*(1.0/sampleCount);
		double R = Math.sqrt(varR);
		return R;
	}
	
	private String getFilename(String fileName, int cycle) {
		if (fileName != null && autoConverge) {
			if (fileName.indexOf(File.pathSeparator) > 0) {
				fileName = fileName.substring(fileName.lastIndexOf(File.pathSeparator) + 1);
				if (fileName.startsWith("cycle") && fileName.indexOf("_") > 0) {
					fileName = fileName.substring(fileName.indexOf("_") + 1);
				}
				fileName = tempDirInput.get() + File.pathSeparator + "cycle" + cycle + "_" + fileName;
			}
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
		
		autoConverge = maxRInput.get() > 1.0;
	}

	private boolean importModels() throws IOException, SAXException, ParserConfigurationException, XMLParserException {
		multiStateFile = multiStateFileInput.get().getPath();
		if (multiStateFile == null || multiStateFile.equals("[[none]]")) {
			multiStateFile = xml1Input.get().getPath() + ".state.multi"; 
		}
		if (!new File(multiStateFile).exists()) {
			throw new IllegalArgumentException("Could not find state file " + multiStateFile);
		}
		// import models
		boolean isResuming = false;
		if (xml2Input.get() == null || xml2Input.get().getName().equals("[[none]]")) {
			// if not specified, assume we are resuming
			model2 = getModelFromFile(xml1Input.get());
			multiStateOut = new PrintStream(xml1Input.get().getAbsolutePath() + ".state.multi.tmp");
			isResuming = true;
		} else {
			model2 = getModelFromFile(xml2Input.get());
			multiStateOut = new PrintStream(xml2Input.get().getAbsolutePath() + ".state.multi");
		}
		return isResuming;
	}

	private void close(boolean isResuming, int cycle) throws IOException {
		fin.close();

		for (Logger logger : loggers) {
			logger.close();
		}
		
		if (isResuming) {
			Files.move(
					new File(xml1Input.get().getAbsolutePath() + ".state.multi.tmp").toPath(), 
					new File(xml1Input.get().getAbsolutePath() + ".state.multi").toPath(),
					java.nio.file.StandardCopyOption.REPLACE_EXISTING);
		}
		
		if (autoConverge) {
			for (Logger logger : loggers) {
				if (!logger.isLoggingToStdout()) {
					String from1 = getFilename(logger.fileNameInput.get(), cycle-1);
					String from2 = getFilename(logger.fileNameInput.get(), cycle);
					String to = logger.fileNameInput.get();
					to = to.substring(to.lastIndexOf(File.pathSeparator) + 1);
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
        fin = new BufferedReader(new FileReader(multiStateFile));
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

        CoreRunnable(int from, int to) {
    		this.from = from;
    		this.to = to;
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
        			model1.state.fromXML(xml);
        			expander.updateState(model1, model2);
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
    
	private void processThreaded(int n)  throws IOException, InterruptedException {
        countDown = new CountDownLatch(nrOfThreads);
        // kick off the threads
        int from = 0;
        int delta = Math.max(1, n / nrOfThreads); 
        int to = delta;
        for (int i = 0; i < nrOfThreads; i++) {
            CoreRunnable coreRunnable = new CoreRunnable(from, to);
            exec.execute(coreRunnable);
            from = to;
            to += delta;
            if (to > n) {
            	to = n;
            }
        }
        countDown.await();
	}

	private void processUnThreaded(int n)  throws IOException, SAXException, ParserConfigurationException, XMLParserException {
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
			model1.state.fromXML(xml);
			updateState(model1, model2);

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
        BufferedReader fin = new BufferedReader(new FileReader(multiStateFile));
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
		if (sampleNr % 10 == 0) {
			if (sampleNr % 50 == 0) {
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
