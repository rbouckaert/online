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
	public final Input<Double> maxRInput = new Input<>("maxR", "maximum acceptable value of Gelman Rubin statistic."
			+ "The chain keeps afterburning (with chainLength steps) till all items in trace log converge. "
			+ "Set to less than 1 to stop after first cycle.", 1.05);
    
    
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
		
		
		boolean isResuming = importModels();

		int n = getStateCount();
		int burnIn = skipBurnin(isResuming, n);
		
		
		initLoggers(0);
		
		
		sampleNr = 0;
		if (nrOfThreads == 1) {
			processUnThreaded(n - burnIn);
		} else {
			// split work over k threads, with work for (n - burnIn)/k states each
			processThreaded(n - burnIn);
		}
		
		close(isResuming);
		
		Long end = System.currentTimeMillis();
		Log.info("Done in " + (end-start)/1000 + " seconds with " + nrOfThreads + " threads");
		System.exit(0);
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
			String fileName = logger.fileNameInput.get();
			if (fileName != null) {
				if (fileName.indexOf(File.pathSeparator) > 0) {
					fileName = fileName.substring(fileName.lastIndexOf(File.pathSeparator) + 1);
					fileName = tempDirInput.get() + File.pathSeparator + "cycle" + cycle + "_" + fileName;
					logger.fileNameInput.setValue(fileName, logger);
				}
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

	private void close(boolean isResuming) throws IOException {
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

		exec.shutdownNow();
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
