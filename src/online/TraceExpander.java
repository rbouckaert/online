package online;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.app.util.Application;
import beast.core.Description;
import beast.core.Input;
import beast.core.Logger;
import beast.core.State;
import beast.core.StateNode;
import beast.core.util.Log;
import beast.util.Randomizer;
import beast.util.XMLParserException;

@Description("Create tree and trace files extending an input multiple-state file with different set of taxa")
public class TraceExpander extends BaseStateExpander {
	final public Input<File> stateInput = new Input<>("state", "state file containing multiple states associated with initial XML file (xml1)", new File("[[none]]"));
    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads","maximum number of threads to use, if less than 1 the number of threads in BeastMCMC is used (default -1)", -1);
	
    private int nrOfThreads;
    private static ExecutorService exec;
    private static CountDownLatch countDown;
	private List<Logger> loggers;
	private Model model2;
	private long sampleNr;
	private BufferedReader fin;
	

	@Override
	public void initAndValidate() {
	}
	
	@Override
	public void run() throws Exception {
		Log.setLevel(Log.Level.debug);
		nrOfThreads = Runtime.getRuntime().availableProcessors();
		if (maxNrOfThreadsInput.get() > 0) {
			nrOfThreads = Math.min(maxNrOfThreadsInput.get(), nrOfThreads);
		}
		if (nrOfThreads > 1) {
		     exec = Executors.newFixedThreadPool(nrOfThreads);
		}

		if (seedInput.get() != null) {
			Randomizer.setSeed(seedInput.get());
		}

		// import models
		model2 = getModelFromFile(xml2Input.get());
		Object o = model2.mcmc.getInput("logger");
		if (o instanceof List<?>) {
			loggers = (List<Logger>) o;
		} else {
			throw new IllegalArgumentException("Expected list of loggers in XML2");
		}
		for (Logger logger : loggers) {
			logger.everyInput.setValue(1, logger);
			logger.initAndValidate();
			logger.init();
		}

		int n = getStateCount();
		sampleNr = 0;
		
		
        BufferedReader fin = new BufferedReader(new FileReader(stateInput.get()));
		if (nrOfThreads == 1) {
			processUnThreaded(n);
		} else {
			// split work over k threads, with work for n/k states each
			processThreaded(n);
		}
		fin.close();

		for (Logger logger : loggers) {
			logger.close();
		}
		Log.debug("Done!");
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
        		Model model1 = getModelFromFile(xml1Input.get());
        		Model model2 = getModelFromFile(xml2Input.get());
                model1.operatorSchedule.setStateFileName(stateInput.get().getAbsolutePath());
                model1.operatorSchedule.restoreFromFile();	
                model2.operatorSchedule.setStateFileName(stateInput.get().getAbsolutePath());
                model2.operatorSchedule.restoreFromFile();	

            	for (int i = from; i < to; i++) {
        			String xml = nextState();
        			model1.state.fromXML(xml);
        			updateState(model1, model2);
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
        model1.operatorSchedule.setStateFileName(stateInput.get().getAbsolutePath());
        model1.operatorSchedule.restoreFromFile();	
        model2.operatorSchedule.setStateFileName(stateInput.get().getAbsolutePath());
        model2.operatorSchedule.restoreFromFile();	

        for (int i = 0; i < n; i++) {
			// get state from file
			String xml = nextState();
			model1.state.fromXML(xml);
			updateState(model1, model2);
			logState(model2.state);
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
        BufferedReader fin = new BufferedReader(new FileReader(stateInput.get()));
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
		
		for (Logger logger : loggers) {
			logger.log(sampleNr);
		}
		sampleNr++;
	} // logState
		
	public static void main(String[] args) throws Exception {
		new Application(new TraceExpander(), "Trace Expander", args);
	}

} // class TraceExpander
