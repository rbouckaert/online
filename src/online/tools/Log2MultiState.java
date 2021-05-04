package online.tools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import beast.app.treeannotator.TreeAnnotator;
import beast.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beast.app.util.Application;
import beast.app.util.LogFile;
import beast.app.util.OutFile;
import beast.app.util.TreeFile;
import beast.app.util.XMLFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.Runnable;
import beast.core.State;
import beast.evolution.tree.Tree;
import beast.util.LogAnalyser;
import beast.util.XMLParser;

@Description("Convert tree and trace log into multi-state file that can be used with TraceExpander")
public class Log2MultiState extends Runnable {
	final public Input<XMLFile> xmlInput = new Input<>("xml","BEAST XML file with loggers", new XMLFile("[[none]]"));
	final public Input<LogFile> logInput = new Input<>("log","Trace log for XML", new LogFile("[[none]]"));
	final public Input<TreeFile> treeInput = new Input<>("tree","Tree log for XML", new TreeFile("[[none]]"));
	final public Input<OutFile> multiStateFileInput = new Input<>("multiStateFile", "output multi state file containing multiple states associated with initial XML file."
			+ "If not specified, use xml+\".state.multi\"", new OutFile("[[none]]"));

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		sanityCheck();

		LogAnalyser trace = new LogAnalyser(logInput.get().getPath(), 0, true, false);
		MemoryFriendlyTreeSet treeSet = new TreeAnnotator().new MemoryFriendlyTreeSet(treeInput.get().getPath(), 0);
		treeSet.reset();
		
		XMLParser parser = new XMLParser();
		Runnable run = parser.parseFile(xmlInput.get());
		if (! (run instanceof MCMC)) {
			throw new IllegalArgumentException("Expected MCMC analysis in xml file");
		}
		MCMC mcmc = (MCMC) run;
		State state = mcmc.startStateInput.get();

		StringBuilder multiStates = new StringBuilder();
		for (int i = 0; i < trace.getTrace(0).length; i++) {
			Tree tree = treeSet.next();
			addState(state, multiStates, trace, i, tree);
		}
		
		output(multiStates);
	}

	
	private void addState(State state, StringBuilder multiStates, LogAnalyser trace, int i, Tree tree) {
		// TODO: do the hard work here
	}

	
	
	private void output(StringBuilder multiStates) throws IOException {
		String multiStateInputFile = multiStateFileInput.get().getPath();
		if (multiStateInputFile == null || multiStateInputFile.equals("[[none]]")) {
			multiStateInputFile = xmlInput.get().getPath() + ".state.multi"; 
		}
        FileWriter outfile = new FileWriter(new File(multiStateInputFile));
        outfile.write(multiStates.toString());
        outfile.close();
	}


	private void sanityCheck() {
		if (xmlInput.get() != null && !xmlInput.get().exists()) {
			throw new IllegalArgumentException("Could not find XML file " + xmlInput.get().getName());
		}
		if (logInput.get() != null && !logInput.get().exists()) {
			throw new IllegalArgumentException("Could not find trace log file " + logInput.get().getName());
		}
		if (treeInput.get() != null && !treeInput.get().exists()) {
			throw new IllegalArgumentException("Could not find tree log file " + treeInput.get().getName());
		}
		if (multiStateFileInput.get() != null && !(multiStateFileInput.get().getName().equals("[[none]]")) &&
				!multiStateFileInput.get().exists()) {
			throw new IllegalArgumentException("Could not find tree log file " + treeInput.get().getName());
		}
	}

	public static void main(String[] args) throws Exception {
		new Application(new MultiState2Log(), "Log to Multi State", args);
	}
}
