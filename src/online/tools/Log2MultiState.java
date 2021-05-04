package online.tools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

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
import beast.core.StateNode;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
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

	
	private void addState(State state, StringBuilder buf, LogAnalyser trace, int i, Tree tree) {
		buf.append("<itsabeastystatewerein version='2.0' sample='0'>\n");
		for (StateNode node : state.stateNodeInput.get()) {
			if (node instanceof Tree) {
				((Tree) node).assignFrom(tree);
			} else if (node instanceof RealParameter) {
				getRealParameterValues((RealParameter) node, trace, i);
			} else if (node instanceof IntegerParameter) {
				getIntegerParameterValues((IntegerParameter) node, trace, i);
			} else if (node instanceof BooleanParameter) {
				getBooleanParameterValues((BooleanParameter) node, trace, i);
			} else {
				Log.warning("Don't know how to retrieve state nodes of type " + node.getClass().getName());
			}
			
			buf.append(node.toXML());
			buf.append("\n");
		}
		buf.append("</itsabeastystatewerein>\n");
	}

	
	
	private void getRealParameterValues(RealParameter param, LogAnalyser trace, int i) {
		String id = param.getID();
		List<String> labels = trace.getLabels();
		if (param.getDimension() == 1) {
			int k = indexOf(id, labels);
			if (k < 0) {
				throw new IllegalArgumentException("cannot find " + id + " in trace file");
			}
			double value = trace.getTrace(k)[i];
			param.setValue(value);
			return;
		} else {
			int k = indexOf(id+"1", labels);
			if (k < 0 && id.lastIndexOf(".") > 0) {
				id = id.substring(id.lastIndexOf("."));
				k = indexOf(id+"1", labels);
				if (k < 0) {
					throw new IllegalArgumentException("cannot find " + id + " in trace file");
				}
			}
			for (int j = 0; j < param.getDimension(); j++) {
				double value = trace.getTrace(k+j)[i];
				param.setValue(j, value);
			}
		}
	}

	private int indexOf(String id, List<String> labels) {
		for (int i = 0; i < labels.size(); i++) {
			if (labels.get(i).equals(id)) {
				return i;
			}
		}
		if (id.lastIndexOf(".") > 0) {
			id = id.substring(id.lastIndexOf("."));
			return indexOf(id, labels);
		}
		return -1;
	}

	private void getIntegerParameterValues(IntegerParameter param, LogAnalyser trace, int i) {
		String id = param.getID();
		List<String> labels = trace.getLabels();
		if (param.getDimension() == 1) {
			int k = indexOf(id, labels);
			if (k < 0) {
				throw new IllegalArgumentException("cannot find " + id + " in trace file");
			}
			int value = (int) Math.round(trace.getTrace(k)[i]);
			param.setValue(value);
			return;
		} else {
			int k = indexOf(id+"1", labels);
			if (k < 0 && id.lastIndexOf(".") > 0) {
				id = id.substring(id.lastIndexOf("."));
				k = indexOf(id+"1", labels);
				if (k < 0) {
					throw new IllegalArgumentException("cannot find " + id + " in trace file");
				}
			}
			for (int j = 0; j < param.getDimension(); j++) {
				int value = (int) Math.round(trace.getTrace(k+j)[i]);
				param.setValue(j, value);
			}
		}
	}

	private void getBooleanParameterValues(BooleanParameter param, LogAnalyser trace, int i) {
		String id = param.getID();
		List<String> labels = trace.getLabels();
		if (param.getDimension() == 1) {
			int k = indexOf(id, labels);
			if (k < 0) {
				throw new IllegalArgumentException("cannot find " + id + " in trace file");
			}
			boolean value = trace.getTrace(k)[i] != 0;
			param.setValue(value);
			return;
		} else {
			int k = indexOf(id+"1", labels);
			if (k < 0 && id.lastIndexOf(".") > 0) {
				id = id.substring(id.lastIndexOf("."));
				k = indexOf(id+"1", labels);
				if (k < 0) {
					throw new IllegalArgumentException("cannot find " + id + " in trace file");
				}
			}
			for (int j = 0; j < param.getDimension(); j++) {
				boolean value = trace.getTrace(k+j)[i] != 0;
				param.setValue(j, value);
			}
		}
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
