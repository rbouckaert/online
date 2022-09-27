package online.tools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.MemoryFriendlyTreeSet;
import beastfx.app.tools.Application;
import beastfx.app.util.LogFile;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import beastfx.app.util.XMLFile;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.MCMC;
import beast.base.inference.Runnable;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beastfx.app.tools.LogAnalyser;
import beast.base.parser.XMLParser;

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
			if (tree == null) {
				Log.warning("\nWARNING: Log file (" + trace.getTrace(0).length + ") and tree files (" + i + ") not the same length");
				Log.warning("WARNING: so there is a possible mismatch between trace and tree file entries.");
				Log.warning("WARNING: Use logcombiner with the `-resample` option to synchronise trace and tree logs.\n");
				break;
			}
			addState(state, multiStates, trace, i, tree);
		}
		
		output(multiStates);
		Log.warning("Done");
	}

	
	private void addState(State state, StringBuilder buf, LogAnalyser trace, int i, Tree tree) {
		buf.append("<itsabeastystatewerein version='2.0' sample='0'>\n");
		for (StateNode node : state.stateNodeInput.get()) {
			if (node instanceof Tree) {
				((Tree) node).assignFromWithoutID(tree);
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
			int k = indexOf(id, labels);		
			if (k < 0) {
				throw new IllegalArgumentException("cannot find " + id + " in trace file");
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
				return i + 1;
			}
		}
		for (int i = 0; i < labels.size(); i++) {
			String label = labels.get(i);
			if (label.startsWith(id) && isAllDigits(label.substring(id.length())))  {
				return i + 1;
			}
		}
		if (id.lastIndexOf(".") > 0) {
			String idWithoutPartitionInfo = id.substring(0, id.lastIndexOf("."));
			return indexOf(idWithoutPartitionInfo, labels);
		}
		return -1;
	}

	private boolean isAllDigits(String substring) {
		for (int i = 0; i < substring.length(); i++) {
			char c = substring.charAt(i);
			if (!Character.isDigit(c) && c != '.') {
				return false;
			}
		}
		return true;
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
			int k = indexOf(id, labels);
			if (k < 0) {
				throw new IllegalArgumentException("cannot find " + id + " in trace file");
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
			int k = indexOf(id, labels);
			if (k < 0) {
				throw new IllegalArgumentException("cannot find " + id + " in trace file");
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
		Log.warning("writing to file " + multiStateInputFile);
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
	}

	public static void main(String[] args) throws Exception {
		new Application(new Log2MultiState(), "Log to Multi State", args);
	}
}
