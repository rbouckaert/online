package online;

import java.io.File;
import java.io.IOException;
import java.util.List;

import beast.app.util.Application;
import beast.app.util.XMLFile;
import beast.core.*;
import beast.core.util.Log;
import beast.util.Randomizer;

@Description("Create state file extending an input state file with different set of taxa")
public class StateExpander extends BaseStateExpander {
	final public Input<File> stateFileInput = new Input<>("stateFile", "state file associated with initial XML file (xml1). "
			+ "If not specified, use xml1+\".state\"", new File("[[none]]"));
	final public Input<XMLFile> xml2Input = new Input<>("xml2", "BEAST XML file with expanded state", new XMLFile("[[none]]"));
	

	@Override
	public void initAndValidate() {
	}
	
	@Override
	public void run() throws Exception {
		String stateFile = stateFileInput.get().getPath();
		if (stateFile == null || stateFile.equals("[[none]]")) {
			stateFile = xml1Input.get().getAbsolutePath() + ".state";
		}
		if (!new File(stateFile).exists()) {
			throw new IllegalArgumentException("Could not find state file " + stateFile);
		}
		
		
		// Log.setLevel(Log.Level.debug);
		if (seedInput.get() != null) {
			Randomizer.setSeed(seedInput.get());
		}

		// import models
		Model model1 = getModelFromFile(xml1Input.get());
		Model model2 = getModelFromFile(xml2Input.get());

		// get state from file
		model1.state.setStateFileName(stateFile);
		model1.state.restoreFromFile();
        model1.operatorSchedule.setStateFileName(stateFile);
        model1.operatorSchedule.restoreFromFile();
        model2.operatorSchedule.setStateFileName(stateFile);
        model2.operatorSchedule.restoreFromFile();
		
		List<String> additions = step1UpdateState(model1, model2);
		if (chainLengthInput.get() > 0) {
			step2OptimiseState(model2, additions);
		}
		
		exportStateFile(model2.state, model1.operatorSchedule);
		Log.debug("Done!");
	}

	private void exportStateFile(State state, OperatorSchedule operatorSchedule) throws IOException {
		String stateFileName = xml2Input.get().getAbsolutePath() + ".state";
		Log.debug("Writing state file to " + stateFileName);
		
		state.setStateFileName(stateFileName);
		state.storeToFile(0);
		
		// also add operator schedule with optimised settings from run of model 1
		operatorSchedule.setStateFileName(stateFileName);
        operatorSchedule.storeToFile();
	} // exportStateFile
		
	public static void main(String[] args) throws Exception {
		new Application(new StateExpander(), "State Expander", args);
	}

}
