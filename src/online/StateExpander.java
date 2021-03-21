package online;

import java.io.File;
import java.io.IOException;

import beast.app.util.Application;
import beast.core.*;
import beast.core.util.Log;
import beast.util.Randomizer;


// TODO: source from tracelog & tree file + write new tree file (& tracelog) 
// TODO: take rates in account

@Description("Create state file extending an input state file with different set of taxa")
public class StateExpander extends BaseStateExpander {
	final public Input<File> stateInput = new Input<>("state", "state file associated with initial XML file (xml1)", new File("[[none]]"));
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		Log.setLevel(Log.Level.debug);
		if (seedInput.get() != null) {
			Randomizer.setSeed(seedInput.get());
		}

		// import models
		Model model1 = getModelFromFile(xml1Input.get());
		Model model2 = getModelFromFile(xml2Input.get());

		// get state from file
		model1.state.setStateFileName(stateInput.get().getAbsolutePath());
		model1.state.restoreFromFile();
        model1.operatorSchedule.setStateFileName(stateInput.get().getAbsolutePath());
        model1.operatorSchedule.restoreFromFile();
        model2.operatorSchedule.setStateFileName(stateInput.get().getAbsolutePath());
        model2.operatorSchedule.restoreFromFile();
		
		updateState(model1, model2);
		
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
