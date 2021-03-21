package online;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import beast.core.Description;
import beast.core.Input;
import beast.core.State;

@Description("State that stores a history of states in a multi-state file")
public class StorableState extends State {
	final public Input<Boolean> storeMultiStateInput = new Input<>("storeMultiState", "if true, stores multi-state file (containing all state files being stored)", true);
	
	
	PrintStream out;
	
	@Override
	public void setStateFileName(final String fileName) {
		super.setStateFileName(fileName);
		
		if (!storeMultiStateInput.get()) {
			return;
		}
		
        if (fileName == null) {
        	throw new IllegalArgumentException("Expected state file to be specified");
        }
        String multiStateFile = fileName + ".multi";
        try {
        	out = new PrintStream(new File(multiStateFile));
        } catch (IOException e) {
        	throw new RuntimeException("Could not open file " + multiStateFile + " for writing: " + e.getMessage());
		}
    }

	@Override
    public void storeToFile(final long sample) {
    	super.storeToFile(sample);

    	if (!storeMultiStateInput.get()) {
			return;
		}
    	
        try {
            out.print(toXML(sample));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


}
