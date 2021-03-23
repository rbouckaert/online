package online;



import beast.core.*;

import beast.core.MCMC;

@Description("Perform MCMC on a partition of the tree -- this assumes that "
		+ "no screen loggin or file logging is required")
public class PartitionMCMC extends MCMC {

	/** restore state from version serialised in XML **/
	public void initState(String xml) {
		state.fromXML(xml);
	}


	// following methods suppress any logging
	@Override
	public void log(long sampleNr) {
	}
	
	@Override
	public void close() {
	}

}
