package online;



import java.io.IOException;
import java.util.Collections;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.core.*;
import beast.core.util.Log;

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

	@Override
	public void run() throws IOException, SAXException, ParserConfigurationException {
        // set up state (again). Other beastObjects may have manipulated the
        // StateNodes, e.g. set up bounds or dimensions
        state.initAndValidate();

        burnIn = burnInInput.get();
        chainLength = chainLengthInput.get();
        state.setEverythingDirty(true);
        posterior = posteriorInput.get();

        burnIn = 0;
        oldLogLikelihood = state.robustlyCalcPosterior(posterior);

        state.storeCalculationNodes();

        
        // do the sampling
        logAlpha = 0;
        debugFlag = Boolean.valueOf(System.getProperty("beast.debug"));

        if (Double.isInfinite(oldLogLikelihood) || Double.isNaN(oldLogLikelihood)) {
            reportLogLikelihoods(posterior, "");
            throw new RuntimeException("Could not find a proper state to initialise. Perhaps try another seed.\nSee http://www.beast2.org/2018/07/04/fatal-errors.html for other possible solutions.");
        }

        loggers = loggersInput.get();
        doLoop();

    } // run;
   

    boolean infinityEncountered = false;
    protected void doLoop() throws IOException {
        int corrections = 0;
        final boolean isStochastic = posterior.isStochastic();
                
        if (burnIn > 0) {
        		Log.warning.println("Please wait while BEAST takes " + burnIn + " pre-burnin samples");
        }
        for (long sampleNr = -burnIn; sampleNr <= chainLength; sampleNr++) {
            final Operator operator = propagateState(sampleNr);

            if (debugFlag && sampleNr % 3 == 0 || sampleNr % 10000 == 0) {
                // check that the posterior is correctly calculated at every third
                // sample, as long as we are in debug mode
            	final double originalLogP = isStochastic ? posterior.getNonStochasticLogP() : oldLogLikelihood;
                final double logLikelihood = isStochastic ? state.robustlyCalcNonStochasticPosterior(posterior) : state.robustlyCalcPosterior(posterior);
                if (Math.abs(logLikelihood - originalLogP) > 1e-6) {
                    reportLogLikelihoods(posterior, "");
                    Log.err.println("At sample " + sampleNr + "\nLikelihood incorrectly calculated: " + originalLogP + " != " + logLikelihood
                    		+ "(" + (originalLogP - logLikelihood) + ")"
                            + " Operator: " + operator.getName());
                }
                if (sampleNr > NR_OF_DEBUG_SAMPLES * 3) {
                    // switch off debug mode once a sufficient large sample is checked
                    debugFlag = false;
	                if (Math.abs(logLikelihood - originalLogP) > 1e-6) {
                        // incorrect calculation outside debug period.
                        // This happens infrequently enough that it should repair itself after a robust posterior calculation
                        corrections++;
                        if (corrections > 100) {
                            // after 100 repairs, there must be something seriously wrong with the implementation
                        	Log.err.println("Too many corrections. There is something seriously wrong that cannot be corrected");
                            state.storeToFile(sampleNr);
                            operatorSchedule.storeToFile();
                            System.exit(1);
                        }
                        oldLogLikelihood = state.robustlyCalcPosterior(posterior);;
                    }
                } else {
	                if (Math.abs(logLikelihood - originalLogP) > 1e-6) {
                        // halt due to incorrect posterior during initial debug period
                        state.storeToFile(sampleNr);
                        operatorSchedule.storeToFile();
                        System.exit(1);
                    }
                }
            } else {
                if (sampleNr >= 0) {
                	operator.optimize(logAlpha);
                }
            }
            callUserFunction(sampleNr);
        }
        if (corrections > 0) {
        	Log.err.println("\n\nNB: " + corrections + " posterior calculation corrections were required. This analysis may not be valid!\n\n");
        }
    }
}

