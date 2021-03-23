package online;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.core.*;
import beast.core.util.Evaluator;
import beast.core.util.Log;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.io.IOException;
import java.util.Collections;

import beast.core.MCMC;

@Description("Perform MCMC on a partition of the tree -- this assumes that "
		+ "1) tree operators restore state nodes, "
		+ "instead of the state restoring state nodes (as the default MCMC does). "
		+ "2) no screen loggin or file logging is required")
public class PartitionMCMC extends MCMC {
	private Tree tree;

    /**
     * Set this to true to enable detailed MCMC debugging information
     * to be displayed.
     */
    public static boolean printDebugInfo = true;

	PartitionMCMC(Tree tree) {
		this.tree = tree;
	}

	@Override
	public void run() throws IOException, SAXException, ParserConfigurationException {
		// set up state (again). Other beastObjects may have manipulated the
		// StateNodes, e.g. set up bounds or dimensions
		state.initAndValidate();
		// also, initialise state with the file name to store and set-up whether
		// to resume from file
		state.setStateFileName(stateFileName);
		operatorSchedule.setStateFileName(stateFileName);

		burnIn = burnInInput.get();
		chainLength = chainLengthInput.get();
		state.setEverythingDirty(true);
		posterior = posteriorInput.get();

		oldLogLikelihood = state.robustlyCalcPosterior(posterior);

		state.storeCalculationNodes();

		// do the sampling
		logAlpha = 0;
		debugFlag = Boolean.valueOf(System.getProperty("beast.debug"));

		if (Double.isInfinite(oldLogLikelihood) || Double.isNaN(oldLogLikelihood)) {
			reportLogLikelihoods(posterior, "");
			throw new RuntimeException(
					"Could not find a proper state to initialise. Perhaps try another seed.\nSee http://www.beast2.org/2018/07/04/fatal-errors.html for other possible solutions.");
		}

		loggers = loggersInput.get();

		// put the loggers logging to stdout at the bottom of the logger list so
		// that screen output is tidier.
		Collections.sort(loggers, (o1, o2) -> {
			if (o1.isLoggingToStdout()) {
				return o2.isLoggingToStdout() ? 0 : 1;
			} else {
				return o2.isLoggingToStdout() ? -1 : 0;
			}
		});

		doLoop();

		Log.debug("[logP=" + oldLogLikelihood + "] " + tree.getRoot().toNewick());
	} // run;

	protected void doLoop() throws IOException {
        final boolean isStochastic = posterior.isStochastic();
        oldLogLikelihood = state.robustlyCalcPosterior(posterior);

        for (long sampleNr = -burnIn; sampleNr <= chainLength; sampleNr++) {

			final Operator operator = propagateState(sampleNr);
			if (sampleNr >= 0) {
				operator.optimize(logAlpha);
			}
			callUserFunction(sampleNr);

			if (posterior.getCurrentLogP() == Double.POSITIVE_INFINITY) {
				throw new RuntimeException(
						"Encountered a positive infinite posterior. This is a sign there may be numeric instability in the model.");
			}
			
			
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
            }
		}
	}

	@Override
	public void log(long sampleNr) {
		// suppress output
	}

	@Override
	public void close() {
		// do nothing
	}

	@Override
	protected Operator propagateState(final long sampleNr) {
		state.store(sampleNr);
		final Operator operator = operatorSchedule.selectOperator();

		final Distribution evaluatorDistribution = operator.getEvaluatorDistribution();
		Evaluator evaluator = null;

		if (evaluatorDistribution != null) {
			evaluator = new Evaluator() {
				@Override
				public double evaluate() {
					double logP = 0.0;

					state.storeCalculationNodes();
					state.checkCalculationNodesDirtiness();

					try {
						logP = evaluatorDistribution.calculateLogP();
					} catch (Exception e) {
						e.printStackTrace();
						System.exit(1);
					}

					// state.restore();
					state.store(sampleNr);

					return logP;
				}
			};
		}
		String stateXML = state.toXML(0);
		double logHastingsRatio = operator.proposal(evaluator);

	
        if (logHastingsRatio != Double.NEGATIVE_INFINITY) {

            if (operator.requiresStateInitialisation()) {
                state.storeCalculationNodes();
                state.checkCalculationNodesDirtiness();
            }

            newLogLikelihood = posterior.calculateLogP();
            if (newLogLikelihood == Double.POSITIVE_INFINITY) {
            	newLogLikelihood = Double.NEGATIVE_INFINITY;
            	logHastingsRatio = Double.NEGATIVE_INFINITY;
            }

            logAlpha = newLogLikelihood - oldLogLikelihood + logHastingsRatio; //CHECK HASTINGS
            if (printDebugInfo) System.err.print("\n" + logAlpha + " " + newLogLikelihood + " " + oldLogLikelihood);

            if (logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
                // accept
                oldLogLikelihood = newLogLikelihood;
                state.acceptCalculationNodes();

                if (sampleNr >= 0) {
                    operator.accept();
                }
                if (printDebugInfo) System.err.print(" accept");
            } else {
                // reject
                if (sampleNr >= 0) {
                    operator.reject(newLogLikelihood == Double.NEGATIVE_INFINITY ? -1 : 0);
                }
                //state.restore();
                state.restoreCalculationNodes();
                if (printDebugInfo) System.err.print(" reject");
				String stateXML2 = state.toXML(0);
				if (!stateXML.equals(stateXML2)) {
					int h = 4;
					h++;
				}
            }
            state.setEverythingDirty(false);
        } else {
            // operation failed
            if (sampleNr >= 0) {
                operator.reject(-2);
            }
            //state.restore();
            if (!operator.requiresStateInitialisation()) {
                state.setEverythingDirty(false);
                state.restoreCalculationNodes();
            }
            if (printDebugInfo) System.err.print(" direct reject");
			String stateXML2 = state.toXML(0);
			if (!stateXML.equals(stateXML2)) {
				int h = 4;
				h++;
			}
        }
        log(sampleNr);
        return operator;
	}
}
