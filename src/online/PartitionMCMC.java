package online;




import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.parser.XMLParser;
import beast.base.parser.XMLParserException;
import beast.base.parser.XMLProducer;
import beastbooster.operators.MultiStepOperatorScheduleForSingleTree;
import online.operators.AfterburnOperatorSchedule;
import online.operators.ExchangeOnPartition;
import online.operators.RandomWalkOnParition;
import online.operators.RateScaleOnPartition;
import online.operators.TreePartition;
import online.operators.UniformOnPartition;

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
		
		if (state != startStateInput.get()) {
			assignState(startStateInput.get(), state);
		}
        // set up state (again). Other beastObjects may have manipulated the
        // StateNodes, e.g. set up bounds or dimensions
        state.initAndValidate();

        burnIn = burnInInput.get();
        chainLength = chainLengthInput.get();
		if (operatorSchedule instanceof AfterburnOperatorSchedule) {
			((AfterburnOperatorSchedule)operatorSchedule).reset((long)(chainLengthProportion * chainLength));
		}
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

		if (state != startStateInput.get()) {
			assignState(state, startStateInput.get());
		}

    } // run;
   

    private void assignState(State stateSource, State stateTarget) {
		List<StateNode> s1 = stateSource.stateNodeInput.get();
		List<StateNode> s2 = stateTarget.stateNodeInput.get();

		for (int i = 0; i < s1.size(); i++) {
			StateNode sn1 = s1.get(i);
			StateNode sn2 = s2.get(i);
			sn2.assignFrom(sn1);
		}		
	}


	boolean infinityEncountered = false;
    protected void doLoop() throws IOException {
        int corrections = 0;
        final boolean isStochastic = posterior.isStochastic();
                
// debugFlag = true;        
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
//	                if (Math.abs(logLikelihood - originalLogP) > 1e-6) {
//                        // incorrect calculation outside debug period.
//                        // This happens infrequently enough that it should repair itself after a robust posterior calculation
//                        corrections++;
//                        if (corrections > 100) {
//                            // after 100 repairs, there must be something seriously wrong with the implementation
//                        	Log.err.println("Too many corrections. There is something seriously wrong that cannot be corrected");
//                            state.storeToFile(sampleNr);
//                            operatorSchedule.storeToFile();
//                            System.exit(1);
//                        }
//                        oldLogLikelihood = state.robustlyCalcPosterior(posterior);
//                    }
//                } else {
//	                if (Math.abs(logLikelihood - originalLogP) > 1e-6) {
//                        // halt due to incorrect posterior during initial debug period
//                        state.storeToFile(sampleNr);
//                        operatorSchedule.storeToFile();
//                        System.exit(1);
//                    }
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


    double chainLengthProportion = 0.5;;
	public void setProportion(Double chainLengthProportion) {
		this.chainLengthProportion = chainLengthProportion;		
	}
    


	static public PartitionMCMC newMCMC(Model model, TreePartition partition, Long chainLength, String definitions) {
		List<Operator> operators = new ArrayList<>();
		if (partition != null) {
			// add partition operators
			ExchangeOnPartition op1 = new ExchangeOnPartition(model.tree, partition, 1.0);
			op1.setID("ExchangeOnPartition");
			UniformOnPartition op2 = new UniformOnPartition(model.tree, partition, 3.0);
			op2.setID("UniformOnPartition");
			operators.add(op1);
			operators.add(op2);
			
			// add RateScale or RandomWalk operator if required for clock parameters
			for (Parameter<?> p : model.parameters) {
				if (Util.isClockModelParameter(p)) {
					if (p instanceof RealParameter) {
						// for parameters representing rate per branch
						RateScaleOnPartition op3 = new RateScaleOnPartition(partition, (RealParameter) p, 1.0);
						op3.setID("RateScaleOnPartition");
						operators.add(op3);
					} else {
						// for parameters representing category (relaxed clock) or indicator (random clock) per branch
						RandomWalkOnParition op3 = new RandomWalkOnParition(partition, (IntegerParameter) p, 1.0);
						op3.setID("RandomWalkOnParition");
						operators.add(op3);
					}
				}
			}
		} else {
			operators.addAll(model.mcmc.operatorsInput.get());
		}

		
		Logger screenlog = new Logger();
		screenlog.initByName("log", model.mcmc.posteriorInput.get(), "logEvery", (int)(long) chainLength);
		
		MultiStepOperatorScheduleForSingleTree subschedule = new MultiStepOperatorScheduleForSingleTree();

		AfterburnOperatorSchedule operatorSchedule = new AfterburnOperatorSchedule();
		operatorSchedule.initByName("subschedule",subschedule);

		MCMC mcmc = new MCMC(); 		
		mcmc.initByName(
				"distribution", model.mcmc.posteriorInput.get(),
				"state", model.mcmc.startStateInput.get(),
				"chainLength", chainLength,
				"operator", operators,
				"logger", screenlog,
				"operatorschedule", operatorSchedule
		);
		subschedule.initByName("operator", model.mcmc.operatorsInput.get());

		
		XMLProducer producer = new XMLProducer();
		String xml = producer.toRawXML(mcmc);
		xml = xml.replaceAll("'beast.base.inference.MCMC'", "'online.PartitionMCMC'");

        xml = xml.replaceAll("\\bbeast.base.evolution.likelihood.ThreadedTreeLikelihood\\b", "beastbooster.likelihood.DuckThreadedTreeLikelihood");
        xml = xml.replaceAll("\\bbeast.base.evolution.likelihood.TreeLikelihood\\b", "beastbooster.likelihood.DuckTreeLikelihood");
		
//		try {
//				PrintStream out = new PrintStream(new File("/tmp/beast.xml"));
//				out.println(xml);
//				out.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
	
		Map<String, String> parserDefinitions = Util.getParserDefinitions(definitions);
		XMLParser parser = new XMLParser(parserDefinitions, null, false);
		try {
			mcmc = (MCMC) parser.parseBareFragment(xml, true);			
		} catch (XMLParserException e) {
			throw new RuntimeException(e);
		}
		return (PartitionMCMC) mcmc;
	}

}

