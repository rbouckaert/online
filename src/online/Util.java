package online;

import java.util.HashMap;
import java.util.Map;

import beast.base.core.BEASTInterface;
import beast.base.inference.parameter.Parameter;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.parser.XMLParserException;

public class Util {

	/**
		 * build up MCMC object including appropriate PartitionOperators
		 * @return
		 * @throws XMLParserException 
		 */
	//	private MCMC newMCMC(Model model, List<String> additions) throws XMLParserException { //List<Operator> operators, Tree tree) {
	//		List<Operator> operators = new ArrayList<>();
	//
	//		TreePartition partition = determinePartition(model, additions);		
	//		if (partition.size() != 0) {
	//			// add partition operators
	//			ExchangeOnPartition op1 = new ExchangeOnPartition(model.tree, partition, 1.0);
	//			op1.setID("ExchangeOnPartition");
	//			UniformOnPartition op2 = new UniformOnPartition(model.tree, partition, 3.0);
	//			op2.setID("UniformOnPartition");
	//			operators.add(op1);
	//			operators.add(op2);
	//			
	//			// add RateScale or RandomWalk operator if required for clock parameters
	//			for (Parameter<?> p : model.parameters) {
	//				if (isClockModelParameter(p))
	//					if (p instanceof RealParameter) {
	//						// for parameters representing rate per branch
	//						RateScaleOnPartition op3 = new RateScaleOnPartition(partition, (RealParameter) p, 1.0);
	//						op3.setID("RateScaleOnPartition");
	//						operators.add(op3);
	//					} else {
	//						// for parameters representing category (relaxed clock) or indicator (random clock) per branch
	//						RandomWalkOnParition op3 = new RandomWalkOnParition(partition, (IntegerParameter) p, 1.0);
	//						op3.setID("RandomWalkOnParition");
	//						operators.add(op3);
	//					}
	//				}
	//		} else {
	//			// set proportion partitioned operators to zero
	//			proportionPartitionedInput.setValue(0.0, this);
	//			operators.addAll(model.mcmc.operatorsInput.get());
	//		}
	//
	//		
	//		Logger screenlog = new Logger();
	//		screenlog.initByName("log", model.mcmc.posteriorInput.get(), "logEvery", (int)(long) chainLengthInput.get());
	//		
	//		MultiStepOperatorScheduleForSingleTree subschedule = new MultiStepOperatorScheduleForSingleTree();
	//
	//		AfterburnOperatorSchedule operatorSchedule = new AfterburnOperatorSchedule();
	//		operatorSchedule.initByName("subschedule",subschedule);
	//
	//		MCMC mcmc = new MCMC(); 		
	//		mcmc.initByName(
	//				"distribution", model.mcmc.posteriorInput.get(),
	//				"state", model.mcmc.startStateInput.get(),
	//				"chainLength", chainLengthInput.get(),
	//				"operator", operators,
	//				"logger", screenlog,
	//				"operatorschedule", operatorSchedule
	//		);
	//		subschedule.initByName("operator", model.mcmc.operatorsInput.get());
	//
	//		
	//		XMLProducer producer = new XMLProducer();
	//		String xml = producer.toRawXML(mcmc);
	//		xml = xml.replaceAll("'beast.base.inference.MCMC'", "'online.PartitionMCMC'");
	//
	//        xml = xml.replaceAll("\\bbeast.base.evolution.likelihood.ThreadedTreeLikelihood\\b", "beastbooster.likelihood.DuckThreadedTreeLikelihood");
	//        xml = xml.replaceAll("\\bbeast.base.evolution.likelihood.TreeLikelihood\\b", "beastbooster.likelihood.DuckTreeLikelihood");
	//		
	////	try {
	////			PrintStream out = new PrintStream(new File("/tmp/beast.xml"));
	////			out.println(xml);
	////			out.close();
	////	} catch (IOException e) {
	////		e.printStackTrace();
	////	}
	//	
	//		Map<String, String> parserDefinitions = getParserDefinitions();
	//		XMLParser parser = new XMLParser(parserDefinitions, null, false);
	//		mcmc = (MCMC) parser.parseBareFragment(xml, true);
	//		
	//		return mcmc;
	//	}
	
		
		static public boolean isClockModelParameter(Parameter<?> p) {
			for (BEASTInterface o : ((BEASTInterface)p).getOutputs()) {
				if (o instanceof BranchRateModel) {
					return true;
				}
			}
			return false;
		}

		static public Map<String, String> getParserDefinitions(String definitions) {
	        Map<String, String> parserDefinitions = new HashMap<>();
	        String [] strs = definitions.split("=",-1);
	        for (int eqIdx = 0; eqIdx<strs.length-1; eqIdx++) {
	            int lastCommaIdx = strs[eqIdx].lastIndexOf(",");

	            if (lastCommaIdx != -1 && eqIdx == 0)
	                throw new IllegalArgumentException("Argument to 'definitions' is not well-formed: expecting comma-separated name=value pairs");

	            String name = strs[eqIdx].substring(lastCommaIdx+1);

	            lastCommaIdx = strs[eqIdx+1].lastIndexOf(",");
	            String value;
	            if (eqIdx+1 == strs.length-1) {
	                value = strs[eqIdx+1];
	            } else {
	                if (lastCommaIdx == -1)
	                    throw new IllegalArgumentException("Argument to 'definitions' is not well-formed: expecting comma-separated name=value pairs");

	                value = strs[eqIdx+1].substring(0, lastCommaIdx);
	            }
	            parserDefinitions.put(name, value);
			}
			return parserDefinitions;
		}
}
