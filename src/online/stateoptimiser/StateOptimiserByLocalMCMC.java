package online.stateoptimiser;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.State;
import beast.core.StateNode;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.Parameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import beast.util.XMLProducer;
import beastbooster.operators.MultiStepOperatorScheduleForSingleTree;
import online.Model;
import online.PartitionMCMC;
import online.Util;
import online.operators.AfterburnOperatorSchedule;
import online.operators.ExchangeOnPartition;
import online.operators.RandomWalkOnParition;
import online.operators.RateScaleOnPartition;
import online.operators.TreePartition;
import online.operators.UniformOnPartition;

@Description("Optimises state by runnin MCMC on nodes and parameters in the partition only")
public class StateOptimiserByLocalMCMC extends BEASTObject implements StateOptimiser {
	final public Input<Long> chainLengthInput = new Input<>("chainLength", "Length of the MCMC chain used after placement of taxa", 1000L);
	final public Input<String> definitionsInput = new Input<>("definitions","comma separated list of definitions used in the XML (like the -D option for BEAST)", "");

	private MCMC mcmc = null;
	
	public StateOptimiserByLocalMCMC() {}
	
	public StateOptimiserByLocalMCMC(Long chainLength, String definitions) {
		initByName("chainLength", chainLength, "definitions", definitions);
	}
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void optimise(Model model, List<String> additions) {
		if (additions.size() == 0) {
			// nothing to do
			return;
		}
		TreePartition partition = determinePartition(model, additions);
		if (mcmc == null) {
			mcmc = newMCMC(model, partition);
		}
		
		try {
			((PartitionMCMC)mcmc).initState(model.state.toXML(0));
			((PartitionMCMC)mcmc).setProportion(1.0);
			
			mcmc.run();
			
			State state = mcmc.startStateInput.get();
			State other = model.state;
			for (int i = 0; i < state.getNrOfStateNodes(); i++) {
				StateNode s1 = other.getStateNode(i);
				StateNode s2 = state.getStateNode(i);
				s1.assignFrom(s2);
			}
			mcmc.run();
		} catch (IOException | SAXException | ParserConfigurationException e) {
			throw new RuntimeException(e);
		}
	}

	private MCMC newMCMC(Model model, TreePartition partition) {
		List<Operator> operators = new ArrayList<>();
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

		
		Logger screenlog = new Logger();
		screenlog.initByName("log", model.mcmc.posteriorInput.get(), "logEvery", (int)(long) chainLengthInput.get());
		
		MultiStepOperatorScheduleForSingleTree subschedule = new MultiStepOperatorScheduleForSingleTree();

		AfterburnOperatorSchedule operatorSchedule = new AfterburnOperatorSchedule();
		operatorSchedule.initByName("subschedule",subschedule);

		MCMC mcmc = new MCMC(); 		
		mcmc.initByName(
				"distribution", model.mcmc.posteriorInput.get(),
				"state", model.mcmc.startStateInput.get(),
				"chainLength", chainLengthInput.get(),
				"operator", operators,
				"logger", screenlog,
				"operatorschedule", operatorSchedule
		);
		subschedule.initByName("operator", model.mcmc.operatorsInput.get());

		
		XMLProducer producer = new XMLProducer();
		String xml = producer.toRawXML(mcmc);
		xml = xml.replaceAll("'beast.core.MCMC'", "'online.PartitionMCMC'");

        xml = xml.replaceAll("\\bbeast.evolution.likelihood.ThreadedTreeLikelihood\\b", "beastbooster.likelihood.DuckThreadedTreeLikelihood");
        xml = xml.replaceAll("\\bbeast.evolution.likelihood.TreeLikelihood\\b", "beastbooster.likelihood.DuckTreeLikelihood");
		
//		try {
//				PrintStream out = new PrintStream(new File("/tmp/beast.xml"));
//				out.println(xml);
//				out.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
	
		Map<String, String> parserDefinitions = Util.getParserDefinitions(definitionsInput.get());
		XMLParser parser = new XMLParser(parserDefinitions, null, false);
		try {
			mcmc = (MCMC) parser.parseBareFragment(xml, true);
		} catch (XMLParserException e) {
			throw new RuntimeException(e);
		}
		return mcmc;
	}



	protected TreePartition determinePartition(Model model, List<String> additions) {
		Set<Integer> values = new HashSet<>();
		for (String taxonName : additions) {
			int nodeNr = indexOf(taxonName, model.tree.getTaxaNames());
			Node newTaxon = model.tree.getNode(nodeNr);
			Node parent = newTaxon.getParent();
			addToPartition(parent, values);
			addToPartition(parent.getLeft(), values);
			addToPartition(parent.getRight(), values);
			
			if (!parent.isRoot()) {
				Node gp = parent.getParent();
				addToPartition(gp, values);
				addToPartition(gp.getLeft(), values);
				addToPartition(gp.getRight(), values);
			}
		}
		
		IntegerParameter index = new IntegerParameter(values.toArray(new Integer[]{}));
		TreePartition partition = new TreePartition(model.tree, index);
		return partition;
	}


	private int indexOf(String taxonName, String[] taxaNames) {
		for (int i = 0; i < taxaNames.length; i++) {
			if (taxonName.equals(taxaNames[i])) {
				return i;
			}
		}
		throw new IllegalArgumentException("Taxon " + taxonName + " not found in tree");
	}
	
	protected void addToPartition(Node node, Set<Integer> values) {
		if (node.isLeaf()) {
			return;
		}
		values.add(node.getNr());
		if (!node.getLeft().isLeaf()) {
			values.add(node.getLeft().getNr());
		}
		if (node.getRight().isLeaf()) {
			values.add(node.getRight().getNr());
		}
	}
}
