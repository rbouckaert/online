package online;



import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.app.util.XMLFile;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.OperatorSchedule;
import beast.core.Runnable;
import beast.core.State;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.Parameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import online.stateoptimiser.StateOptimiser;
import online.stateoptimiser.StateOptimiserByLocalMCMC;
import online.treeexpander.BinarySearchExpander;

// take rates in account in estimated parameters
// take group sizes in account in estimated parameters
// automatically determines chain length, based on Gelman Rubin or other statistic

@Description("Base class for create a new state extending an input state with different set of taxa")
public class BaseStateExpander extends beast.core.Runnable {
	final public Input<XMLFile> xml1Input = new Input<>("xml1", "BEAST XML file with initial state", new XMLFile("[[none]]"));
	
	final public Input<Long> chainLengthInput = new Input<>("chainLength", "Length of the MCMC chain used after placement of taxa", 1000L);
	final public Input<Double> proportionPartitionedInput = new Input<>("proportionPartitioned", "proportion of MCMC chain only using operators local to the part of the tree that changed. "
			+ "If no additional sequences are deteceted this is set to 0.", 0.75);
	
	final public Input<Long> seedInput = new Input<>("seed", "Specify a random number generator seed");
	final public Input<String> definitionsInput = new Input<>("definitions","comma separated list of definitions used in the XML (like the -D option for BEAST)", "");

	public BaseStateExpander() {
	}
	public BaseStateExpander(Long chainLength) {
		chainLengthInput.setValue(chainLength, this);
	}


	@Override
	public void initAndValidate() {
	}


	@Override
	public void run() throws Exception {
	}

	/* initialise part of model2 that is not already initialised by that of model1 (and remove parts of the state/tree if necessary) */
	public List<String> step1UpdateState(Model model1, Model model2) throws IOException, SAXException, ParserConfigurationException, XMLParserException {
		
		List<String> exclusions = new ArrayList<>();
		List<String> additions = new ArrayList<>();
		determineInclusionsExclusions(model1, model2, exclusions, additions);

		// copy pare of the state that model1 and model2 have in common
		copyCommonStateNodes(model1, model2, additions.size() - exclusions.size());

		
		// remove taxa from model1 that are not in model2
		removeExclusions(model1.tree.getRoot(), exclusions);
		
		// initialise tree of model2 with taxa from model1
		// position additional taxa at locations with high support
		BinarySearchExpander expander = new BinarySearchExpander();
		expander.expandTree(model1, model2, additions, hasGroupSizes);
		
		
		Log.debug(model2.tree.getRoot().toNewick());
		return additions;

	} // updateState
		

	public void step2OptimiseState(Model model2, List<String> additions) throws IOException, SAXException, ParserConfigurationException, XMLParserException {
		afterBurner(model2, additions, proportionPartitionedInput.get());		
	}
	
	public void step3RunMCMC(Model model2) throws IOException, SAXException, ParserConfigurationException {
		if (model2.mcmc2 == null) {
			PartitionMCMC mcmc = PartitionMCMC.newMCMC(model2, null, chainLengthInput.get(), definitionsInput.get());
			mcmc.setProportion(0.0);
			model2.mcmc2 = mcmc;
			//model2.state = mcmc.startStateInput.get();
		}
		model2.mcmc2.run();
	}

	/** run short MCMC chain on subset of nodes around newTaxon 
	 * @throws ParserConfigurationException 
	 * @throws SAXException 
	 * @throws IOException 
	 * @throws XMLParserException **/
	
	//MCMC mcmc = null;
	StateOptimiser optimiser;
	protected void afterBurner(Model model, List<String> additions, double proportionPartitioned) throws IOException, SAXException, ParserConfigurationException, XMLParserException {
		if (optimiser == null) {
			optimiser = new StateOptimiserByLocalMCMC(chainLengthInput.get(), definitionsInput.get());
		}
		
		optimiser.optimise(model, additions);
//		if (mcmc == null) {
//			mcmc = newMCMC(model, additions);
//		}
//		
//		((PartitionMCMC)mcmc).initState(model.state.toXML(0));
//		((PartitionMCMC)mcmc).setProportion(proportionPartitioned);
//		
//		mcmc.run();
//		
////		try {
////			PrintStream o = new PrintStream(new File("/tmp/beast1.xml"));
////			XMLProducer p = new XMLProducer();
////			o.println(p.toRawXML(mcmc));
////			o.close();
////		} catch (Exception e) {
////			// ignore
////		}
////		
////		try {
////			PrintStream o = new PrintStream(new File("/tmp/beast2.xml"));
////			XMLProducer p = new XMLProducer();
////			o.println(p.toRawXML(model.mcmc));
////			o.close();
////		} catch (Exception e) {
////			// ignore
////		}
//		
//		
//		
//		State state = mcmc.startStateInput.get();
//		State other = model.state;
//		for (int i = 0; i < state.getNrOfStateNodes(); i++) {
//			StateNode s1 = other.getStateNode(i);
//			StateNode s2 = state.getStateNode(i);
//			s1.assignFrom(s2);
//		}
	}
	
	protected Node removeExclusions(Node node, List<String> taxaToExclude) {
		if (node.isLeaf()) {
			if (Collections.binarySearch(taxaToExclude, node.getID()) < 0) {
				return node;
			} else {
				return null;
			}
		} else {
			Node left_ = node.getLeft(); 
			Node right_ = node.getRight(); 
			left_ = removeExclusions(left_, taxaToExclude);
			right_ = removeExclusions(right_, taxaToExclude);
			if (left_ == null && right_ == null) {
				return null;
			}
			if (left_ == null) {
				return right_;
			}
			if (right_ == null) {
				return left_;
			}
			node.removeAllChildren(false);
			node.addChild(left_);
			node.addChild(right_);
			return node;
		}
	} // removeExclusions

	
	protected void determineInclusionsExclusions(Model model1, Model model2, List<String> exclusions,
				List<String> additions) {
		String [] taxa1 = model1.tree.getTaxaNames();
		Arrays.sort(taxa1);
		String [] taxa2 = model2.tree.getTaxaNames();
		Arrays.sort(taxa2);
		int i = 0, j = 0;
		while (i < taxa1.length && j < taxa2.length) {
			int k = taxa1[i].compareTo(taxa2[j]);
			if (k == 0) {
				i++;
				j++;
			} else if (k < 0) {
				exclusions.add(taxa1[i]);
				i++;
			} else { // k > 0
				additions.add(taxa2[j]);
				j++;
			}
		}
		while (i < taxa1.length) {
			exclusions.add(taxa1[i]);
			i++;
		}
		while (j < taxa2.length) {
			additions.add(taxa2[j]);
			j++;			
		}
	} // determineInclusionsExclusions


	boolean hasGroupSizes = false;
	// copy all state-nodes unless they are a tree, or they are parameters with different dimensions (like rates)
	protected void copyCommonStateNodes(Model model1, Model model2, int deltaTaxaCount) {
		State state1 = model1.state; 
		State state2 = model2.state;
		List<StateNode> s1 = state1.stateNodeInput.get();
		List<StateNode> s2 = state2.stateNodeInput.get();

		for (int i = 0; i < s1.size(); i++) {
			StateNode sn1 = s1.get(i);
			StateNode sn2 = s2.get(i);
			if (!(sn1 instanceof Tree)) {
				if (sn1 instanceof Parameter<?>) {
					if (((Parameter<?>)sn1).getDimension() == ((Parameter<?>)sn2).getDimension()) {
						if (sn1 instanceof IntegerParameter && sn1.getID().startsWith("bGroupSizes")) {
							hasGroupSizes = true;
						}
					}
				}
			}
		}		
		if (hasGroupSizes) {
			for (StateNodeInitialiser init : model2.mcmc.initialisersInput.get()) {
				init.initStateNodes();
			}
			// initialises internal node count in tree
			model2.state.setEverythingDirty(true);
			model2.state.storeCalculationNodes();
			model2.state.checkCalculationNodesDirtiness();
			model2.posterior.calculateLogP();
			model2.state.acceptCalculationNodes();
		}
		
		Map<String, StateNode> snMap2 = new HashMap<>();
		for (int i = 0; i < s2.size(); i++) {
			snMap2.put(s2.get(i).getID(), s2.get(i));
		}		
		
		for (int i = 0; i < s1.size(); i++) {
			StateNode sn1 = s1.get(i);
			if (snMap2.containsKey(sn1.getID())) {
				StateNode sn2 = snMap2.get(sn1.getID());
				// sanity check
				if (!sn1.getID().equals(sn2.getID()) || sn1.getClass() != sn2.getClass()) {
					throw new IllegalArgumentException("Different states found");
				}
				
				// copy under some conditions
				if (!(sn1 instanceof Tree)) {
					if (sn1 instanceof Parameter<?>) {
						if (((Parameter<?>)sn1).getDimension() == ((Parameter<?>)sn2).getDimension()) {
							if (sn1 instanceof IntegerParameter && sn1.getID().startsWith("bGroupSizes")) {
								sn2.assignFrom(sn1);
								IntegerParameter p = (IntegerParameter) sn2;
								int k = 0;
								for (int j = 0; j < deltaTaxaCount; j++) {
									p.setValue(k, p.getValue(k) + 1);
									k += 1;
									if (k == p.getDimension()) {
										k = 0;
									}
								}
							} else {
								sn2.assignFrom(sn1);
							}
						} else {
							model1.parameters.add((Parameter<?>)sn1);
							model2.parameters.add((Parameter<?>)sn2);
							// get meta data values from tree, copied when tree1 was copied to tree2
							Parameter<?> p1 = (Parameter<?>) sn1;
							Parameter<?> p2 = (Parameter<?>) sn2;
							if (p1 instanceof RealParameter) {
								for (int j = 0; j < p1.getDimension(); j++) {
									((RealParameter)p2).setValue((Double)p1.getValue(j));
								}
							} else {
								for (int j = 0; j < p1.getDimension(); j++) {
									((IntegerParameter)p2).setValue((Integer)p1.getValue(j));
								}
							}
						}
					} else {
						sn2.assignFrom(sn1);
					}
				}
			}
		}
	} // copyCommonStateNodes

	
	protected Model getModelFromFile(XMLFile xmlFile) throws SAXException, IOException, ParserConfigurationException, XMLParserException {
		if (xmlFile == null || xmlFile.getName().equals("[[none]]")) {
			throw new IllegalArgumentException("XML file not specified");
		}
        Map<String, String> parserDefinitions = Util.getParserDefinitions(definitionsInput.get());

		XMLParser parser = new XMLParser(parserDefinitions);
		Runnable runnable = parser.parseFile(xmlFile);
		Model model = new Model();
		model.state = (State) runnable.getInputValue("state");
		for (StateNode sn : model.state.stateNodeInput.get()) {
			if (sn instanceof Tree) {
				if (model.tree == null) {
					model.tree = (Tree) sn;
				} else {
					// TODO: expand to multiple trees
					throw new IllegalArgumentException("More than one tree found in the analysis. Can only handle a single tree");
				}
			}
		}
		for (StateNode sn : model.state.stateNodeInput.get()) {
			if (sn instanceof Parameter && ((Parameter<?>)sn).getDimension() > model.tree.getLeafNodeCount() && Util.isClockModelParameter((Parameter<?>)sn)) {
				Node [] nodes = model.tree.getNodesAsArray();
				for (int j = 0; j < ((Parameter<?>)sn).getDimension(); j++) {
					nodes[j].setMetaData(sn.getID(), ((Parameter<?>)sn).getValue(j));
				}				
			}
		}		

		
		model.mcmc = (MCMC) runnable;
		model.posterior = (Distribution) runnable.getInputValue("distribution");
		model.operatorSchedule = (OperatorSchedule) runnable.getInputValue("operatorschedule");
		if (model.operatorSchedule == null && runnable instanceof MCMC) {
			model.operatorSchedule = ((MCMC) runnable).getOperatorSchedule();
		}
		return model;
	} // getModelFromFile
	

}
