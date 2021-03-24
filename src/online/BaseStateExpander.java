package online;



import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.app.util.XMLFile;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.OperatorSchedule;
import beast.core.Runnable;
import beast.core.State;
import beast.core.StateNode;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.Parameter;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import beast.util.XMLProducer;
import online.operators.AfterburnOperatorSchedule;
import online.operators.ExchangeOnPartition;
import online.operators.TreePartition;
import online.operators.UniformOnPartition;

//TODO: take rates in account in estimated parameters
//TODO: take group sizes in account in estimated parameters
//TODO: figure out a way to automatically determine chain length, perhaps based on the Gelman Rubin statistic

@Description("Base class for create a new state extending an input state with different set of taxa")
public class BaseStateExpander extends beast.core.Runnable {
	final public Input<XMLFile> xml1Input = new Input<>("xml1", "BEAST XML file with initial state", new XMLFile("[[none]]"));
	
	final public Input<Long> chainLengthInput = new Input<>("chainLength", "Length of the MCMC chain used after placement of taxa", 1000L);
	final public Input<Double> proportionPartitionedInput = new Input<>("proportionPartitioned", "proportion of MCMC chain only using operators local to the part of the tree that changed. "
			+ "If no additional sequences are deteceted this is set to 0.", 0.75);
	
	final public Input<Long> seedInput = new Input<>("seed", "Specify a random number generator seed");

	Map<String, Integer> map;
	Node internalNode;

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

	
	protected void updateState(Model model1, Model model2) throws IOException, SAXException, ParserConfigurationException, XMLParserException {
		// copy pare of the state that model1 and model2 have in common
		copyCommonStateNodes(model1, model2);
		
		List<String> exclusions = new ArrayList<>();
		List<String> additions = new ArrayList<>();
		determineInclusionsExclusions(model1, model2, exclusions, additions);
		
		// remove taxa from model1 that are not in model2
		removeExclusions(model1.tree.getRoot(), exclusions);
		
		// initialise tree of model2 with taxa from model1
		int leafNodeCount = model2.tree.getLeafNodeCount();
		initialiseTree(model1, model2);
		
		// position additional taxa at locations with high support
		for (String taxon : additions) {
			addTaxon(model2, taxon, leafNodeCount);
			positionAdditions(model2, taxon);
		}
		
		afterBurner(model2, additions);
		
		Log.debug(model2.tree.getRoot().toNewick());

	} // updateState
		

	protected void initialiseTree(Model model1, Model model2) {
		Tree tree2 = model2.tree;
		
		// map is used to check the numbering of taxa is not upset -- leaf nodes 0...n-1, internal nodes n...2n-2
		map = new HashMap<>();
		for (int i = 0; i < tree2.getLeafNodeCount(); i++) {
			map.put(tree2.getNode(i).getID(), i);
		}
		
		for (Node node : tree2.getNodesAsArray()) {
			node.removeAllChildren(false);
			node.setParent(null);
			// node.setID(null);
		}
		
        final Tree otherTree = (Tree) model1.tree;
        Node root = tree2.getRoot();
        root.assignFrom(tree2.getNodesAsArray(), otherTree.getRoot());
        root.setParent(null);
    }

	protected void addTaxon(Model model2, String taxon, int leafNodeCount) {
		Tree tree2 = model2.tree;


		Node child = new Node();
		child.setID(taxon);
		child.setNr(map.get(taxon));
		child.setHeight(0.0);
		if (tree2.hasDateTrait()) {
			TraitSet traitSet = tree2.getDateTrait();
			String pattern = traitSet.getTraitName();
	        if (pattern.equals(TraitSet.DATE_TRAIT) ||
	        		pattern.equals(TraitSet.AGE_TRAIT) ||
	                pattern.equals(TraitSet.DATE_FORWARD_TRAIT) ||
	                pattern.equals(TraitSet.DATE_BACKWARD_TRAIT)) {
	        	child.setMetaData(pattern, traitSet.getValue(taxon));
	        }
		}
		
		Node newRoot = new Node();
		newRoot.addChild(model2.tree.getRoot());
		newRoot.addChild(child);
		newRoot.setHeight(model2.tree.getRoot().getHeight() * (model2.tree.getNodeCount()+1.0)/model2.tree.getNodeCount());
		setRoot(model2.tree, newRoot);

		renumberInternal(tree2.getRoot(), tree2.getNodesAsArray(), map, new int[]{leafNodeCount});
		// tree2.initAndValidate();
	} // initialiseTree

	/** changes root of tree in such a way that the original root
	 * node remains root node, swapping nodes if required
	 */
	protected void setRoot(Tree tree, Node newRoot) {
		if (tree.getRoot() == newRoot) {
			return;
		}
		// swap root nodes
		Node oldRoot = tree.getRoot();
		List<Node> children = new ArrayList<>();
		children.addAll(oldRoot.getChildren());

		List<Node> children2 = new ArrayList<>();
		children2.addAll(newRoot.getChildren());
		Node parent = oldRoot.getParent();		
		
		oldRoot.removeAllChildren(false);
		parent.removeChild(oldRoot);
		newRoot.removeAllChildren(false);

		for (Node child : children) {
			if (child != newRoot) { 
				newRoot.addChild(child);
			} else {
				throw new IllegalArgumentException("Don't know how to handle this");
			}
		}
		if (parent != newRoot) {
			parent.addChild(newRoot);
		}
		for (Node child : children2) {
			if (child != oldRoot) { 
				oldRoot.addChild(child);
			} else {
				oldRoot.addChild(newRoot);
			}
		}
		oldRoot.setParent(null);
		
		double tmp = newRoot.getHeight();
		newRoot.setHeight(oldRoot.getHeight());
		oldRoot.setHeight(tmp);
	}

	protected int renumberInternal(Node node, Node [] nodes, Map<String, Integer> map, int[] nr) {
		for (Node child : node.getChildren()) {
			renumberInternal(child, nodes, map, nr);
		}
		if (!node.isLeaf()) {
			node.setNr(nr[0]);
			nr[0]++;
		} else { // node is leaf
			if (node.getNr() != map.get(node.getID())) {
				Log.debug(node.getID() + " " + node.getNr()+" => " + map.get(node.getID()));
				node.setNr(map.get(node.getID()));
			}
		}
		nodes[node.getNr()] = node;
		return nr[0];
	}

	protected void positionAdditions(Model model2, String taxon) throws IOException, SAXException, ParserConfigurationException, XMLParserException {
		// adding a single taxon
		State state = model2.state;
		Distribution posterior = model2.posterior;
		
		// calc logP when new taxon is outgroup
        state.setEverythingDirty(true);
        state.storeCalculationNodes();
        state.checkCalculationNodesDirtiness();
    	double logP = posterior.calculateLogP();
		state.acceptCalculationNodes();
Log.debug("[" + logP + "] " + model2.tree.getRoot().toNewick());		

		// move node that attaches halfway left and right
		int nodeNr = map.get(taxon);
		Node newTaxon = model2.tree.getNode(nodeNr);
		Node root = model2.tree.getRoot();
		internalNode = root;
		// Node newTaxon = model2.tree.getNode(model2.tree.getLeafNodeCount());
		tryLeftRight(newTaxon,
				root.getLeft() == newTaxon ? root.getRight() : root.getLeft(),
				state, posterior, model2.tree, logP);
	} // addAdditions

	/** run short MCMC chain on subset of nodes aroun newTaxon 
	 * @throws ParserConfigurationException 
	 * @throws SAXException 
	 * @throws IOException 
	 * @throws XMLParserException **/
	
	MCMC mcmc = null;
	protected void afterBurner(Model model, List<String> additions) throws IOException, SAXException, ParserConfigurationException, XMLParserException {
		
		if (mcmc == null) {
			mcmc = newMCMC(model, additions);
		}
		
		((PartitionMCMC)mcmc).initState(model.state.toXML(0));
		((PartitionMCMC)mcmc).setProportion(proportionPartitionedInput.get());
		
		mcmc.run();
		
//		try {
//			PrintStream o = new PrintStream(new File("/tmp/beast1.xml"));
//			XMLProducer p = new XMLProducer();
//			o.println(p.toRawXML(mcmc));
//			o.close();
//		} catch (Exception e) {
//			// ignore
//		}
//		
//		try {
//			PrintStream o = new PrintStream(new File("/tmp/beast2.xml"));
//			XMLProducer p = new XMLProducer();
//			o.println(p.toRawXML(model.mcmc));
//			o.close();
//		} catch (Exception e) {
//			// ignore
//		}
		
		
		
		State state = mcmc.startStateInput.get();
		State other = model.state;
		for (int i = 0; i < state.getNrOfStateNodes(); i++) {
			StateNode s1 = other.getStateNode(i);
			StateNode s2 = state.getStateNode(i);
			s1.assignFrom(s2);
		}
	}
	
	/**
	 * build up MCMC object including appropriate PartitionOperators
	 * @return
	 * @throws XMLParserException 
	 */
	private MCMC newMCMC(Model model, List<String> additions) throws XMLParserException { //List<Operator> operators, Tree tree) {
		List<Operator> operators = new ArrayList<>();

		TreePartition partition = determinePartition(model, additions);		
		if (partition.size() != 0) {
			// add partition operators
			ExchangeOnPartition op1 = new ExchangeOnPartition(model.tree, partition, 1.0);
			op1.setID("ExchangeOnPartition");
			UniformOnPartition op2 = new UniformOnPartition(model.tree, partition, 3.0);
			op2.setID("UniformOnPartition");
			operators.add(op1);
			operators.add(op2);

			// TODO: add RateScale operator if required
		} else {
			// set proportion partitioned operators to zero
			proportionPartitionedInput.setValue(0.0, this);
		}
		operators.addAll(model.mcmc.operatorsInput.get());

		
		Logger screenlog = new Logger();
		screenlog.initByName("log", model.mcmc.posteriorInput.get(), "logEvery", (int)(long) chainLengthInput.get());

		MCMC mcmc = new MCMC(); 
		
		mcmc.initByName(
				"distribution", model.mcmc.posteriorInput.get(),
				"state", model.mcmc.startStateInput.get(),
				"chainLength", chainLengthInput.get(),
				"operator", operators,
				"logger", screenlog,
				"operatorschedule", new AfterburnOperatorSchedule()
		);

		
		XMLProducer producer = new XMLProducer();
		String xml = producer.toRawXML(mcmc);
		xml = xml.replaceAll("'beast.core.MCMC'", "'online.PartitionMCMC'");
		
//	try {
//			PrintStream out = new PrintStream(new File("/tmp/beast.xml"));
//			out.println(xml);
//			out.close();
//	} catch (IOException e) {
//		e.printStackTrace();
//	}
		
		XMLParser parser = new XMLParser();
		mcmc = (MCMC) parser.parseBareFragment(xml, true);
		
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

	protected void tryLeftRight(Node newTaxon, Node child, State state, Distribution posterior, Tree tree, double logP) {
		double originalHeigt = internalNode.getHeight();
		
		Node left = child.getLeft();
		Node right = child.getRight();
		double logPleft = tryBranch(newTaxon, left, state, posterior, tree);
		double logPright = tryBranch(newTaxon, right, state, posterior, tree);
		if (logPleft < logP && logPright < logP) {
			// restore internalNode above child
			child = internalNode.getParent();
			positionOnBranch(newTaxon, child, tree);
			newTaxon.getParent().setHeight(originalHeigt);
			return;
		}
		if (logPleft < logPright) {
			if (!right.isLeaf()) {
				tryLeftRight(newTaxon, right, state, posterior, tree, logPright);
			}
			return;
		}
		positionOnBranch(newTaxon, left, tree);
		if (!left.isLeaf()) {
			tryLeftRight(newTaxon, left, state, posterior, tree, logPleft);
		}
	}

	protected double tryBranch(Node newTaxon, Node node, State state, Distribution posterior, Tree tree) {

        state.storeCalculationNodes();
        positionOnBranch(newTaxon, node, tree);
        
//        state.store(-1);
//        state.setEverythingDirty(true);
        state.checkCalculationNodesDirtiness();
    	double logP = posterior.calculateLogP();
		state.acceptCalculationNodes();
Log.debug("[" + logP + "] " + tree.getRoot().toNewick());		
		return logP;
	}

	private void positionOnBranch(Node newTaxon, Node node, Tree tree) {
		// remove attachments of internalNode
		Node newRoot = null;
		if (internalNode.isRoot()) {
			newRoot = internalNode.getLeft() == newTaxon ? internalNode.getRight() : internalNode.getLeft();
			if (newRoot.getParent() != null) {
				newRoot.getParent().removeChild(newRoot);
			}
			newRoot.setParent(null);
		} else {
			internalNode.getParent().removeChild(internalNode);
			internalNode.getParent().addChild(internalNode.getLeft() == newTaxon ? internalNode.getRight() : internalNode.getLeft());
		}
		internalNode.removeChild(internalNode.getLeft() == newTaxon ? internalNode.getRight() : internalNode.getLeft());
		
		// add internalNode above node, halfway along the branch
		Node parent = node.getParent();
		if (parent == null) {
			newRoot = internalNode;
			// setRoot(tree, internalNode);
		} else {
			parent.removeChild(node);
			parent.addChild(internalNode);
			internalNode.setHeight((parent.getHeight() + node.getHeight())/2);		
		}
		internalNode.addChild(node);
		if (newRoot != null) {
			setRoot(tree, newRoot);
			internalNode = newRoot;
		}
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


	// copy all state-nodes unless they are a tree, or they are parameters with different dimensions (like rates)
	protected void copyCommonStateNodes(Model model1, Model model2) {
		State state1 = model1.state; 
		State state2 = model2.state;
		List<StateNode> s1 = state1.stateNodeInput.get();
		List<StateNode> s2 = state2.stateNodeInput.get();
		for (int i = 0; i < s1.size(); i++) {
			StateNode sn1 = s1.get(i);
			StateNode sn2 = s2.get(i);
			// sanity check
			if (!sn1.getID().equals(sn2.getID()) || sn1.getClass() != sn2.getClass()) {
				throw new IllegalArgumentException("Different states found");
			}

			// copy under some conditions
			if (!(sn1 instanceof Tree)) {
				if (sn1 instanceof Parameter<?>) {
					if (((Parameter<?>)sn1).getDimension() == ((Parameter<?>)sn2).getDimension()) {
						sn2.assignFrom(sn1);
					} else {
						model1.parameters.add((Parameter<?>)sn1);
						model2.parameters.add((Parameter<?>)sn2);
					}
				} else {
					sn2.assignFrom(sn1);
				}
			}
		}
	} // checkStatesAreCompatible

	
	protected Model getModelFromFile(XMLFile xmlFile) throws SAXException, IOException, ParserConfigurationException, XMLParserException {
		if (xmlFile == null || xmlFile.getName().equals("[[none]]")) {
			throw new IllegalArgumentException("XML file not specified");
		}
		XMLParser parser = new XMLParser();
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
		model.mcmc = (MCMC) runnable;
		model.posterior = (Distribution) runnable.getInputValue("distribution");
		model.operatorSchedule = (OperatorSchedule) runnable.getInputValue("operatorschedule");
		if (model.operatorSchedule == null && runnable instanceof MCMC) {
			model.operatorSchedule = ((MCMC) runnable).getOperatorSchedule();
		}
		return model;
	} // getModelFromFile


}
