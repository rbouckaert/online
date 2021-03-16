package online;

import java.io.File;
import java.io.IOException;
import java.util.*;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.app.util.Application;
import beast.app.util.XMLFile;
import beast.core.*;
import beast.core.Runnable;
import beast.core.parameter.Parameter;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.XMLParser;
import beast.util.XMLParserException;

@Description("Create state file extending an input state file with different set of taxa")
public class StateExpander extends Runnable {
	final public Input<XMLFile> xml1Input = new Input<>("xml1", "BEAST XML file with initial state", new XMLFile("[[none]]"));
	final public Input<File> stateInput = new Input<>("state", "state file associated with initial XML file (xml1)", new File("[[none]]"));
	final public Input<XMLFile> xml2Input = new Input<>("xml2", "BEAST XML file with expanded state", new XMLFile("[[none]]"));

	Map<String, Integer> map;

	/** container of bits relevant to updating the state of models **/
	class Model {
		State state;
		Tree tree;
		List<Parameter<?>> parameters;
		Distribution posterior;
		OperatorSchedule operatorSchedule;
		
		Model() {
			parameters = new ArrayList<>();
		}
	}
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		// import models
		Model model1 = getModelFromFile(xml1Input.get());
		Model model2 = getModelFromFile(xml2Input.get());

		// get state from file
		model1.state.setStateFileName(stateInput.get().getAbsolutePath());
		model1.state.restoreFromFile();
        model1.operatorSchedule.setStateFileName(stateInput.get().getAbsolutePath());
        model1.operatorSchedule.restoreFromFile();
		
		updateState(model1, model2);
		
		exportStateFile(model2.state, model1.operatorSchedule);
		Log.warning("Done!");
	}

	private void exportStateFile(State state, OperatorSchedule operatorSchedule) throws IOException {
		String stateFileName = xml2Input.get().getAbsolutePath() + ".state";
		Log.warning("Writing state file to " + stateFileName);
		
		state.setStateFileName(stateFileName);
		state.storeToFile(0);
		
		// also add operator schedule with optimised settings from run of model 1
		operatorSchedule.setStateFileName(stateFileName);
        operatorSchedule.storeToFile();
	} // exportStateFile
	

	private void updateState(Model model1, Model model2) {
		// copy pare of the state that model1 and model2 have in common
		copyCommonStateNodes(model1, model2);
		
		List<String> exclusions = new ArrayList<>();
		List<String> additions = new ArrayList<>();
		determineInclusionsExclusions(model1, model2, exclusions, additions);
		if (additions.size() != 1) {
			throw new IllegalArgumentException("Adding more than 1 taxon is not implemented yet");
			// TODO: implement adding more than 1 taxon
		}
		
		// remove taxa from model1 that are not in model2
		removeExclusions(model1.tree.getRoot(), exclusions);
		
		// initialise tree of model2 with taxa from model1
		initialiseTree(model1, model2, additions);
		
		// position additional taxa at locations with high support
		positionAdditions(model2, additions);
	} // updateState
		

	private void initialiseTree(Model model1, Model model2, List<String> additions) {
		Tree tree2 = model2.tree;
		
		// map is used to check the numbering of taxa is not upset -- leaf nodes 0...n-1, internal nodes n...2n-2
		map = new HashMap<>();
		for (int i = 0; i < tree2.getLeafNodeCount(); i++) {
			map.put(tree2.getNode(i).getID(), i);
		}
		
		for (Node node : tree2.getNodesAsArray()) {
			node.removeAllChildren(false);
			node.setParent(null);
			node.setID(null);
		}
		
		tree2.assignFrom(model1.tree);
		
		for (int i = 0; i < additions.size(); i++) {
			String taxon = additions.get(i);

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
			model2.tree.setRoot(newRoot);
		}
		renumberInternal(tree2.getRoot(), map, new int[]{model1.tree.getLeafNodeCount() + additions.size()});
		tree2.initAndValidate();
	} // initialiseTree

	private int renumberInternal(Node node, Map<String, Integer> map, int[] nr) {
		for (Node child : node.getChildren()) {
			renumberInternal(child, map, nr);
		}
		if (!node.isLeaf()) {
			node.setNr(nr[0]);
			nr[0]++;
		} else { // node is leaf
			if (node.getNr() != map.get(node.getID())) {
				System.out.println(node.getID() + " " + node.getNr()+" => " + map.get(node.getID()));
				node.setNr(map.get(node.getID()));
			}
		}
		return nr[0];
	}

	
	private void positionAdditions(Model model2, List<String> additions) {
		// adding a single taxon
		String taxon = additions.get(0);
		State state = model2.state;
		Distribution posterior = model2.posterior;
		
		// calc logP when new taxon is outgroup
        state.storeCalculationNodes();
        state.checkCalculationNodesDirtiness();
    	double logP = posterior.calculateLogP();
		state.acceptCalculationNodes();
Log.warning("[" + logP + "] " + model2.tree.getRoot().toNewick());		

		// move node that attaches halfway left and right
		Node root = model2.tree.getRoot();
		int nodeNr = map.get(taxon);
		Node newTaxon = model2.tree.getNode(nodeNr);
		// Node newTaxon = model2.tree.getNode(model2.tree.getLeafNodeCount());
		tryLeftRight(root, newTaxon,
				root.getLeft() == newTaxon ? root.getRight() : root.getLeft(),
				state, posterior, model2.tree, logP);
		
		
		
		
		Log.warning(model2.tree.getRoot().toNewick());

	} // addAdditions

	private void tryLeftRight(Node internalNode, Node newTaxon, Node child, State state, Distribution posterior, Tree tree, double logP) {
		double originalHeigt = internalNode.getHeight();
		
		Node left = child.getLeft();
		Node right = child.getRight();
		double logPleft = tryBranch(internalNode, newTaxon, left, state, posterior, tree);
		double logPright = tryBranch(internalNode, newTaxon, right, state, posterior, tree);
		if (logPleft < logP && logPright < logP) {
			// restore internalNode above child
			positionOnBranch(internalNode, newTaxon, child, tree);
			internalNode.setHeight(originalHeigt);
			return;
		}
		if (logPleft < logPright) {
			if (!right.isLeaf()) {
				tryLeftRight(internalNode, newTaxon, right, state, posterior, tree, logPright);
			}
			return;
		}
		positionOnBranch(internalNode, newTaxon, child.getLeft(), tree);
		if (!left.isLeaf()) {
			tryLeftRight(internalNode, newTaxon, left, state, posterior, tree, logPleft);
		}
	}

	private double tryBranch(Node internalNode, Node newTaxon, Node node, State state, Distribution posterior, Tree tree) {

        state.storeCalculationNodes();
        positionOnBranch(internalNode, newTaxon, node, tree);
        
//        state.store(-1);
//        state.setEverythingDirty(true);
        state.checkCalculationNodesDirtiness();
    	double logP = posterior.calculateLogP();
		state.acceptCalculationNodes();
Log.warning("[" + logP + "] " + tree.getRoot().toNewick());		
		return logP;
	}

	private void positionOnBranch(Node internalNode, Node newTaxon, Node node, Tree tree) {
		// remove attachments of internalNode
		if (internalNode.isRoot()) {
			Node newRoot = internalNode.getLeft() == newTaxon ? internalNode.getRight() : internalNode.getLeft();
			tree.setRoot(newRoot);
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
			tree.setRoot(internalNode);
		} else {
			parent.removeChild(node);
			parent.addChild(internalNode);
			internalNode.setHeight((parent.getHeight() + node.getHeight())/2);		
		}
		internalNode.addChild(node);
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

	
	private void determineInclusionsExclusions(Model model1, Model model2, List<String> exclusions,
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


	// copy all state-nodes unless they are a tree, or they are parameters with differenjava -cp AARS.jar:/home/rbou019/.beast/2.6/BEAST/lib/beast.jar beast.app.beastapp.BeastMain -java ClassI_II_protozyme.xmlt dimensions (like rates)
	private void copyCommonStateNodes(Model model1, Model model2) {
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

	
	private Model getModelFromFile(XMLFile xmlFile) throws SAXException, IOException, ParserConfigurationException, XMLParserException {
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
		model.posterior = (Distribution) runnable.getInputValue("distribution");
		
		// TODO: remove next 2 debug lines
		model.posterior = ((CompoundDistribution) model.posterior).pDistributions.get().get(1);
		model.posterior = ((CompoundDistribution) model.posterior).pDistributions.get().get(0);

		model.operatorSchedule = (OperatorSchedule) runnable.getInputValue("operatorschedule");
		if (model.operatorSchedule == null && runnable instanceof MCMC) {
			model.operatorSchedule = ((MCMC) runnable).getOperatorSchedule();
		}
		return model;
	} // getModelFromFile

	
	public static void main(String[] args) throws Exception {
		new Application(new StateExpander(), "State Expander", args);
	}

}
