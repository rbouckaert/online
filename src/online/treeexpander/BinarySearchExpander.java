package online.treeexpander;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.State;
import beast.core.StateNode;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.Parameter;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.XMLParserException;
import online.Model;

@Description("Tree Expander that adds new taxa in by binary search of "
		+ "posterior fit starting from the root.")
public class BinarySearchExpander implements TreeExpander {
	private Map<String, Integer> map;
	private Node internalNode;
	
	private boolean useLikelihoodNotPosterior;

	
	@Override
	public void expandTree(Model model1, Model model2, List<String> additions, boolean useLikelihoodNotPosterior) {
		this.useLikelihoodNotPosterior = useLikelihoodNotPosterior;
		int leafNodeCount = model2.tree.getLeafNodeCount();
		initialiseTree(model1, model2);
		
		Log.info.print("Adding " + additions.size() + " taxa:");
		int k = 0;
		for (String taxon : additions) {
			if (addTaxon(model1, model2, taxon, leafNodeCount)) {
				try {
					positionAdditions(model2, taxon);
				} catch (IOException | SAXException | ParserConfigurationException | XMLParserException e) {
					e.printStackTrace();
					throw new RuntimeException(e);
				}
			}
			k++;
			if (k%10 == 0) {
				Log.info.print("|");
			} else {
				Log.info.print(".");
			}
		}
		Log.info.println("Done");
		
	}
	
	public void initialiseTree(Model model1, Model model2) {
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
        if (tree2.hasDateTrait()) {
        	// TODO: robustify, taking all leaf dates in account?
        	TraitSet trait = tree2.getDateTrait();
//        	// calc average height difference from trait
        	Node node = otherTree.getNode(0);
        	double delta = Math.abs(trait.getDate(trait.getValue(node.getID())) - trait.getDate(node.getHeight()));
        	if (delta > 0) {
        		shiftNodes(tree2.getRoot(), delta);
        	}
        }
        root.setParent(null);
    }

	private void shiftNodes(Node node, double delta) {
		node.setHeight(node.getHeight() + delta);
		for (Node child : node.getChildren()) {
			shiftNodes(child, delta);
		}
	}
	
	
	/** returns false if new taxon is older than root, so should not be propagated down the tree **/
	protected boolean addTaxon(Model model1, Model model2, String taxon, int leafNodeCount) {
		boolean newTaxonIsYoungerThanRoot = true;
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
		double h = model2.tree.getRoot().getHeight() * (model2.tree.getNodeCount()+1.0)/model2.tree.getNodeCount();
		if (h < child.getHeight()) {
			// child is a new taxon older than the root
			h = child.getHeight() + 0.1;
			newTaxonIsYoungerThanRoot = false;
		}
		newRoot.setHeight(h);
		setRoot(model2.tree, newRoot);

		for (int i = 0; i < model1.parameters.size(); i++) {
			((StateNode)model1.parameters.get(i)).assignFrom((StateNode)model2.parameters.get(i));
		}
		
		// renumber nodes and update clock model parameters
		renumberInternal(tree2.getRoot(), tree2.getNodesAsArray(), map, new int[]{leafNodeCount}, model1.parameters, model2.parameters);
		
		// set meta data for child and newRoot
		setupMetaData(child.getNr(), model2.parameters);
		setupMetaData(newRoot.getNr(), model2.parameters);
		setupMetaData(newRoot.getLeft().getNr(), model2.parameters);
		setupMetaData(newRoot.getRight().getNr(), model2.parameters);

		return newTaxonIsYoungerThanRoot;
	} // addTaxon

	private void setupMetaData(int i, List<Parameter<?>> metaData) {
		for (int k = 0; k < metaData.size(); k++) {
			Parameter<?> p = metaData.get(k);
			if (i < p.getDimension()) {
				if (p instanceof RealParameter) {
					// TODO: robustify this 
					((RealParameter) p).setValue(i, 1.0);
				} else if (p instanceof BooleanParameter) {
					// TODO: robustify this?
					((BooleanParameter) p).setValue(i, false);
				} else {
					// TODO: robustify this? 
					IntegerParameter ip = (IntegerParameter) p;
					ip.setValue(i, (ip.getUpper() - ip.getLower())/2);
				}
			}
		}
	}
	
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

	protected int renumberInternal(Node node, Node [] nodes, Map<String, Integer> map, int[] nr,
			List<Parameter<?>> originalMetaData, List<Parameter<?>> metaData) {
		for (Node child : node.getChildren()) {
			renumberInternal(child, nodes, map, nr, originalMetaData, metaData);
		}
		if (!node.isLeaf()) {
			int i = node.getNr();
			int j = nr[0];
			updateMetaData(i, j, originalMetaData, metaData);
			
			// update node number
			node.setNr(nr[0]);
			nr[0]++;			
		} else { // node is leaf
			int i = node.getNr();
			if (node.getID() == null || map == null || map.get(node.getID()) == null) {
				Log.warning("WARNING: programmer error -- Taxon found " + node.getID() + " that should have been removed");
				Log.warning("WARNING: Expect a crash.");
			}
			if (i != map.get(node.getID())) {
				Log.debug(node.getID() + " " + i +" => " + map.get(node.getID()));
				node.setNr(map.get(node.getID()));
				int j = map.get(node.getID());
				updateMetaData(i, j, originalMetaData, metaData);
			}
		}
		nodes[node.getNr()] = node;
		return nr[0];
	}

	private void updateMetaData(int i, int j, List<Parameter<?>> originalMetaData, List<Parameter<?>> metaData) {
		if (i != j) {
			for (int k = 0; k < metaData.size(); k++) {
				Parameter<?> p = metaData.get(k);
				if (j < p.getDimension()) {
					Parameter<?> p0 = originalMetaData.get(k);
					if (p instanceof RealParameter) {
						((RealParameter) p).setValue(j, (Double)(p0.getValue(i)));
					} else if (p instanceof BooleanParameter) {
						((BooleanParameter) p).setValue(i, (Boolean)(p0.getValue(i)));
					} else {
						((IntegerParameter) p).setValue(j, (Integer)(p0.getValue(i)));
					}
				}
			}
		}		
	}
	
	private void positionAdditions(Model model2, String taxon) throws IOException, SAXException, ParserConfigurationException, XMLParserException {
		// adding a single taxon
		State state = model2.state;
		Distribution posterior = model2.posterior;
		if (useLikelihoodNotPosterior) {
			posterior = ((CompoundDistribution) posterior).pDistributions.get().get(1);
		}
		
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
			if (!right.isLeaf() && right.getHeight() > newTaxon.getHeight()) {
				tryLeftRight(newTaxon, right, state, posterior, tree, logPright);
			}
			return;
		}
		positionOnBranch(newTaxon, left, tree);
		if (!left.isLeaf() && left.getHeight() > newTaxon.getHeight()) {
			tryLeftRight(newTaxon, left, state, posterior, tree, logPleft);
		}
	}
		
	protected double tryBranch(Node newTaxon, Node node, State state, Distribution posterior, Tree tree) {
	
	    state.storeCalculationNodes();
	    
	    boolean success = positionOnBranch(newTaxon, node, tree);
	    if (!success) {
	    	return Double.NEGATIVE_INFINITY;
	    }
	    
	//    state.store(-1);
	//    state.setEverythingDirty(true);
	    state.checkCalculationNodesDirtiness();
		double logP = posterior.calculateLogP();
		state.acceptCalculationNodes();
	Log.debug("[" + logP + "] " + tree.getRoot().toNewick());		
		return logP;
	}
	
	/**
	 * @param newTaxon
	 * @param node
	 * @param tree
	 * @return if node could successfully be placed on branch
	 */
	private boolean positionOnBranch(Node newTaxon, Node node, Tree tree) {
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
			
			double newHeight = (parent.getHeight() + node.getHeight())/2;
			if (newHeight < newTaxon.getHeight()) {
				newHeight = newTaxon.getHeight();
				if (newHeight > parent.getHeight()) {
					return false;
				}
			}
			internalNode.setHeight(newHeight);		
		}
		internalNode.addChild(node);
		if (newRoot != null) {
			setRoot(tree, newRoot);
			internalNode = newRoot;
		}
		return true;
	}
}
