package online.operators;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;

@Description("Specifies connected subset of nodes in a tree")
public class TreePartition extends BEASTObject {
	public Input<TreeInterface> treeInput = new Input<>("tree", "beast tree for which partition is specified", Validate.REQUIRED);
	public Input<IntegerParameter> partitionInput = new Input<>("partition", "node numbers in the tree specifying the partition", Validate.REQUIRED);
	
	Integer [] partition;
	TreeInterface tree;
	
	public TreePartition() {
	}
	
	public TreePartition(TreeInterface tree, IntegerParameter partition) {
		initByName("tree", tree, "partition", partition);
	}

	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		partition = partitionInput.get().getValues();
	}
	
	public void update() {
		partition = partitionInput.get().getValues();
	}
	
	public int getRandomNode() {
		int result = partition[Randomizer.nextInt(partition.length)];
		return result;
	}
	
	public int getRandomNodeNotRoot() {
		for (int i = 0; i < 1000; i++) {
			int nodeNr = getRandomNode();
			if (!tree.getNode(nodeNr).isRoot()) {
				return nodeNr;
			}
		}
		return -1;
	}

	public int size() {
		return partition.length;
	}
	
}
