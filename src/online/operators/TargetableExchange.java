package online.operators;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

/*
 * KNOWN BUGS: WIDE operator cannot be used on trees with 4 or less tips!
 */

@Description("Nearest neighbour interchange but with the restriction that node height must remain consistent.")
public class TargetableExchange extends TreeOperator implements TargetableOperator {


	private int target = -1;
	
	@Override
	public void setTarget(int target) {
		this.target = target;
	}
	
	@Override
	public double proposal() {
		if (target >= 0) {
			return proposal(target);
		} else {
			throw new RuntimeException("Target is not set, use Exchange instead");
		}
	}
	
	@Override
	public void initAndValidate() {
		
	}

	/**
	 * override this for proposals,
	 *
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal
	 *         should not be accepted *
	 */
	@Override
	public double proposal(int target) {
		final Tree tree = treeInput.get(this);

		double logHastingsRatio = 0;
		logHastingsRatio = narrow(tree, target);

		return logHastingsRatio;
	}

	private int isg(final Node n) {
		return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
	}

	private int sisg(final Node n) {
		return n.isLeaf() ? 0 : isg(n);
	}

	/**
	 * WARNING: Assumes strictly bifurcating beast.tree.
	 */
	public double narrow(final Tree tree, int target) {

		final int internalNodes = tree.getInternalNodeCount();
		if (internalNodes <= 1) {
			return Double.NEGATIVE_INFINITY;
		}

		Node grandParent = tree.getNode(target);
        if (grandParent.isLeaf())
            return Double.NEGATIVE_INFINITY;

		Node parentIndex = grandParent.getLeft();
		Node uncle = grandParent.getRight();
		if (parentIndex.getHeight() < uncle.getHeight()) {
			parentIndex = grandParent.getRight();
			uncle = grandParent.getLeft();
		}

		if (parentIndex.isLeaf()) {
			// tree with dated tips
			return Double.NEGATIVE_INFINITY;
		}

		int validGP = 0;
		{
			for (int i = internalNodes + 1; i < 1 + 2 * internalNodes; ++i) {
				validGP += isg(tree.getNode(i));
			}
		}

		final int c2 = sisg(parentIndex) + sisg(uncle);

		final Node i = (Randomizer.nextBoolean() ? parentIndex.getLeft() : parentIndex.getRight());
		exchangeNodes(i, uncle, parentIndex, grandParent);

		final int validGPafter = validGP - c2 + sisg(parentIndex) + sisg(uncle);

		return Math.log((float) validGP / validGPafter);
	}

	
	/* exchange sub-trees whose root are i and j */
	protected void exchangeNodes(Node i, Node j, Node p, Node jP) {
		// precondition p -> i & jP -> j
		replace(p, i, j);
		replace(jP, j, i);
		// postcondition p -> j & p -> i
	}

	@Override
	public boolean canHandleLeafTargets() {
		return false;
	}


	@Override
	public String getName() {
    	return "TargetableExchange(" + getID() + ")"; 
	}
}
