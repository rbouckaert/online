package online.operators;


import beast.core.Description;
import beast.evolution.operators.BactrianNodeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

@Description("Node operator that proposes node heights in traversal order")
public class TargetableBactrianNodeOperator extends BactrianNodeOperator implements TargetableOperator {
	
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
			throw new RuntimeException("Target is not set, use BactrianNodeOperator insetad");
		}
	}
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
	}
	
    @Override
    public double proposal(int target) {
        Tree tree = treeInput.get(this);

        // select target node

        // Abort if no non-root internal nodes
        if (tree.getInternalNodeCount()==1)
            return Double.NEGATIVE_INFINITY;
        
        Node node = tree.getNode(target);
        if (node.isLeaf() || node.isRoot())
            return Double.NEGATIVE_INFINITY;
        
        double upper = node.getParent().getHeight();
        double lower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
        
        double scale = kernelDistribution.getScaler(0, Double.NaN, scaleFactor);

        // transform value
        double value = node.getHeight();
        double y = (upper - value) / (value - lower);
        y *= scale;
        double newValue = (upper + lower * y) / (y + 1.0);
        
        if (newValue < lower || newValue > upper) {
        	throw new RuntimeException("programmer error: new value proposed outside range");
        }
        
        node.setHeight(newValue);

        double logHR = Math.log(scale) + 2.0 * Math.log((newValue - lower)/(value - lower));
        return logHR;
    }

	@Override
	public boolean canHandleLeafTargets() {		
		return false;
	}

	@Override
	public boolean canHandleRootTargets() {		
		return false;
	}
 
	@Override
	public String getName() {
    	return "TargetableBactrianNodeOperator(" + getID() + ")"; 
	}
}
