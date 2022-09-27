package beast.evolution.tree.coalescent;

import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.core.Input.Validate;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.distribution.ParametricDistribution;

@Description("Distribution on the difference between logs of rates on branches and parent rates")
public class CorrelatedClockDistribution extends Distribution {
    final public Input<TreeInterface> treeInput = new Input<>("tree", 
    		"phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);
    final public Input<BranchRateModel.Base> clockModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.", Validate.REQUIRED);
    final public Input<ParametricDistribution> distInput = new Input<>("distr", 
    		"Distribution used to calculate prior, e.g. normal, beta, gamma.", Validate.REQUIRED);

    private ParametricDistribution dist;
    private TreeInterface tree;
    private BranchRateModel.Base clockModel;
    private double [] logRates = null;
    
    @Override
    public void initAndValidate() {
        dist = distInput.get();
        tree = treeInput.get();
        clockModel = clockModelInput.get();
		logRates = new double[tree.getNodeCount()];
        calculateLogP();
    }

    
    
    @Override
    public double calculateLogP() {
    	Node [] nodes = tree.getNodesAsArray();
    	for (int i = 0; i < nodes.length-1; i++) {
    		logRates[i] = clockModel.getRateForBranch(nodes[i]);
    	}
    	
    	logP = 0;
    	for (int i = 0; i < nodes.length-1; i++) {
    		final Node parent =nodes[i].getParent();
    		if (!parent.isRoot()) {
    			logP += dist.logDensity(Math.log(logRates[parent.getNr()]) - Math.log(logRates[i]));
    		}
    	}
    	return logP;
    }
    
	@Override
	public List<String> getArguments() {
		return null;	
	}

	@Override
	public List<String> getConditions() {
		return null;
	}

	@Override
	public void sample(State state, Random random) {
	}

}
