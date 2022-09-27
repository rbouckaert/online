package online.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.util.Randomizer;

@Description("A random walk operator that selects a random dimension of the integer parameter from a partition "
		+ "and perturbs the value a random amount within +/- windowSize.")
public class RandomWalkOnParition extends Operator implements PartitionOperator {
	final public Input<TreePartition> partitionInput = new Input<>("partition", "specifies part of the tree to be operated on");
    final public Input<IntegerParameter> parameterInput = new Input<>("parameter", "the parameter that is operated on", Validate.REQUIRED);
    final public Input<Integer> windowSizeInput =
            new Input<Integer>("windowSize", "the size of the window both up and down when using uniform interval OR standard deviation when using Gaussian", Input.Validate.REQUIRED);

    private TreePartition partition;
    private int windowSize;
	
    public RandomWalkOnParition() {}
    public RandomWalkOnParition(TreePartition partition, IntegerParameter parameter, double weight) {
    	this(partition, parameter, 1, weight);
    }
    public RandomWalkOnParition(TreePartition partition, IntegerParameter parameter, int windowSize, double weight) {
    	initByName("partition", partition, "parameter", parameter, "windowSize", windowSize, "weight", weight);
    }
 
    
    @Override
	public void initAndValidate() {
        windowSize = windowSizeInput.get();
        partition = partitionInput.get();
    }

    /**
     * override this for proposals,
     * returns log of hastingRatio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        final IntegerParameter param = parameterInput.get();
        param.startEditing(this);

        final int i = param.getDimension() %2 == 0 ?
        				partition.getRandomNodeNotRoot() :
        				partition.getRandomNode();
        
        final int value = param.getValue(i);
        final int newValue = value + Randomizer.nextInt(2 * windowSize + 1) - windowSize;

        if (newValue < param.getLower() || newValue > param.getUpper()) {
            // invalid move, can be rejected immediately
            return Double.NEGATIVE_INFINITY;
        }
        if (newValue == value) {
            // this saves calculating the posterior
            return Double.NEGATIVE_INFINITY;
        }

        param.setValue(i, newValue);

        return 0.0;
    }

    @Override
    public void optimize(final double logAlpha) {
        // nothing to optimise
    }

}
