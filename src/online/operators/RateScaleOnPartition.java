package online.operators;

import java.text.DecimalFormat;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;

@Description("Scale rate parameter for relaxed clock, but only change values in a given tree partition")
public class RateScaleOnPartition extends Operator implements PartitionOperator {
	final public Input<TreePartition> partitionInput = new Input<>("partition", "specifies part of the tree to be operated on");
    public final Input<RealParameter> parameterInput = new Input<>("parameter", "the parameter that is operated on", Validate.REQUIRED);
    public final Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: larger means more bold proposals", 0.75);
    final public Input<Double> scaleUpperLimit = new Input<>("upper", "Upper Limit of scale factor", 1.0 - 1e-8);
    final public Input<Double> scaleLowerLimit = new Input<>("lower", "Lower limit of scale factor", 1e-8);

    TreePartition partition;
    RealParameter parameter;
    double scaleFactor;
    private double upper, lower;
    TreeInterface tree;
	
    public RateScaleOnPartition() {}
    public RateScaleOnPartition(TreePartition partition, RealParameter parameter, double weight) {
    	this(partition, parameter, 0.75, weight);
    }
    public RateScaleOnPartition(TreePartition partition, RealParameter parameter, double scaleFactor, double weight) {
    	initByName("partition", partition, "parameter", parameter, "scaleFactor", scaleFactor, "weight", weight);
    }
    
	@Override
	public void initAndValidate() {
		partition = partitionInput.get();
		parameter = parameterInput.get();
		scaleFactor = scaleFactorInput.get();
		tree = partition.treeInput.get();
        upper = scaleUpperLimit.get();
        lower = scaleLowerLimit.get();
	}

	@Override
	public double proposal() {
		// find suitable rate to scale
        int i = parameter.getDimension() %2 == 0 ?
				partition.getRandomNodeNotRoot() :
				partition.getRandomNode();
				
		// do the scaling
		double scaleOne = getScaler();
        final double newValue = scaleOne * parameter.getValue(i);
        double hastingsRatio = - Math.log(scaleOne);
        if( outsideBounds(newValue, parameter) ) {
            return Double.NEGATIVE_INFINITY;
        }
        parameter.setValue(i, newValue);
		return hastingsRatio;
	}
	
    protected boolean outsideBounds(final double value, final RealParameter param) {
        final Double l = param.getLower();
        final Double h = param.getUpper();

        return (value < l || value > h);
    }

    protected double getScaler() {
        return (scaleFactor + (Randomizer.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
    }


    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        double delta = calcDelta(logAlpha);
        delta += Math.log(1.0 / scaleFactor - 1.0);
        setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = Math.max(Math.min(value, upper), lower);
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double sf = Math.pow(scaleFactor, ratio);

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }
}
