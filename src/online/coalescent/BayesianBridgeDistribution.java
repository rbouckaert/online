package online.coalescent;

import org.apache.commons.math.distribution.Distribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.distribution.ParametricDistribution;

@Description("Flexible, heavy-tailed, Bayesian bridge distribution")
public class BayesianBridgeDistribution extends ParametricDistribution {
	final public Input<Double> alphaInput = new Input<>("alpha", "ranges from 0 to 1 and changes the shape, where smaller alpha places more mass near zero", 0.25);
	final public Input<Double> scaleInput = new Input<>("scale", ">0 is termed the global scale", 1.0);
	
	private double alpha, scale;
	
	@Override
	public void initAndValidate() {
		alpha = alphaInput.get();
		scale = scaleInput.get();
		if (alpha < 0 || alpha > 1) {
			Log.warning("WARNING: BayesianBridgeDistribution alpah should be in [0,1)");
		}
	}
	
	@Override
	public double logDensity(double x) {
		final double logP = -Math.pow(Math.abs(x/scale), alpha);
		return logP;
	}

	@Override
	public Distribution getDistribution() {
		throw new RuntimeException("not implemented yet");
	}

}
