package beast.evolution.tree.coalescent;

import org.apache.commons.math.distribution.Distribution;

import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;
import beast.math.distributions.ParametricDistribution;

@Description("Parametric distribution that is high within epsilon from 0 and low otherwise")
public class SlapAndSpikeDistribution extends ParametricDistribution {
	final public Input<Double> epsilonInput = new Input<>("epsilon", "value deemed close enough to zero to warant high density", 1e-5);
	final public Input<Double> deltaLogPInput = new Input<>("deltaLogP", "difference in logP between near zero and away from zero", Math.log(2));
	
	private double epsilon, deltaLogP;
	
	@Override
	public void initAndValidate() {
		epsilon = epsilonInput.get();
		deltaLogP = deltaLogPInput.get();
		if (deltaLogP < 0) {
			Log.warning("WARNING: SlapAndSpikeDistribution delta log p =" + deltaLogP + " < 0 so small values are penalised!");
		}
	}
	
	@Override
	public double logDensity(double x) {
		if (Math.abs(x) < epsilon) {
			return deltaLogP;
		}
		return 0;
	}

	@Override
	public Distribution getDistribution() {
		throw new RuntimeException("not implemented yet");
	}

}
