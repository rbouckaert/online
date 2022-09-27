package online.coalescent;

import org.apache.commons.math.distribution.Distribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.distribution.ParametricDistribution;

@Description("scale mixture of normals from https://arxiv.org/pdf/2105.07119.pdf")
public class NormalMixtureDistribution extends ParametricDistribution {
	final public Input<Double> siInput = new Input<>("xi", "slab width that bounds  the  variance  of  increments  to xi^2. "
			+ "A slab width of 2 asserts that there is at most 5% probability for a rate to be greaterthan 50 times its parent rate", 2.0);
	final public Input<Double> lambdaInput = new Input<>("lambda", ">0 ???", 1.0);
	final public Input<Double> muInput = new Input<>("mu", ">0 ???", 1.0);

	
	
	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub

	}

	@Override
	public Distribution getDistribution() {
		// TODO Auto-generated method stub
		return null;
	}

}
