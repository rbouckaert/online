package online.stateoptimiser;

import java.util.List;

import beast.core.Description;
import online.Model;

@Description("Attempt at improving the state of a model to make it closer to the posterior support")
public interface StateOptimiser {
	
	/**
	 * @param model contains the tree and parameters to be optimised
	 * @param additions represents set of nodes added during tree expansion
	 */
	public void optimise(Model model, List<String> additions);

}
