package online.treeexpander;

import java.util.List;

import beast.base.core.Description;
import online.Model;

@Description("Class for updating a tree in a state based on another tree")
public interface TreeExpander {
	
	/* The tree in model2 will have taxa added from the list of additions
	 * The tree from model1 has a subset of taxa of model2, and can be used
	 * as guide tree for adding taxa
	 * If the useLikelihoodNotPosterior flag is set, the likelihood needs to
	 * be used as guide instead of the posterior (which may include a difficult
	 * to calculate prior for the increasing number of nodes in the tree
	 */
	public void expandTree(Model model1, Model model2, List<String> additions,
			boolean useLikelihoodNotPosterior);

}
