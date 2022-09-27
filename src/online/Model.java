package online;

import java.util.ArrayList;
import java.util.List;

import beast.base.inference.Distribution;
import beast.base.inference.MCMC;
import beast.base.inference.OperatorSchedule;
import beast.base.inference.State;
import beast.base.inference.parameter.Parameter;
import beast.base.evolution.tree.Tree;

/** container of bits relevant to updating the state of models **/
public class Model {
	public State state;
	public Tree tree;
	public List<Parameter<?>> parameters;
	public Distribution posterior;
	public OperatorSchedule operatorSchedule;
	public MCMC mcmc;
	public PartitionMCMC mcmc2;
	
	Model() {
		parameters = new ArrayList<>();
	}
}