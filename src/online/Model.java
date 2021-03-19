package online;

import java.util.ArrayList;
import java.util.List;

import beast.core.Distribution;
import beast.core.MCMC;
import beast.core.OperatorSchedule;
import beast.core.State;
import beast.core.parameter.Parameter;
import beast.evolution.tree.Tree;

/** container of bits relevant to updating the state of models **/
class Model {
	State state;
	Tree tree;
	List<Parameter<?>> parameters;
	Distribution posterior;
	OperatorSchedule operatorSchedule;
	MCMC mcmc;
	
	Model() {
		parameters = new ArrayList<>();
	}
}