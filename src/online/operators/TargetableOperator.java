package online.operators;

import beast.core.Description;

@Description("Operator that can be set up to target a particular node in tree/dimension in parameter")
public interface TargetableOperator {

	/** sets target dimension **/
	abstract public void setTarget(int target);

	/** 
	 * Propose new state for target node or dimension
	 * @param target dimension to propose for
	 * @return log Hasting's ratio
	 */
	abstract public double proposal(int target);


	/** tell which nodes in the tree can be handled **/
	
	abstract boolean canHandleLeafTargets();
	
	default boolean canHandleRootTargets() {
		return false;
	};
	
	default public boolean canHandleInternlTargets() {
		return true;
	}
	
	/** override for meta-operators that are not always Targetable **/
	default boolean isTargetable() {
		return true;
	}
}
