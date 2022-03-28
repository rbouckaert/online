package online.operators;



import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.OperatorSchedule;
import beast.core.StateNode;
import beast.core.parameter.Parameter;
import beast.core.util.Log;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.operators.Exchange;
import beast.evolution.operators.SubtreeSlide;
import beast.evolution.operators.TreeOperator;
import beast.evolution.operators.Uniform;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;

@Description("Operator schedule that recognises MultiStepOperators and selects them for the desired "
		+ "number of steps consecutively")
public class MultiStepOperatorScheduleForSingleTree extends OperatorSchedule {
	final public Input<Integer> proposalsPerNodeInput = new Input<>("proposalsPerNode", "number of proposals done for a node before moving on to the next node", 3);
	final public Input<List<GenericTreeLikelihood>> potentialTargetsInput = new Input<>("potentialTarget", "likelihoods affected by the node proposal that potenitally can be targetable", new ArrayList<>());
	final public Input<List<Targetable>> targetsInput = new Input<>("target", "likelihoods affected by the node proposal. "
			+ "If not specified, traverse the model graph in search of targetable objects.", new ArrayList<>());
	// public Input<List<TargetableOperator>> operatorInput = new Input<>("operator", "operator the CompoundOperator chooses from with probability proportional to its weight", new ArrayList<>(), Validate.REQUIRED);
	final public Input<Boolean> fullTraverseInput = new Input<>("fullTraverse", "whether to visit every node once (false), or on every node visit (true)" , true);
	final public Input<Boolean> includeLeafsInput = new Input<>("includeLeafs", "whether to visit leaf nodes (true) or nor (false)" , false);
	
	final public Input<String> weightsInput = new Input<>("weights", "space separated list of doubles. If specified will be used as state node dimensions");
	final public Input<Boolean> autoWeightInput = new Input<>("autoweight", "whether to automatically reweight operators, or keep given weights" , true);
	

    // private List<TargetableOperator> operators = new ArrayList<>();
	private int proposalsPerNode;
	private List<Targetable> targets;
	private int [] order;
	private int target;
	private Tree tree;
	private boolean fullTraverse, includeLeafs;

    /** last operator used -- record the last choice for parameter tuning **/
    // public TargetableOperator lastOperator;


	// number of attempts for choosing a random operator, e.g. 
	// NNI may be rejected for leaf nodes, so give up if such
    // operators are chosen after MAX_ATTEMPTS times these are
    // selected
	private final static int MAX_ATTEMPTS = 10;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		fullTraverse = fullTraverseInput.get();
		includeLeafs = includeLeafsInput.get();
		
		proposalsPerNode = proposalsPerNodeInput.get();
		targets = targetsInput.get();
		if (targets.size() == 0) {
			discoverTargets(this);
		}
		operators = new ArrayList<>();// operatorsInput.get();
		
		processPotentialTargetables();
		
		processOperators();
		for (Operator p : operatorsInput.get()) {
			addOperator(p);
		}

		for (Operator p : operatorsInput.get()) {
			addOperator(p);
		}
	}
    
	private void discoverTargets(BEASTInterface o2) {
		for (BEASTInterface o : o2.getOutputs()) {
			if (o instanceof MCMC) {
				discoverTargetsDown(o);
			} else {
				discoverTargets(o);
			}
		}
	}

	private void discoverTargetsDown(BEASTInterface o) {
		for (BEASTInterface o2 : o.listActiveBEASTObjects()) {
			if (o2 instanceof Targetable) {
				if (!targets.contains((Targetable) o2)) {
					targets.add((Targetable) o2);
				}
			} else {
				discoverTargetsDown(o2);
			}
		}
	}


	private void setUpWeights() {
	   	// collect stateNodes
		List<StateNode> stateNodes = new ArrayList<>();		
        for (final Operator op : operators) {
        	List<StateNode> nodes = op.listStateNodes();
        	for (StateNode sn : nodes) {
        		if (!stateNodes.contains(sn)) {
        			stateNodes.add(sn);
        		}
        	}
        }
        
        // initialise dimensions
        double [] stateNodeDimensions = new double[stateNodes.size()];
        int totalDimension = 0;
        if (weightsInput.get() != null) {
        	String [] strs = weightsInput.get().trim().split("\\s+");
        	if (strs.length != stateNodes.size()) {
        		throw new IllegalArgumentException("expected " + stateNodes.size() + " weights, but got " + strs.length);
        	}
            totalDimension = 0;
        	for (int i = 0; i < strs.length; i++) {
        		stateNodeDimensions[i] = Double.parseDouble(strs[i]);
        		Log.warning(stateNodes.get(i).getID() + ": " + stateNodeDimensions[i]);
            	totalDimension += stateNodeDimensions[i];
        	}
        } else {
            totalDimension = 0;
            for (int i = 0; i < stateNodeDimensions.length; i++) {
            	StateNode sn = stateNodes.get(i);
            	if (sn instanceof Parameter) {
            		stateNodeDimensions[i] = ((Parameter) sn).getDimension();
            	} else {
            		stateNodeDimensions[i] = 2 * ((TreeInterface) sn).getNodeCount();
            	}
            	totalDimension += stateNodeDimensions[i];
            }
        }
        
        
        double [] stateNodeAcceptCounts = new double[stateNodes.size()];
        Arrays.fill(stateNodeAcceptCounts, 1.0);
        Log.warning("#parameters = " + totalDimension);
        
        // initialise stateNodeWeights
        double [] stateNodeWeights = new double[stateNodes.size()];
        System.arraycopy(stateNodeDimensions, 0, stateNodeWeights, 0, stateNodeWeights.length);
        
        // initialise stateNodeOperators
        List[] stateNodeOperators = new List[stateNodes.size()];
        for (int i = 0; i < stateNodeOperators.length; i++) {
        	stateNodeOperators[i] = new ArrayList<>();
        }
        for (final Operator op : operators) {
        	List<StateNode> nodes = op.listStateNodes();
        	for (StateNode sn : nodes) {
            	int i = stateNodes.indexOf(sn);
            	stateNodeOperators[i].add(op);
        	}
        }
        
        // determine total sum of dimensions of parameters operators operate on
        int [] operatorDimensions = new int[operators.size()];
        for (int i = 0; i < operators.size(); i++) {
        	for (StateNode sn : operators.get(i).listStateNodes()) {
            	if (sn instanceof Parameter) {
            		operatorDimensions[i] += ((Parameter) sn).getDimension();
            	} else {
            		operatorDimensions[i] += 2 * ((TreeInterface) sn).getNodeCount();
            	}        		 
        	}
        }
        
        // initialise operatorsWeights
        double [][] operatorWeights = new double[stateNodes.size()][];
        for (int i = 0; i < stateNodes.size(); i++) {
        	List<Operator> sOperators = stateNodeOperators[i];
        	operatorWeights[i] = new double[sOperators.size()];
        	double dim = stateNodeDimensions[i];
        	for (int j = 0; j < sOperators.size(); j++) {
        		Operator operator = sOperators.get(j);
        		int k = operators.indexOf(operator);
        		operatorWeights[i][j] = operator.getWeight() * dim / operatorDimensions[k];
        	}
        	operatorWeights[i] = Randomizer.getNormalized(operatorWeights[i]);
        	for (int j = 1; j < operatorWeights[i].length; j++) {
        		operatorWeights[i][j] += operatorWeights[i][j-1]; 
        	}
        
        }

        // calculate normalised weights
        double [] normalizedWeights = new double[operators.size()];
        for (int i = 0; i < operatorWeights.length; i++) {
        	for (int j = 0; j < operatorWeights[i].length; j++) {
        		int k = operators.indexOf(stateNodeOperators[i].get(j));
        		double w = (j == 0 ? operatorWeights[i][0] : operatorWeights[i][j] - operatorWeights[i][j-1]); 
        		normalizedWeights[k] += stateNodeDimensions[i] * w;
        	}
        }
        normalizedWeights = Randomizer.getNormalized(normalizedWeights);
        System.out.println(Arrays.toString(normalizedWeights));
        setNormalizedWeights(normalizedWeights);

        // calc cumulative probabilities
        double [] cumulativeProbs = new double[normalizedWeights.length];
        cumulativeProbs[0] = normalizedWeights[0];
        for (int i = 1; i < operators.size(); i++) {
            cumulativeProbs[i] = normalizedWeights[i] + cumulativeProbs[i - 1];
        }
        setCumulativeProbs(cumulativeProbs);
	}


	/*
	 * Replace potential targetables with actual targetables in model
	 */
	private void processPotentialTargetables() {
		for (GenericTreeLikelihood likelihood : potentialTargetsInput.get()) {
			// TODO Auto-generated method stub			
		}
	}



	@Override
	public void addOperator(Operator p) {
		if (p.getClass() == Uniform.class) {
			Operator bp = new TargetableBactrianNodeOperator();
			p = initialiseOperator(p, bp);
		} else if (p.getClass() == SubtreeSlide.class) {
			Operator bp = new TargetableSubTreeSlide();
			p = initialiseOperator(p, bp);
		} else if (p.getClass() == Exchange.class && ((Exchange)p).isNarrowInput.get()) {
			Operator bp = new TargetableExchange();
			p = initialiseOperator(p, bp);
		}
		super.addOperator(p);
		processOperators();
	}
	
	private Operator initialiseOperator(Operator p, Operator bp) {
		Log.warning("replacing " + p.getID() + " with " + bp.getClass().getSimpleName());

		List<Object> os = new ArrayList<>();
		Set<String> inputNames = new LinkedHashSet<>();
		for (Input<?> input : p.listInputs()) {
			inputNames.add(input.getName());
		}
		
		for (Input<?> input : bp.listInputs()) {
			if (inputNames.contains(input.getName())) {
				Object value = p.getInputValue(input.getName());
				if (value != null && !(value instanceof List && ((List<?>)value).size() == 0)) {
				    os.add(input.getName());
				    os.add(value);
				}	
			}
		}
		bp.initByName(os.toArray());
		bp.setID(p.getID());
		return bp;
	}
	
	@Override
	protected void addOperators(Collection<Operator> ops) {
		super.addOperators(ops);
		processOperators();
	}
	
	
	private void processOperators() {
		for (Operator operator : operators) {
			if (operator instanceof TreeOperator) {
				if (tree != null && tree != ((TreeOperator) operator).treeInput.get()) {
					throw new IllegalArgumentException("This operator schedule can only handle 1 tree but found another " + tree.getID() + " " + ((TreeOperator) operator).treeInput.get().getID());
				}
				tree = ((TreeOperator) operator).treeInput.get();
				if (fullTraverse) {
					if (includeLeafs) {
						order = new int[tree.getNodeCount() + tree.getInternalNodeCount() * 2];
					} else {
						order = new int[tree.getInternalNodeCount() * 3];
					}
				} else {
					if (includeLeafs) {
						order = new int[tree.getNodeCount()];
					} else {
						order = new int[tree.getInternalNodeCount()];
					}
				}
			}
		}

		currentStep = 0;
	}
	
	

	/** establish post-order traversal on internal nodes **/
	private void traverse(Node node, int[] is) {
		if (node.isLeaf()) {
			if (includeLeafs) {
				order[is[0]++] = node.getNr();
			}
		} else {
			if (fullTraverse) {
				order[is[0]++] = node.getNr();
			}
			if (Randomizer.nextBoolean()) {
				traverse(node.getLeft(), is);
				order[is[0]++] = node.getNr();
				traverse(node.getRight(), is);
			} else {
				traverse(node.getRight(), is);
				order[is[0]++] = node.getNr();
				traverse(node.getLeft(), is);
			}
			if (fullTraverse) {
				order[is[0]++] = node.getNr();
			}
		}		
	}

	@Override
	protected void reweightOperators() {
		if (autoWeightInput.get()) {
			setUpWeights();
		} else {
			super.reweightOperators();
		}
	}
	
    private int currentStep;	

    
    @Override
	public Operator selectOperator() {
    	
        while (true) {
    		Operator operator = super.selectOperator();
    		
        	if (operator instanceof TargetableOperator && ((TargetableOperator) operator).isTargetable()) {
        		TargetableOperator targetableOperator = (TargetableOperator) operator;
        		int attemptTargetable = 0;

        		while (true) {
					if (currentStep == stepCount()) {
						target = order[order.length - 1];
						// set DuckTreeLikelihood targets, if target is not a leaf
						if (!tree.getNode(target).isLeaf()) {
							for (Targetable t : targets) {
								t.setTarget(target);
							}
						}
						((TargetableOperator) operator).setTarget(target);
						currentStep = 0;
						return operator;
					}
					if (currentStep == 0) {
						// first step, determine current post-order
						traverse(tree.getRoot(), new int[1]);
					}
	
					// set target nodes in Targetables (only if target node changes)
					if (currentStep % proposalsPerNode == 0) {
						target = order[currentStep / proposalsPerNode];
						if (!tree.getNode(target).isLeaf()) {
							for (Targetable t : targets) {
								t.setTarget(target);
							}
						}
					}
					((TargetableOperator) operator).setTarget(target);
	                Node node = tree.getNode(target);
					if ((node.isLeaf() && targetableOperator.canHandleLeafTargets()) ||
						(node.isRoot() && targetableOperator.canHandleRootTargets()) ||
						(!node.isLeaf() && !node.isRoot() && targetableOperator.canHandleInternlTargets())) {
						currentStep++;
						return operator;
					}
					attemptTargetable++;
					if (attemptTargetable >= MAX_ATTEMPTS) {
						currentStep++;
						if (currentStep >= stepCount()) {
							currentStep = 0;
						}
						attemptTargetable = 0;
					}
					do {
						operator = super.selectOperator();
					} while (!(operator instanceof TargetableOperator));
	        	}
			} else {
				return operator;
			}
		}
		
//		if (multiStepOperator != null) {
//			currentStep++;
//			multiStepOperator.setStepNr(currentStep);
//			if (currentStep < stepCount) {
//				return (Operator) multiStepOperator;
//			} else {
//				multiStepOperator = null;
//			}
//		}
//		
//		Operator operator = super.selectOperator();
//		if (operator instanceof MultiStepOperator) {
//			multiStepOperator = (MultiStepOperator) operator;
//			currentStep = 0;
//			multiStepOperator.setStepNr(currentStep);
//			stepCount = multiStepOperator.stepCount();
//		}
// 		return operator;
    }

	public int stepCount() {
		return proposalsPerNode * order.length;
	}

}
