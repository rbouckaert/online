package online;

import java.io.File;
import java.io.IOException;
import java.util.*;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.app.util.Application;
import beast.app.util.XMLFile;
import beast.core.*;
import beast.core.Runnable;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.Parameter;
import beast.core.util.Evaluator;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.XMLParser;
import beast.util.XMLParserException;

@Description("Create state file extending an input state file with different set of taxa")
public class StateExpander extends Runnable {
	final public Input<XMLFile> xml1Input = new Input<>("xml1", "BEAST XML file with initial state", new XMLFile("[[none]]"));
	final public Input<File> stateInput = new Input<>("state", "state file associated with initial XML file (xml1)", new File("[[none]]"));
	final public Input<XMLFile> xml2Input = new Input<>("xml2", "BEAST XML file with expanded state", new XMLFile("[[none]]"));
	
	final public Input<Long> chainLengthInput = new Input<>("chainLength", "Length of the MCMC chain used after placement of taxa", 1000L);
	final public Input<Long> seedInput = new Input<>("seed", "Specify a random number generator seed");


	Map<String, Integer> map;
	Tree tree;

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
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		Log.setLevel(Log.Level.debug);
		if (seedInput.get() != null) {
			Randomizer.setSeed(seedInput.get());
		}

		// import models
		Model model1 = getModelFromFile(xml1Input.get());
		Model model2 = getModelFromFile(xml2Input.get());

		// get state from file
		model1.state.setStateFileName(stateInput.get().getAbsolutePath());
		model1.state.restoreFromFile();
        model1.operatorSchedule.setStateFileName(stateInput.get().getAbsolutePath());
        model1.operatorSchedule.restoreFromFile();
		
		updateState(model1, model2);
		
		exportStateFile(model2.state, model1.operatorSchedule);
		Log.debug("Done!");
	}

	private void exportStateFile(State state, OperatorSchedule operatorSchedule) throws IOException {
		String stateFileName = xml2Input.get().getAbsolutePath() + ".state";
		Log.debug("Writing state file to " + stateFileName);
		
		state.setStateFileName(stateFileName);
		state.storeToFile(0);
		
		// also add operator schedule with optimised settings from run of model 1
		operatorSchedule.setStateFileName(stateFileName);
        operatorSchedule.storeToFile();
	} // exportStateFile
	

	private void updateState(Model model1, Model model2) throws IOException, SAXException, ParserConfigurationException {
		// copy pare of the state that model1 and model2 have in common
		copyCommonStateNodes(model1, model2);
		
		List<String> exclusions = new ArrayList<>();
		List<String> additions = new ArrayList<>();
		determineInclusionsExclusions(model1, model2, exclusions, additions);
		
		// remove taxa from model1 that are not in model2
		removeExclusions(model1.tree.getRoot(), exclusions);
		
		// initialise tree of model2 with taxa from model1
		int leafNodeCount = model2.tree.getLeafNodeCount();
		initialiseTree(model1, model2);
		
		// position additional taxa at locations with high support
		for (String taxon : additions) {
			addTaxon(model2, taxon, leafNodeCount);
			positionAdditions(model2, taxon);
		}
	} // updateState
		

	private void initialiseTree(Model model1, Model model2) {
		Tree tree2 = model2.tree;
		
		// map is used to check the numbering of taxa is not upset -- leaf nodes 0...n-1, internal nodes n...2n-2
		map = new HashMap<>();
		for (int i = 0; i < tree2.getLeafNodeCount(); i++) {
			map.put(tree2.getNode(i).getID(), i);
		}
		
		for (Node node : tree2.getNodesAsArray()) {
			node.removeAllChildren(false);
			node.setParent(null);
			// node.setID(null);
		}
		
        final Tree otherTree = (Tree) model1.tree;
        Node root = tree2.getRoot();
        root.assignFrom(tree2.getNodesAsArray(), otherTree.getRoot());
        root.setParent(null);
    }

	
	private void addTaxon(Model model2, String taxon, int leafNodeCount) {
		Tree tree2 = model2.tree;


		Node child = new Node();
		child.setID(taxon);
		child.setNr(map.get(taxon));
		child.setHeight(0.0);
		if (tree2.hasDateTrait()) {
			TraitSet traitSet = tree2.getDateTrait();
			String pattern = traitSet.getTraitName();
	        if (pattern.equals(TraitSet.DATE_TRAIT) ||
	        		pattern.equals(TraitSet.AGE_TRAIT) ||
	                pattern.equals(TraitSet.DATE_FORWARD_TRAIT) ||
	                pattern.equals(TraitSet.DATE_BACKWARD_TRAIT)) {
	        	child.setMetaData(pattern, traitSet.getValue(taxon));
	        }
		}
		
		Node newRoot = new Node();
		newRoot.addChild(model2.tree.getRoot());
		newRoot.addChild(child);
		newRoot.setHeight(model2.tree.getRoot().getHeight() * (model2.tree.getNodeCount()+1.0)/model2.tree.getNodeCount());
		setRoot(model2.tree, newRoot);

		renumberInternal(tree2.getRoot(), tree2.getNodesAsArray(), map, new int[]{leafNodeCount});
		// tree2.initAndValidate();
	} // initialiseTree

	/** changes root of tree in such a way that the original root
	 * node remains root node, swapping nodes if required
	 */
	private void setRoot(Tree tree, Node newRoot) {
		if (tree.getRoot() == newRoot) {
			return;
		}
		// swap root nodes
		Node oldRoot = tree.getRoot();
		List<Node> children = new ArrayList<>();
		children.addAll(oldRoot.getChildren());

		List<Node> children2 = new ArrayList<>();
		children2.addAll(newRoot.getChildren());
		Node parent = oldRoot.getParent();		
		
		oldRoot.removeAllChildren(false);
		parent.removeChild(oldRoot);
		newRoot.removeAllChildren(false);

		for (Node child : children) {
			if (child != newRoot) { 
				newRoot.addChild(child);
			} else {
				throw new IllegalArgumentException("Don't know how to handle this");
			}
		}
		if (parent != newRoot) {
			parent.addChild(newRoot);
		}
		for (Node child : children2) {
			if (child != oldRoot) { 
				oldRoot.addChild(child);
			} else {
				oldRoot.addChild(newRoot);
			}
		}
		oldRoot.setParent(null);
		
		double tmp = newRoot.getHeight();
		newRoot.setHeight(oldRoot.getHeight());
		oldRoot.setHeight(tmp);
	}

	private int renumberInternal(Node node, Node [] nodes, Map<String, Integer> map, int[] nr) {
		for (Node child : node.getChildren()) {
			renumberInternal(child, nodes, map, nr);
		}
		if (!node.isLeaf()) {
			node.setNr(nr[0]);
			nr[0]++;
		} else { // node is leaf
			if (node.getNr() != map.get(node.getID())) {
				Log.debug(node.getID() + " " + node.getNr()+" => " + map.get(node.getID()));
				node.setNr(map.get(node.getID()));
			}
		}
		nodes[node.getNr()] = node;
		return nr[0];
	}

	Node internalNode;
	
	private void positionAdditions(Model model2, String taxon) throws IOException, SAXException, ParserConfigurationException {
		// adding a single taxon
		State state = model2.state;
		Distribution posterior = model2.posterior;
		
		// calc logP when new taxon is outgroup
        state.storeCalculationNodes();
        state.checkCalculationNodesDirtiness();
    	double logP = posterior.calculateLogP();
		state.acceptCalculationNodes();
Log.debug("[" + logP + "] " + model2.tree.getRoot().toNewick());		

		// move node that attaches halfway left and right
		int nodeNr = map.get(taxon);
		Node newTaxon = model2.tree.getNode(nodeNr);
		Node root = model2.tree.getRoot();
		internalNode = root;
		// Node newTaxon = model2.tree.getNode(model2.tree.getLeafNodeCount());
		tryLeftRight(newTaxon,
				root.getLeft() == newTaxon ? root.getRight() : root.getLeft(),
				state, posterior, model2.tree, logP);
		
		
		runAfterBurn(model2, newTaxon);
		
		Log.debug(model2.tree.getRoot().toNewick());

	} // addAdditions

	/** run short MCMC chain on subset of nodes aroun newTaxon 
	 * @throws ParserConfigurationException 
	 * @throws SAXException 
	 * @throws IOException **/
	private void runAfterBurn(Model model, Node newTaxon) throws IOException, SAXException, ParserConfigurationException {
		
		tree = model.tree;
		TreePartition partition = determinePartition(model, newTaxon);
		// add partition operators
		ExchangeOnPartition op1 = new ExchangeOnPartition(model.tree, partition, 1.0);
		op1.setID("ExchangeOnPartition");
		UniformOnPartition op2 = new UniformOnPartition(model.tree, partition, 3.0);
		op2.setID("UniformOnPartition");
		List<Operator> operators = new ArrayList<>();
		operators.add(op1);
		operators.add(op2);
		
		MCMC mcmc = newMCMC(model.mcmc, operators);
		mcmc.run();
		
	}

	private TreePartition determinePartition(Model model, Node newTaxon) {
		Set<Integer> values = new HashSet<>();
		Node parent = newTaxon.getParent();
		addToPartition(parent, values);
		addToPartition(parent.getLeft(), values);
		addToPartition(parent.getRight(), values);
		
		if (!parent.isRoot()) {
			Node gp = parent.getParent();
			addToPartition(gp, values);
			addToPartition(gp.getLeft(), values);
			addToPartition(gp.getRight(), values);
		}
		
		IntegerParameter index = new IntegerParameter(values.toArray(new Integer[]{}));
		TreePartition partition = new TreePartition(model.tree, index);
		return partition;
	}

	private void addToPartition(Node node, Set<Integer> values) {
		if (node.isLeaf()) {
			return;
		}
		values.add(node.getNr());
		if (!node.getLeft().isLeaf()) {
			values.add(node.getLeft().getNr());
		}
		if (node.getRight().isLeaf()) {
			values.add(node.getRight().getNr());
		}
	}

	private void tryLeftRight(Node newTaxon, Node child, State state, Distribution posterior, Tree tree, double logP) {
		double originalHeigt = internalNode.getHeight();
		
		Node left = child.getLeft();
		Node right = child.getRight();
		double logPleft = tryBranch(newTaxon, left, state, posterior, tree);
		double logPright = tryBranch(newTaxon, right, state, posterior, tree);
		if (logPleft < logP && logPright < logP) {
			// restore internalNode above child
			child = internalNode.getParent();
			positionOnBranch(newTaxon, child, tree);
			newTaxon.getParent().setHeight(originalHeigt);
			return;
		}
		if (logPleft < logPright) {
			if (!right.isLeaf()) {
				tryLeftRight(newTaxon, right, state, posterior, tree, logPright);
			}
			return;
		}
		positionOnBranch(newTaxon, left, tree);
		if (!left.isLeaf()) {
			tryLeftRight(newTaxon, left, state, posterior, tree, logPleft);
		}
	}

	private double tryBranch(Node newTaxon, Node node, State state, Distribution posterior, Tree tree) {

        state.storeCalculationNodes();
        positionOnBranch(newTaxon, node, tree);
        
//        state.store(-1);
//        state.setEverythingDirty(true);
        state.checkCalculationNodesDirtiness();
    	double logP = posterior.calculateLogP();
		state.acceptCalculationNodes();
Log.debug("[" + logP + "] " + tree.getRoot().toNewick());		
		return logP;
	}

	private void positionOnBranch(Node newTaxon, Node node, Tree tree) {
		// remove attachments of internalNode
		Node newRoot = null;
		if (internalNode.isRoot()) {
			newRoot = internalNode.getLeft() == newTaxon ? internalNode.getRight() : internalNode.getLeft();
			if (newRoot.getParent() != null) {
				newRoot.getParent().removeChild(newRoot);
			}
			newRoot.setParent(null);
		} else {
			internalNode.getParent().removeChild(internalNode);
			internalNode.getParent().addChild(internalNode.getLeft() == newTaxon ? internalNode.getRight() : internalNode.getLeft());
		}
		internalNode.removeChild(internalNode.getLeft() == newTaxon ? internalNode.getRight() : internalNode.getLeft());
		
		// add internalNode above node, halfway along the branch
		Node parent = node.getParent();
		if (parent == null) {
			newRoot = internalNode;
			// setRoot(tree, internalNode);
		} else {
			parent.removeChild(node);
			parent.addChild(internalNode);
			internalNode.setHeight((parent.getHeight() + node.getHeight())/2);		
		}
		internalNode.addChild(node);
		if (newRoot != null) {
			setRoot(tree, newRoot);
			internalNode = newRoot;
		}
	}

	protected Node removeExclusions(Node node, List<String> taxaToExclude) {
		if (node.isLeaf()) {
			if (Collections.binarySearch(taxaToExclude, node.getID()) < 0) {
				return node;
			} else {
				return null;
			}
		} else {
			Node left_ = node.getLeft(); 
			Node right_ = node.getRight(); 
			left_ = removeExclusions(left_, taxaToExclude);
			right_ = removeExclusions(right_, taxaToExclude);
			if (left_ == null && right_ == null) {
				return null;
			}
			if (left_ == null) {
				return right_;
			}
			if (right_ == null) {
				return left_;
			}
			node.removeAllChildren(false);
			node.addChild(left_);
			node.addChild(right_);
			return node;
		}
	} // removeExclusions

	
	private void determineInclusionsExclusions(Model model1, Model model2, List<String> exclusions,
				List<String> additions) {
		String [] taxa1 = model1.tree.getTaxaNames();
		Arrays.sort(taxa1);
		String [] taxa2 = model2.tree.getTaxaNames();
		Arrays.sort(taxa2);
		int i = 0, j = 0;
		while (i < taxa1.length && j < taxa2.length) {
			int k = taxa1[i].compareTo(taxa2[j]);
			if (k == 0) {
				i++;
				j++;
			} else if (k < 0) {
				exclusions.add(taxa1[i]);
				i++;
			} else { // k > 0
				additions.add(taxa2[j]);
				j++;
			}
		}
		while (i < taxa1.length) {
			exclusions.add(taxa1[i]);
			i++;
		}
		while (j < taxa2.length) {
			additions.add(taxa2[j]);
			j++;			
		}
	} // determineInclusionsExclusions


	// copy all state-nodes unless they are a tree, or they are parameters with differenjava -cp AARS.jar:/home/rbou019/.beast/2.6/BEAST/lib/beast.jar beast.app.beastapp.BeastMain -java ClassI_II_protozyme.xmlt dimensions (like rates)
	private void copyCommonStateNodes(Model model1, Model model2) {
		State state1 = model1.state; 
		State state2 = model2.state;
		List<StateNode> s1 = state1.stateNodeInput.get();
		List<StateNode> s2 = state2.stateNodeInput.get();
		for (int i = 0; i < s1.size(); i++) {
			StateNode sn1 = s1.get(i);
			StateNode sn2 = s2.get(i);
			// sanity check
			if (!sn1.getID().equals(sn2.getID()) || sn1.getClass() != sn2.getClass()) {
				throw new IllegalArgumentException("Different states found");
			}

			// copy under some conditions
			if (!(sn1 instanceof Tree)) {
				if (sn1 instanceof Parameter<?>) {
					if (((Parameter<?>)sn1).getDimension() == ((Parameter<?>)sn2).getDimension()) {
						sn2.assignFrom(sn1);
					} else {
						model1.parameters.add((Parameter<?>)sn1);
						model2.parameters.add((Parameter<?>)sn2);
					}
				} else {
					sn2.assignFrom(sn1);
				}
			}
		}
		
	} // checkStatesAreCompatible

	
	private Model getModelFromFile(XMLFile xmlFile) throws SAXException, IOException, ParserConfigurationException, XMLParserException {
		if (xmlFile == null || xmlFile.getName().equals("[[none]]")) {
			throw new IllegalArgumentException("XML file not specified");
		}
		XMLParser parser = new XMLParser();
		Runnable runnable = parser.parseFile(xmlFile);
		Model model = new Model();
		model.state = (State) runnable.getInputValue("state");
		for (StateNode sn : model.state.stateNodeInput.get()) {
			if (sn instanceof Tree) {
				if (model.tree == null) {
					model.tree = (Tree) sn;
				} else {
					// TODO: expand to multiple trees
					throw new IllegalArgumentException("More than one tree found in the analysis. Can only handle a single tree");
				}
			}
		}
		model.mcmc = (MCMC) runnable;
		model.posterior = (Distribution) runnable.getInputValue("distribution");
		model.operatorSchedule = (OperatorSchedule) runnable.getInputValue("operatorschedule");
		if (model.operatorSchedule == null && runnable instanceof MCMC) {
			model.operatorSchedule = ((MCMC) runnable).getOperatorSchedule();
		}
		return model;
	} // getModelFromFile

	
	/**
	 * build up MCMC object with appropriate tree operators
	 * TODO: need rate operators as well?
	 * @return
	 */
	private MCMC newMCMC(MCMC other, List<Operator> operators) {
		
		Logger screenlog = new Logger();
		screenlog.initByName("log", other.posteriorInput.get(), "logEvery", (int)(long) chainLengthInput.get());

		MCMC mcmc = new MCMC() {
			   @Override
			    public void run() throws IOException, SAXException, ParserConfigurationException {				   
			        // set up state (again). Other beastObjects may have manipulated the
			        // StateNodes, e.g. set up bounds or dimensions
			        state.initAndValidate();
			        // also, initialise state with the file name to store and set-up whether to resume from file
			        state.setStateFileName(stateFileName);
			        operatorSchedule.setStateFileName(stateFileName);

			        burnIn = burnInInput.get();
			        chainLength = chainLengthInput.get();
			        state.setEverythingDirty(true);
			        posterior = posteriorInput.get();

			        oldLogLikelihood = state.robustlyCalcPosterior(posterior);

			        state.storeCalculationNodes();

			        
			        // do the sampling
			        logAlpha = 0;
			        debugFlag = Boolean.valueOf(System.getProperty("beast.debug"));

			        if (Double.isInfinite(oldLogLikelihood) || Double.isNaN(oldLogLikelihood)) {
			            reportLogLikelihoods(posterior, "");
			            throw new RuntimeException("Could not find a proper state to initialise. Perhaps try another seed.\nSee http://www.beast2.org/2018/07/04/fatal-errors.html for other possible solutions.");
			        }

			        loggers = loggersInput.get();

			        // put the loggers logging to stdout at the bottom of the logger list so that screen output is tidier.
			        Collections.sort(loggers, (o1, o2) -> {
			            if (o1.isLoggingToStdout()) {
			                return o2.isLoggingToStdout() ? 0 : 1;
			            } else {
			                return o2.isLoggingToStdout() ? -1 : 0;
			            }
			        });

			        doLoop();
			        
			        Log.debug("[logP=" + oldLogLikelihood + "] " + tree.getRoot().toNewick());
			    } // run;
			   
			    protected void doLoop() throws IOException {
			        for (long sampleNr = -burnIn; sampleNr <= chainLength; sampleNr++) {
			        				        	
			            final Operator operator = propagateState(sampleNr);
		                if (sampleNr >= 0) {
		                	operator.optimize(logAlpha);
		                }
			            callUserFunction(sampleNr);
			            
			            if (posterior.getCurrentLogP() == Double.POSITIVE_INFINITY) {
			            	throw new RuntimeException("Encountered a positive infinite posterior. This is a sign there may be numeric instability in the model.");
			            }
			        }
			    }

			   @Override
			   public void log(long sampleNr) {
				   // suppress output
			   }
			   
			   @Override
			   public void close() {
				   // do nothing
			   }
			   
			   @Override
			   protected Operator propagateState(final long sampleNr) {
			        state.store(sampleNr);
			        final Operator operator = operatorSchedule.selectOperator();

			        if (printDebugInfo) System.err.print("\n" + sampleNr + " " + operator.getName()+ ":");

			        final Distribution evaluatorDistribution = operator.getEvaluatorDistribution();
			        Evaluator evaluator = null;

			        if (evaluatorDistribution != null) {
			            evaluator = new Evaluator() {
			                @Override
			                public double evaluate() {
			                    double logP = 0.0;

			                    state.storeCalculationNodes();
			                    state.checkCalculationNodesDirtiness();

			                    try {
			                        logP = evaluatorDistribution.calculateLogP();
			                    } catch (Exception e) {
			                        e.printStackTrace();
			                        System.exit(1);
			                    }

			                    state.restore();
			                    state.store(sampleNr);

			                    return logP;
			                }
			            };
			        }
			        double logHastingsRatio = operator.proposal(evaluator);

			        if (logHastingsRatio != Double.NEGATIVE_INFINITY) {

			            if (operator.requiresStateInitialisation()) {
			                state.storeCalculationNodes();
			                state.checkCalculationNodesDirtiness();
			            }

			            newLogLikelihood = posterior.calculateLogP();
			            if (newLogLikelihood == Double.POSITIVE_INFINITY) {
			            	newLogLikelihood = Double.NEGATIVE_INFINITY;
			            	logHastingsRatio = Double.NEGATIVE_INFINITY;
			            }

			            logAlpha = newLogLikelihood - oldLogLikelihood + logHastingsRatio; //CHECK HASTINGS
			            if (printDebugInfo) System.err.print(logAlpha + " " + newLogLikelihood + " " + oldLogLikelihood);

			            if (logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
			                // accept
			                oldLogLikelihood = newLogLikelihood;
			                state.acceptCalculationNodes();

			                if (sampleNr >= 0) {
			                    operator.accept();
			                }
			                if (printDebugInfo) System.err.print(" accept");
			            } else {
			                // reject
			                if (sampleNr >= 0) {
			                    operator.reject(newLogLikelihood == Double.NEGATIVE_INFINITY ? -1 : 0);
			                }
			                // state.restore();
			                state.restoreCalculationNodes();
			                if (printDebugInfo) System.err.print(" reject");
			            }
			            state.setEverythingDirty(false);
			        } else {
			            // operation failed
			            if (sampleNr >= 0) {
			                operator.reject(-2);
			            }
			            // state.restore();
			            if (!operator.requiresStateInitialisation()) {
			                state.setEverythingDirty(false);
			                state.restoreCalculationNodes();
			            }
			            if (printDebugInfo) System.err.print(" direct reject");
			        }
			        log(sampleNr);
			        return operator;
			    }
		};
		
		mcmc.initByName(
				"distribution", other.posteriorInput.get(), 
				"chainLength", chainLengthInput.get(),
				"operator", operators,
				"logger", screenlog,
				"operatorschedule", new OperatorSchedule()
		);
		return mcmc;
	}

	
	public static void main(String[] args) throws Exception {
		new Application(new StateExpander(), "State Expander", args);
	}

}
