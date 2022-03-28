package online.stateoptimiser;

import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.State;
import beast.core.StateNode;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;
import online.Model;
import online.PartitionMCMC;
import online.operators.TreePartition;

@Description("Optimises state by runnin MCMC on nodes and parameters in the partition only")
public class StateOptimiserByLocalMCMC extends BEASTObject implements StateOptimiser {
	final public Input<Long> chainLengthInput = new Input<>("chainLength",
			"Length of the MCMC chain used after placement of taxa", 1000L);
	final public Input<String> definitionsInput = new Input<>("definitions",
			"comma separated list of definitions used in the XML (like the -D option for BEAST)", "");

	private MCMC mcmc = null;

	public StateOptimiserByLocalMCMC() {
	}

	public StateOptimiserByLocalMCMC(Long chainLength, String definitions) {
		initByName("chainLength", chainLength, "definitions", definitions);
	}

	@Override
	public void initAndValidate() {
	}

	@Override
	public void optimise(Model model, List<String> additions) {
		if (additions.size() == 0) {
			// nothing to do
			return;
		}
		TreePartition partition = determinePartition(model, additions);
		if (mcmc == null) {
			mcmc = PartitionMCMC.newMCMC(model, partition, chainLengthInput.get(), definitionsInput.get());
		}

		try {
			((PartitionMCMC) mcmc).initState(model.state.toXML(0));
			((PartitionMCMC) mcmc).setProportion(1.0);

			mcmc.run();

			State state = mcmc.startStateInput.get();
			State other = model.state;
			for (int i = 0; i < state.getNrOfStateNodes(); i++) {
				StateNode s1 = other.getStateNode(i);
				StateNode s2 = state.getStateNode(i);
				s1.assignFrom(s2);
			}
		} catch (IOException | SAXException | ParserConfigurationException e) {
			throw new RuntimeException(e);
		}
	}

	protected TreePartition determinePartition(Model model, List<String> additions) {
		Set<Integer> values = new HashSet<>();
		for (String taxonName : additions) {
			int nodeNr = indexOf(taxonName, model.tree.getTaxaNames());
			Node newTaxon = model.tree.getNode(nodeNr);
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
		}

		IntegerParameter index = new IntegerParameter(values.toArray(new Integer[] {}));
		TreePartition partition = new TreePartition(model.tree, index);
		return partition;
	}

	private int indexOf(String taxonName, String[] taxaNames) {
		for (int i = 0; i < taxaNames.length; i++) {
			if (taxonName.equals(taxaNames[i])) {
				return i;
			}
		}
		throw new IllegalArgumentException("Taxon " + taxonName + " not found in tree");
	}

	protected void addToPartition(Node node, Set<Integer> values) {
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
}
