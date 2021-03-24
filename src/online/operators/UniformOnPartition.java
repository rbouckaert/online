/*
* File Uniform.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
/*
 * UniformOperator.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package online.operators;

import java.text.DecimalFormat;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;


@Description("Randomly selects true internal tree node (i.e. not the root) and move node height uniformly in interval " +
        "restricted by the nodes parent and children.")
public class UniformOnPartition extends TreeOperator implements PartitionOperator {
	final public Input<TreePartition> partitionInput = new Input<>("partition", "specifies part of the tree to be operated on");
	TreePartition partition;
	
	// empty constructor to facilitate construction by XML + initAndValidate
    public UniformOnPartition() {
    }

    public UniformOnPartition(Tree tree) {
        try {
            initByName(treeInput.getName(), tree);
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            throw new RuntimeException("Failed to construct Uniform Tree Operator.");
        }
    }

    public UniformOnPartition(TreeInterface tree, TreePartition partition, double weight) {
        initByName("tree", tree, "partition", partition, "weight", weight);
	}

    double scaleFactor = 0.75;
    
	@Override
    public void initAndValidate() {
    	partition = partitionInput.get();
    }

    /**
     * change the parameter and return the hastings ratio.
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
    	int i = partition.getRandomNode();
    	
        Node node = treeInput.get().getNode(i);
        int attempt = 0;
        while (node.isLeaf()) {
        	i = partition.getRandomNode();
        	attempt++;
        	if (attempt == 100) {
//                this.node = null;
        		return Double.NEGATIVE_INFINITY;
        	}
        }

//        this.node = node;
//        this.originalHeight = node.getHeight();
        
        if (node.isRoot()) {
            final Node root = node;
//            final double lower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
//            final double upper = lower * 2;
//            final double newValue = (Randomizer.nextDouble() * (upper - lower)) + lower;
//            node.setHeight(newValue);
//
//            return 0.0;

            double scale = (scaleFactor + (Randomizer.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
            final double newHeight = root.getHeight() * scale;

            if (newHeight < Math.max(root.getLeft().getHeight(), root.getRight().getHeight())) {
//                this.node = null;
                return Double.NEGATIVE_INFINITY;
            }
            root.setHeight(newHeight);
            return -Math.log(scale);
        }
        
    	
		// Abort if no non-root internal nodes
        if (i < 0) {
//            this.node = null;
            return Double.NEGATIVE_INFINITY;
        }
        
        final double upper = node.getParent().getHeight();
        final double lower = Math.max(node.getLeft().getHeight(), node.getRight().getHeight());
        final double newValue = (Randomizer.nextDouble() * (upper - lower)) + lower;
        node.setHeight(newValue);

        return 0.0;
    }

    
//    Node node;
//    double originalHeight;
    
    @Override
    public void reject(int reason) {
    	super.reject(reason);
    	
//    	if (node != null) {
//    		node.setHeight(originalHeight);
//    	}
    }
    
    @Override
    public void optimize(final double logAlpha) {
        //if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(1.0 / scaleFactor - 1.0);
            setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
        //}
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = value;
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double sf = Math.pow(scaleFactor, ratio);

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

}
