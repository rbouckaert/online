package beast.evolution.tree.coalescent;

import java.util.ArrayList;
import java.util.List;

import beast.core.*;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.BayesianSkyline;
import beast.evolution.tree.coalescent.IntervalType;

@Description("Bayesian skyline plot that integrates out population sizes under an inverse gamma prior")
public class IBSP extends BayesianSkyline {
    final public Input<RealParameter> populationShapeInput = new Input<>("populationShape", "Shape of the inverse gamma prior distribution on population sizes.", Validate.REQUIRED);
    final public Input<RealParameter> populationMeanInput = new Input<>("populationMean", "Mean of the inverse gamma prior distribution on population sizes.", Validate.REQUIRED);
    public Input<Double> ploidyInput = new Input<>("ploidy", "Ploidy (copy number) for the gene, typically a whole number or half (default is 2) "
    		+ "autosomal nuclear: 2, X: 1.5, Y: 0.5, mitrochondrial: 0.5.", 2.0);

    private RealParameter populationShape;
    private RealParameter populationMean;
    
	 // alpha: alpha parameter of inverse Gamma prior on pop sizes
	 // beta: ditto but beta parameters
	 // ploidy: copy number of gene
    private double alpha, beta, ploidy;
    
    public IBSP() {
		popSizeParamInput.setRule(Validate.FORBIDDEN);
	}

    @Override
    public void initAndValidate() {
    	super.initAndValidate();
    	
    	populationShape = populationShapeInput.get();
    	populationMean = populationMeanInput.get();
    	ploidy = ploidyInput.get();
    	if (ploidy <= 0) {
    		throw new IllegalArgumentException("ploidy should be a positive number, not " + ploidy);
    	}			
    }
	
	@Override
	public double calculateLogP() {
        if (!m_bIsPrepared) {
            prepare();
        }

        
        alpha = populationShape.getValue();
        beta = populationMean.getValue() * (alpha - 1.0);

        logP = 0.0;

        int groupIndex = 0;
        Integer [] groupSizes = this.groupSizes.getValues();
        int subIndex = 0;

        List<Integer> lineageCounts = new ArrayList<>();
        List<Double> intervalSizes = new ArrayList<>();
        
        for (int j = 0; j < intervals.getIntervalCount(); j++) {
            lineageCounts.add(intervals.getLineageCount(j));
            intervalSizes.add(intervals.getInterval(j));
            if (intervals.getIntervalType(j) == IntervalType.COALESCENT) {
                subIndex += 1;
            }
            if (subIndex >= groupSizes[groupIndex]) {
            	logP += analyticalLogP(lineageCounts, groupSizes[groupIndex], intervalSizes);
            	
                groupIndex += 1;
                subIndex = 0;
                lineageCounts.clear();
                intervalSizes.clear();
            }
        }
        return logP;
    }
	
	/**
	 * Analytically integrates out population sizes on an epoch in a tree
	 * @param lineageCount: number of lineages at bottom of the epoch
	 * @param eventCount: number of coalescent events in epoch (this excludes tip being sampled)
	 * @param intervalSizes: array of interval sizes
	 * @return
	 */
    private double analyticalLogP(
    		List<Integer> lineageCounts, 
    		int eventCounts,
    		List<Double> intervalSizes 
    		) {

        double partialGamma = 0.0;
        // contributions of intervals
        for (int i = 0; i < lineageCounts.size(); i++) {
        	partialGamma += intervalSizes.get(i) * lineageCounts.get(i) * (lineageCounts.get(i) - 1.0) / 2.0;
        }

        double logGammaRatio = 0.0;
        for (int i = 0; i < eventCounts; i++) {
            logGammaRatio += Math.log(alpha + i);
        }

        final double logP = 
        		- (alpha + eventCounts) * Math.log(beta + partialGamma / ploidy) 
        		+ alpha * Math.log(beta) 
        		- eventCounts * Math.log(ploidy) 
        		+ logGammaRatio;

        return logP;
    }
    
    

}
