package online.math;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import beast.app.util.Application;
import beast.app.util.LogFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.core.Runnable;
import beast.util.LogAnalyser;

@Description("Various ways to compare two distributions")
public class DistributionComparator extends Runnable {

	
	public enum ConvergenceCriterion {GR, SplitR, KS, mean, KDE, none}
	final public static String convergenceCriterionDescription = "Criterion for testig convergence:"
			+ "none for always accepting equality, "
    		+ "GR for Gelman-Rubin statistic, "
    		+ "SplitR for use split-R estimate of Gelman-Rubin statistic, "
    		+ "KS for Kolmogorov Smirnov test at p=5% "
    		+ "KDE for kernel density estimate "
    		+ "mean for checking difference of means with stdev=(2*error1+2*error2) ";
	final public Input<List<LogFile>> traceInput = new Input<>("log", "two or more trace files to compare", new ArrayList<>());
	final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trace logs to used as burn-in (and will be ignored)", 10);
    final public Input<ConvergenceCriterion> criterionInput = new Input<>("criterion", convergenceCriterionDescription, ConvergenceCriterion.SplitR, ConvergenceCriterion.values());

	
    private boolean verbose = false;
	final static String space = "                                                ";

	@Override
	public void initAndValidate() {
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}
	
	@Override
	public void run() throws Exception {
		verbose = true;
		int burnInPercentage = burnInPercentageInput.get();
		List<LogFile> traceFiles = traceInput.get();
		LogAnalyser [] trace = new LogAnalyser[traceFiles.size()];
		for (int i = 0; i < trace.length; i++) {
			trace[i] = new LogAnalyser(traceFiles.get(i).getAbsolutePath(), burnInPercentage, true, false);
		}
		
		ConvergenceCriterion criterion = criterionInput.get();
		
		for (int i = 0; i < trace.length - 1; i++) {
			Log.info(traceFiles.get(i).getName() + "--" + traceFiles.get(i+1).getName());
			calcStats(trace[i], trace[i+1], criterion);
		}
	}
	

	public double calcStats(String fileName1, String fileName2, ConvergenceCriterion criterion) throws IOException {
		LogAnalyser trace1 = new LogAnalyser(fileName1, 0, true, false);
		LogAnalyser trace2 = new LogAnalyser(fileName2, 0, true, false);
		return calcStats(trace1, trace2, criterion);
	}	

	public double calcStats(LogAnalyser log1, LogAnalyser log2, ConvergenceCriterion criterion) throws IOException {
		Double [][] trace1 = log1.getTraces();
		Double [][] trace2 = log2.getTraces();

		if (criterion == ConvergenceCriterion.mean) {
			log1.calcStats();
			log2.calcStats();
		}
		
		double maxStat = Double.MIN_VALUE;
		double minStat = Double.MAX_VALUE;
		for (int i = 0; i < log1.getLabels().size(); i++) {
			double stat = maxStat;
			switch (criterion) {
			case GR:
				stat = calcGRStat(trace1[i+1], trace2[i+1]);
				break;
			case SplitR:
				stat =	calcSplitGRStat(trace1[i+1], trace2[i+1]);
				break;
			case KS:
				stat = calsKSStat(trace1[i+1], trace2[i+1]);
				break;
			case mean:
				stat = calcMeanStat(log1.getMean(i+1), log2.getMean(i+1), 
						log1.getStdError(i+1), log2.getStdError(i+1));
				break;
			case KDE:
				stat = calsKDEStat(trace1[i+1], trace2[i+1]);
				break;
			default:
			}
			maxStat = Math.max(maxStat, stat);
			minStat = Math.min(minStat, stat);
			if (verbose) {
				String label = log1.getLabels().get(i);
				Log.info(label + (label.length() < space.length() ? space.substring(label.length()) : " ") + " " + stat);
				
			}
		}
		switch (criterion) {
			case KS:
				return minStat;
			default:
				return maxStat;
		}
	}

	final static private int RANGE = 100;
	private double calsKDEStat(Double[] trace1, Double[] trace2) {
		int n = trace1.length + trace2.length;
		double [] joint = new double[n];
		for (int i = 0; i < trace1.length; i++) {
			joint[i] = trace1[i];
		}
		for (int i = 0; i < trace2.length; i++) {
			joint[trace1.length + i] = trace2[i];
		}
		
		double mean = 0, sumsq = 0, min = joint[0], max = joint[0];
		for (double d : joint) {
			mean += d;
			sumsq += d*d;
			min = Math.min(min, d);
			max = Math.max(max, d);
		}
		mean /= n;
		double stdev0 = Math.sqrt((sumsq * sumsq - n * mean * mean) / (n-1));
		double stdev = 0.9 * stdev0 / Math.pow(n, 0.2);
		
		// pre-calculate contribution to plots
		double delta = (max - min) / RANGE;
		int tail = (int)(2*stdev / delta + 0.5);
		tail = 20;
		stdev = Math.sqrt(10.0*delta);
		
		double [] kernel = new double[tail*2 + 1];
		kernel[tail] = N(0.0, stdev);
		for (int i = 0; i < tail; i++) {
			kernel[i] = N((tail-i)*delta, stdev);
			kernel[tail*2-i] = kernel[i];
		}
		
		
		// create plots
		double [] plot1 = new double[RANGE];
		double [] plot2 = new double[RANGE];
		createPlot(trace1, plot1, min, kernel, delta, tail);
		createPlot(trace2, plot2, min, kernel, delta, tail);
		
		// calculate difference in plots
		double diff12  = 0, diff21 = 0;
		for (int i = 0; i < RANGE; i++) {
			double diff = plot1[i] - plot2[i];
			if (diff > 0) {
				diff12 += diff;
			} else {
				diff21 -= diff;
			}
		}
		// diff12 should be equal to diff21
		
		return diff12 + diff21;
	}

	private void createPlot(Double[] trace1, double[] plot1, double min, double[] kernel, double delta, int tail) {
		for (double d : trace1) {
			int centre = (int)((d-min) / delta + 0.5);
			int lower = Math.max(0,  centre-tail);
			int upper = Math.min(RANGE-1, centre + tail);
			for (int i = lower; i < upper; i++) {
				plot1[i] += kernel[i-centre+tail];
			}
		}
		
		// normalise
		double sum = 0;
		for (double d : plot1) {
			sum += d;
		}
		for (int i = 0; i < plot1.length; i++) {
			plot1[i] /= sum;
		}
	}

	private double N(double x, double stdev) {
		double f = 1.0/Math.sqrt(2.0*Math.PI*stdev*stdev) * Math.exp(-x*x/(2.0*stdev*stdev));
		return f;
	}

	private double calcMeanStat(double mean1, double mean2, double stdErr1, double stdErr2) {
		double diff = Math.abs(mean1 - mean2);
		double stdev =  stdErr1 * 2 + stdErr2 * 2; 
		return diff/stdev;
	}

	private double calsKSStat(Double[] trace1, Double[] trace2) {
		double [] x0 = toDouble(trace1);
		double [] y0 = toDouble(trace2);
		
		KolmogorovSmirnovTest test = new KolmogorovSmirnovTest();
		double p = test.kolmogorovSmirnovTest(x0, y0);
		return p;
	}

	private double[] toDouble(Double[] x) {
		double [] x0 = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			x0[i] = x[i];
		}
		return x0;
	}

	/** original Gelman Rubin statistic for 2 chains **/	
	private double calcGRStat(Double[] trace1, Double[] trace2) {
		double sampleCount = trace1.length;
		if (trace2.length != sampleCount) {
			throw new IllegalArgumentException("Expected traces of the same length");
		}
		
		// calc means and squared means
		double mean1 = 0, mean2 = 0, sumsq1 = 0, sumsq2 = 0;
		for (Double d : trace1) {
			mean1 += d;
			sumsq1 += d * d;
		}
		mean1 /= sampleCount;
		for (Double d : trace2) {
			mean2 += d;
			sumsq2 += d * d;
		}
		mean2 /= sampleCount;

		// calculate variances for both chains
		double var1 = (sumsq1 - mean1 * mean1 * sampleCount)/(sampleCount - 1);
		double var2 = (sumsq2 - mean2 * mean2 * sampleCount)/(sampleCount - 1);
		
		// average variance for this item
		double fW = (var1 + var2) / 2;
		if (fW == 0) {
			return 1;
		}

		// sum to get totals
		double totalMean = (mean1 + mean2) / 2;
		double totalSq = mean1*mean1 + mean2*mean2;
		
		// variance for joint
		double fB = (totalSq - totalMean * totalMean * 2);
		
		
		double varR = ((sampleCount - 1.0)/sampleCount) + (fB/fW)*(1.0/sampleCount);
		double R = Math.sqrt(varR);
		return R;
	}

	/** Split-R, following 
	 * Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A. and Rubin, D.B.. 
	 * Bayesian data analysis. CRC press. 2013.
	 */
	private double calcSplitGRStat(Double[] trace1, Double[] trace2) {
		if (trace2.length != trace1.length) {
			throw new IllegalArgumentException("Expected traces of the same length");
		}
		int sampleCount = trace1.length/2;
		int sampleCountb = trace1.length - sampleCount;
		
		// calc means and squared means
		double mean1a = 0, mean2a = 0, sumsq1a = 0, sumsq2a = 0;
		double mean1b = 0, mean2b = 0, sumsq1b = 0, sumsq2b = 0;
		for (int i = 0; i < sampleCount; i++) {
			final double d = trace1[i];
			mean1a += d;
			sumsq1a += d * d;
			final double d2 = trace2[i];
			mean2a += d2;
			sumsq2a += d2 * d2;
		}
		mean1a /= sampleCount;
		mean2a /= sampleCount;
		for (int i = sampleCount; i < trace1.length; i++) {
			final double d = trace1[i];
			mean1b += d;
			sumsq1b += d * d;
			final double d2 = trace2[i];
			mean2b += d2;
			sumsq2b += d2 * d2;
		}
		mean1b /= sampleCountb;
		mean2b /= sampleCountb;
		

		// calculate variances for both chains
		double var1a = (sumsq1a - mean1a * mean1a * sampleCount)/(sampleCount - 1);
		double var2a = (sumsq2a - mean2a * mean2a * sampleCount)/(sampleCount - 1);
		double var1b = (sumsq1b - mean1b * mean1b * sampleCountb)/(sampleCountb - 1);
		double var2b = (sumsq2b - mean2b * mean2b * sampleCountb)/(sampleCountb - 1);
		
		// average variance for this item
		double fW = (var1a + var2a + var1b + var2b) / 4;
		if (fW == 0) {
			return 1;
		}

		// sum to get totals
		double totalMean = (mean1a + mean2a + mean1b + mean2b) / 4;
		double totalSq = mean1a*mean1a + mean2a*mean2a + mean1b*mean1b + mean2b*mean2b;
		
		// variance for joint
		double fB = (totalSq - totalMean * totalMean * 4)/3.0;
		
		
		double varR = ((sampleCount - 1.0)/sampleCount) + (fB/fW)*(1.0/sampleCount);
		double R = Math.sqrt(varR);
		return R;
	}



	public static void main(String[] args) throws Exception {
		new Application(new DistributionComparator(), "Distribution Comparator", args);
	}
	
}
