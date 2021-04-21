package test.online;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;

import beast.app.util.Application;
import beast.app.util.OutFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.util.Log;
import beast.util.Randomizer;

@Description("Generates trait set with taxa sampled back in time")
public class SampleTimeGenerator extends Runnable {
	final public Input<Integer> taxonCountInput = new Input<>("taxa", "number of taxa to generate sample times for", 100);
	final public Input<OutFile> outputInput = new Input<>("out", "output file, or stdout if not specified",
			new OutFile("[[none]]"));


	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		int n = taxonCountInput.get();
		
		Double [] times = new Double[n];
		for (int i = 0; i < n-1; i++) {
			times[i] = Randomizer.nextDouble();
		}
		times[n-1] = 0.0;
		
		Arrays.sort(times, new Comparator<Double>() {
			@Override
			public int compare(Double o1, Double o2) {
				if (o1 > o2) {
					return -1;
				} else if (o1 < o2) {
					return 1;
				}
				return 0;
			}
		});
		
		// open file for writing
		PrintStream out = System.out;
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
			Log.warning("Writing to file " + outputInput.get().getPath());
			out = new PrintStream(outputInput.get());
		}

		out.println("<trait id=\"datetrait\" spec=\"beast.evolution.tree.TraitSet\" traitname=\"date-backward\">");
		for (int i = 0; i < n; i++) {
			out.println("t"+i+" = " + times[i] + (i<n-1?",":""));
		}
		out.println("                <taxa id=\"TaxonSet.dna\" spec=\"TaxonSet\" alignment=\"@dna\"/>"); 
		out.println("</trait>");
		
		if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
        	out.close();
        } else {
        	Thread.sleep(250); // make sure stdout is finished before printing "Done"
        }
        Log.warning.println("Done.");
	}

	public static void main(String[] args) throws Exception {
		new Application(new SampleTimeGenerator(), "Sample Time Generator", args);
	}
}
