package online.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import beast.app.util.Application;
import beast.app.util.XMLFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Runnable;
import beast.core.State;
import beast.util.XMLParser;

@Description("Convert multi-state file produced by StorableState to tree and trace logs")
public class MultiState2Log extends Runnable {
	final public Input<XMLFile> xmlInput = new Input<>("xml","BEAST XML file with loggers", new XMLFile("[[none]]"));
	final public Input<File> multiStateFileInput = new Input<>("multiStateFile", "state file containing multiple states associated with initial XML file."
			+ "If not specified, use xml+\".state.multi\"", new File("[[none]]"));


	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		if (xmlInput.get() != null && !xmlInput.get().exists()) {
			throw new IllegalArgumentException("Could not find XML file " + xmlInput.get().getName());
		}
		
		String multiStateInputFile = multiStateFileInput.get().getPath();
		if (multiStateInputFile == null || multiStateInputFile.equals("[[none]]")) {
			multiStateInputFile = xmlInput.get().getPath() + ".state.multi"; 
		}
		if (!new File(multiStateInputFile).exists()) {
			throw new IllegalArgumentException("Could not find multi state file " + multiStateInputFile);
		}

		XMLParser parser = new XMLParser();
		Runnable run = parser.parseFile(xmlInput.get());
		if (! (run instanceof MCMC)) {
			throw new IllegalArgumentException("Expected MCMC analysis in xml file");
		}
		MCMC mcmc = (MCMC) run;
		State state = mcmc.startStateInput.get();
		for (Logger logger : mcmc.loggersInput.get()) {
			logger.init();
		}
        
		BufferedReader fin = new BufferedReader(new FileReader(multiStateInputFile));
		String stateXML = null;
		long sampleNr = 0;
		do {
			stateXML = nextState(fin);
			state.fromXML(stateXML);
			state.robustlyCalcPosterior(mcmc.posteriorInput.get());
			for (Logger logger : mcmc.loggersInput.get()) {
				logger.log(sampleNr);
			}
			sampleNr++;
		} while (stateXML != null);

		for (Logger logger : mcmc.loggersInput.get()) {
			logger.close();
		}

        fin.close();
	}

	private String nextState(BufferedReader fin) throws IOException {
		StringBuilder b = new StringBuilder();
		while (fin.ready()) {
			String str = fin.readLine();
			b.append(str);
			if (str.startsWith("</itsabeastystatewerein>")) {
				return b.toString();
			}
		}
		return null;
	}

	public static void main(String[] args) throws Exception {
		new Application(new MultiState2Log(), "Multi State to Log", args);
	}

}
