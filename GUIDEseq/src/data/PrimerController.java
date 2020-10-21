package data;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

public class PrimerController {
	ArrayList<SamplePrimer> a = new ArrayList<SamplePrimer>();

	public PrimerController(File file) {
		try {
			Scanner s = new Scanner(file);
			while(s.hasNextLine()) {
				String line = s.nextLine();
				SamplePrimer sp = SamplePrimer.parse(line);
				a.add(sp);
			}
			s.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	public SamplePrimer getSamplePrimer(File f) {
		String name = f.getName();
		for(SamplePrimer sp: a) {
			if(sp.sampleNameMatches(name)) {
				return sp;
			}
		}
		return null;
	}
	
	public String getPrimer(File f) {
		SamplePrimer sp = getSamplePrimer(f);
		if(sp!=null) {
			return sp.getPrimer();
		}
		return null;
	}

	public File getRefSeq(File f) {
		SamplePrimer sp = getSamplePrimer(f);
		if(sp!=null) {
			return sp.getRef();
		}
		return null;
	}
	public String getChr(File f) {
		SamplePrimer sp = getSamplePrimer(f);
		if(sp!=null) {
			return sp.getChr();
		}
		return null;
	}
	public boolean isLB(File f) {
		SamplePrimer sp = getSamplePrimer(f);
		if(sp!=null) {
			return sp.isLB();
		}
		System.err.println("cannot find file "+f.getAbsolutePath());
		return false;
	}

}
