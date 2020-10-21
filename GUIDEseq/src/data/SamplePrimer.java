package data;

import java.io.File;

public class SamplePrimer {
	private String sample, primer, chr;
	private File ref;
	private boolean LB;
	public SamplePrimer(String sample, String primer, File ref, String chr, boolean LB) {
		this.sample = sample;
		this.primer = primer;
		this.ref = ref;
		this.chr = chr;
		this.LB = LB;
	}
	public static SamplePrimer parse(String line) {
		String[] parts = line.split("\t");
		String sample = parts[0];
		String primer = parts[1].replace("TCAGACGTGTGCTCTTCCGATCT", "");
		File ref = new File(parts[2]);
		String chr = parts[3];
		boolean LB = true;
		if(parts[4].contentEquals("LB")) {
			LB = true;
		}
		else if(parts[4].contentEquals("RB")) { 
			LB = false;
		}
		else {
			System.err.println("LB or RB should be filled in: "+line);
			System.exit(0);
		}
		return new SamplePrimer(sample, primer, ref, chr, LB);
	}
	public boolean sampleNameMatches(String str) {
		if(sample.contentEquals(str)) {
			return true;
		}
		return false;
	}
	public String getPrimer() {
		// TODO Auto-generated method stub
		return primer;
	}
	public File getRef() {
		return ref;
	}
	public String getChr() {
		return this.chr;
	}
	public boolean isLB() {
		return this.LB;
	}

}
