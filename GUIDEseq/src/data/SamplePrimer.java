package data;

import java.io.File;

public class SamplePrimer {
	private String sample, primer, chr;
	private File ref;
	public SamplePrimer(String sample, String primer, File ref, String chr) {
		this.sample = sample;
		this.primer = primer;
		this.ref = ref;
		this.chr = chr;
	}
	public static SamplePrimer parse(String line) {
		String[] parts = line.split("\t");
		String sample = parts[0];
		String primer = parts[1].replace("TCAGACGTGTGCTCTTCCGATCT", "");
		File ref = new File(parts[2]);
		String chr = parts[3];
		return new SamplePrimer(sample, primer, ref, chr);
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

}
