package data;

import java.io.File;

/**
 *@param sample (private String)
 *@param primer (private String)
 *@param chr (private String)
 *@param ref (private File)
 *@param LB (private boolean)
 *@param DNAsample (private String)
 *@param genotype (private String)
 *@param ecotype (private String)
 *@param P5 (private String)
 *@param RBisForward (private boolean)
 *@param RBpos (private int)
 *@param LBisForward (private boolean)
 *@param LBpos (private int)
 *@param file (private File)
 *@param run (private String)
 *@param UMI (private boolean)
 */
public class SamplePrimer {
	private String sample; 
	private String primer; 
	private String chr;
	private File ref;
	private boolean LB;
	private String DNAsample;
	private String genotype;
	private String ecotype;
	private String P5;
	private boolean RBisForward;
	private int RBpos;
	private boolean LBisForward;
	private int LBpos;
	private File file;
	private String run;
	private boolean UMI;
	
	public SamplePrimer(String sample, String primer, File ref, String chr, boolean LB, String DNAsample, String genotype, String ecotype, String P5, String run, boolean uMI) {
		this.sample = sample;
		this.primer = primer;
		this.ref = ref;
		this.chr = chr;
		this.LB = LB;
		this.DNAsample = DNAsample;
		this.genotype = genotype;
		this.ecotype = ecotype;
		this.P5 = P5;
		this.run = run;
		this.UMI = uMI;
	}
	public static SamplePrimer parse(String line) {
		String[] parts = line.split("\t");
		String sample = parts[0];
		String primer = parts[1].replace("TCAGACGTGTGCTCTTCCGATCT", "");
		File ref = new File(parts[2]);
		String chr = parts[3];
		String P5 = parts[5];
		String dna = parts[6];
		String genotype = parts[7];
		String ecotype = parts[8];
		String run = parts[9];
		
		boolean UMI = false;
		//if(parts.length>=17) {
		//	UMI = Boolean.parseBoolean(parts[16]);
		//}
		
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
		return new SamplePrimer(sample, primer, ref, chr, LB, dna, genotype, ecotype, P5, run, UMI);
	}
	public boolean sampleNameMatches(String str) {
		String sampleBam = sample+".sorted.bam";
		String umiSampleBam = sample+".umi.sorted.bam";
		if(sample.contentEquals(str) || sampleBam.contentEquals(str)
				|| umiSampleBam.contentEquals(str)) {
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
	/**
	 * @return the value of String "chr" of the current instance
	 */
	public String getChr() {
		return this.chr;
	}
	public boolean isLB() {
		return this.LB;
	}
	public String toString() {
		String LB = "LB";
		if(!this.LB) {
			LB = "RB";
		}
		return sample+"\t"+LB+"\t"+this.getChr()+"\t"+getPrimer()+"\t"+DNAsample+"\t"+genotype+"\t"+ecotype;
	}
	/**
	 * Returns true the P5 primer has been specified, or false if the P7 primer has been specified. Otherwise exit program.
	 * @return true or false 
	 */
	public boolean isFirstOfPairFlag() {
		if(P5.contentEquals("P5")) {
			return true;
		}
		else if(P5.contentEquals("P7")) {
			return false;
		}
		else {
			System.err.println("Either P5 or P7 has to be specified");
			System.exit(0);
		}
		return false;
	}
	/**
	 * Sets RBpos from indexRB
	 * @param pos
	 * @param rBisForward
	 */
	public void setTDNARBPos(int pos, boolean rBisForward) {
		this.RBisForward = rBisForward;
		this.RBpos = pos;
	}
	/** 
	 * @return RBpos
	 */
	public int getTDNARBPos() {
		return RBpos;
	}
	/**
	 * Sets LBpos from indexLB
	 * @param pos
	 * @param lBisForward
	 */
	public void setTDNALBPos(int pos, boolean lBisForward) {
		this.LBisForward = lBisForward;
		this.LBpos = pos;
	}
	/** 
	 * @return LBpos
	 */
	public int getTDNALBPos() {
		return LBpos;
	}
	/**
	 * 
	 * @return this.LBisForward
	 */
	public boolean getTDNALBisForward() {
		return this.LBisForward;
	}
	public boolean getTDNARBisForward() {
		return this.RBisForward;
	}
	public static int getMaxMismatches() {
		return 1;
	}
	public String getSample() {
		return sample;
	}
	public void setFile(File f) {
		this.file = f;
	}
	public File getFile() {
		return file;
	}
	public String getPrimerPart(int number) {
		if(this.primer!=null) {
			primer = this.primer.toUpperCase();
		}
		if(primer.length()>number) {
			return primer.substring(0, number);
		}
		return primer;
	}
	public long getMinSupport() {
		return 1;
	}
	public String getDNAsample() {
		return this.DNAsample;
	}
	
	public String getEcotype() {
		return this.ecotype;
	}
	
	public String getSampleString() {
		return this.getSample()+"\t"+this.genotype+"\t"+this.DNAsample;
	}
	public static String getSampleStringHeader() {
		return "Sample\tGenotype\tDNA";
	}
	public int getMinQuality() {
		return 20;
	}
	public String getRun() {
		return run;
	}
	public boolean hasUMI() {
		return UMI;
	}
}
