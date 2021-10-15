package data;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;

public class SATag {
	private String contig;
	private int position;
	private boolean forward;
	private Cigar cigar;
	private int mapq;
	private int notSure;
	
	public String getContig() {
		return contig;
	}
	public int getPosition() {
		return position;
	}
	public boolean isForward() {
		return forward;
	}
	public Cigar getCigar() {
		return cigar;
	}
	public int getMapq() {
		return mapq;
	}
	public int getNotSure() {
		return notSure;
	}
	
	public SATag(String c, int pos, boolean fw, Cigar cig, int mapq, int i) {
		this.contig = c;
		this.position = pos;
		this.forward = fw;
		this.cigar = cig;
		this.mapq = mapq;
		this.notSure = i;
	}
	public static SATag parseSATag(String s) {
		if(s==null) {
			return null;
		}
		String[] parts = s.split(",");
		Cigar c = Utils.parseCigar(parts[3]);
		int pos = Integer.parseInt(parts[1]);
		boolean forward = parts[2].contentEquals("+")?true:false;
		int mapq = Integer.parseInt(parts[4]);
		int other = Integer.parseInt(parts[5]);
		String contig = parts[0];
		return new SATag(contig, pos, forward, c, mapq, other);
	}
	public boolean contigMatches(String s) {
		return this.getContig().contentEquals(s);
	}
	public int getCigarLength() {
		return cigar.getCigarElements().size();
	}
	public String toString() {
		String s = ",";
		return contig+s+position+s+forward+s+cigar;
	}
	public boolean equals(SATag s) {
		if(s.getContig().contentEquals(contig)) {
			if(s.getPosition() == this.position) {
				if(s.isForward()==this.forward) {
					if(this.cigar.toString().contentEquals(s.getCigar().toString())) {
						return true;
					}
				}
			}
		}
		return false;
	}
	public int getPositionEnd() {
		if(this.cigar.getCigarElement(0).getOperator()==CigarOperator.M) {
			int pos = this.position+cigar.getCigarElement(0).getLength()-1;
			return pos;
		}
		else if(this.cigar.getCigarElement(1).getOperator()==CigarOperator.M) {
			int pos = this.position+cigar.getCigarElement(1).getLength()-1;
			return pos;
		}
		System.err.println("hier");
		return -1;
	}
}
