package data;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SAMRecordWrap{
	private static final int ANCHORLENGTH = Translocation.ANCHORLENGTH;
	private static final int MINMAPPINGQUALITY = Translocation.MINMAPPINGQUALITY;
	private Position pos = null;
	public static final String testName = TranslocationController.testName;
	SAMRecord s;
	private ArrayList<SATag> satags = null;
	public SAMRecordWrap(SAMRecord s) {
		this.s = s;
		fillSATags();
;	}
	private void fillSATags() {
		String SATagString = this.getSATag();
		if(SATagString == null) {
			return;
		}
		String[] satagsString = SATagString.split(";");
		if(satagsString != null) {
			satags = new ArrayList<SATag>();
			for(String s: satagsString) {
				SATag sa = SATag.parseSATag(s);
				satags.add(sa);
			}
		}
	}
	/**Returns true if the SAMRecord is considered to be an anchor read
	 * 
	 * @param sam
	 * @param anchorLength
	 * @param countPartial
	 * @return
	 */
	public boolean isAnchor(SamplePrimer sp) {
		if(!s.getContig().contentEquals(sp.getChr())) {
			//only get the opposite reads
			if(s.getFirstOfPairFlag() != sp.isFirstOfPairFlag()) {
				int matchNucleotides = getMatchLength();
				if(s.getMappingQuality()>MINMAPPINGQUALITY && matchNucleotides>=150) {
					return true;
				}
			}
		}
		return false;
	}
	public boolean isPartialAnchor(int anchorLength, SamplePrimer sp) {
		if(!s.getContig().contentEquals(sp.getChr())) {
			//only get the opposite reads
			if(s.getFirstOfPairFlag() != sp.isFirstOfPairFlag()) {
				int matchNucleotides = getMatchLength();
				if(s.getMappingQuality()>MINMAPPINGQUALITY && matchNucleotides>=anchorLength) {
					return true;
				}
			}
		}
		return false;
	}
	private int getMatchLength() {
		if (s.getCigarLength()==1) {
			return s.getCigar().getLastCigarElement().getLength();
		}
		if (s.getCigarLength()==2) {
			CigarElement ce1 = s.getCigar().getCigarElement(0);
			CigarElement ce2 = s.getCigar().getCigarElement(1);
			if ((s.getReadNegativeStrandFlag()==false) && (ce1.getOperator()==CigarOperator.M)) {
				return ce1.getLength();
			}
			if ((s.getReadNegativeStrandFlag()==true) && (ce2.getOperator()==CigarOperator.M)) {
				return ce2.getLength();
			}
		}//the situation below can occur when there is a small deletion or insertion (of 1 bp) at the end of the read, especially with longer reads.
		//but it can also be that there is a filler with a piece of genomic sequence in it, causing the detection of a longer insertion or deletion, but calling the filler as genomic, which is not desirable
		if (s.getCigarLength()==4) { 
			CigarElement ce1 = s.getCigar().getCigarElement(0);
			CigarElement ce2 = s.getCigar().getCigarElement(1);
			CigarElement ce3 = s.getCigar().getCigarElement(2);
			CigarElement ce4 = s.getCigar().getCigarElement(3);
			if ((s.getReadNegativeStrandFlag()==false) && ((ce2.getOperator()==CigarOperator.I) || (ce2.getOperator()==CigarOperator.D))&& (ce1.getOperator()==CigarOperator.M) && (ce3.getOperator()==CigarOperator.M)) {
				if (ce2.getLength()==1){
					return ce1.getLength()+ce3.getLength();
				}
				if (ce2.getLength()>1){
					return ce1.getLength();
				}
			}
			if ((s.getReadNegativeStrandFlag()==true) && ((ce3.getOperator()==CigarOperator.I) || (ce3.getOperator()==CigarOperator.D))&& (ce2.getOperator()==CigarOperator.M) && (ce4.getOperator()==CigarOperator.M)) {
				if (ce3.getLength()==1){
					return ce2.getLength()+ce4.getLength();
				}
				if (ce3.getLength()>1){
					return ce4.getLength();
				}
			}
		}
		return 0;
	}
	public String getContig() {
		return s.getContig();
	}
	public int getAlignmentStart() {
		return s.getAlignmentStart();
	}
	public int getAlignmentEnd() {
		return s.getAlignmentEnd();
	}
	public boolean getFirstOfPairFlag() {
		return s.getFirstOfPairFlag();
	}
	public boolean isSecondaryOrSupplementary() {
		return s.isSecondaryOrSupplementary();
	}
	public Object getAttribute(String string) {
		return s.getAttribute(string);
	}
	public String getReadString() {
		return s.getReadString();
	}
	public boolean getReadNegativeStrandFlag() {
		return s.getReadNegativeStrandFlag();
	}
	public boolean hasAttribute(String string) {
		return s.hasAttribute(string);
	}
	public boolean getContigSATagIsContig(String chr) {
		SATag tag = satags.get(0);
		if(tag!=null) {
			return tag.contigMatches(chr);
		}
		System.err.println("requested getContigSATagIsContig from read without SA");
		return false;
	}
	public String getSATag() {
		if(s.hasAttribute("SA")) {
			return (String)s.getAttribute("SA");
		}
		return null;
	}
	public boolean isForwardSATag() {
		SATag tag = satags.get(0);
		if(tag!=null) {
			return tag.isForward();
		}
		System.err.println("requested isForwardSATag from read without SA");
		return false;
	}
	public boolean isForwardSecondSATag() {
		SATag tag = satags.get(1);
		if(tag!=null) {
			return tag.isForward();
		}
		System.err.println("requested isForwardSecondSATag from read without SA");
		return false;
	}
	public SATag getGenomicSATag(String tdna) {
		if (!getContig().equals(tdna) && s.getCigarLength()==2) {
			return new SATag(this.getContig(), this.getAlignmentStart(), !this.getReadNegativeStrandFlag(),
					this.getCigar(), this.getMappingQuality(),0);
		}
		else {
			SATag one = satags.get(0);
			if (!one.contigMatches(tdna) && one.getCigarLength()==2) {
				return one;
			}
			else {
				if(satags.size()>1) {
					SATag two = satags.get(1);
					if(!two.contigMatches(tdna) && two.getCigarLength()==2) {
						return two;
					}
				}
			}
		}
		return null;
		
	}
	public String getGenomicPart(String tdna) {
		SATag genomic = getGenomicSATag(tdna);
		if(genomic!=null) {
			int length = 0;
			if(!genomic.isForward()) {
				length = genomic.getCigar().getFirstCigarElement().getLength();
			}
			else {
				length = genomic.getCigar().getLastCigarElement().getLength();
			}
			return  this.getReadString().substring(this.getReadString().length()-length, this.getReadString().length());
		}
		return null;
	}
	//TODO: remove this method
	public int getSALength() {
		return satags.size()*6;
	}
	public String getReadName() {
		return s.getReadName();
	}
	public Cigar getCigar() {
		return s.getCigar();
	}
	public int getCigarLength() {
		return s.getCigarLength();
	}
	public int getSACigarLength() {
		SATag tag = satags.get(0);
		if(tag!=null) {
			return tag.getSACigarLength();
		}
		System.err.println("requested getSACigarLength from read without SA");
		return -1;
	}
	public int getSASecondCigarLength() {
		SATag tag = satags.get(1);
		if(tag!=null) {
			return tag.getSACigarLength();
		}
		System.err.println("requested getSASecondCigarLength from read without SA");
		return -1;
	}
	public static int countChar(String str, char c)
	{
	    int count = 0;
	    for(int i=0; i < str.length(); i++)
	    {    if(str.charAt(i) == c)
	            count++;
	    }
	    return count;
	}
	public String getCigarString() {
		return s.getCigarString();
	}
	public String getContigSATag() {
		SATag tag = satags.get(0);
		if(tag!=null) {
			return tag.getContig();
		}
		System.err.println("requested getContigSATag from read without SA");
		return null;
	}
	public String getContigSecondSATag() {
		SATag tag = satags.get(1);
		if(tag!=null) {
			return tag.getContig();
		}
		System.err.println("requested getContigSecondSATag from read without SA");
		return null;
	}
	public boolean getContigSecondSATagIsContig(String chr) {
		SATag tag = satags.get(1);
		if(tag!=null) {
			return tag.contigMatches(chr);
		}
		System.err.println("requested getContigSecondSATagIsContig from read without SA");
		return false;
	}
	public int getPosSATagEnd() {
		SATag tag = satags.get(0);
		if(tag!=null) {
			return tag.getPositionEnd();
		}
		System.err.println("requested getPosSATagEnd from read without SA");
		return -1;
	}
	public int getPosSATagEnd2() {
		SATag tag = satags.get(1);
		if(tag!=null) {
			return tag.getPositionEnd();
		}
		System.err.println("requested getPosSATagEnd2 from read without SA");
		return -1;
	}
	/**
	 * Get the position on the chromosome of the first alignment on the SA tag. 
	 * @param s
	 * @return a signed decimal integer with the chromosomal position
	 */
	public int getPosSATag() {
		SATag tag = satags.get(0);
		if(tag!=null) {
			return tag.getPosition();
		}
		System.err.println("requested getPosSATag from read without SA");
		return -1;
	}

	/**
	 * Get the position on the chromosome of the second alignment on the SA tag. 
	 * @param s
	 * @return a signed decimal integer with the chromosomal position
	 */
	public int getPosSATag2() {
		SATag tag = satags.get(1);
		if(tag!=null) {
			return tag.getPosition();
		}
		System.err.println("requested getPosSATag2 from read without SA");
		return -1;
	}
	
	public boolean getMateNegativeStrandFlag() {
		return s.getMateNegativeStrandFlag();
	}
	public int getMappingQuality() {
		return s.getMappingQuality();
	}
	public boolean isDuplicate() {
		return s.getDuplicateReadFlag();
	}
	public boolean isSecondaryAlignment() {
		return s.isSecondaryAlignment();
	}
	public void reverseComplement() {
		s.reverseComplement();
	}
	/**
	 * If the read aligns to the plasmid, but there is a secondary alignment,
	 * and the secondary alignment is not to the plasmid, then get the location in that secondary alignment.
	 * and if the secondary alignment is to the plasmid, but the tertiary alignment is not, then get the location on that tertiary alignment.
	 * If instead the read does not align to the plasmid, then get the end of the primary alignment. 
	 * @param s
	 * @param sp
	 * @return a chromosomal position integer. 
	 */
	public Position getPosition2(SamplePrimer sp) {
		if(this.pos!=null) {
			return pos;
		}
		SATag genomic = this.getGenomicSATag(sp.getChr());
		if(genomic ==null) {
			return null;
		}
		if(genomic.isForward()) {
			pos = new Position(genomic.getContig(),genomic.getPosition());
		}
		else {
			pos = new Position(genomic.getContig(), genomic.getPositionEnd());
		}
		return pos;
	}
	public String getAlignmentLocation() {
		return this.getContig()+":"+this.getAlignmentStart()+"-"+this.getAlignmentEnd()+" ("+this.getReadNegativeStrandFlag()+")";
	}
	/**
	 * Return true if either the NM (edit distance) of "src" hasn't been set, or if the number of mismatches is 0.
	 * Return false if the number of mismatches exceeds the limit. 
	 * @param srec
	 * @return true or false
	 */
	public static boolean getNMis0(SAMRecord srec) {
		if(srec.hasAttribute("NM")) {
			int tag = (Integer) srec.getAttribute("NM");
			return tag <= 0;
		}
		return true;
	}
	/**The read that is amplified using specific primers located in the T-DNA
	 * the read does not have to be primarily located in the T-DNA itself
	 * 
	 * @param sp
	 * @return
	 */
	public boolean isStartRead(SamplePrimer sp) {
		if(s.getFirstOfPairFlag()==sp.isFirstOfPairFlag()) {
			return true;
		}
		return false;
	}
}
