package data;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SAMRecordWrap{
	private static final int ANCHORLENGTH = Translocation.ANCHORLENGTH;
	private static final int MINMAPPINGQUALITY = Translocation.MINMAPPINGQUALITY;
	SAMRecord s;
	public SAMRecordWrap(SAMRecord s) {
		this.s = s;
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
		String contig = getContigSATag();
		if(contig!=null) {
			return contig.contentEquals(chr);
		}
		return false;
	}
	public String getSATag() {
		if(s.hasAttribute("SA")) {
			return (String)s.getAttribute("SA");
		}
		return null;
	}
	public boolean isForwardSATag() {
		String SATag = getSATag();
		if(SATag != null) {
			String signString = SATag.split(",|;")[2];
			if (signString.equals("+")){
					return true;
			}
			else {
				return false;
			}
		}
		System.err.println("requested isForwardSATag from read without SA");
		return false;
	}
	public boolean isForwardSecondSATag() {
		String SATag = getSATag();
		if(SATag != null) {
			String signString = SATag.split(",|;")[8];
			if (signString.equals("+")){
					return true;
			}
			else {
				return false;
			}
		}
		System.err.println("requested isForwardSecondSATag from read without SA");
		return false;
	}
	public int getSALength() {
		String SATag = getSATag();
		if(SATag != null) {
			String[] SAList = SATag.split(",|;");
			int SALength = SAList.length;
			return SALength;
		}
		System.err.println("requested getSALength from read without SA");
		return -1;
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
		String SATag = getSATag();
		if(SATag != null) {
			String cigarString = SATag.split(",|;")[3];
			int countS = countChar(cigarString, 'S');
			int countM = countChar(cigarString, 'M');
			int totalCount = countS+countM;
			return totalCount;
		}
		System.err.println("requested getSACigarLength from read without SA");
		return -1;
	}
	public int getSASecondCigarLength() {
		String SATag = getSATag();
		if(SATag != null) {
			String cigarString = SATag.split(",|;")[9];
			int countS = countChar(cigarString, 'S');
			int countM = countChar(cigarString, 'M');
			int totalCount = countS+countM;
			return totalCount;
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
	private String getContigSATag() {
		String SATag = getSATag();
		if(SATag != null) {
			return SATag.split(",|;")[0];
		}
		System.err.println("requested getContigSATag from read without SA");
		return null;
	}
	private String getContigSecondSATag() {
		String SATag = getSATag();
		if(SATag != null) {
			String[] split = SATag.split(",|;");
			if(split.length>6) {
				return SATag.split(",|;")[6];
			}
			return null;
		}
		System.err.println("requested getContigSecondSATag from read without SA");
		return null;
	}
	public boolean getContigSecondSATagIsContig(String chr) {
		String contig = getContigSecondSATag();
		if(contig!=null) {
			return contig.contentEquals(chr);
		}
		return false;
	}
	public int getPosSATagEnd() {
		String SATag = getSATag();
		if(SATag != null) {
			String intString = SATag.split(",|;")[1];
			int pos = Integer.parseInt(intString);
			String cigarString = SATag.split(",|;")[3];
			int indexFirstM = cigarString.indexOf("M");
			int indexFirstS = cigarString.indexOf("S");
			if (indexFirstS < indexFirstM) {
				int posSM = Integer.parseInt(cigarString.substring(indexFirstS+1, indexFirstM));
				int position = pos+posSM-1;
				return position;
			}
			if (indexFirstM < indexFirstS) {
				int posSM = Integer.parseInt(cigarString.substring(0, indexFirstM));
				int position = pos+posSM-1;
				return position;
			}
			int position =-1;
			return position;
		}
		System.err.println("requested getPosSATagEnd from read without SA");
		return -1;
	}
	public int getPosSATagEnd2() {
		String SATag = getSATag();
		if(SATag != null) {
			String intString = SATag.split(",|;")[7];
			int pos = Integer.parseInt(intString);
			String cigarString = SATag.split(",|;")[9];
			int indexFirstM = cigarString.indexOf("M");
			int indexFirstS = cigarString.indexOf("S");
			//int indexFirstH = cigarString.indexOf("H");
			if (indexFirstS < indexFirstM) {
				int posSM = Integer.parseInt(cigarString.substring(indexFirstS+1, indexFirstM));
				int position = pos+posSM-1;
				return position;
			}
			if (indexFirstM < indexFirstS) {
				int posSM = Integer.parseInt(cigarString.substring(0, indexFirstM));
				int position = pos+posSM-1;
				return position;
			}
			int position =-1;
			return position;
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
		String SATag = getSATag();
		if(SATag != null) {
			String intString = SATag.split(",|;")[1];
			return Integer.parseInt(intString);
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
		String SATag = getSATag();
		if(SATag != null) {
			String intString = SATag.split(",|;")[7];
			return Integer.parseInt(intString);
		}
		System.err.println("requested getPosSATag2 from read without SA");
		return -1;
	}
	public String getMateReferenceName() {
		return s.getMateReferenceName();
	}
	public boolean getMateNegativeStrandFlag() {
		return s.getMateNegativeStrandFlag();
	}
	public int getMappingQuality() {
		return s.getMappingQuality();
	}
	public boolean getDuplicateReadFlag() {
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
	public int getPosition2(SamplePrimer sp) {
		if (!getContig().equals(sp.getChr()) && s.getCigarLength()==2) {
			if (s.getReadNegativeStrandFlag()==true){
				return s.getAlignmentEnd();
			}
			if (s.getReadNegativeStrandFlag()==false){
				return s.getAlignmentStart();
			}
		}
		else {
			if (!getContigSATagIsContig(sp.getChr()) && getSACigarLength()==2) {
				if (isForwardSATag()) {
					return getPosSATag();
				}
				else if (!isForwardSATag()) {
					return getPosSATagEnd();
				}
			}
			else {
				if ((getSALength()>6) && !getContigSecondSATagIsContig(sp.getChr()) && getSASecondCigarLength() == 2) {
					if (isForwardSecondSATag()) {
						return getPosSATag2();
					}
					else if (!isForwardSecondSATag()) {
						return getPosSATagEnd2();
					}
				}
			}
		}
		return -1;
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
}