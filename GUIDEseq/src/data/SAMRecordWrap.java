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
	private Position pos = null;
	public static final String testName = TranslocationController.testName;
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
	public String getGenomicPart(String tdna) {
		int length = 0;
		if (!getContig().equals(tdna) && s.getCigarLength()==2) {
			if (s.getReadNegativeStrandFlag()==true){
				pos = new Position(s.getContig(),s.getAlignmentEnd());
				length = this.getCigar().getFirstCigarElement().getLength();
			}
			if (s.getReadNegativeStrandFlag()==false){
				pos = new Position(s.getContig(),s.getAlignmentStart());
				length = this.getCigar().getLastCigarElement().getLength();
			}
			return  this.getReadString().substring(this.getReadString().length()-length, this.getReadString().length());
		}
		else {
			if (!getContigSATagIsContig(tdna) && getSACigarLength()==2) {
				if (isForwardSATag()) {
					length = getSATagLastElement();
				}
				else if (!isForwardSATag()) {
					length = getSATagFirstElement();
				}
				return  this.getReadString().substring(this.getReadString().length()-length, this.getReadString().length());
			}
			else {
				if ((getSALength()>6) && !getContigSecondSATagIsContig(tdna) && getSASecondCigarLength() == 2) {
					if (isForwardSecondSATag()) {
						length = getSATag2LastElement();
					}
					else if (!isForwardSecondSATag()) {
						length = getSATag2FirstElement();
					}
					return  this.getReadString().substring(this.getReadString().length()-length, this.getReadString().length());
				}
			}
		}
		return null;
	}
	/** Find the SATag that is closest to the position that is given
	 *  Generally this is the genomic location of the junction
	 *  it will return the SATag if that is located within 20bp of the given location
	 *  
	 *  if only one SATag was found on the asked chromosome that one is returned
	 *  
	 * @param chr the contig that needs to match
	 * @param location
	 * @return
	 */
	private String getSATagLocation(String chr, int location) {
		String[] SATags = this.getSATags();
		if(SATags.length==1) {
			return SATags[0];
		}
		int countCorrectChrs = 0;
		String lastCorrectChrSA = null;
		for(String satag: SATags) {
			String contig = getContig(satag);
			if(contig.contentEquals(chr)) {
				countCorrectChrs++;
				lastCorrectChrSA = satag;
				int position = getPosition(satag);
				//recalculate position for reverse reads
				if(!isForward(satag)) {
					Cigar c = Utils.parseCigar(getCigar(satag));
					if(c.getFirstCigarElement().getOperator()==CigarOperator.M) {
						//off by 1
						position = position+c.getFirstCigarElement().getLength()-1;
					}
				}
				//distance of 20 is perhaps a bit arbitrary
				if(Math.abs(position-location)<20) {
					return satag;
				}
			}
		}
		//we only had one option, so return this one anyway
		if(countCorrectChrs==1) {
			return lastCorrectChrSA;
		}
		return null;
	}
	/**Get the correct part of the read, based on the cigar string and the mapping orientation
	 * if first element is M, take that sequence
	 * otherwise if last element is M, take that sequenc 
	 * or null if that was not found
	 * @param c
	 * @param isForward
	 * @return
	 */
	private String getReadString(Cigar c, boolean isForward) {
		if(c.getFirstCigarElement().getOperator()==CigarOperator.M) {
			String seq = null;
			if(isForward) {
				seq = this.getReadString().substring(0, c.getFirstCigarElement().getLength());
			}
			else {
				seq = this.getReadString().substring(c.getReadLength()-c.getFirstCigarElement().getLength());
			}
			return seq;
		}
		else if(c.getLastCigarElement().getOperator()==CigarOperator.M) {
			String seq = null;
			if(isForward) {
				seq = this.getReadString().substring(c.getReadLength()-c.getLastCigarElement().getLength());
			}
			else {
				seq = this.getReadString().substring(0, c.getLastCigarElement().getLength());
			}
			return seq;
		}
		return null;
	}
	private static String getCigar(String satag) {
		return satag.split(",")[3];
	}
	private static String getContig(String satag) {
		return satag.split(",")[0];
	}
	private static int getPosition(String satag) {
		String pos = satag.split(",")[1]; 
		return Integer.parseInt(pos);
	}
	private static boolean isForward(String satag) {
		String fw = satag.split(",")[2];
		return fw.contentEquals("+");
	}
	private String[] getSATags() {
		String SATag = this.getSATag();
		return SATag.split(";");
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
	public String getContigSATag() {
		String SATag = getSATag();
		if(SATag != null) {
			return SATag.split(",|;")[0];
		}
		System.err.println("requested getContigSATag from read without SA");
		return null;
	}
	public String getContigSecondSATag() {
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

	public int getSATagFirstElement() {
		String SATag = getSATag();
		if(SATag != null) {
			String SACigar = SATag.split(",|;")[3];
			String[] CigarElements = SACigar.split("S|D|M|I|H");
			int SALastInt = Integer.valueOf(CigarElements[0]);
			return SALastInt;
		}
		System.err.println("requested SA element from read without SA");
		return -1;
	}
	public int getSATagLastElement() {
		String SATag = getSATag();
		if(SATag != null) {
			String SACigar = SATag.split(",|;")[3];
			String[] CigarElements = SACigar.split("S|D|M|I|H");
			int SALength = CigarElements.length;
			int SALastInt = Integer.valueOf(CigarElements[SALength-1]);
			return SALastInt;
		}
		System.err.println("requested SA element from read without SA");
		return -1;
	}
	public int getSATag2FirstElement() {
		String SATag = getSATag();
		if(SATag != null) {
			String SACigar = SATag.split(",|;")[9];
			String[] CigarElements = SACigar.split("S|D|M|I|H");
			int SALastInt = Integer.valueOf(CigarElements[0]);
			return SALastInt;
		}
		System.err.println("requested SA element from read without SA");
		return -1;
	}
	public int getSATag2LastElement() {
		String SATag = getSATag();
		if(SATag != null) {
			String SACigar = SATag.split(",|;")[9];
			String[] CigarElements = SACigar.split("S|D|M|I|H");
			int SALength = CigarElements.length;
			int SALastInt = Integer.valueOf(CigarElements[SALength-1]);
			return SALastInt;
		}
		System.err.println("requested SA element from read without SA");
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
		
		if (!getContig().equals(sp.getChr()) && s.getCigarLength()==2) {
			if (s.getReadNegativeStrandFlag()==true){
				pos = new Position(s.getContig(),s.getAlignmentEnd());
			}
			if (s.getReadNegativeStrandFlag()==false){
				pos = new Position(s.getContig(),s.getAlignmentStart());
			}
			return pos;
		}
		else {
			if (!getContigSATagIsContig(sp.getChr()) && getSACigarLength()==2) {
				if (isForwardSATag()) {
					pos = new Position(this.getContigSATag(),getPosSATag());
				}
				else if (!isForwardSATag()) {
					pos = new Position(this.getContigSATag(),getPosSATagEnd());
				}
				return pos;
			}
			//Sometimes the primary alignment and the second alignment is both on the T-DNA plasmid if there is a large templated filler matching the plasmid.
			//in that case the tertiary alignment is on the genome, so it is there where you have to find the genomic position.
			else {
				if ((getSALength()>6) && !getContigSecondSATagIsContig(sp.getChr()) && getSASecondCigarLength() == 2) {
					if (isForwardSecondSATag()) {
						pos = new Position(this.getContigSecondSATag(),getPosSATag2()); 
					}
					else if (!isForwardSecondSATag()) {
						pos = new Position(this.getContigSecondSATag(),getPosSATagEnd2());
					}
					return pos;
				}
			}
		}
		return null;
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
