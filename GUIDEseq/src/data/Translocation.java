package data;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.stream.Collectors;

import dnaanalysis.InsertionSolverTwoSides;
import dnaanalysis.RandomInsertionSolverTwoSides;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

public class Translocation {
	private static final String NOGENOMIC = "NOGENOMIC";
	private static final String NOTRANSLOCATION = "NOTRANSLOCATION";
	private static final int refSize = 500;
	private static final int minSizeInsertionSolver = 6;
	public static final int ANCHORLENGTH = 50;
	public static final int MINMAPPINGQUALITY = 50;
	private static final String adapterSeq = "AGATCGGAAGAGCG"; //temporarily only 14 bp of the adapter. Later methods below should take the part they need from the full adapter sequence
	private ArrayList<SAMRecordWrap> sams;
	private HashMap<String, Integer> names;
	private String filler = "";
	private String hom = "";
	private String jType;
	private String error;
	private String ref;
	private String realPositionCounter;
	private SamplePrimer sp;
	private String warning;
	private InsertionSolverTwoSides is;
	private Consensus TDNAcons = null;
	private Consensus TDNAfullcons = null;
	private Consensus TDNAgencons = null;
	private Consensus junctionMinimumcons;
	private StringBuffer splitReads = new StringBuffer();
	private String sampleMatchingDNA = "";
	private String sampleNonMatchingDNA = "";
	private boolean multipleEvents = false;
	private int countLBWeird;
	public static final String testName = "A00379:349:HM7WFDSXY:4:2642:1108:31344";
	
	
	/**
	 * Constructor to make an instance of the Translocation class. Populate "sams" with reads (after some filtering and reverse complementing), and "names" with names of those reads.
	 * @param s
	 * @param sp
	 */
	public Translocation(SAMRecordWrap s, SamplePrimer sp) {
		sams = new ArrayList<SAMRecordWrap>();
		names = new HashMap<String, Integer>();
		this.sp = sp;
		addSam(s);
	}
	
	/**
	 * If the read aligns primarily to the plasmid, and the primary and secondary alignments in the SA tag are both the plasmid, return false.
	 * In the other cases, if the first position in the SA tag is not the plasmid, if the read is reverse and if being the first read corresponds to whether the P5 adapter is used in the T-DNA side, take the reverse complement of the CIGAR.
	 * Also do this if the second secondary alignment is genomic.
	 * Then add samrecords to "sams"
	 * and add samrecord names to "names"
	 * @param s
	 * @return true or false
	 */
	public boolean addSam(SAMRecordWrap s) {
		if((s.getContig().equals(sp.getChr()) && s.getContigSATagIsContig(sp.getChr()) && s.getContigSecondSATagIsContig(sp.getChr())) || 
				s.isSecondaryAlignment()){
			if(s.getReadName().contentEquals(testName)) { 
			System.err.println("The read "+s.getReadName()+" with CIGAR "+s.getCigarString()+" has not been added");}
			return false;
		}
		boolean takeRC = false;
		// only take RC of anchors if necessary, because this already happened for T-DNA reads.
		if((s.getFirstOfPairFlag()!=sp.isFirstOfPairFlag()) && (s.getReadNegativeStrandFlag()==false)) {
			takeRC = true;
		}
		if(s.getReadName().contentEquals(testName)) {
			System.err.println("added "+s.getReadString());
			System.err.println("added "+s.getCigarString());
		}
		if(takeRC == true) {
			s.reverseComplement();
		}
		sams.add(s);
		names.put(s.getReadName(),0);
		return true;
	}

	
	/**
	 * @return number of reads. This should be the sum of anchors and T-DNA reads.
	 */
	public int getNrSupportingReads() {
		return sams.size();
	}
	
	/**An anchor is defined as the read which does not contain the initial integration
	 * read. In our T-DNA experiments P7 is the T-DNA read.
	 * 
	 * Count #reads that perfectly map to the genome. Make a separate count of the anchors that align at least the length of "anchorlength" to the genome.
	 * NM = 0
	 * 151M (dependendent on read length) 
	 * MAPQ>50
	 * 
	 * @return
	 */
	public int getNrAnchors() {
		int count = 0;
		for(SAMRecordWrap sam: sams) {
			if(sam.isAnchor(sp)) {
				count++;
			}
		}
		return count;
	}
	
	/**
	 * @return StringBuffer "sb" containing headers.
	 */
	public static String getHeader() {
		StringBuffer sb  = new StringBuffer();
		String s = "\t";
		sb.append("Plasmid").append(s);
		sb.append("Ecotype").append(s);
		sb.append("Run").append(s);
		sb.append("DNANonMatching").append(s);
		sb.append("DNAMatching").append(s);
		sb.append("multipleEvents").append(s);
		sb.append("umis").append(s);
		sb.append("tailSize").append(s);
		sb.append("NrTailAnchors").append(s);
		sb.append("NrSupportingReads").append(s);
		sb.append("NrAnchors").append(s);
		sb.append("NrAnchorsIncludingPartial").append(s);
		sb.append("countLBWeird").append(s);
		sb.append("Chr").append(s);
		sb.append("Position").append(s);
		sb.append("RealPosition").append(s);
		sb.append("realPositionCounter").append(s);
		sb.append("DistanceToLBRB").append(s);
		sb.append("LB/RB").append(s);
		sb.append("IGVPos").append(s);
		sb.append("isForward").append(s);
		sb.append("CigarString").append(s);
		sb.append("TDNAPartMostRepeated").append(s);
		sb.append("JunctionMostRepeated").append(s);
		sb.append("GenomicPartMostRepeated").append(s);
		sb.append("Type").append(s);
		sb.append("MinimumJunctionConsensus").append(s);
		sb.append("NrSupportingMinimumJunction").append(s);
		sb.append("getFractionHighestContributorToMinimumJunction").append(s);
		sb.append("Homology").append(s);
		sb.append("HomologyLength").append(s);
		sb.append("Filler").append(s);
		sb.append("FillerLength").append(s);
		sb.append("FillerIsTemplated").append(s);
		sb.append("getLargestMatchString").append(s);
		sb.append("getSubS").append(s);
		sb.append("getSubS2").append(s);
		sb.append("getType").append(s);
		sb.append("getLengthS").append(s);
		sb.append("getPosS").append(s);
		sb.append("getFirstHit").append(s);
		sb.append("getFirstPos").append(s);
		sb.append("TotalLength").append(s);
		sb.append("RefSequence").append(s);
		sb.append("error").append(s);
		sb.append("warning").append(s);
		sb.append("isOK");
		
		
		return sb.toString();
	}
	public String toString(boolean debug) {
		//long start = System.currentTimeMillis();
		StringBuffer sb  = new StringBuffer();
		String s = "\t";
		sb.append(getPlasmid()).append(s);
		sb.append(getEcotype()).append(s);
		sb.append(getRun()).append(s);
		sb.append(sampleNonMatchingDNA).append(s);
		sb.append(sampleMatchingDNA).append(s);
		sb.append(this.multipleEvents).append(s);
		sb.append(this.getBarcodes()).append(s);
		sb.append(getTailSize()).append(s);
		sb.append(getTailAnchors()).append(s);
		sb.append(getNrSupportingReads()).append(s);
		sb.append(getNrAnchors()).append(s);
		sb.append(getNrPartialAnchors()).append(s);
		sb.append(countLBWeird).append(s);
		sb.append(getContigMate()).append(s);
		sb.append(getPosition()).append(s);
		sb.append(getRealPosition()).append(s);
		sb.append(realPositionCounter).append(s);
		sb.append(getDistanceToLBRB()).append(s);
		sb.append(sp.isLB()? "LB":"RB").append(s);
		sb.append(getIGVPos()).append(s);
		sb.append(isForwardText()).append(s);
		sb.append(getCigarString()).append(s);
		sb.append(getTDNASequenceMostRepeated()).append(s);
		sb.append(getTranslocationSequenceMostRepeated()).append(s);
		sb.append(getGenomicSequenceMostRepeated()).append(s);
		sb.append(this.jType).append(s);
		sb.append(getMinimumJunctionConsensus()).append(s);
		sb.append(getNrHighestContributorToMinimumJunction()).append(s);
		sb.append(getFractionHighestContributorToMinimumJunction()).append(s);
		sb.append(getHomology()).append(s);
		sb.append(getHomology().length()).append(s);
		sb.append(getFiller()).append(s);
		sb.append(getFiller().length()).append(s);
		if(!debug) {
			sb.append(getFillerIsTemplated(-100,100,6)).append(s);
		}
		else {
			sb.append("debug").append(s);
		}
		if(is!=null) {
			sb.append(is.getLargestMatchString()).append(s);
			sb.append(is.getSubS()).append(s);
			sb.append(is.getSubS2()).append(s);
			sb.append(is.getType()).append(s);
			sb.append(is.getLengthS()).append(s);
			sb.append(is.getPosS()).append(s);
			sb.append(is.getFirstHit()).append(s);
			sb.append(is.getFirstPos()).append(s);
		}
		else {
			sb.append("").append(s);
			sb.append("").append(s);
			sb.append("").append(s);
			sb.append("").append(s);
			sb.append("").append(s);
			sb.append("").append(s);
			sb.append("").append(s);
			sb.append("").append(s);
		}
		sb.append(-1*getHomology().length()+getFiller().length()).append(s);
		sb.append(ref).append(s);
		sb.append(error).append(s);
		sb.append(warning).append(s);
		sb.append(isOK());
		
		
		String ret = sb.toString();

		return ret;
	}
	public int getNrPartialAnchors() {
		int count = 0;
		for(SAMRecordWrap sam: sams) {
			if(sam.isPartialAnchor(ANCHORLENGTH, sp)) {
				count++;
			}
		}
		return count;
	}

	/**
	 * @return the length between the translocation position, and the furthest base from the furthest read that does not primarily match to the plasmid
	 */
	private int getTailSize() {
		int distance = 0;
		int curPosition = this.getRealPosition();
		for(SAMRecordWrap s: sams) {
			if(!s.getContig().contentEquals(sp.getChr())) {
				int disStart = Math.abs(curPosition-s.getAlignmentStart());
				int disEnd = Math.abs(curPosition-s.getAlignmentEnd());
				int maxDis = Math.max(disStart, disEnd);
				if(maxDis>distance) {
					distance = maxDis;
				}
			}
		}
		return distance;
	}
	private int getTailAnchors() {
		//int sizeOfTailAnchors = 150;
		return getNrAnchors();
	}
	
	
	private String getPlasmid() {
		return sp.getChr();
	}
	private String getRun() {
		return sp.getRun();
	}
	private String getEcotype() {
		return sp.getEcotype();
	}
	
	private String getBarcodes() {
		HashMap<String, Integer> bc = new HashMap<String, Integer>();
		for(SAMRecordWrap s: sams) {
			//only one of the two
			if(s.getFirstOfPairFlag() && !s.isSecondaryOrSupplementary()) {
				String umi = (String) s.getAttribute("RX");
				if(umi!=null) {
					if(bc.containsKey(umi)) {
						bc.put(umi,bc.get(umi)+1);
					}
					else {
						bc.put(umi,1);
					}
				}
			}
		}
		StringBuffer sb = new StringBuffer();
		for(String key: bc.keySet()) {
			if(sb.length()>0) {
				sb.append("|");
			}
			sb.append(key).append(":").append(bc.get(key));
		}
		return sb.toString();
	}
	private String getTDNASequenceMostRepeated() {
		if (this.TDNAcons.size() != 0){
		return this.TDNAcons.getMostRepeatedString();
		}
		else {
			return "cannot get most repeated";
		}
	}
	
	/** 
	 * @param sr
	 * @param replaceWithSpaces
	 * @param sp
	 * @return the part that matches with the T-DNA, excluding the genomic part or the filler.
	 */
	private static String getTDNAPart(SAMRecordWrap sr, SamplePrimer sp) {
		//get the T-DNA reads
		if (sr.getFirstOfPairFlag()== sp.isFirstOfPairFlag())  {
			CigarElement ce1 = sr.getCigar().getCigarElement(0);
			CigarElement ce2 = sr.getCigar().getCigarElement(1);
			String SATag = (String) sr.getAttribute("SA");
    		String[] SAList = SATag.split(",|;");
    		int SALength = SAList.length;
    		String SACigar = SATag.split(",|;")[3];
    		int indexFirstM = SACigar.indexOf("M");
			int indexFirstS = SACigar.indexOf("S");
			if(sr.getReadName().contentEquals(testName)) {
				System.out.println("getTDNApart - checkpoint1");
			}
			//read mapped on T-DNA, so take sequence from here
			if ((sr.getContig().equals(sp.getChr())==true) && (sr.getCigarLength()==2)) {
				if (!sr.getReadNegativeStrandFlag()) {
					if (ce1.getOperator().equals(CigarOperator.M)) {
						int mLength = ce1.getLength();
						String tdnaPart = sr.getReadString().substring(0, mLength);
						return tdnaPart;
					}
				}
				else if (sr.getReadNegativeStrandFlag()==true) {
					if (ce2.getOperator().equals(CigarOperator.M)) {
						int mLength = ce2.getLength();
						String tdnaPart = sr.getReadString().substring(0, mLength);
						return tdnaPart;
					}
				}
			}
			//read not mapped on T-DNA
			else {
				//System.out.println(sr.getReadName() + " "+sr.getSATag()+" "+sr.getCigarString());
				if (sr.getContigSATagIsContig(sp.getChr())  && sr.getSACigarLength()==2) {
					if (!sr.isForwardSATag()) {
						if (indexFirstM > indexFirstS) {
							int mLength = Integer.parseInt(SACigar.substring(indexFirstS+1, indexFirstM));
							String tdnaPart = sr.getReadString().substring(0, mLength);
							//System.out.println("taking "+0+"-"+mLength);
							//System.out.println(tdnaPart);
							return tdnaPart;
						}
					}
					if (sr.isForwardSATag()) {
						if (indexFirstM < indexFirstS) {
							int mLength = Integer.parseInt(SACigar.substring(0, indexFirstM));
							String tdnaPart = sr.getReadString().substring(0, mLength);
							//System.out.println("taking "+0+"-"+mLength);
							//System.out.println(tdnaPart);
							return tdnaPart;
						}
					}
					//System.out.println("not taking anything when read not mapped on T-DNA");
				}
				else {
					if ((SALength>6) && sr.getContigSecondSATagIsContig(sp.getChr()) && sr.getSASecondCigarLength()==2 ) {
						String SACigar2 = SATag.split(",|;")[9];
						int indexFirstM2 = SACigar2.indexOf("M");
						int indexFirstS2 = SACigar2.indexOf("S");
						if (!sr.isForwardSecondSATag()) {
							if (indexFirstS2 < indexFirstM2) {
								int mLength = Integer.parseInt(SACigar2.substring(indexFirstS2+1, indexFirstM2));
								String tdnaPart = sr.getReadString().substring(0, mLength);
								return tdnaPart;
							}
						}
						else if (sr.isForwardSecondSATag()) {
							if (indexFirstM2 < indexFirstS2) {
								int mLength = Integer.parseInt(SACigar2.substring(0, indexFirstM2));
								String tdnaPart = sr.getReadString().substring(0, mLength);
								return tdnaPart;
							}
						}
					}
				}
			}
		}
		if(sr.getReadName().contentEquals(testName)) { 
			System.err.println("The T-DNA part of read "+sr.getReadName()+" with CIGAR "+sr.getCigarString()+" will be considered null");
			}
		return null;
	}
	private static String getGenomicPart(SAMRecordWrap sr, SamplePrimer sp) {
		if (sr.getFirstOfPairFlag()== sp.isFirstOfPairFlag())  {
			CigarElement ce1 = sr.getCigar().getCigarElement(0);
			CigarElement ce2 = sr.getCigar().getCigarElement(1);
			String SATag = (String) sr.getAttribute("SA");
    		String[] SAList = SATag.split(",|;");
    		int SALength = SAList.length;
    		String SACigar = SATag.split(",|;")[3];
    		int indexFirstM = SACigar.indexOf("M");
			int indexFirstS = SACigar.indexOf("S");
			int indexLastS = SACigar.lastIndexOf("S");
			//don't allow info beyond 150bp, but allow shorter.
			int readlength = 150;
			if (sr.getReadString().length() <150) {
				readlength = sr.getReadString().length();
			}
			if(sr.getReadName().contentEquals(testName)) {
				System.out.println("getgenomicpart - checkpoint1");
			}
			if (sr.getContig().equals(sp.getChr())==true) {
				if (sr.getContigSATagIsContig(sp.getChr()) && (SALength==12) && !sr.getContigSecondSATagIsContig(sp.getChr())) {
					if(sr.getReadName().contentEquals(testName)) {
						System.out.println("getgenomicpart - checkpoint2");
					}
					String SACigar2 = SATag.split(",|;")[9];
		    		int indexFirstM2 = SACigar2.indexOf("M");
					int indexFirstS2 = SACigar2.indexOf("S");
					int indexLastS2 = SACigar2.lastIndexOf("S");
					if ((indexFirstM2 < indexFirstS2) && (indexFirstS2==indexLastS2)) {
						int mLength = Integer.parseInt(SACigar2.substring(0, indexFirstM2));
						if (sr.getReadString().length()-mLength < readlength){
							String genomicPart = sr.getReadString().substring(sr.getReadString().length()-mLength, readlength);
							return genomicPart;
						}
					}
					if ((indexFirstM2 > indexFirstS2) && (indexFirstS2==indexLastS2)) {
						if(sr.getReadName().contentEquals(testName)) {
							System.out.println("getgenomicpart - checkpoint3");
						}
						int mLength = Integer.parseInt(SACigar2.substring(indexFirstS2+1, indexFirstM2));
						if (sr.getReadString().length()-mLength < readlength){
							String genomicPart = sr.getReadString().substring(sr.getReadString().length()-mLength, readlength);
							return genomicPart;
						}
					}
				}
				if (!sr.getContigSATagIsContig(sp.getChr())){
					if (indexFirstS==indexLastS) {
						if (indexFirstM < indexFirstS) {
							if(sr.getReadName().contentEquals(testName)) {
								System.out.println("getgenomicpart - checkpoint4");
							}
							int mLength = Integer.parseInt(SACigar.substring(0, indexFirstM));
							if (sr.getReadString().length()-mLength < readlength){
								String genomicPart = sr.getReadString().substring(sr.getReadString().length()-mLength, readlength);
								return genomicPart;
							}
						}
						if (indexFirstS < indexFirstM) {
							if(sr.getReadName().contentEquals(testName)) {
								System.out.println("getgenomicpart - checkpoint5");
							}
							int mLength = Integer.parseInt(SACigar.substring(indexFirstS+1, indexFirstM));
							if (sr.getReadString().length()-mLength < readlength){
								String genomicPart = sr.getReadString().substring(sr.getReadString().length()-mLength, readlength);
								return genomicPart;
							}
						}
					}
					else if (sr.getReadString().contains(adapterSeq) && (Utils.cigarStringFollowsSMS(SACigar))) {
						if(sr.getReadName().contentEquals(testName)) {
							System.out.println("getgenomicpart - checkpoint6");
						}
						int adapterLength = sr.getReadString().length()-sr.getReadString().indexOf(adapterSeq);
						int s1Length = Integer.parseInt(SACigar.substring(0, indexFirstS));
						int s2Length = Integer.parseInt(SACigar.substring(indexFirstM+1, indexLastS));
						int mLength = Integer.parseInt(SACigar.substring(indexFirstS+1, indexFirstM));
						if (sr.isForwardSATag() && adapterLength==s2Length) {
							String genomicPart = sr.getReadString().substring(s1Length, s1Length+mLength);
							return genomicPart;
						}
						if (!sr.isForwardSATag() && adapterLength==s1Length) {
							String genomicPart = sr.getReadString().substring(s1Length, s1Length+mLength);
							return genomicPart;
						}
					}
				}
			}
			if (!sr.getContig().equals(sp.getChr())) {
				if(sr.getReadName().contentEquals(testName)) {
					System.out.println("getgenomicpart - checkpoint7");
				}
				if (sr.getCigarLength()==2) {
					if(sr.getReadName().contentEquals(testName)) {
						System.out.println("getgenomicpart - checkpoint8");
					}
					if ((ce1.getOperator().equals(CigarOperator.M) && (sr.getReadNegativeStrandFlag()==true))) {
						if(sr.getReadName().contentEquals(testName)) {
							System.out.println("getgenomicpart - checkpoint9a");
						}
						int mLength = ce1.getLength();
						if (sr.getReadString().length()-mLength < readlength){
							String genomicPart = sr.getReadString().substring(sr.getReadString().length()-mLength, readlength);
							if(sr.getReadName().contentEquals(testName)) {
								System.out.println("getgenomicpart - checkpoint9b -: " +genomicPart);
							}
							return genomicPart;
							
							
						}
					}
					if ((ce2.getOperator().equals(CigarOperator.M) && (sr.getReadNegativeStrandFlag()==false))) {
						if(sr.getReadName().contentEquals(testName)) {
							System.out.println("getgenomicpart - checkpoint10");
						}
						int mLength = ce2.getLength();
						if (sr.getReadString().length()-mLength < readlength){
							String genomicPart = sr.getReadString().substring(sr.getReadString().length()-mLength, readlength);
							return genomicPart;
						}
					}
				}
				//if there is no adapter, then the filler likely has a longer genomic piece than the genomic end itself, and then the sequence is not useable.
				if ((sr.getCigarLength()==3) && (sr.getReadString().contains(adapterSeq))) {
					if(sr.getReadName().contentEquals(testName)) {
						System.out.println("getgenomicpart - checkpoint11");
					}
					int adapterLength = sr.getReadString().length()-sr.getReadString().indexOf(adapterSeq);
					CigarElement ce3 = sr.getCigar().getCigarElement(2);
					if ((ce2.getOperator().equals(CigarOperator.M) && (sr.getReadNegativeStrandFlag()==true)) && (ce1.getLength()==adapterLength)) {
						String genomicPart = sr.getReadString().substring(ce1.getLength(), ce1.getLength()+ce2.getLength());
						return genomicPart;		
					}
					if ((ce2.getOperator().equals(CigarOperator.M) && (sr.getReadNegativeStrandFlag()==false)) && (ce3.getLength()==adapterLength)) {
						String genomicPart = sr.getReadString().substring(ce1.getLength(), ce1.getLength()+ce2.getLength());
						return genomicPart;	
					}
				}
			}
		}
		if(sr.getReadName().contentEquals(testName)) { 
			System.err.println("The genomic part of read "+sr.getReadName()+" with CIGAR "+sr.getCigarString()+" will be considered null");
		}
		return null;
	}
	private boolean getFillerIsTemplated(int start, int end, int maxTries) {
		if(this.isOK() && this.getFiller().length()>=minSizeInsertionSolver){
			//there was a bug here if the size was too small
			//get a piece left and right of the TDNA
			//TODO get the correct piece
			String left = this.getTDNASequenceRelative(start, end);
			//this value now contains the number as well
			int adjustmentLeft = start;//Integer.parseInt(left.length());
			
			String right = this.getGenomicSequenceRelative(start,end); 
			int adjustmentRight = start;//Integer.parseInt(-1);
			/*
			System.out.println("left:"+left);
			System.out.println("right:"+left);
			System.out.println("filler:"+getFiller());
			System.out.println(this.getSams().size());
			System.out.println(this.isOK());
			System.out.println(this.error);
			System.out.println(this.warning);
			System.out.println(this.getContigMate());
			*/
			
			is = new InsertionSolverTwoSides(left, right,this.getFiller(),"test");
			is.setAdjustedPositionLeft(adjustmentLeft);		
			is.setAdjustedPositionRight(adjustmentRight);
			is.search(true, true);
			is.setMaxTriesSolved(maxTries);
			is.setMinimumMatch(minSizeInsertionSolver, false);
			is.solveInsertion();
			//now determine if this is random or not
			//one peculiar thing is if the flanks overlap it is not quite fair anymore
			//if(is.getLargestMatchString().contentEquals("tgctctagccaatacgcaa")) {
			//	System.out.println(left);
			//	System.exit(0);
			//}
			
			RandomInsertionSolverTwoSides ris = new RandomInsertionSolverTwoSides(left,right, getFiller());
			boolean isFlankInsert = ris.isNonRandomInsert(0.9, is.getLargestMatch());
			return isFlankInsert;
			
		}
		return false;
	}
	/**
	 * 
	 * @param start
	 * @param end
	 * @return the sequence in the reference genome between "start" and "end"
	 */
	private String getGenomicSequenceRelative(int start, int end) {
		ReferenceSequenceFile rsf = ReferenceSequenceFileFactory.getReferenceSequenceFile(sp.getRef());
	    //SAMRecordIterator r = sr.iterator();
	    String chr = this.getContigMate();
	    int startPos = this.getRealPosition()+start;
	    int endPos = this.getRealPosition()+end;
	    ReferenceSequence rs = rsf.getSequence(chr);
	    //System.out.println(this.error);
	    //System.out.println(chr+":"+getRealPosition());
	    //System.out.println(startPos);
	    //System.out.println(endPos);
	    if(startPos<0) {
	    	startPos = 0;
	    }
	    //another safety
	    if(endPos>=rs.getBaseString().length()) {
	    	endPos = rs.getBaseString().length();
	    }
	    String part = rs.getBaseString().substring(startPos,endPos);
	    if(this.isForward()) {
	    	part = Utils.reverseComplement(part);
	    }
	    return part;
	}
	private String getTDNASequenceRelative(int start, int end) {
		ReferenceSequenceFile rsf = ReferenceSequenceFileFactory.getReferenceSequenceFile(sp.getRef());
	    //SAMRecordIterator r = sr.iterator();
	    String chr = sp.getChr();
	    
	    ReferenceSequence rs = rsf.getSequence(chr);
	    String tdna = this.getTDNASequenceMostRepeated();
	    
	    ConsensusInt con = new ConsensusInt();
	    int forwardReads = 0;
	    int reverseReads = 0;
	    
	    for(SAMRecordWrap sam: sams) {
	    	if(sam.getReadString().startsWith(tdna)) {// && sam.getContig().equals(sp.getChr())) {
	    		if(sam.getContig().contentEquals(sp.getChr())) {
		    		if(!sam.getReadNegativeStrandFlag()) {
			    		int pos = sam.getAlignmentEnd();
			    		con.add(pos);
			    		forwardReads++;
		    		}
		    		else {
		    			int pos = sam.getAlignmentStart();
		    			con.add(pos);
		    			reverseReads++;
		    		}
	    		}
	    		if(sam.getContigSATagIsContig(sp.getChr())) {
	    			if(sam.isForwardSATag()) {
    					int pos = sam.getPosSATagEnd();
    					con.add(pos);
    					forwardReads++;
    				}
    				else {
    					int pos = sam.getPosSATag();
    					con.add(pos);
    					reverseReads++;
    				}
	    		}
	    	}
	    }
	    int pos = con.getMostRepeatedInt();
	    if(pos<0) {
	    	return null;
	    }
	    int startPos = pos+start;
		int endPos = pos+end;
		String seq = rs.getBaseString().substring(startPos, endPos);
		if(reverseReads>forwardReads) {
			seq = Utils.reverseComplement(seq);
		}
		return seq;
	}
	private int getDistanceToLBRB() {
		if(sp.isLB()) {
			return this.getDistanceToLB();
		}
		return this.getDistanceToRB();
	}
	public boolean isOK() {
		return !this.hasErrors();
	}
	private String isForwardText() {
		if(this.isForward()) {
			return "forward";
		}
		return "reverse";
	}
	private String getHomology() {
		return this.hom;
	}
	
	
	/**
	 * @return this.filler
	 */
	private String getFiller() {
		return this.filler;
	}

	
	private int getDistanceToRB() {
		int RBPos = sp.getTDNARBPos();
		ArrayList<Integer> posA = new ArrayList<Integer>();
		for(SAMRecordWrap sam: sams) {
			if(sam.getFirstOfPairFlag()==false) {
				if ((sam.getContig().equals(sp.getChr())) && (sam.getCigarLength()==2)){
					if ((sp.getTDNARBisForward()==false) && (sam.getReadNegativeStrandFlag()==false)) { //A00379:349:HM7WFDSXY:4:1610:23384:24330-4 & A00379:349:HM7WFDSXY:4:2158:11451:16313 (BL29_LZB2_RB_1_30044710) & A00379:349:HM7WFDSXY:4:1218:1796:29528-3 (BL29_LZB2_RB_1_27525504)
						int distance = (sam.getAlignmentEnd())-RBPos;
						posA.add(distance);
					}
					if (sp.getTDNARBisForward() && (sam.getReadNegativeStrandFlag()==true) ) { //A00379:349:HM7WFDSXY:4:1443:4155:21120-3 (BL25_LZB1_RB_2_7137561)
						int distance = (RBPos-sam.getAlignmentStart());
						posA.add(distance);
					}
				}
				else {
					if (sam.getContigSATagIsContig(sp.getChr()) && (sam.getSACigarLength()==2)){ 
						if ((sp.getTDNARBisForward()==false) && sam.isForwardSATag() ) { //A00379:349:HM7WFDSXY:4:1258:25102:13244-2 (BL29_LZB2_RB_1_30044710) & A00379:349:HM7WFDSXY:4:1347:6406:9392-2 (BL29_LZB2_RB_1_27525504)
							int distance = sam.getPosSATagEnd()-RBPos;
							posA.add(distance);
						}
						if (sp.getTDNARBisForward() && !sam.isForwardSATag()) {//A00379:349:HM7WFDSXY:4:1250:15031:28902-1 (BL25_LZB1_RB_2_7137561)
							int distance = RBPos-sam.getPosSATag();
							posA.add(distance);
						}
					}
					else {
						if (sam.getSALength()>6 && sam.getContigSecondSATagIsContig(sp.getChr()) && sam.getSASecondCigarLength()==2){
							if ((sp.getTDNARBisForward()==false) && sam.isForwardSecondSATag()) { //
								int distance = sam.getPosSATagEnd2()-RBPos;
								posA.add(distance);
							}
							if ((sp.getTDNARBisForward()==true) && !sam.isForwardSecondSATag()) {//M02948:174:000000000-JBDYN:1:1101:17357:17713 (LZ64-1-R-4-17861618)
								int distance = RBPos-sam.getPosSATag2();
								posA.add(distance);
							}
						}
					}
				}
			}
		}
		if(posA.size()>0) {
			return ConsensusInt.consensusInt(posA);
		}
		return Integer.MIN_VALUE;
	}
	private int getDistanceToLB() {
		int LBPos = sp.getTDNALBPos();
		ArrayList<Integer> posA = new ArrayList<Integer>();
		for(SAMRecordWrap sam: sams) {
			if(sam.getFirstOfPairFlag()==false) {
				if ((sam.getContig().equals(sp.getChr())) && (sam.getCigarLength()==2)){
					if (sp.getTDNARBisForward() && !sam.getReadNegativeStrandFlag()) { 
						int distance = (sam.getAlignmentEnd())-LBPos;
						posA.add(distance);
					}
					if (!sp.getTDNARBisForward() && sam.getReadNegativeStrandFlag() ) { 
						int distance = (LBPos-sam.getAlignmentStart());
						posA.add(distance);
					}
				}
				else {
					if (sam.getContigSATagIsContig(sp.getChr()) && sam.getSACigarLength()==2){ 
						if (sp.getTDNARBisForward() && sam.isForwardSATag()) { 
							int distance = sam.getPosSATagEnd()-LBPos;
							posA.add(distance);
						}
						if (!sp.getTDNARBisForward() && !sam.isForwardSATag()) {
							int distance = LBPos-sam.getPosSATag();
							posA.add(distance);
						}
					}
					else {
						if (sam.getSALength()>6 && sam.getContigSecondSATagIsContig(sp.getChr()) && sam.getSASecondCigarLength()==2){
							if (sp.getTDNARBisForward() && sam.isForwardSecondSATag()) { 
								int distance = sam.getPosSATagEnd2()-LBPos;
								posA.add(distance);
							}
							if (!sp.getTDNARBisForward() && !sam.isForwardSecondSATag()) {
								int distance = LBPos-sam.getPosSATag2();
								posA.add(distance);
							}
						}
					}
				}
			}
		}
		if(posA.size()>0) {
			return ConsensusInt.consensusInt(posA);
		}
		return Integer.MIN_VALUE;
	}

	
	private String getGenomicSequenceMostRepeated() {
		if(TDNAgencons == null) {
			System.err.println("Can't call getGenomicSequenceMostRepeated before splitTDNAReads is called");
		}
		if (this.TDNAgencons.size() != 0){
		return this.TDNAgencons.getMostRepeatedString();
		}
		else {
			return NOTRANSLOCATION;
		}
	}
	
	/**
	 * @return fraction of full T-DNA reads that contribute to the consensus
	 */
	private double getFractionHighestContributorToMinimumJunction(){
		if(junctionMinimumcons == null) {
			System.err.println("Can't call getFractionHighestContributorToMinimumJunction before splitTDNAReads is called");
		}
		if (this.junctionMinimumcons.size() != 0){
			return junctionMinimumcons.getMostRepeatedStringFraction();
		}
		else {
			return 0;
		}
	}
	
	private double getNrHighestContributorToMinimumJunction(){
		if(junctionMinimumcons == null) {
			System.err.println("Can't call getNrHighestContributorToMinimumJunction before splitTDNAReads is called");
		}
		if (this.junctionMinimumcons.size() != 0){
			return junctionMinimumcons.getMostRepeatedStringNr();
		}
		else {
			return 0;
		}
	}
	
	public void getMinimumJunction() {
		String junction = this.getTranslocationSequenceMostRepeated();
		String genomic =  this.getGenomicSequenceMostRepeated();
		String tdna = this.getTDNASequenceMostRepeated();
		int startTDNA = junction.indexOf(tdna);
		int endTDNA = startTDNA+tdna.length();
		int startGenomic = junction.indexOf(genomic);
		int minimumJunction = startGenomic+20;
		if (minimumJunction < (endTDNA)) {
			this.addError("error: homology between T-DNA and genome >20 bp");
		}
		if(startTDNA <0 || startGenomic <0) {
			this.addError("error: either TDNA or genomic not found");
			jType= "Not determined";
		}
		else {
			if(endTDNA>=startGenomic) {
				hom = junction.substring(startGenomic, endTDNA).toUpperCase();
				jType= "NON-FILLER";
			}
			else {
				filler = junction.substring(endTDNA, startGenomic).toUpperCase();
				jType= "FILLER";
			}
		}
		ArrayList<SAMRecordWrap> tdnareads = new ArrayList<SAMRecordWrap>();
		this.junctionMinimumcons = new Consensus();
		for(SAMRecordWrap sr: sams) {
			if(sr.getFirstOfPairFlag()==sp.isFirstOfPairFlag()) {
				tdnareads.add(sr);
			}
		}
		for(SAMRecordWrap sr: tdnareads) {
			int readlength = minimumJunction;
			if (sr.getReadString().length() >= minimumJunction) {
				String readSequence = sr.getReadString().substring(0, readlength);
				junctionMinimumcons.add(readSequence);
			}
			else {
				removeSam(sr);
				if(sr.getReadName().equals(testName)) { 
					System.err.println("The read "+testName+" has been removed because it is smaller than the minimum junction");
					}
			}
		}
	}
	private String getMinimumJunctionConsensus() {
		if (junctionMinimumcons!=null && this.junctionMinimumcons.size() != 0){
			return this.junctionMinimumcons.getConsensusString();
		}
		else {
			return NOTRANSLOCATION;
		}
	}
	
	public String getTranslocationSequenceMostRepeated() {
		if (TDNAfullcons!=null && this.TDNAfullcons.size() != 0){
			return this.TDNAfullcons.getMostRepeatedString();
		}
		else {
			return NOTRANSLOCATION;
		}
	}

	private String getCigarString() {
		ArrayList<String> cigars = new ArrayList<String>();
		for(SAMRecordWrap s: sams) {
			cigars.add(s.getCigarString());
		}
		Consensus c = new Consensus(cigars);
		return c.getMostRepeatedString();
	}
	/**
	 * @return combination of the chromosome name that the mate aligns to, and the consensus position
	 */
	public String getIGVPos() {
		return this.getContigMate()+":"+this.getRealPosition();
	}
	/**
	 * if the number of supporting reads is greater than 0, for s in sams, if the mate is not primarily aligned to the plasmid, return the name of the chromosome the mate is primarily aligned to.
	 * @return the name of the chromosome the mate aligns to, or null
	 */
	public String getContigMate() {
		if(getNrSupportingReads()>0) {
			for(SAMRecordWrap s: sams) {
				if(!s.getMateReferenceName().equals(sp.getChr())) {
					return s.getMateReferenceName();
				}
			}
			for(SAMRecordWrap s: sams) {
				return s.getMateReferenceName();
			}
		}
		return null;
	}
	public boolean isForward() {
		if(getNrSupportingReads()>0) {
			return sams.get(0).getMateNegativeStrandFlag()==false;
		}
		return false;
	}
	/**
	 * get the consensus position without adding debug info
	 * @return consensus position or mate position
	 */
	public int getRealPosition() {
		return getRealPosition(false);
	}
	/**
	 * Get the position as usual from the T-DNA reads.
	 * From the anchors, also get the position. Even though the anchors usually do not match at the position, this information is written in the anchor anyway. But not always.
	 * When it is written in the anchor, this position is the same for all anchors. If for some reason the position is not given correctly, then it will be different for different anchors.
	 * Then the most occurring position will be the real position, and this will be returned. 
	 * The reason why this (calling the real position) is necessary is that sometimes with large fillers you have a bigger alignment in the filler than outside (for the T-DNA read). This causes the actual position to shift.
	 * @param debug
	 * @return corrected position.
	 */
	public int getRealPosition(boolean debug) {
		ArrayList<Integer> positions = new ArrayList<Integer>();
		for(SAMRecordWrap s:sams) {
			if(s.isStartRead(sp)) {
				int position=s.getPosition2(sp);
				if (position!=-1) {
					positions.add(position);
				}
			}
		}
		if(positions.size()>0) {
			//because the read that is not the one connected to the T-DNA
			//can also determine the consensus this might be a bit tricky
			int pos = ConsensusInt.consensusInt(positions);
			int counter = 0;
			int total = 0;
			for(Integer i: positions) {
				if(i == pos) {
					counter++;
				}
				total++;
			}
			this.realPositionCounter = counter+" of "+total;
			return pos;
		}
		this.addError("No positions found");
		return Integer.MAX_VALUE;
	}
	/**
	 * @return sams
	 */
	public ArrayList<SAMRecordWrap> getSams(){
		return sams;
	}
	/**
	 * @return chromosomal position integer of the first record in "sams"
	 */
	public int getPosition() {
		//String consensus = getCigarString();
		return sams.get(0).getPosition2(sp);
		//found a bug here!
		//return sams.get(0).getMateAlignmentStart();
	}
	
	public void splitTDNAReads() {
		ArrayList<SAMRecordWrap> tdnareads = new ArrayList<SAMRecordWrap>();
		ArrayList<String> badreads = new ArrayList<String>();
		ArrayList<SAMRecordWrap> genomicReadsAll = new ArrayList<SAMRecordWrap>();
		
		splitReads.append(sp.getSample()+"\r\n"+this.getContigMate()+":"+this.getRealPosition()+" forward: "+this.isForward()+" primer "+sp.getPrimer());
		splitReads.append("\t"+this.getNrAnchors());
		splitReads.append("\r\n\n");
		
		for(SAMRecordWrap sr: sams) {
			if(sr.isAnchor(sp)) {
				splitReads.append("anchor:\t"+sr.getReadName()+" "+sr.getAlignmentLocation()+"\n");
			}
			else if(sr.isPartialAnchor(ANCHORLENGTH, sp)) {
				splitReads.append("part. anchor:\t"+sr.getReadName()+" "+sr.getAlignmentLocation()+"\n");
			}
		}

		for(SAMRecordWrap sr: sams) {
			if(sr.getFirstOfPairFlag()==sp.isFirstOfPairFlag())  {
				tdnareads.add(sr);
			}
		}
		this.TDNAcons = new Consensus();
		this.TDNAfullcons = new Consensus();
		this.TDNAgencons = new Consensus();
		splitReads.append("\r\n\nTDNA reads - only the T-DNA part\r\n");
		for(SAMRecordWrap sr: tdnareads) {
			String readPart = Translocation.getTDNAPart(sr, sp);
			if (readPart != null) {
				TDNAcons.add(readPart);
				String cigar = Utils.getString(sr.getCigarString(),24);
				splitReads.append(cigar+readPart+" "+sr.getReadName()).append("\r\n");
				if(sr.getReadName().equals(testName)) { 
					System.out.println("The T-DNA part of read "+testName+" with cigar "+cigar+" has been added");
					}
			}
			if (readPart == null) {
				badreads.add(sr.getReadName());
				if(sr.getReadName().equals(testName)) { 
					System.err.println("The read "+testName+" will be removed because the T-DNA part is null");
					}
			}
		}
		splitReads.append(Utils.getString("consensus TDNA parts", 24)).append(TDNAcons.getConsensusString()).append("\r\n");
		splitReads.append(Utils.getString("most frequent TDNA part", 24)).append(TDNAcons.getMostRepeatedString()).append(" no. reads matching con: "+TDNAcons.getMostRepeatedStringNr()).append(" fraction reads matching con: " +TDNAcons.getMostRepeatedStringFraction());

		splitReads.append("\r\n\nTDNA reads - full reads\r\n");
		for(SAMRecordWrap sr: tdnareads) {
			if (!badreads.contains(sr.getReadName())) {
			//limit to 150bp, because different platforms were used.
			int readlength = 150;
				if (sr.getReadString().length() <150) {
					readlength = sr.getReadString().length();
				}
				String readSequence = sr.getReadString().substring(0, readlength);
				if (readSequence != null) {
					TDNAfullcons.add(readSequence);
					String cigar = Utils.getString(sr.getCigarString(),24);
					splitReads.append(cigar+readSequence+" "+sr.getReadName()).append("\r\n");
					if(sr.getReadName().equals(testName)) { 
						System.out.println("The trimmed full read "+testName+" with cigar "+cigar+" has been added");
						}
				}
				if (readSequence == null) {
					badreads.add(sr.getReadName());
					if(sr.getReadName().equals(testName)) { 
						System.err.println("The read "+testName+" will be removed because the full T-DNA read is null");
					}
				}
			
			}
		}
		splitReads.append(Utils.getString("consensus junction", 24)).append(TDNAfullcons.getConsensusString()).append("\r\n");
		splitReads.append(Utils.getString("most frequent junction", 24)).append(TDNAfullcons.getMostRepeatedString()).append(" no. reads matching con: "+TDNAfullcons.getMostRepeatedStringNr()).append(" fraction reads matching con: " +TDNAfullcons.getMostRepeatedStringFraction());
		
		splitReads.append("\r\n\nTDNA reads - only the genomic part\r\n");
		for(SAMRecordWrap sr: tdnareads) {
			if (!badreads.contains(sr.getReadName())) {
				String readPart = Translocation.getGenomicPart(sr, sp);
				if ((readPart != null) && (readPart.length()>=15)) {
					String minimumPart = readPart.substring(0, 15);
					TDNAgencons.add(minimumPart);
					String cigar = Utils.getString(sr.getCigarString(),24);
					splitReads.append(cigar+minimumPart+" "+sr.getReadName()).append("\r\n");
					if(sr.getReadName().equals(testName)) { 
						System.out.println("The trimmed genomic part of read "+testName+" with cigar "+cigar+" has been added");
						}
				}
				else {
					badreads.add(sr.getReadName());
					if(sr.getReadName().equals(testName)) { 
						System.err.println("The read "+testName+" will be removed because the genomic part is null or shorter than 15 bp");
						}
				}
			}
		}
		splitReads.append(Utils.getString("consensus genomic parts", 24)).append(TDNAgencons.getConsensusString()).append("\r\n");
		splitReads.append(Utils.getString("minimum consensus genomic parts", 24)).append(TDNAgencons.getConsensusStringMinimum()).append("\r\n");
		splitReads.append(Utils.getString("most frequent genomic", 24)).append(TDNAgencons.getMostRepeatedString()).append(" no. reads matching con: "+TDNAgencons.getMostRepeatedStringNr()).append(" fraction reads matching con: " +TDNAgencons.getMostRepeatedStringFraction());
		splitReads.append("\r\nend of TDNA\r\n");
		
		//finally remove the bad reads altogether
		this.removeSams(badreads);
		
		splitReads.append("GENOMIC reads\r\n");
		for(SAMRecordWrap sr: sams) {
			if (sr.getFirstOfPairFlag()!=sp.isFirstOfPairFlag()){
				genomicReadsAll.add(sr);
				String cigar = Utils.getString(sr.getCigarString(),24);
				splitReads.append(cigar+sr.getReadString()+" "+sr.getReadName()+"\r\n");
				if(sr.getReadName().equals(testName)) { 
					System.out.println("The anchor with name "+testName+" has been added");
					}
			}
		}
		splitReads.append("end of GENOMIC reads\r\n========\r\n");
	}
		
	/** Method for removing multiple sam records.
	 *  be very careful when applying this method as this might alter the outcomes.
	 * 
	 * @param badreads ArrayList of names that are to be removed
	 */
	private void removeSams(ArrayList<String> badreads) {
		for(String s: badreads) {
			removeSam(s);
		}
	}
	/** method to remove a read name based on the name.
	 *  Be very careful with applying this method as it will alter the outcomes of the translocaion
	 * 
	 * @param readName the name of the read to be removed
	 * @return true if the sam read has been removed
	 */
	private boolean removeSam(String readName) {
		boolean removed = false;
		for(int i=sams.size()-1;i>=0;i--) {
			if(sams.get(i).getReadName().contentEquals(readName)) {
				sams.remove(i);
				removed = true;
			}
		}
		names.remove(readName);
		return removed;
	}
	/**
	 * If "names" contains the mapping for the readName key.
	 * @param readName
	 * @return true or false
	 */
	public boolean containsRecord(String readName) {
		if(names.containsKey(readName)) {
			return true;
		}
		return false;
		
	}
	/**
	 * Add warning message "string" to object "warning"
	 * @param string
	 */
	private void addWarning(String string) {
		if(warning == null) {
			warning = string;
		}
		else {
			if(warning.length()>0) {
				warning+=":";
			}
			warning += string;
		}
	}
	private void addError(String string) {
		if(error == null) {
			error = string;
		}
		else {
			if(error.length()>0) {
				error+=":";
			}
			error += string;
		}
	}
	private boolean hasErrors() {
		if(this.getNrPartialAnchors() == 0) {
			return true;
		}
		if(this.getGenomicSequenceMostRepeated().equals(NOGENOMIC)) {
			return true;
		}
		if(this.getTranslocationSequenceMostRepeated().equals(NOTRANSLOCATION)) {
			return true;
		}
		if(!this.getTranslocationSequenceMostRepeated().contains(this.getTDNASequenceMostRepeated())) {
			return true;
		}
		if(this.error!= null && error.length()>0) {
			return true;
		}
		return false;
	}
	/**Adds the reference sequence based on the calculated junction position from realPosition
	 * takes into account the orientation of the event
	 * 
	 * @param rsf
	 */
	public void addRefSequence(ReferenceSequenceFile rsf) {
		int start = -1;
		int end = -1;
		String genomicSequenceMostRepeated = getGenomicSequenceMostRepeated().toUpperCase();
		if (!genomicSequenceMostRepeated.contentEquals(NOGENOMIC)) {
			//these local variables are to prevent additional calls to these expensive methods
			//reuse them and don't overwrite
			int realPosition = this.getRealPosition();
			String contigMate = this.getContigMate();
			ReferenceSequence seq = rsf.getSequence(contigMate);
			
			if(this.isForward()) {
				start = realPosition-refSize;
				//safety
				if(start<1) {
					start = 1;
				}
				end = realPosition;
				if (end <= seq.length()) {
					this.ref = rsf.getSubsequenceAt(contigMate, start, end).getBaseString();
					//take the reverse complement
					//bug, the ref also needs to be updated
					ref = Utils.reverseComplement(this.ref);
					if(ref.startsWith(genomicSequenceMostRepeated)==false) {
						int getMismatches = getMismatches(ref,genomicSequenceMostRepeated);
						if(getMismatches <= SamplePrimer.getMaxMismatches()) {
							addWarning("Ref sequence deviates from sequence in reads: "+getMismatches+" mismatches");
						}
						else {
							this.addError("Ref sequence deviates from sequence in reads: "+getMismatches+" mismatches");
						}
					}
				}
				else if (end > seq.length()) {
					this.addError("Probably trying to get a coordinate from the wrong chromosome");
				}
			}
			else {
				start = realPosition;
				end = realPosition+refSize;
				if (start <= rsf.getSequence(contigMate).length()) {
					if(end>seq.length()) {
						end = seq.length();
					}
					this.ref = rsf.getSubsequenceAt(contigMate, start, end).getBaseString();
					String compRef = this.ref;
					if(compRef.startsWith(genomicSequenceMostRepeated)==false) {
						int getMismatches = getMismatches(compRef,genomicSequenceMostRepeated);
						if(getMismatches <= SamplePrimer.getMaxMismatches()) {
							addWarning("Ref sequence from "+this.getContigMate() +" deviates from sequence in reads: "+getMismatches+" mismatches");
						}
						else {
							this.addError("Ref sequence from "+this.getContigMate() +" deviates from sequence in reads: "+getMismatches+" mismatches. compRef= "+compRef+", genomic= " +this.getGenomicSequenceMostRepeated());
						}
					}
				}
				else if (start > seq.length()) {
					this.addError("Trying to get a coordinate from the wrong chromosome");
				}
			}
		}
		else {
			this.addError("Cannot get refseq because no genomic seq found");
		}
	}
	
	private int getMismatches(String ref, String g) {
		int len = Math.min(ref.length(), g.length());
		int mismatches = 0;
		for(int i=0;i<len;i++) {
			if(ref.charAt(i)!=g.charAt(i)) {
				mismatches++;
			}
		}
		return mismatches;
	}
	
	public void printDebug() {
		String tl = this.toString();
		for(SAMRecordWrap sr: sams) {
			tl+="\n\t";
			tl+=sr.getReadString()+" "+sr.getContig()+":"+sr.getAlignmentStart()+" "+sr.getCigarString() +" "+sr.getFirstOfPairFlag()+" "+sr.getMappingQuality();// +" "+sr.getAttributes().toString();
		}
		System.out.println(tl);
		System.out.println("TDNA SEQ:");
		System.out.println(this.getTDNASequenceMostRepeated());
		System.out.println("GENOMIC SEQ:");
		System.out.println(this.getGenomicSequenceMostRepeated());
		System.out.println(getRealPosition(true));
		System.out.println(this.getTranslocationSequenceMostRepeated());
		
	} 
	
	public String getReads() {
		return splitReads.toString();
	}
	public void addFoundOtherSampleMatchingDNA(String sample) {
		if(sampleMatchingDNA.length()!=0) {
			this.sampleMatchingDNA+=",";
		}
		this.sampleMatchingDNA+=sample;
	}
	public void addFoundOtherSampleNonMatchingDNA(String sample) {
		if(sampleNonMatchingDNA.length()!=0) {
			this.sampleNonMatchingDNA+=",";
		}
		this.sampleNonMatchingDNA+=sample;
	}
	public void setMultipleEvents() {
		this.multipleEvents = true;
	}
	public boolean removeSam(SAMRecordWrap srec) {
		return this.removeSam(srec.getReadName());
	}

	public String matchEnd(String string) {
	    //Create a return string that is blank.
	    String ret = "";
	    //Create the maximum amount that we are going 
	    //to go both backward and forward in the string
	    int maxInt = (string.length() / 2) + 1;
	    //We start at iteration 0 for the first check and go up to the maxInt.
	    for (int i = 0; i < maxInt; i++)
	    {
	      //Create a temporary string will change that we will check for the value.
	      String temp = string.substring(0, i);
	      //If the temporary string is at the beginning and end, we set our return
	      //String equal to it.
	      if (string.endsWith(temp))
	      {
	        ret = temp;
	      }
	    }
	    //What we end up with is the return string that has the most characters
	    //At the the end.
	    return ret;
	  }
}

