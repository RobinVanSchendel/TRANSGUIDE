package data;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import dnaanalysis.InsertionSolverTwoSides;
import dnaanalysis.RandomInsertionSolverTwoSides;
import htsjdk.samtools.*;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

public class Translocation {
	private static final String NOGENOMIC = "NOGENOMIC";
	private static final String NOTRANSLOCATION = "NOTRANSLOCATION";
	private static final int refSize = 500;
	private static final int minSizeInsertionSolver = 6;
	private static final int ANCHORLENGTH = 50;
	private static final int MINMAPPINGQUALITY = 50;
	public ArrayList<SAMRecord> sams;
	public HashMap<String, Integer> names;
	private String filler = "";
	private String hom = "";
	private String error;
	private String ref;
	private String realPositionCounter;
	private SamplePrimer sp;
	private String warning;
	private InsertionSolverTwoSides is;
	private Consensus TDNAcons;
	private Consensus genomeCons;
	private Consensus genomePrimerCons;
	private ConsensusInt genomeConsInt;
	private Consensus junctionAll;
	private String sampleMatchingDNA = "";
	private String sampleNonMatchingDNA = "";
	private boolean multipleEvents = false;
	private int countLBWeird;

	public Translocation(SAMRecord s, SamplePrimer sp) {
		sams = new ArrayList<SAMRecord>();
		names = new HashMap<String, Integer>();
		this.sp = sp;
		addSam(s);
	}
	public boolean addSam(SAMRecord s) {
		//do not add sams that have a primary and secondary alignment in contig one
		if(s.getContig().equals(sp.getChr()) && getContigSATagIsContig(s,sp.getChr())) {
			return false;
		}
		
		byte[] quals = s.getBaseQualities();
		String qualString = s.getBaseQualityString();
		byte minQual = Translocation.getMinQuality(s);
		int i=0;
		for(byte b: quals) {
			//System.out.println(b+" "+qualString.charAt(i));
			i++;
		}
		//System.out.println(minQual);
		if(minQual<20) {
			//return false;
		}
		//System.exit(0);
		//orientation should be correct, otherwise don't add
		//remove this filter for now
		//if(s.getContig().equals(SAMReader.str) && s.getReadNegativeStrandFlag() == SAMReader.forwardRB) {
		//	return;
		//}
		//all others can be added
		boolean takeRC = false;
		if(!s.getContig().contentEquals(sp.getChr()) && s.getReadNegativeStrandFlag()!=true 
				&& s.getFirstOfPairFlag()!=sp.isFirstOfPairFlag()) {
			takeRC = true;
		}
		else if(!s.getContig().contentEquals(sp.getChr()) && s.getReadNegativeStrandFlag()==true 
				&& s.getFirstOfPairFlag()==sp.isFirstOfPairFlag()) {
			takeRC = true;
		}
		if(takeRC) {
			s.setAttribute("OC", s.getCigarString());
			s = Translocation.reverseComplementSAMRecord(s);
		}
		
		sams.add(s);
		names.put(s.getReadName(),0);
		return true;
	}
	private static byte getMinQuality(SAMRecord s) {
		byte min = Byte.MAX_VALUE;
		for(byte b: s.getBaseQualities()) {
			if(b<min) {
				min = b;
			}
		}
		return min;
	}
	private static SAMRecord reverseComplementSAMRecord(SAMRecord sr) {
		SAMRecord srec = sr.deepCopy();
		String cigar = srec.getCigarString();
		boolean currentNegativeStrandFlag = srec.getReadNegativeStrandFlag();
		
		srec.reverseComplement();
		//bug so reverse it myself
		Cigar tempCigar = srec.getCigar();
		
		//bug so reverse it myself
		if(srec.getReadNegativeStrandFlag()==currentNegativeStrandFlag) {
			srec.setReadNegativeStrandFlag(!currentNegativeStrandFlag);
		}
		
		if(tempCigar.numCigarElements()>1 && cigar.equals(srec.getCigarString())) {
			Cigar rev = new Cigar();
			//reverse the cigar
			for(int i = tempCigar.numCigarElements()-1;i>=0;i--) {
				rev.add(tempCigar.getCigarElement(i));
			}
			srec.setCigar(rev);
		}
		return srec;
	}
	public int getNrSupportingReads() {
		return sams.size();
	}
	/**normal way to get the number of anchor
	 * 
	 * @param addPartial
	 * @return
	 */
	public int getNrAnchors(boolean addPartial) {
		return getNrAnchors(addPartial, ANCHORLENGTH);
	}
	/**An anchor is defined as the read which does not contain the initial integration
	 * read. In our T-DNA experiments P7 is the T-DNA read.
	 * 
	 * Count #reads that perfectly map to the genome.
	 * NM = 0
	 * 151M (dependendent on read length) 
	 * MAPQ>50
	 * 
	 * @return
	 */
	public int getNrAnchors(boolean addPartial, int anchorLength) {
		int count = 0;
		for(SAMRecord sam: sams) {
			if(!sam.getContig().contentEquals(sp.getChr())) {
				//only get the opposite reads
				if(sam.getFirstOfPairFlag() != sp.isFirstOfPairFlag()) {
					//System.out.println("Still here first "+sam.getReadName() + " "+Translocation.getNMis0(sam) +" "+ sam.getAttribute("NM"));
					//no mismatches
					int matchNucleotides = Translocation.getMatchLength(sam);
					//System.out.println(this.getPosition());
					//System.out.println(matchNucleotides+" "+sam.getMappingQuality());
					//System.out.println(sam.getCigarString());
					if(sam.getMappingQuality()>MINMAPPINGQUALITY && matchNucleotides>=anchorLength) {
						if(this.getIGVPos().contentEquals("4:5085682") && anchorLength==150) {
							System.out.println(sam.toString());
							System.out.println(sam.getDuplicateReadFlag()); 
						}
						
						count++;
					}
				}
				//only if not counted already
				else if(addPartial) {
					int matchNucleotides = Translocation.getMatchLength(sam);
					if(sam.getMappingQuality()>MINMAPPINGQUALITY && matchNucleotides>=anchorLength) {
						count++;
					}
				}
			}
		}
		return count;
	}
	//get the the anchor lenght of this read, which is the first or last cigar, which has to be M and then the size
	private static int getMatchLength(SAMRecord sam) {
		CigarElement last = sam.getCigar().getLastCigarElement();
		//require a match at the end, otherwise it is not an anchor
		//all reads are now turned such that the cigar far right is the last one
		if(last.getOperator() == CigarOperator.M){
			return last.getLength();
		}
		/*
		else {
			System.err.println("What do we have here");
			System.err.println(sam.getCigarString());
			System.err.println(last);
			System.err.println(sam.toString());
		}
		*/
		return 0;
	}
	public int getSizePrimary() {
		int count = 0;
		for(SAMRecord s: sams) {
			if(!s.isSecondaryAlignment()) {
				count++;
			}
		}
		return count;
	}
	public int getSizeSecondary() {
		int count = 0;
		for(SAMRecord s: sams) {
			if(s.isSecondaryAlignment()) {
				count++;
			}
		}
		return count;
	}
	public static String getHeader() {
		StringBuffer sb  = new StringBuffer();
		String s = "\t";
		sb.append("DNANonMatching").append(s);
		sb.append("DNAMatching").append(s);
		sb.append("multipleEvents").append(s);
		sb.append("umis").append(s);
		sb.append("tailSize").append(s);
		sb.append("NrTailAnchors").append(s);
		sb.append("NrSupportingReads").append(s);
		sb.append("NrAnchors").append(s);
		sb.append("NrAnchorsIncludingPartial").append(s);
		sb.append("NrSupportJunction").append(s);
		sb.append("NrNotSupportJunction").append(s);
		sb.append("NrSupportJunctionHQ").append(s);
		sb.append("NrNotSupportJunctionGQ").append(s);
		sb.append("countLBWeird").append(s);
		sb.append("NrPrimaryAlignment").append(s);
		sb.append("NrSecondaryAlignment").append(s);
		sb.append("Chr").append(s);
		sb.append("Position").append(s);
		sb.append("RealPosition").append(s);
		sb.append("realPositionCounter").append(s);
		sb.append("DistanceToLBRB").append(s);
		sb.append("LB/RB").append(s);
		sb.append("IGVPos").append(s);
		sb.append("isForward").append(s);
		sb.append("CigarString").append(s);
		sb.append("JunctionSequence").append(s);
		sb.append("TDNASequence").append(s);
		sb.append("GenomicSequence").append(s);
		sb.append("getTDNASequenceNew").append(s);
		sb.append("getGenomicSequenceNew").append(s);
		sb.append("getGenomicSequenceIntNew").append(s);
		sb.append("getGenomicSequenceIntMinMaxNew").append(s);
		sb.append("getGenomicSequenceIntMedian").append(s);
		sb.append("getGenomicSequenceIntMean").append(s);
		sb.append("getGenomicSequenceIntSD").append(s);
		sb.append("getMostRepeatedPreGenome").append(s);
		sb.append("getMostRepeatedConsensus").append(s);
		sb.append("Type").append(s);
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
		sb.append("isOK").append(s);
		sb.append("getTDNASequenceDiffReads").append(s);
		sb.append("getTDNASequenceDiffReadsHighestContributor");
		
		return sb.toString();
	}
	public String toString(boolean debug) {
		//long start = System.currentTimeMillis();
		StringBuffer sb  = new StringBuffer();
		String s = "\t";
		sb.append(sampleNonMatchingDNA).append(s);
		sb.append(sampleMatchingDNA).append(s);
		sb.append(this.multipleEvents).append(s);
		sb.append(this.getBarcodes()).append(s);
		sb.append(getTailSize()).append(s);
		sb.append(getTailAnchors()).append(s);
		sb.append(getNrSupportingReads()).append(s);
		sb.append(getNrAnchors(false)).append(s);
		sb.append(getNrAnchors(true)).append(s);
		sb.append(getNrSupportJunction(false)).append(s);
		sb.append(this.getNrNotSupportJunction(false)).append(s);
		sb.append(getNrSupportJunction(true)).append(s);
		sb.append(this.getNrNotSupportJunction(true)).append(s);
		sb.append(countLBWeird).append(s);
		sb.append(getSizePrimary()).append(s);
		sb.append(getSizeSecondary()).append(s);
		sb.append(getContigMate()).append(s);
		sb.append(getPosition()).append(s);
		sb.append(getRealPosition()).append(s);
		sb.append(realPositionCounter).append(s);
		sb.append(getDistanceToLBRB()).append(s);
		sb.append(sp.isLB()? "LB":"RB").append(s);
		sb.append(getIGVPos()).append(s);
		sb.append(isForwardText()).append(s);
		sb.append(getCigarString()).append(s);
		sb.append(getTranslocationSequence()).append(s);
		sb.append(getTDNASequence(false)).append(s);
		sb.append(getGenomicSequence()).append(s);
		sb.append(getTDNASequenceNew()).append(s);
		sb.append(getGenomicSequenceNew()).append(s);
		sb.append(getGenomicSequenceIntNew()).append(s);
		sb.append(getGenomicSequenceIntMinMaxNew()).append(s);
		sb.append(this.genomeConsInt.getMedian()).append(s);
		sb.append(this.genomeConsInt.getMean()).append(s);
		sb.append(this.genomeConsInt.getSD()).append(s);
		sb.append(junctionAll.getMostRepeatedString()).append(s);
		sb.append(junctionAll.getConsensusString()).append(s);
		sb.append(getType()).append(s);
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
		sb.append(isOK()).append(s);
		sb.append(getTDNASequenceDiffReads()).append(s);
		sb.append(getTDNASequenceDiffReadsHighestContributor());
		//sb.append(getFillerSequenceDiffReads());
		//sb.append(getFillerSequenceDiffReadsHighestContributor()).append(s);
		
		//sb.append("\n");

		//for(SAMRecord sam: sams) {
			//sb.append("\t"+sam.getSAMString());
		//}
		
		String ret = sb.toString();
		//long end = System.currentTimeMillis();
		//long duration = end-start;
		/*
		System.out.println(duration +" nrSams: "+this.sams.size());
		if(this.sams.size()==7441) {
			for(SAMRecord sam: sams) {
				System.out.println(sam.getReadString());
				System.out.println(sam.getContig()+":"+sam.getAlignmentStart()+"-"+sam.getAlignmentEnd());
			}
			System.out.println(sp.getTDNALBPos());
			System.out.println(sp.getTDNARBPos());
			System.exit(0);
		}
		*/
		return ret;
	}
	private int getTailSize() {
		int distance = 0;
		int curPosition = this.getRealPosition();
		for(SAMRecord s: sams) {
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
		int sizeOfTailAnchors = 150;
		return getNrAnchors(false,sizeOfTailAnchors);
	}
	private String getBarcodes() {
		HashMap<String, Integer> bc = new HashMap<String, Integer>();
		for(SAMRecord s: sams) {
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
	private String getGenomicSequenceIntMinMaxNew() {
		if(this.genomeConsInt.getMin()== this.genomeConsInt.getMax()) {
			return this.genomeConsInt.getMin()+"";
		}
		return this.genomeConsInt.getMin()+"-"+this.genomeConsInt.getMax();
	}
	private String getGenomicSequenceIntNew() {
		return this.genomeConsInt.getMostRepeatedInt()+" found "+this.genomeConsInt.getMostRepeatedIntNr()+" from "+this.genomeConsInt.size();
	}
	private String getGenomicSequenceNew() {
		return this.genomeCons.getConsensusString();
	}
	private String getTDNASequenceNew() {
		return this.TDNAcons.getConsensusString();
	}
	private String getFillerSequenceDiffReads() {
		HashMap<String, Integer> seqsDiff = new HashMap<String, Integer>();
		String tdna = this.getTDNASequence();
		//System.out.println("getFillerSequenceDiffReads");
		//System.out.println(tdna);
		ArrayList<String> al = new ArrayList<String>();
		for(SAMRecord sr: sams) {
			if(!sr.getContig().contentEquals(sp.getChr())) {
				String part = getSoftClippedPart(sr);
				CigarElement ce = sr.getCigar().getFirstCigarElement();
				//take the RC here so that the filler is correctly aligned
				//if(ce.getOperator() == CigarOperator.S) {
				//	String part = Utils.reverseComplement(part);
				
				//al.add(Utils.reverseComplement(part));
				//System.out.println("\t"+sr.getCigarString()+"\t"+sr.getReadString());
//				System.out.println(sr.getReadName());
	//			System.out.println(part);
			}
		}
		for(String s: al) {
			//System.out.println(s);
		}
		//Integer.parseInt("s");
		//System.exit(0);
		
		return null;
	}
	private static String getSoftClippedPart(SAMRecord sr) {
		if(sr.getCigarLength()==2) {
			CigarElement ce = sr.getCigar().getCigarElement(0);
			CigarElement ce2 = sr.getCigar().getCigarElement(1);
			if(ce.getOperator() == CigarOperator.M && ce2.getOperator() == CigarOperator.S) {
				return sr.getReadString().substring(ce.getLength());
			}
			else if(ce.getOperator() == CigarOperator.S && ce2.getOperator() == CigarOperator.M) {
				return sr.getReadString().substring(0, ce.getLength());
			}
		}
		return null;
	}
	private static String getMatchedPart(SAMRecord sr, boolean replaceWithSpaces) {
		if(sr.getCigarLength()==2) {
			CigarElement ce = sr.getCigar().getCigarElement(0);
			CigarElement ce2 = sr.getCigar().getCigarElement(1);
			if(ce.getOperator() == CigarOperator.S && ce2.getOperator() == CigarOperator.M) {
				String sPart = sr.getReadString().substring(0, ce.getLength());
				if(replaceWithSpaces) {
					return Translocation.getSpaces(sPart.length()) +  sr.getReadString().substring(ce.getLength());
				}
				return sr.getReadString().substring(ce.getLength());
			}
			else if(ce.getOperator() == CigarOperator.M && ce2.getOperator() == CigarOperator.S) {
				String mPart = sr.getReadString().substring(ce.getLength());
				if(replaceWithSpaces) {
					return sr.getReadString().substring(0,ce.getLength()) + Translocation.getSpaces(mPart.length()) ;
				}
				return sr.getReadString().substring(0,ce.getLength());
			}
			//hard clipped can be returned as is
			else if(ce.getOperator()==CigarOperator.H || ce2.getOperator()==CigarOperator.H) {
				return sr.getReadString();
			}
			System.err.println("weird reads..."+sr.getCigarString());
		}
		else if(sr.getCigarLength()==1) {
			CigarElement ce = sr.getCigar().getCigarElement(0);
			if(ce.getOperator() == CigarOperator.M) {
				return sr.getReadString();
			}
		}
		return null;
	}
	private static String getSpaces(int length) {
		StringBuffer sb = new StringBuffer(length);
		for(int i=0;i<length;i++) {
			sb.append(" ");
		}
		return sb.toString();
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
	    String tdna = this.getTDNASequence();
	    for(SAMRecord sam: sams) {
	    	if(sam.getReadString().startsWith(tdna) && sam.getContig().equals(sp.getChr())) {
	    		if(!sam.getReadNegativeStrandFlag()) {
		    		int pos = sam.getAlignmentEnd();
		    		//System.out.println(pos);
		    		int startPos = pos+start;
		    		int endPos = pos+end;
		    		String seq = rs.getBaseString().substring(startPos, endPos);
		    		//System.out.println(sam.getCigarString());
		    		//System.out.println(sam.getReadString());
		    		//System.out.println(sam.getReadNegativeStrandFlag());
		    		//System.out.println(seq);
		    		return seq;
	    		}
	    		//reverse strand, junction is at start of alignment
	    		else {
	    			int pos = sam.getAlignmentStart();
		    		//System.out.println(pos);
		    		int startPos = pos+start;
		    		int endPos = pos+end;
		    		String seq = rs.getBaseString().substring(startPos, endPos);
		    		//need the rc
		    		seq = Utils.reverseComplement(seq);
		    		//System.out.println(sam.getCigarString());
		    		//System.out.println(sam.getReadString());
		    		//System.out.println(sam.getReadNegativeStrandFlag());
		    		//System.out.println(seq);
		    		return seq;
	    		}
	    	}
	    }
	    
	    
		return null;
	}
	private int getDistanceToLBRB() {
		if(sp.isLB()) {
			return this.getDistanceToLB();
		}
		return this.getDistanceToRB();
	}
	
	private int getDistanceToRB() {
		String seq = this.getTDNASequence();
		int RBPos = sp.getTDNARBPos();
		ArrayList<Integer> posA = new ArrayList<Integer>();
		for(SAMRecord sam: sams) {
			if(sam.getReadString().startsWith(seq) && sam.getContig().equals(sp.getChr())) {
				if(!sp.getTDNARBisForward() && !sam.getReadNegativeStrandFlag()) {
					int distance = sam.getAlignmentEnd()-RBPos;
					posA.add(distance);
				}
				else if(sp.getTDNARBisForward() && !sam.getReadNegativeStrandFlag()) {
					int distance = RBPos-sam.getAlignmentStart();
					posA.add(distance);
				}
				else {
					//Probably better to know if sample is directed at LB or RB
					return Integer.MAX_VALUE;
					/*
					System.err.println(sp.getTDNARBisForward());
					System.err.println(sp.getPrimer());
					System.err.println(sam.getReadNegativeStrandFlag());
					System.err.println(sam.getReadString());
					System.err.println("RB: Some logical erro in my thinking...");
					System.exit(0);
					*/
				}
			}
		}
		if(posA.size()>0) {
			return consensusInt(posA);
		}
		return Integer.MIN_VALUE;
	}
	private String getTDNASequence() {
		return this.getTDNASequence(false);
	}
	public boolean isOK() {
		return !this.hasErrors();
	}
	public String isLBorRB() {
		int disLB = this.getDistanceToLB();
		int disRB = this.getDistanceToRB();
		int disLBAbs = Math.abs(disLB);
		int disRBAbs = Math.abs(disRB);
		if(disLBAbs<disRBAbs) {
			return "LB";
		}
		return "RB";
	}
	private int getNrSupportJunction(boolean qual) {
		if(this.hasErrors()) {
			return -1;
		}
		String junctionMin = this.getMinimalJunction();
		int count = 0;
		for(SAMRecord s: sams) {
			if(qual) {
				int minQual = Translocation.getMinQuality(s);
				if(minQual<sp.getMinQuality()) {
					continue;
				}
			}
			//not enough to match only TDNA
			if(this.getTDNASequence().contentEquals(s.getReadString())) {
				continue;
			}
			//only the ones which are supposed to be equal to the junction
			if(s.getFirstOfPairFlag()==sp.isFirstOfPairFlag()) {
				String seq = s.getReadString();
				if(seq.contains(junctionMin) || junctionMin.contains(seq)) {
					count++;
				}
			}
		}
		return count;
	}
	private int getNrNotSupportJunction(boolean qual) {
		if(this.hasErrors()) {
			return -1;
		}
		String junctionMin = this.getMinimalJunction();
		String genome =  this.getGenomicSequence();
		int count = 0;
		if(this.getRealPosition()==8019692) {
			System.out.println(junctionMin);
			System.out.println("Comparing the following seqs");
		}
		countLBWeird=0;
		for(SAMRecord s: sams) {
			if(qual) {
				int minQual = Translocation.getMinQuality(s);
				if(minQual<sp.getMinQuality()) {
					continue;
				}
			}
			String seq = s.getReadString();
			//stuff that matches genome is ok
			if(genome.startsWith(s.getReadString()) || s.getReadString().startsWith(genome)) {
				continue;
			}
			//only the ones which are supposed to be equal to the junction
			if(s.getFirstOfPairFlag()==sp.isFirstOfPairFlag()) {
				if(!seq.contains(junctionMin) && !junctionMin.contains(seq)) {
					if(this.getRealPosition()==8019692) {
						System.out.println(Translocation.getMinQuality(s)+"\t"+s.getCigarString()+"\t"+seq+"\t"+s.toString()+"\t"+s.getMateReferenceName()+":"+s.getMateAlignmentStart() );
					}
					count++;
				}
				else if(!seq.startsWith(sp.getPrimer()) && s.getCigar().getFirstCigarElement().getOperator()!=CigarOperator.H) {
					if(this.getRealPosition()==8019692) {
						System.out.println("x"+Translocation.getMinQuality(s)+"\t"+s.getCigarString()+"\t"+seq);
					}
					count++;
				}
				if(qual && seq.startsWith("GATAATTCAATTCGGCGTTAATTCAGTACATTAAAAACGTCCGCAATGTGTTACTAGATCGACCGGCATGCAAGCT")) {
					this.countLBWeird++;
				}
			}
		}
		if(this.getRealPosition()==8019692) {
			//System.exit(0);
		}
		return count;
	}

	private String isForwardText() {
		if(this.isForward()) {
			return "forward";
		}
		return "reverse";
	}
	private String getMinimalJunction() {
		return getMinimalJunction(false);
	}
	/**get only the junction, 30bp of TDNA and 30bp of genome (plus filler)
	 * 
	 * @return
	 */
	private String getMinimalJunction(boolean debug) {
		int size = 30;
		int locJunction = -1;
		String tdnaOrig = this.getTDNASequence();
		String junction = this.getTranslocationSequence(debug);
		String genomicOrig = this.getGenomicSequence();
		//might overshoot otherwise
		String tdna = Utils.longestCommonSubstring(tdnaOrig, junction);
		String genomic = Utils.longestCommonSubstring(genomicOrig, junction);
		//System.out.println("genomic: "+genomic);
		//System.out.println("genomicOrig: "+genomicOrig);
		//System.out.println(genomic.indexOf(genomicOrig));
		//junction has to start with tdna
		if(junction.indexOf(tdna)!=0) {
			addError("TDNA not ok");
		}
		if(genomicOrig.indexOf(genomic)!=0) {
			locJunction = findGenomicPart(genomicOrig, junction);
			if(debug) {
				System.out.println("getMinimalJunctionDEBUG");
				System.out.println(genomicOrig);
				System.out.println(genomic);
				System.out.println(junction);
				System.out.println("locJunction: "+locJunction);
			}
			if(locJunction==-1) {
				addError("genomic not ok");
			}
			//System.out.println(genomicOrig);
			//System.out.println(genomic);
		}
		
		if(genomic.length()>size) {
			genomic = genomic.substring(0, size);
		}
		int start = 0;
		int end = junction.indexOf(genomic)+genomic.length();
		//overwrite with value calculated earlier
		if(locJunction!=-1) {
			end = junction.length()-locJunction;
		}
		if(end<start) {
			return junction+"INSERT_TOO_LONG?";
		}
		return junction.substring(start, end);
	}
	private int findGenomicPart(String genomicOrig, String junction) {
		//speed up?
		int location = junction.indexOf(genomicOrig);
		if(location!=-1) {
			return location;
		}
		
		int max = Math.min(genomicOrig.length(), junction.length());
		int maxHit = -1;
		//BUG: should be a minumum of 5 nucleotides
		for(int i=5;i<=max;i++) {
			String seqJ = junction.substring(junction.length()-i);
			if(genomicOrig.startsWith(seqJ)) {
				maxHit = i;
			}
		}
		return maxHit;
	}
	private String getHomology() {
		return this.hom;
	}
	private String getFiller() {
		return this.filler;
	}
	private String getType() {
		if(this.hasErrors()) {
			return "Probably not correct";
		}
		if(this.getNrAnchors(false)==0) {
			return "Probably not correct";
		}
		String junction = this.getTranslocationSequence();
		String genomic = this.getGenomicSequence();
		String tdnaLCS = Utils.longestCommonSubstring(this.getTDNASequence(), junction);
		int location = findGenomicPart(genomic, junction);
		String genomicLCS = Utils.longestCommonSubstring(genomic, junction);
		//System.out.println(location);
		//System.out.println(genomicLCS);
		//System.out.println(junction);
		//System.out.println(getGenomicSequence());
		//System.out.println(tdnaLCS);
		
		int startTDNA = junction.indexOf(tdnaLCS);
		int endTDNA = startTDNA+tdnaLCS.length();
		int startGenomic = -1;
		if(!genomic.startsWith(genomicLCS) && location != -1) {
			startGenomic = junction.length()-location;
		}
		else {
			startGenomic = junction.indexOf(genomicLCS);
		}
		if(endTDNA <0 || startGenomic <0) {
			return "Probably not correct";
		}
		if(endTDNA>=startGenomic) {
			hom = this.getTranslocationSequence().substring(startGenomic, endTDNA);
			return "NON-FILLER";
		}
		else {
			filler = this.getTranslocationSequence().substring(endTDNA, startGenomic);
			return "FILLER";
		}
	}
	private int getDistanceToLB() {
		String seq = this.getTDNASequence();
		int LBPos = sp.getTDNALBPos();
		ArrayList<Integer> posA = new ArrayList<Integer>();
		for(SAMRecord sam: sams) {
			if(sam.getReadString().startsWith(seq) && sam.getContig().equals(sp.getChr())) {
				//System.out.println(sam.getReadNegativeStrandFlag());
				//System.out.println(seq);
				
				if(sp.getTDNALBisForward() && !sam.getReadNegativeStrandFlag()) {
					int distance = sam.getAlignmentEnd()-LBPos;
					posA.add(distance);
				}
				else if(!sp.getTDNALBisForward() && !sam.getReadNegativeStrandFlag()) {
					int distance = LBPos-sam.getAlignmentStart();
					posA.add(distance);
				}
				else {
					//probably this is a read at the RB side
					return Integer.MAX_VALUE;
					/*
					System.err.println(sp.getTDNALBisForward());
					System.err.println(sp.getPrimer());
					System.err.println(sam.getReadNegativeStrandFlag());
					System.err.println(sam.getReadString());
					
					System.err.println("LB: Some logical erro in my thinking...");
					System.exit(0);
					*/
				}
			}
		}
		if(posA.size()>0) {
			return consensusInt(posA);
		}
		return Integer.MIN_VALUE;
	}
	private String getGenomicSequence() {
		return getGenomicSequence(false);
	}

	private String getGenomicSequence(boolean debug) {
		//search for the sequences that give the genomic sequence
		ArrayList<String> seqs = new ArrayList<String>();
	
		for(SAMRecord srec: sams) {
			if(!srec.isSecondaryAlignment() && srec.getFirstOfPairFlag() == sp.isFirstOfPairFlag()
					&& cigarStringFollowsMSH(srec.getCigarString()) 
					&& !srec.getContig().equals(sp.getChr())
					&& getNMis0(srec)) {
				if(debug) {
					System.out.println(srec.getCigarString()+" "+srec.getDuplicateReadFlag() + " "+srec.getAlignmentStart()+"-"+srec.getAlignmentEnd() + " "+ srec.getMappingQuality()+ " "+srec.getReadName());
				}
				//add the correct part of the seq
				CigarElement ce = srec.getCigar().getCigarElement(0);
				CigarElement ce2 = srec.getCigar().getCigarElement(0);
				if((ce.getOperator() == CigarOperator.S || ce.getOperator() == CigarOperator.H) &&
						ce2.getOperator() == CigarOperator.M){
					if(debug) {
						System.out.println("SEQ1\t"+srec.getReadString().substring(ce.getLength()));
					}
					String part = "";
					if(ce.getOperator() == CigarOperator.H) {
						part = srec.getReadString();
					}
					else {
						part = srec.getReadString().substring(ce.getLength());
					}
					int end = 150;
					//System.out.println(srec.getCigarString());
					//trim to 150
					if(end>part.length()) {
						end = part.length();
					}
					part = part.substring(0, end);
					seqs.add(part);
				}
				if(debug) {
					System.out.println(" adding "+srec.getContig());
					System.out.println(srec.getReadString());
				}
			}
			else if(debug) {
				/*
				System.out.println("not survived filer");
				System.out.println(srec.getCigarString());
				System.out.println(!srec.isSecondaryAlignment());
				System.out.println(srec.getFirstOfPairFlag() == sp.isFirstOfPairFlag());
				System.out.println(!srec.getContig().equals(sp.getChr()));
				System.out.println(getNMis0(srec));
				System.out.println(srec.getReadString());
				*/
			}
		}
		if(debug) {
			System.out.println("size genomic "+seqs.size() +" number of sam records "+sams.size());
			for(String s: seqs) {
				System.out.println("\t\t"+s);
			}
		}
		//sometimes this contains one read or so, so don't rely on this then
		//5 is an arbitrary number
		if(seqs.size()>=5) {
			Consensus c = new Consensus(seqs);
			String seq = c.getMostRepeatedString();
			//System.out.println("CONSENSUS");
			//System.out.println(seq);
			//System.out.println("CONSENSUS");
			return seq;
		}
		
		//some location have only the second read mapping to the genome
		int maxGenomeSizePart = 50;
		for(SAMRecord srec: sams) {
			if(debug) {
				System.out.println("SAMRecord Iteration");
				System.out.println(srec.getCigarString());
				System.out.println(!srec.isSecondaryAlignment());
				System.out.println(srec.getFirstOfPairFlag() != sp.isFirstOfPairFlag());
				System.out.println(srec.getContig());
				System.out.println(srec.getFirstOfPairFlag());
			}
			if(!srec.isSecondaryAlignment() && srec.getFirstOfPairFlag() != sp.isFirstOfPairFlag()
					&& cigarStringFollowsMSH(srec.getCigarString()) 
					&& !srec.getContig().equals(sp.getChr())) {
				if(debug) {
					System.out.println("sec:");
					System.out.println(srec.getCigarString());
					System.out.println(srec.getReadString());
				}
				CigarElement ce = srec.getCigar().getCigarElement(0);
				if(ce.getOperator() == CigarOperator.S) {
					String seq = srec.getReadString().substring(ce.getLength());
					if(seq.length()>maxGenomeSizePart) {
						seq = seq.substring(0,maxGenomeSizePart);
					}
					if(debug) {
						System.out.println("seq2 :"+seq);
					}
					seqs.add(seq);
				}
				else {
					String seq = srec.getReadString().substring(0,ce.getLength());
					if(debug) {
						System.out.println("seq1 :\t"+seq);
					}
					if(seq.length()>maxGenomeSizePart) {
						seq = seq.substring(seq.length()-maxGenomeSizePart);
					}
					seq = Utils.reverseComplement(seq);
					seqs.add(seq);
					if(debug) {
						System.out.println("seq1rc:\t"+seq);
					}
					
				}
				
			}
		}
		if(debug) {
			System.out.println("size genomic "+seqs.size());
			for(String seq: seqs) {
				System.out.println("\t"+seq);
			}
		}
		if(seqs.size()>0) {
			Consensus c = new Consensus(seqs);
			String seq = c.getMostRepeatedString();
			return seq;
		}
		//give up
		return NOGENOMIC;
	}
	private String getTDNASequenceDiffReads(){
		HashMap<String, Integer> seqsDiff = new HashMap<String, Integer>();
		//simplest case if second in pair and primary alignment
		
		for(SAMRecord sam: sams) {
			if(!sam.isSecondaryAlignment() && sam.getContig().equals(sp.getChr())) {
				//up unto the M
				CigarElement ce = sam.getCigar().getCigarElement(0);
				//always correct?
				if(ce.getOperator() != CigarOperator.M && sam.getCigar().numCigarElements() == 2) {
					//System.err.println("bug still here");
					//System.out.println(sam.getCigarString());
					//System.out.println(sam.getReadString());
				}
				else {
					int pos = ce.getLength();
					String seq = sam.getReadString().substring(0, pos); 
					if(!seqsDiff.containsKey(seq)) {
						seqsDiff.put(seq, 0);
					}
					seqsDiff.put(seq,seqsDiff.get(seq)+1);
				}
			}
			if(sam.isSecondaryAlignment() && sam.getContig().equals(sp.getChr())) {
				String seq = sam.getReadString();
				if(!seqsDiff.containsKey(seq)) {
					seqsDiff.put(seq, 0);
				}
				seqsDiff.put(seq,seqsDiff.get(seq)+1);
			}
		}
		double total = 0;
		ArrayList<Integer> al = new ArrayList<Integer>();
		for(String key: seqsDiff.keySet()) {
			total+= seqsDiff.get(key);
			al.add(seqsDiff.get(key));
		}
		//sort, highest first
		Collections.sort(al, Collections.reverseOrder());
		String highestTwo = "";
		if(al.size()>=2) {
			String highest = ""+al.get(0)/total;
			String secondHighest = "|"+al.get(1)/total;
			highestTwo = highest+secondHighest;
			String one = null;
			String two = null;
			for(String key: seqsDiff.keySet()) {
				if(seqsDiff.get(key)==al.get(0)) {
					one = key;
				}
				if(seqsDiff.get(key)==al.get(1)) {
					two = key;
				}
			}
			int sizeDiff = Math.abs(one.length()-two.length());
			highestTwo+="|"+sizeDiff;
			
		}
		return seqsDiff.size()+"|"+total+"|"+highestTwo;
	}
	private double getTDNASequenceDiffReadsHighestContributor(){
		HashMap<String, Integer> seqsDiff = new HashMap<String, Integer>();
		//simplest case if second in pair and primary alignment
		
		for(SAMRecord sam: sams) {
			if(!sam.isSecondaryAlignment() && sam.getContig().equals(sp.getChr())) {
				//up unto the M
				CigarElement ce = sam.getCigar().getCigarElement(0);
				//always correct?
				if(ce.getOperator() != CigarOperator.M && sam.getCigar().numCigarElements() == 2) {
					//System.err.println("bug still here");
					//System.out.println(sam.getCigarString());
					//System.out.println(sam.getReadString());
				}
				else {
					int pos = ce.getLength();
					String seq = sam.getReadString().substring(0, pos); 
					if(!seqsDiff.containsKey(seq)) {
						seqsDiff.put(seq, 0);
					}
					seqsDiff.put(seq,seqsDiff.get(seq)+1);
				}
			}
			//most go here!
			if(sam.isSecondaryAlignment() && sam.getContig().equals(sp.getChr())) {
				String seq = sam.getReadString();
				if(!seqsDiff.containsKey(seq)) {
					seqsDiff.put(seq, 0);
				}
				seqsDiff.put(seq,seqsDiff.get(seq)+1);
			}
		}
		double total = 0;
		int max = 0;
		for(String key: seqsDiff.keySet()) {
			total+= seqsDiff.get(key);
			if(seqsDiff.get(key)>max) {
				max = seqsDiff.get(key);
			}
		}
		if(total>0) {
			return max/total;
		}
		return -1;
	}
	private String getTDNASequence(boolean debug) {
		ArrayList<String> seqs = new ArrayList<String>();
		HashMap<String, Integer> seqsDiff = new HashMap<String, Integer>();
		//simplest case if second in pair and primary alignment
		if(debug) {
			System.out.println("####TDNA adding:");
		}
		for(SAMRecord sam: sams) {
			if(!sam.isSecondaryAlignment() && sam.getContig().equals(sp.getChr())) {
				//up unto the M
				CigarElement ce = sam.getCigar().getCigarElement(0);
				//only do it if this is not the case
				if(ce.getOperator() == CigarOperator.M) {
					int pos = ce.getLength();
					String seq = sam.getReadString().substring(0, pos); 
					seqs.add(seq);
					if(!seqsDiff.containsKey(seq)) {
						seqsDiff.put(seq, 0);
					}
					seqsDiff.put(seq,seqsDiff.get(seq)+1);
					if(debug) {
						System.out.println("putting1 "+seq+" "+sam.getCigarString() +" "+sam.getReadString()+" "+sam.getContig()+" "+sam.getReadName());						
					}
				}
			}
			//most go here!
			if(sam.isSecondaryAlignment() && sam.getContig().equals(sp.getChr())) {
				String seq = sam.getReadString();
				seqs.add(seq);
				if(!seqsDiff.containsKey(seq)) {
					seqsDiff.put(seq, 0);
				}
				seqsDiff.put(seq,seqsDiff.get(seq)+1);
				if(debug) {
					System.out.println("putting2 "+seq+" "+sam.getCigarString()  );						
				}
				/*
				if(this.getPosition()==2923296) {
					for(SAMRecord sam2: sams) {
						System.out.println(sam2.isSecondaryAlignment()+"\t"+sam2.getReadString());
						System.out.println(sam2.getCigarString());
						System.out.println(sam2.getReadName());
						System.out.println(sam2.getFirstOfPairFlag());
						System.out.println(sam2.getMappingQuality());
					}
					System.out.println(this.getNrAnchors(false));
					//System.exit(0);
				}
				*/
			}
		}
		if(debug) {
			System.out.println("####TDNAdiff");
			for(String s: seqsDiff.keySet()) {
				System.out.println(seqsDiff.get(s)+"\t"+s);
			}
			System.out.println("####TDNAdiff end");
		}
		if(seqs.size()>0) {
			Consensus c = new Consensus(seqs);
			String seq = c.getMostRepeatedString();
			return seq;
		}
		return "unknown";
	}
	public String getTranslocationSequence() {
		return getTranslocationSequence(false);
	}

	public String getTranslocationSequence(boolean debug) {
		ArrayList<String> seqs = new ArrayList<>();
		if(debug) {
			System.out.println("getTranslocationSequence".toUpperCase());
		}
		for(SAMRecord s:sams) {
			if(!s.isSecondaryAlignment() && cigarStringFollowsMSH(s.getCigarString())) {
				if(s.getReadString().startsWith(sp.getPrimer())) {
					seqs.add(s.getReadString());
					if(debug) {
						System.out.println("\t"+s.getReadString());
						System.out.println("\t"+s.getContig());
						System.out.println("\t"+s.getCigarString());
						System.out.println("\t"+s.getReadName());
					}
				}
				if(debug) {
					//System.out.println(s.getReadString());
					//System.out.println(s.getFirstOfPairFlag());
					//System.out.println(s.getContig());
					//System.out.println(s.getCigarString());
					//System.out.println(s.getReadString().startsWith(sp.getPrimer()));
				}
			}
			if(debug) {
				//System.out.println("HIER!\t"+s.getReadString());
				//System.out.println(s.getCigarString());
				//System.out.println(!s.isSecondaryAlignment());
			}
		}
		if(debug) {
			for(String s: seqs) {
				System.out.println("\t"+s.length()+"\t"+s);
			}
			System.out.println("currently have: "+seqs.size());
		}
		if(seqs.size()>0) {
			Consensus c = new Consensus(seqs);
			String seq = c.getMostRepeatedString();
			if(debug) {
				System.out.println("CONSEN\t"+seq);
			}
			//only sequences with primer included are now added, so simply return
			return seq;
		}
		if(debug) {
			System.out.println("END OF getTranslocationSequence".toUpperCase());
		}
		return NOTRANSLOCATION;
	}
	private static boolean cigarStringFollowsMSH(String cigarString) {
		Pattern p = Pattern.compile("\\d*[MSH]\\d*[MSH]");
		Matcher m = p.matcher(cigarString);
		boolean b = m.matches();
		return b;
	}

	private String getCigarString() {
		ArrayList<String> cigars = new ArrayList<String>();
		for(SAMRecord s: sams) {
			cigars.add(s.getCigarString());
		}
		Consensus c = new Consensus(cigars);
		return c.getMostRepeatedString();
	}
	//TODO this consensus string likely only works if we align
	//TODO the sequences correctly
	//20190506 my solution is to put only substrings of the seq in the list
	//this way I can avoid the problem of having tiling seqs
	private static Consensus consensusString(ArrayList<String> list) {
		Consensus c = new Consensus(list);
		return c;
	}
	private static Integer consensusInt(ArrayList<Integer> list) {
		Integer mostRepeatedWord 
	    = list.stream()
	          .collect(Collectors.groupingBy(w -> w, Collectors.counting()))
	          .entrySet()
	          .stream()
	          .max(Comparator.comparing(Entry::getValue))
	          .get()
	          .getKey();
		return mostRepeatedWord;
	}

	public String getIGVPos() {
		return this.getContigMate()+":"+this.getRealPosition();
	}

	public String getContigMate() {
		if(getNrSupportingReads()>0) {
			for(SAMRecord s: sams) {
				if(!s.getMateReferenceName().equals(sp.getChr())) {
					return s.getMateReferenceName();
				}
			}
			//System.err.println("Hier gaat ie naar hetzelfde chromosoom");
			for(SAMRecord s: sams) {
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
	public int getRealPosition() {
		return getRealPosition(false);
	}
	public int getRealPosition(boolean debug) {
		ArrayList<Integer> seqs = new ArrayList<Integer>();
		for(SAMRecord s:sams) {
			if(!s.isSecondaryAlignment() && cigarStringFollowsMSH(s.getCigarString()) 
					&& !s.getContig().equals(sp.getChr())) {
				//no mismatches
				if(getNMis0(s)) {
					//maybe get the actual position
					String cigar = s.getCigarString();
					if(s.getAttribute("OC")!=null) {
						cigar = (String) s.getAttribute("OC");
					}
					int indexS = cigar.indexOf("S");
					int indexM = cigar.indexOf("M");
					if(indexS>=0 && indexS<indexM) {
						seqs.add(s.getAlignmentStart());
						if(debug) {
							System.out.println(s.getAlignmentStart()+" "+getNMis0(s));
						}
					}
					else {
						seqs.add(s.getAlignmentEnd());
						if(debug) {
							System.out.println(s.getAlignmentEnd()+" "+getNMis0(s));
						}

					}
					
				}
				else if(debug) {
					System.out.println(s.getAlignmentStart()+" "+getNMis0(s)+" not used");
				}
			}
		}
		if(seqs.size()>0) {
			//because the read that is not the one connected to the T-DNA
			//can also determine the consensus this might be a bit tricky
			int pos = consensusInt(seqs);
			int counter = 0;
			int total = 0;
			for(Integer i: seqs) {
				if(i == pos) {
					counter++;
				}
				total++;
			}
			this.realPositionCounter = counter+" from "+total;
			return pos;
		}
		return sams.get(0).getMateAlignmentStart();
		
	}
	public ArrayList<SAMRecord> getSams(){
		return sams;
	}
	
	public int getPosition() {
		//String consensus = getCigarString();
		return Translocation.getPosition(sams.get(0), sp);
		//found a bug here!
		//return sams.get(0).getMateAlignmentStart();
	}
	public static int getPosition(SAMRecord s, SamplePrimer sp) {
		//info is in the SA tag
		if(s.isSecondaryAlignment() && !getContigSATagIsContig(s, sp.getChr())) {
			int position = getPosSATag(s);
			//System.out.println("SA Tag position "+position );
			return position;
		}
		//rough estimate is in the mate
		return s.getMateAlignmentStart();
	}

	private static boolean getContigSATagIsContig(SAMRecord s, String str) {
		if(s.isSecondaryAlignment()) {
			String SATag = (String) s.getAttribute("SA");
			String contigString = SATag.split(",")[0];
			return contigString.equals(str);
		}
		return false;
	}


	private static int getPosSATag(SAMRecord s) {
		if(s.isSecondaryAlignment()) {
			String SATag = (String) s.getAttribute("SA");
			String intString = SATag.split(",")[1];
			return Integer.parseInt(intString);
		}
		return Integer.MIN_VALUE;
	}
	public boolean containsRecord(String readName) {
		if(names.containsKey(readName)) {
			return true;
		}
		return false;
		
	}
	public static boolean getNMis0(SAMRecord srec) {
		if(srec.hasAttribute("NM")) {
			int tag = (Integer) srec.getAttribute("NM");
			/*
			if(tag>0) {
				String md = (String) srec.getAttribute("MD");
				String[] parts = md.split("[\\^A-Z]");
				//System.out.println(parts[0]);
				int size = Integer.parseInt(parts[0]);
				//TODO change to non hardcoded
				if(size>=150) {
					return true;
				}
			}
			*/
			//System.out.println("Mismatches: "+tag);
			return tag <= SamplePrimer.getMaxMismatches();
			//return tag == 0;
		}
		return true;
	}
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
		if(this.getNrAnchors(false) == 0) {
			return true;
		}
		if(this.getGenomicSequence().equals(NOGENOMIC)) {
			return true;
		}
		if(this.getTranslocationSequence().equals(NOTRANSLOCATION)) {
			return true;
		}
		if(!this.getTranslocationSequence().contains(this.getTDNASequence())) {
			return true;
		}
		if(this.error!= null && error.length()>0) {
			return true;
		}
		return false;
	}
	public void addRefSequence(ReferenceSequenceFile rsf) {
		if(this.isOK()) {
			int start = -1;
			int end = -1;
			if(this.isForward()) {
				start = this.getRealPosition()-refSize;
				end = this.getRealPosition();
				this.ref = rsf.getSubsequenceAt(this.getContigMate(), start, end).getBaseString();
				//take the reverse complement
				ref = Utils.reverseComplement(ref);
			}
			else {
				start = this.getRealPosition();
				end = this.getRealPosition()+refSize;
				if(end>rsf.getSequence(this.getContigMate()).length()) {
					end = rsf.getSequence(this.getContigMate()).length();
				}
				this.ref = rsf.getSubsequenceAt(this.getContigMate(), start, end).getBaseString();
			}
			if(!this.ref.startsWith(this.getGenomicSequence())) {
				int getMismatches = getMismatches(ref,this.getGenomicSequence());
				if(getMismatches <= SamplePrimer.getMaxMismatches()) {
					addWarning("Ref sequence has "+getMismatches+" mismatches with genomic sequence");
				}
				else {
					this.addError("Ref sequence deviates from sequence in reads "+getMismatches+" mismatches");
				}
			}
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
		for(SAMRecord sr: sams) {
			tl+="\n\t";
			tl+=sr.getReadString()+" "+sr.getContig()+":"+sr.getAlignmentStart()+" "+sr.getCigarString() +" "+sr.getFirstOfPairFlag()+" "+sr.getMappingQuality();// +" "+sr.getAttributes().toString();
		}
		System.out.println(tl);
		System.out.println("TDNA SEQ:");
		System.out.println(this.getTDNASequence(true));
		System.out.println("GENOMIC SEQ:");
		System.out.println(this.getGenomicSequence(true));
		System.out.println("getMinimalJunction:");
		System.out.println(getMinimalJunction(true));
		System.out.println(getRealPosition(true));
		
	}
	/**Debug method
	 * 
	 */
	public String getReads() {
		StringBuffer sb = new StringBuffer();
		sb.append(sp.getSample()+"\r\n"+this.getContigMate()+":"+this.getRealPosition()+" forward: "+this.isForward()+" primer "+sp.getPrimer());
		sb.append("\t"+this.getNrAnchors(false));
		
		sb.append("\r\n\t");
		
//		Assembler a = new Assembler();
		ArrayList<SAMRecord> tdnareads = new ArrayList<SAMRecord>();
		ArrayList<SAMRecord> genomicReads = new ArrayList<SAMRecord>();
		ArrayList<SAMRecord> genomicReadsIncPrimer = new ArrayList<SAMRecord>();
		ArrayList<SAMRecord> genomicReadsAll = new ArrayList<SAMRecord>();
		for(SAMRecord sr: sams) {
			if(sr.getContig().contentEquals(sp.getChr())) {
				tdnareads.add(sr);
			}
			else if(!sr.getContig().contentEquals(sp.getChr()) && (sr.getReadString().startsWith(sp.getPrimer()) ||
					sr.getCigar().getFirstCigarElement().getOperator() ==CigarOperator.H)) {
				genomicReadsIncPrimer.add(sr);
			}
			else if(!sr.getContig().contentEquals(sp.getChr())) {
				genomicReads.add(sr);
			}
			else {
				System.err.println("What's this?");
				System.err.println(sr.getReadString());
			}
			if((sr.getFirstOfPairFlag()==sp.isFirstOfPairFlag())
					) {
					genomicReadsAll.add(sr)	;
			}
			//boolean startWithPrimer = sr.getReadString().startsWith(sp.getPrimer());
			//String aligned = getAligned(sr);
			//sb.append("\t"+startWithPrimer);
			//sb.append("\t"+sr.getCigarString());
			//sb.append("\t"+sr.getFirstOfPairFlag());
			//sb.append(Translocation.getString(sr.getContig(),15));
			//sb.append(sr.getReadString());
			//sb.append("\t"+sr.getReadNegativeStrandFlag());
			//sb.append("\t"+sr.getAlignmentStart()+"-"+sr.getAlignmentEnd());
			//sb.append("\t"+sr.getReadName());
			//sb.append("\r\n\t");
			//a.addFragment(new Fragment(sr.getReadString()));
		}
		sb.append("\r\n\nTDNA reads\r\n");
		//Assembler a = new Assembler();
		this.TDNAcons = new Consensus();
		this.genomeCons = new Consensus();
		this.genomeConsInt = new ConsensusInt();
		this.genomePrimerCons = new Consensus();
		this.junctionAll = new Consensus();
		for(SAMRecord sr: tdnareads) {
			String cigar = Translocation.getString(sr.getCigarString(),24);
			//String map = Translocation.getString(""+sr.getFirstOfPairFlag(),4);
			//cigar = map+cigar;
			String readPart = Translocation.getMatchedPart(sr, false);
			TDNAcons.add(readPart);
			sb.append(cigar+readPart).append("\r\n");
		}
		sb.append(Translocation.getString("consensusTDNA", 24)).append(TDNAcons.getConsensusString()).append("\r\n");
		sb.append(Translocation.getString("most frequent TDNA", 24)).append(TDNAcons.getMostRepeatedString()).append(" "+TDNAcons.getMostRepeatedStringNr()).append(" fraction:" +TDNAcons.getMostRepeatedStringFraction());
		sb.append("\r\nend of TDNA\r\n");
		sb.append("Genomic reads that include the primer\r\n");
		for(SAMRecord sr: genomicReadsIncPrimer) {
			String cigar = Translocation.getString(sr.getCigarString(),24);
			//String map = Translocation.getString(""+sr.getFirstOfPairFlag(),4);
			//cigar = map+cigar;
			sb.append(cigar+sr.getReadString()+"\t"+sr.getFirstOfPairFlag()).append("\r\n");
		}
		int count = 0;
		for(SAMRecord sr: genomicReadsIncPrimer) {
			String cigar = Translocation.getString(sr.getCigarString(),20);
			String map = Translocation.getString(""+sr.getAlignmentStart(),10);
			cigar = map+cigar;
			String readPart = Translocation.getMatchedPart(sr, true);
			String readPartNoSpaces = Translocation.getMatchedPart(sr, false);
			//check if first part is a Hard clip, because that will not provide any info here
			if(sr.getCigar().getFirstCigarElement().getOperator()==CigarOperator.H) {
				//System.out.println("Hard "+sr.getReadName());
				for(SAMRecord tempSr: sams) {
					if(tempSr.getContig().contentEquals(sp.getChr()) && tempSr.getReadName().contentEquals(sr.getReadName())) {
						//add this one
						genomePrimerCons.add(tempSr.getReadString());
						genomeCons.add(readPartNoSpaces);
						genomeConsInt.add(tempSr.getReadLength()-sr.getReadLength());
					}
				}
				
			}
			else {
				genomePrimerCons.add(sr.getReadString());
				genomeCons.add(readPartNoSpaces);
				//if(sr.getFirstOfPairFlag()==sp.isFirstOfPairFlag()) {
				CigarElement ce = sr.getCigar().getFirstCigarElement();
				genomeConsInt.add(ce.getLength());
			}
			//}
			//sb.append(cigar+readPart).append("\r\n");
			count++;
		}
		sb.append(Translocation.getString("consGenomeIncPrimer", 24)).append(genomePrimerCons.getConsensusString()).append("\r\n");
		sb.append(Translocation.getString("freqGenomeIncPrimer", 24)).append(genomePrimerCons.getMostRepeatedString()).append(" "+genomePrimerCons.getMostRepeatedStringNr()).append("\r\n");
		sb.append(Translocation.getString("consensusGenome", 30)).append(genomeCons.getConsensusString()).append(" "+genomeCons.getMostRepeatedStringNr());
		sb.append("\n").append(Translocation.getString("consensusGenomeINT", 30)).append(genomeConsInt.getMostRepeatedInt())
		.append(" ").append(genomeConsInt.getMin()+":"+genomeConsInt.getMax());
		sb.append("\r\nend of GENOMIC reads including primer\r\n");
		sb.append("rest of GENOMIC reads\r\n");
		for(SAMRecord sr: genomicReadsAll) {
			if(sr.getContig().contentEquals(sp.getChr())) {
				junctionAll.add(sr.getReadString());
				sb.append(sr.getReadString()+"\r\n");
			}
			else {
				junctionAll.add(Translocation.getSoftClippedPart(sr));
				sb.append(sr.getReadString()+"\r\n");
			}
		}
		//get the filler sequence
		
		
		for(SAMRecord sr: genomicReads) {
			String cigar = Translocation.getString(sr.getCigarString(),20);
			String map = Translocation.getString(""+sr.getMappingQuality(),4);
			cigar = map+cigar;
			//sb.append(cigar+sr.getReadString()).append("\r\n");
		}
		sb.append("end of GENOMIC reads\r\n========\r\n");
		
		//sb.append(this.toString());
		//System.out.println("BEFORE ASSEMBLY1");
		//System.out.println(a.toString());
		//System.out.println("BEFORE ASSEMBLY2");
		/*
		if(this.getRealPosition()==4010692) {
			System.out.println("Found it "+sp.getSample());
			System.out.println(genomeCons.getConsensusString());
			System.out.println(genomeCons.getMostRepeatedString());
			System.out.println(genomePrimerCons.getConsensusString());
			System.out.println(genomePrimerCons.getMostRepeatedString());
			
			System.out.println("##");
			for(SAMRecord sr: genomicReadsIncPrimer) {
				System.out.println(sr.getReadString());
			}
			//System.out.println("====");
			for(SAMRecord sr: sams) {
				//System.out.println(sr.getReadName()+"\t"+sr.getContig()+"\t"+sr.getCigarString()+"\t"+sr.getReadString());
			}
			//System.exit(0);
			//a.assembleAll();
			//sb.append(a.toString());
		}
		*/
		//a.assembleAll();
		//if(a.getFragments().size()>1) {
			//System.out.println("ASSEMBLY>1");
		//}
		//System.out.println(a.toString());
		//System.out.println("AFTER ASSEMBLY");
		//System.out.println(this.getTranslocationSequence(true));
		
		
		return sb.toString();
	}
	private static String getString(String cigarString, int size) {
		StringBuffer s = new StringBuffer(cigarString);
		if(s.length()>size) {
			return s.toString().substring(0, size);
		}
		while(s.length()<size) {
			s.append(" ");
		}
		return s.toString();
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
}
