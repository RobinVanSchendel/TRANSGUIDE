package data;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class Translocation {
	private static final String NOGENOMIC = "NOGENOMIC";
	private static final String NOTRANSLOCATION = "NOTRANSLOCATION";
	private static final int refSize = 500;
	private static final int maxSizeTranslocationSequence = 151;
	public ArrayList<SAMRecord> sams;
	public HashMap<String, Integer> names;
	private String filler = "";
	private String hom = "";
	private String error;
	private String ref;
	private String realPositionCounter;
	private MyOptions options;

	public Translocation(SAMRecord s, MyOptions options) {
		sams = new ArrayList<SAMRecord>();
		names = new HashMap<String, Integer>();
		this.options = options;
		addSam(s);
	}
	public boolean addSam(SAMRecord s) {
		//do not add sams that have a primary and secondary alignment in contig one
		if(s.getContig().equals(options.getChr()) && getContigSATagIsContig(s,options.getChr())) {
			return false;
		}
		//orientation should be correct, otherwise don't add
		//remove this filter for now
		//if(s.getContig().equals(SAMReader.str) && s.getReadNegativeStrandFlag() == SAMReader.forwardRB) {
		//	return;
		//}
		//all others can be added
		
		sams.add(s);
		names.put(s.getReadName(),0);
		return true;
	}
	public int getNrSupportingReads() {
		return sams.size();
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
	public int getNrAnchors(boolean addPartial) {
		int count = 0;
		ArrayList<Integer> uniquePos = new ArrayList<Integer>();
		for(SAMRecord sam: sams) {
			//only get the opposite reads
			if(sam.getFirstOfPairFlag() != options.isFirstOfPairFlag()) {
				//System.out.println("Still here first "+sam.getReadName() + " "+Translocation.getNMis0(sam) +" "+ sam.getAttribute("NM"));
				//no mismatches
				if(Translocation.getNMis0(sam)) {
					//mapping quality of >50 && only one cigar part
					//System.out.println("Still here NM0 "+sam.getReadName());
					if(sam.getMappingQuality()>50) {
						//System.out.println("Still here MQ>50 "+sam.getReadName() +" "+ sam.getCigarString());
						//complete Match
						if(sam.getCigarLength()==1) {
							//cigar has to be a match
							if(sam.getCigarString().endsWith("M")) {
								if(!uniquePos.contains(sam.getAlignmentStart())) {
									//System.out.println("+1");
									count++;
									uniquePos.add(sam.getAlignmentStart());
								}
								//System.out.println(sam.getCigarString());
							}
						}
						//partial match
						else if(addPartial && sam.getCigarLength()==2) {
							int lengthMatch = 0;
							int position = 0;
							if(sam.getCigar().getCigarElement(0).getOperator() == CigarOperator.M) {
								lengthMatch = sam.getCigar().getCigarElement(0).getLength();
								//the start determines uniqueness
								position = sam.getAlignmentStart();
							}
							else if(sam.getCigar().getCigarElement(1).getOperator() == CigarOperator.M) {
								lengthMatch = sam.getCigar().getCigarElement(1).getLength();
								//the end determines uniqueness
								position = sam.getAlignmentEnd();
							}
							//System.out.println("Partial match "+lengthMatch);
							//maybe 70 is a bit arbitrary
							if(lengthMatch>=70) {
								if(!uniquePos.contains(position)) {
									count++;
									uniquePos.add(position);
								}
							}
						}
					}
				}
			}
		}
		return count;
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
		sb.append("NrSupportingReads").append(s);
		sb.append("NrAnchors").append(s);
		sb.append("NrAnchorsIncludingPartial").append(s);
		sb.append("NrSupportJunction").append(s);
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
		sb.append("Type").append(s);
		sb.append("Homology").append(s);
		sb.append("Filler").append(s);
		sb.append("RefSequence").append(s);
		sb.append("error").append(s);
		sb.append("isOK").append(s);
		sb.append("getTDNASequenceDiffReads").append(s);
		sb.append("getTDNASequenceDiffReadsHighestContributor").append(s);
		
		return sb.toString();
	}
	public String toString() {
		StringBuffer sb  = new StringBuffer();
		String s = "\t";
		sb.append(getNrSupportingReads()).append(s);
		sb.append(getNrAnchors(false)).append(s);
		sb.append(getNrAnchors(true)).append(s);
		sb.append(getNrSupportJunction()).append(s);
		sb.append(getSizePrimary()).append(s);
		sb.append(getSizeSecondary()).append(s);
		sb.append(getContigMate()).append(s);
		sb.append(getPosition()).append(s);
		sb.append(getRealPosition()).append(s);
		sb.append(realPositionCounter).append(s);
		sb.append(getDistanceToLBRB()).append(s);
		sb.append(isLBorRB()).append(s);
		sb.append(getIGVPos()).append(s);
		sb.append(isForwardText()).append(s);
		sb.append(getCigarString()).append(s);
		sb.append(getTranslocationSequence()).append(s);
		sb.append(getTDNASequence(false)).append(s);
		sb.append(getGenomicSequence()).append(s);
		sb.append(getType()).append(s);
		sb.append(getHomology()).append(s);
		sb.append(getFiller()).append(s);
		sb.append(ref).append(s);
		sb.append(error).append(s);
		sb.append(isOK()).append(s);
		sb.append(getTDNASequenceDiffReads()).append(s);
		sb.append(getTDNASequenceDiffReadsHighestContributor()).append(s);
		
		//sb.append("\n");

		for(SAMRecord sam: sams) {
			//sb.append("\t"+sam.getSAMString());
		}
		
		return sb.toString();
	}
	private int getDistanceToLBRB() {
		int disLB = this.getDistanceToLB();
		int disRB = this.getDistanceToRB();
		int disLBAbs = Math.abs(disLB);
		int disRBAbs = Math.abs(disRB);
		if(disLBAbs<disRBAbs) {
			return disLB;
		}
		return disRB;
	}
	
	private int getDistanceToRB() {
		String seq = this.getTDNASequence();
		int RBPos = options.getTDNARBPos();
		ArrayList<Integer> posA = new ArrayList<Integer>();
		for(SAMRecord sam: sams) {
			if(sam.getReadString().startsWith(seq) && sam.getContig().equals(options.getChr())) {
				int disStart = Math.abs(RBPos-sam.getAlignmentEnd());
				int disEnd = Math.abs(RBPos-sam.getAlignmentStart());
				if(disStart<disEnd) {
					posA.add(RBPos-sam.getAlignmentEnd());
				}
				else {
					posA.add(RBPos-sam.getAlignmentStart());
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
	private int getNrSupportJunction() {
		if(this.hasErrors()) {
			return -1;
		}
		String junctionMin = this.getMinimalJunction();
		int count = 0;
		for(SAMRecord s: sams) {
			String seq = s.getReadString();
			if(seq.contains(junctionMin)) {
				count++;
			}
			else if(Utils.reverseComplement(seq).contains(junctionMin)) {
				count++;
			}
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
		String junction = this.getTranslocationSequence();
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
		
		if(tdna.length()>size) {
			tdna = tdna.substring(tdna.length()-size);
		}
		if(genomic.length()>size) {
			genomic = genomic.substring(0, size);
		}
		int start = junction.indexOf(tdna);
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
		for(int i=1;i<=max;i++) {
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
		int LBPos = options.getTDNALBPos();
		ArrayList<Integer> posA = new ArrayList<Integer>();
		for(SAMRecord sam: sams) {
			if(sam.getReadString().startsWith(seq) && sam.getContig().equals(options.getChr())) {
				int disStart = Math.abs(LBPos-sam.getAlignmentEnd());
				int disEnd = Math.abs(LBPos-sam.getAlignmentStart());
				if(disStart<disEnd) {
					//System.out.println(LBPos-sam.getAlignmentEnd());
					posA.add(sam.getAlignmentEnd()-LBPos);
				}
				else {
					//System.out.println(LBPos-sam.getAlignmentStart());
					posA.add(sam.getAlignmentStart()-LBPos);
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
			if(!srec.isSecondaryAlignment() && srec.getFirstOfPairFlag() == options.isFirstOfPairFlag()
					&& cigarStringFollowsMSH(srec.getCigarString()) 
					&& !srec.getContig().equals(options.getChr())
					&& getNMis0(srec)) {
				if(debug) {
					System.out.println(srec.getCigarString()+" "+srec.getDuplicateReadFlag() + " "+srec.getAlignmentStart()+"-"+srec.getAlignmentEnd() + " "+ srec.getMappingQuality()+ " "+srec.getReadName());
				}
				//add the correct part of the seq
				CigarElement ce = srec.getCigar().getCigarElement(0);
				if(ce.getOperator() == CigarOperator.S) {
					if(debug) {
						System.out.println("SEQ1\t"+srec.getReadString().substring(ce.getLength()));
					}
					seqs.add(srec.getReadString().substring(ce.getLength()));
				}
				//M, so the S part follows
				else {
					//reverse complement
					String seqPart = srec.getReadString().substring(0, ce.getLength());
					if(debug) {
						System.out.println("seqPart\t"+seqPart);
						System.out.println("SEQ2\t"+Utils.reverseComplement(seqPart));
					}
					seqs.add(Utils.reverseComplement(seqPart));
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
				System.out.println(srec.getFirstOfPairFlag() == options.isFirstOfPairFlag());
				System.out.println(!srec.getContig().equals(options.getChr()));
				System.out.println(getNMis0(srec));
				System.out.println(srec.getReadString());
				*/
			}
		}
		if(debug) {
			System.out.println("size genomic "+seqs.size());
			for(String s: seqs) {
				System.out.println("\t\t"+s);
			}
		}
		//sometimes this contains one read or so, so don't rely on this then
		//5 is an arbitrary number
		if(seqs.size()>=5) {
			String seq = consensusString(seqs);
			return seq;
		}
		
		//some location have only the second read mapping to the genome
		int maxGenomeSizePart = 50;
		for(SAMRecord srec: sams) {
			if(!srec.isSecondaryAlignment() && srec.getFirstOfPairFlag() != options.isFirstOfPairFlag()
					&& cigarStringFollowsMSH(srec.getCigarString()) 
					&& !srec.getContig().equals(options.getChr())) {
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
			String seq = consensusString(seqs);
			return seq;
		}
		//give up
		return NOGENOMIC;
	}
	private String getTDNASequenceDiffReads(){
		HashMap<String, Integer> seqsDiff = new HashMap<String, Integer>();
		//simplest case if second in pair and primary alignment
		
		for(SAMRecord sam: sams) {
			if(!sam.isSecondaryAlignment() && sam.getContig().equals(options.getChr())) {
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
			if(sam.isSecondaryAlignment() && sam.getContig().equals(options.getChr())) {
				String seq = sam.getReadString();
				if(!seqsDiff.containsKey(seq)) {
					seqsDiff.put(seq, 0);
				}
				seqsDiff.put(seq,seqsDiff.get(seq)+1);
			}
		}
		double total = 0;
		for(String key: seqsDiff.keySet()) {
			total+= seqsDiff.get(key);
		}
		return seqsDiff.size()+"|"+total;
	}
	private double getTDNASequenceDiffReadsHighestContributor(){
		HashMap<String, Integer> seqsDiff = new HashMap<String, Integer>();
		//simplest case if second in pair and primary alignment
		
		for(SAMRecord sam: sams) {
			if(!sam.isSecondaryAlignment() && sam.getContig().equals(options.getChr())) {
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
			if(sam.isSecondaryAlignment() && sam.getContig().equals(options.getChr())) {
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
		
		for(SAMRecord sam: sams) {
			if(!sam.isSecondaryAlignment() && sam.getContig().equals(options.getChr())) {
				//up unto the M
				CigarElement ce = sam.getCigar().getCigarElement(0);
				//always correct?
				if(ce.getOperator() != CigarOperator.M && sam.getCigar().numCigarElements() == 2) {
					//System.err.println("bug still here");
					//System.out.println(sam.getCigarString()+" "+sam.getAlignmentStart()+"-"+sam.getAlignmentEnd());
					//System.out.println(sam.getReadString());
					//System.out.println(sam.getFirstOfPairFlag());
					//System.exit(0);
				}
				//only do it if this is not the case
				else {
					int pos = ce.getLength();
					String seq = sam.getReadString().substring(0, pos); 
					seqs.add(seq);
					if(!seqsDiff.containsKey(seq)) {
						seqsDiff.put(seq, 0);
					}
					seqsDiff.put(seq,seqsDiff.get(seq)+1);
				}
			}
			//most go here!
			if(sam.isSecondaryAlignment() && sam.getContig().equals(options.getChr())) {
				String seq = sam.getReadString();
				seqs.add(seq);
				if(!seqsDiff.containsKey(seq)) {
					seqsDiff.put(seq, 0);
				}
				seqsDiff.put(seq,seqsDiff.get(seq)+1);
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
			return consensusString(seqs);
		}
		return "unknown";
	}

	public String getTranslocationSequence() {
		//String consensus = getCigarString();
		
		ArrayList<String> seqs = new ArrayList<>();
		for(SAMRecord s:sams) {
			if(!s.isSecondaryAlignment() && cigarStringFollowsMSH(s.getCigarString())) {
				if(s.getReadString().length()>maxSizeTranslocationSequence) {
					seqs.add(s.getReadString().substring(0, maxSizeTranslocationSequence));
				}
				else {
					seqs.add(s.getReadString());
				}
			}
		}
		if(seqs.size()>0) {
			String seq = consensusString(seqs);
			//check if we need the reverse complement
			for(SAMRecord s: sams) {
				if(s.getReadString().startsWith(seq) && s.getFirstOfPairFlag() == options.isFirstOfPairFlag()) {
					//only take the revComplement if s is on the reverse strand
					if(!s.getContig().equals(options.getChr()) && s.getReadNegativeStrandFlag()) {
						SAMRecord srev = s.deepCopy();
						srev.reverseComplement();
						return srev.getReadString();
					}
					else {
						return seq; 
					}
					
				}
			}
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
		return consensusString(cigars);
	}
	//TODO this consensus string likely only works if we align
	//TODO the sequences correctly
	//20190506 my solution is to put only substrings of the seq in the list
	//this way I can avoid the problem of having tiling seqs
	private static String consensusString(ArrayList<String> list) {
		String mostRepeatedWord 
	    = list.stream()
	          .collect(Collectors.groupingBy(w -> w, Collectors.counting()))
	          .entrySet()
	          .stream()
	          .max(Comparator.comparing(Entry::getValue))
	          .get()
	          .getKey();
		return mostRepeatedWord;
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
				if(!s.getMateReferenceName().equals(options.getChr())) {
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
					&& !s.getContig().equals(options.getChr())) {
				//no mismatches
				if(getNMis0(s)) {
					//maybe get the actual position
					if(s.getCigar().getCigarElement(0).getOperator() == CigarOperator.S) {
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
		return Translocation.getPosition(sams.get(0), options);
		//found a bug here!
		//return sams.get(0).getMateAlignmentStart();
	}
	public static int getPosition(SAMRecord s, MyOptions options) {
		//info is in the SA tag
		if(s.isSecondaryAlignment() && !getContigSATagIsContig(s, options.getChr())) {
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
			return tag == 0;
		}
		return true;
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
				this.ref = rsf.getSubsequenceAt(this.getContigMate(), start, end).getBaseString();
			}
			if(!this.ref.startsWith(this.getGenomicSequence())) {
				this.addError("Ref sequence deviates from sequence in reads");
			}
		}
	}
	public void printDebug() {
		String tl = this.toString();
		for(SAMRecord sr: sams) {
			tl+="\n\t";
			tl+=sr.getReadString()+" "+sr.getAlignmentStart();// +" "+sr.getAttributes().toString();
		}
		System.out.println(tl);
		System.out.println("TDNA SEQ:");
		System.out.println(this.getTDNASequence(true));
		System.out.println("GENOMIC SEQ:");
		System.out.println(this.getGenomicSequence(true));
		System.out.println(getMinimalJunction(true));
		System.out.println(getRealPosition(true));
		
	}
}
