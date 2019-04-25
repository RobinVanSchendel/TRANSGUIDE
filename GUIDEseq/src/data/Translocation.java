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

public class Translocation {
	public ArrayList<SAMRecord> sams;

	public Translocation(SAMRecord s) {
		sams = new ArrayList<SAMRecord>();
		addSam(s);
	}

	public void addSam(SAMRecord s) {
		//do not add sams that have a secondary alignment in contig one
		if(!getContigSATagIsContig(s,SAMReader.str)) {
			sams.add(s);
		}
	}
	public int getSize() {
		return sams.size();
	}
	/**An anchor is defined as the read which does not contain the initial integration
	 * read. In our T-DNA experiments P7 is the T-DNA read.
	 * 
	 * Count #reads that perfectly map to the genome.
	 * NM = 0
	 * 151M (dependendent on read length) 
	 * 
	 * @return
	 */
	public int getNrAnchors() {
		int count = 0;
		for(SAMRecord sam: sams) {
			//only get the opposite reads
			if(sam.getFirstOfPairFlag() != SAMReader.getFirstOfPairFlag) {
				//no mismatches
				if(Translocation.getNMis0(sam)) {
					//mapping quality of >50 && only one cigar part
					if(sam.getMappingQuality()>50 && sam.getCigarLength()==1) {
						//cigar has to be a match
						if(sam.getCigarString().endsWith("M")) {
							count++;
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
	public String toString() {
		StringBuffer sb  = new StringBuffer();
		String s = "\t";
		sb.append(getSize()).append(s);
		sb.append(getNrAnchors()).append(s);
		sb.append(getSizePrimary()).append(s);
		sb.append(getSizeSecondary()).append(s);
		sb.append(getContigMate()).append(s);
		sb.append(getPosition()).append(s);
		sb.append(getRealPosition()).append(s);
		sb.append(getIGVPos()).append(s);
		sb.append(isForward()).append(s);
		sb.append(getCigarString()).append(s);
		sb.append(getTranslocationSequence()).append(s);
		sb.append(getFirstContigSequence()).append(s);
		sb.append(getSecondContigSequence()).append(s);
		sb.append(sams.get(0).getReadName()).append(s);
		//sb.append("\n");

		for(SAMRecord sam: sams) {
			//sb.append("\t"+sam.getSAMString());
		}
		
		return sb.toString();
	}
	private String getSecondContigSequence() {
		//get the second in pair and derive the ones with a cigar sequence
		//has to be the primary alignment
		return "still implement!";
		//return null;
	}

	private String getFirstContigSequence() {
		ArrayList<String> seqs = new ArrayList<String>();
		//simplest case if second in pair and primary alignment
		for(SAMRecord sam: sams) {
			if(!sam.isSecondaryAlignment() && sam.getContig().equals(SAMReader.str)) {
				//up unto the M
				String cigar = sam.getCigarString();
				String position = cigar.substring(0, cigar.indexOf("M"));
				int pos = Integer.parseInt(position);
				seqs.add(sam.getReadString().substring(0, pos));
			}
			if(sam.isSecondaryAlignment() && sam.getContig().equals(SAMReader.str)) {
				seqs.add(sam.getReadString());
			}
		}
		if(seqs.size()>0) {
			return consensusString(seqs);
		}
		return "unknown";
		/*
		for(SAMRecord sam: sams) {
			if(sam.isSecondaryAlignment()) {
				seqs.add(sam.getReadString());
			}
		}
		if(seqs.size()>0) {
			return consensusString(seqs);
		}
		//still need to implement
		String seq = getConsensusSeqPrimary();
		for(SAMRecord sam: sams) {
			if(sam.getReadString().equals(seq) && Translocation.getNMis0(sam)){
				String cigar = sam.getCigarString();
				int indexS = cigar.indexOf("S");
				int indexM = cigar.indexOf("M");
				if(indexS<indexM) {
					String pos = cigar.substring(0, indexS);
					int position = Integer.parseInt(pos);
					return seq.substring(0,position);
				}
				else {
					String pos = cigar.substring(indexM+1, indexS);
					int position = Integer.parseInt(pos);
					return seq.substring(0,position);
				}
			}
		}
		return null;
		*/
	}

	private String getTranslocationSequence() {
		//String consensus = getCigarString();
		ArrayList<String> seqs = new ArrayList<>();
		for(SAMRecord s:sams) {
			if(!s.isSecondaryAlignment() && cigarStringFollowsMSH(s.getCigarString())) {
				seqs.add(s.getReadString());
			}
		}
		if(seqs.size()>0) {
			String seq = consensusString(seqs);
			//check if we need the reverse complement
			for(SAMRecord s: sams) {
				if(s.getReadString().equals(seq) && s.getFirstOfPairFlag() == SAMReader.getFirstOfPairFlag) {
					//only take the revComplement if s is on the reverse strand
					if(!s.getContig().equals(SAMReader.str) && s.getReadNegativeStrandFlag()) {
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
		return "UNKNOWN";
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

	public String getIGVPos() {
		return this.getContigMate()+":"+this.getRealPosition();
	}

	public String getContigMate() {
		if(getSize()>0) {
			for(SAMRecord s: sams) {
				if(!s.getMateReferenceName().equals(SAMReader.str)) {
					return s.getMateReferenceName();
				}
			}
		}
		return null;
	}

	public boolean isForward() {
		if(getSize()>0) {
			return sams.get(0).getMateNegativeStrandFlag()==false;
		}
		return false;
	}
	public int getRealPosition() {
		ArrayList<String> seqs = new ArrayList<>();
		for(SAMRecord s:sams) {
			if(!s.isSecondaryAlignment() && cigarStringFollowsMSH(s.getCigarString()) 
					&& !s.getContig().equals(SAMReader.str)) {
				//no mismatches
				if(getNMis0(s)) {
					seqs.add(s.getReadString());
				}
			}
		}
		if(seqs.size()>0) {
			//because the read that is not the one connected to the T-DNA
			//can also determine the consensus this might be a bit tricky
			String seq = consensusString(seqs);
			for(SAMRecord s: sams) {
				if(s.getReadString().equals(seq)) {
					//only take the revComplement if s is on the reverse strand
					//let the orientation depend on the position of the clip
					if(s.getCigar().getCigarElement(0).getOperator() == CigarOperator.S) {
						return s.getAlignmentStart();
					}
					else {
						return s.getAlignmentEnd();
					}
				}
			}
		}
		return sams.get(0).getMateAlignmentStart();
		
	}
	public int getPosition() {
		//String consensus = getCigarString();
		return sams.get(0).getMateAlignmentStart();
	}
	public static int getPosition(SAMRecord s) {
		//info is in the SA tag
		if(s.isSecondaryAlignment() && !getContigSATagIsContig(s, SAMReader.str)) {
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
		for(SAMRecord s: sams) {
			if(s.getReadName().equals(readName)) {
				return true;
			}
		}
		return false;
	}
	public static boolean getNMis0(SAMRecord srec) {
		if(srec.hasAttribute("NM")) {
			int tag = (Integer) srec.getAttribute("NM");
			return tag == 0;
		}
		return true;
	}
}
