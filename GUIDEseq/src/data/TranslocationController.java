package data;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

public class TranslocationController {
	private HashMap<String, ArrayList<Translocation>> trans = new HashMap<String, ArrayList<Translocation>>();
	public static final int MAXDIST = 5;
	public static final int MAXANCHORDIST = 2000;
	private int primerStart, primerEnd;
	public static final String LB = "GTTTACACCACAATATATCCTGCCA";
	public static final String RB = "GTTTACCCGCCAATATATCCTGTCA";
	private SamplePrimer sp;
	private HashMap<String, Translocation> searchRealPositions = new HashMap<String, Translocation>();
	boolean debug = true;
	public static final String testName = "A00379:349:HM7WFDSXY:4:1274:28962:25254";
	
	public TranslocationController(SamplePrimer sp) {
		this.sp = sp;
	}
	/**
	 * If the position is not faulty,
	 * @param s
	 * @return
	 */
	private Translocation getNearestTranslocation(SAMRecord s) {
		if (Translocation.getPosition2(s, sp)!=-1) {  //filter for faulty positions
			int pos = Translocation.getPosition2(s, sp);
			Translocation nearest = null;
			int minDis = Integer.MAX_VALUE;
			if(trans.get(s.getMateReferenceName()) == null) {
				return null;
			}
			for(Translocation tl: trans.get(s.getMateReferenceName())) {
				//if(tl.isForward() != s.getMateNegativeStrandFlag()) {
					int tempPos = tl.getPosition1();
					int distance = Math.abs(pos-tempPos);
					if(distance<minDis) {
						minDis = distance;
						nearest = tl; 
					}
					else if(minDis>=0 && minDis<MAXDIST && nearest!=null) {
						return nearest;
					}
					//}
			}
			if (minDis<MAXDIST) {
				return nearest;
			}
			return null;
		}
		else {
			return null;
		}
	}
	public Translocation addTranslocation(SAMRecord s, int maxReads) {
		Translocation nearest = getNearestTranslocation(s);
//		int pos = Translocation.getPosition2(s, options);
		if(nearest !=null) {
			if(nearest.getSams().size()<maxReads) {
				nearest.addSam(s);
			}
			return nearest;
		}
		else {
			Translocation tl = new Translocation(s, sp);
			if(s.getReadName().contentEquals(testName)) { 
			System.out.println("checkpoint1");}
			//sometimes the sam is not added due to filtering of secondary alignments
			if(tl.getNrSupportingReads()>0) {
				ArrayList<Translocation> al = trans.get(s.getMateReferenceName());
				if(al==null) {
					al = new ArrayList<Translocation>();
					trans.put(s.getMateReferenceName(), al);
				}
				al.add(tl);
				if(s.getReadName().contentEquals(testName)) { 
					System.out.println("checkpoint2");}
				return tl;
			}
			else {
				//System.out.println(tl.getContigMate()+" has 0 supporting reads");
			}
		}
		return null;
	}
    

	public void printContents(BufferedWriter bw, long minSupport) {
		//sorting doesn't work at the moment
		//trans.sort(comparing(Translocation::getNrSupportingReads, Collections.reverseOrder()));
		//System.out.println("File\t"+Translocation.getHeader());
		try {
			int counter = 0;
			for(String key: trans.keySet()){
				int printNr = 0;
				for(Translocation tl: trans.get(key)) {
					if(tl.getSams().size()>0) {
						if(tl.getNrAnchors(false)>=minSupport) {
							//String output = null;
							String output = tl.toString(debug);
							bw.write(sp.getSampleString()+"\t"+output+"\n");
							counter++;
							printNr++;
						}
					}
				}
				System.out.println("Written "+printNr+" : "+trans.get(key).size()+" events with >= "+minSupport+" support "+sp.getSample());
			}
			System.out.println("Written "+counter+" events with >= "+minSupport+" support");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
	}
	/**
	 * This adds the anchors. They can be at most 2000 bp away from the T-DNA primers. Usually they should not come close to that amount.
	 * Maybe have to change the way to get the position in some cases. 
	 * @param sr
	 */
	public void addMates(SamReader sr) {
		int around = MAXANCHORDIST;
		int duplicates = 0;
		int count = 0;
		for(String key: trans.keySet()){
			for(Translocation tl: trans.get(key)) {
				count++;
				SAMRecordIterator sri = sr.query(tl.getContigMate(), tl.getPosition1()-around, tl.getPosition1()+around, false);
				while(sri.hasNext()) {
					SAMRecord srec = sri.next();
					if(srec.getReadName().contentEquals(testName)) {
					}
					if(srec.getContig()!=null) {
						if(srec.getReadName().contentEquals(testName)) {
						}
						//below duplicates are not included, only anchors (first of pair), and the read has to align primarily to the chromosome, which will remove some anchors from very small fragments, or fillers with very large perfect alignments
						if((isDuplicate(srec)==false) && (tl.containsRecord(srec.getReadName())==true) && (srec.getFirstOfPairFlag()==true) && (srec.getContig().equals(sp.getChr())==false)) { 
							if(srec.getReadName().contentEquals(testName)) {
							}
							boolean added = tl.addSam(srec);
							if(srec.getReadName().contentEquals(testName)) {
								System.err.println("added??? "+added);
							}
						}
						//sometimes hard clipped reads went through our duplicate filter
						if(isDuplicate(srec)) {
							if(srec.getReadName().contentEquals(testName)) {
								System.err.println("removing "+srec.getReadName());
							}
							tl.removeSam(srec);
						}
						if(isDuplicate(srec)) {
							duplicates++;
						}
					}
				}
				sri.close();
				if(count%1000==0) {
					System.out.println("Already processed "+count+" mates "+trans.get(key).size() + " from chr "+key);
				}
				//System.out.println(tl.getReads());
			}
			
			
		}
		System.out.println("Found "+duplicates+" duplicate mates (should be >0)");
	}

	private boolean isDuplicate(SAMRecord srec) {
		//if it has an UMI, that is already taken care of
		//disabled this as it does not seem to work
		//if(sp.hasUMI()) {
		//	return false;
		//}
		return srec.getDuplicateReadFlag();
		
	}
	/**
	 * 
	 * @param rsf
	 */
	public void addRefGenomePart(ReferenceSequenceFile rsf) {
		for(String key: trans.keySet()){
			for(Translocation tl: trans.get(key)) {
				if(tl.isOK()) {
					tl.addRefSequence(rsf);
				}
			}
		}
	}

	public void launchAnalysis(int maxTrans, int maxReadsPerTrans, int maxReads) {
		if(maxTrans>=0) {
			System.err.println("Warning: program will stop after finding "+maxTrans+" translocations");		}
		if(maxReads>=0) {
			System.err.println("Warning: program will stop after "+maxReads+" reads");		
		}
		SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(sp.getFile());
		SamReader sr2 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(sp.getFile());
		if(!sr.hasIndex()) {
			System.out.println("no index available for this bam file. Please run samtools index "+sp.getFile());
			System.exit(0);
		}
	    ReferenceSequenceFile rsf = ReferenceSequenceFileFactory.getReferenceSequenceFile(sp.getRef());
	    //SAMRecordIterator r = sr.iterator();
	    String chr = sp.getChr();
	    ReferenceSequence rs = rsf.getSequence(chr);
	    if(rs == null) {
	    	System.err.println("Chromosome "+chr+" is not part of the reference genome");
	    	System.exit(0);
	    }
	    else {
	    	//System.out.println("Chromosome "+chr+" was found");
	    }
	    String seq = rs.getBaseString().toUpperCase();
	    int index = seq.indexOf(sp.getPrimer());
	    int start = -1;
	    int end = -1;
	    boolean positiveStrand = true;
	    if(index == -1) {
	    	index = seq.indexOf(Utils.reverseComplement(sp.getPrimer()));
	    	positiveStrand = false;
	    }
	    else {
	    	//System.out.println("Primer found on + strand at position: "+index+" ("+options.getPrimer()+")");
	    }
	   
	    if(index == -1) {
	    	System.err.println("Primer could not be found ["+sp.getPrimer()+"]");
	    	System.exit(0);
	    }
		start = index;
    	end = index+sp.getPrimer().length();
    	primerStart = start+1;
    	primerEnd = end;
    	System.out.println("Primerstart:" +primerStart +", Primerend:" +primerEnd +", positivestrand:"+positiveStrand +", primerseq:"+sp.getPrimer());
	    SAMRecordIterator r = sr.iterator () ;
	    
	    int count = 0;
	    int NM0Count = 0;
	    int duplicateFlag = 0;
	    int Ncount = 0;
	 
        while(r.hasNext()) {
        	count++;
        	SAMRecord srec = r.next();
        	if(debug) {
        		//System.err.println("Still here");
        		if(srec.getReadName().contentEquals(testName)) {
        			System.out.println("BEFORE RCing");
        			System.out.println("read name: " +srec.getReadName());
        			System.out.println("CIGAR: " +srec.getCigarString());
        			System.out.println("contig: " +srec.getContig());
        			System.out.println("add info: " +srec.toString());
        			System.out.println("sequence: " +srec.getReadString());
        			System.out.println("alignment start: " +srec.getAlignmentStart());
        			System.out.println("alignment end: " +srec.getAlignmentEnd());
        			System.out.println("Read negative strand flag: " +srec.getReadNegativeStrandFlag());

        			if (srec.getAttribute("SA") != null) {
        				System.out.println("pos SATAg: " +getPosSATag(srec));
        				System.out.println("sign SATag: " +getForwardSATag(srec));
    					System.out.println("pos SATag end: " +getPosSATagEnd(srec));
        			}
        			SAMRecord temp = sr2.queryMate(srec);
        			System.out.println("mate name: " +temp.getReadName());
        			System.out.println("CIGAR: " +temp.getCigarString());
        			System.out.println("contig: " +temp.getContig());
        			System.out.println("add info: " +temp.toString());
        			System.out.println("sequence: " +temp.getReadString());
        			System.out.println("alignment start: " +temp.getAlignmentStart());
        			System.out.println("alignment end: " +temp.getAlignmentEnd());
        			System.out.println("Read negative strand flag: " +temp.getReadNegativeStrandFlag());
        			if (temp.getAttribute("SA") != null) {
        				System.out.println("pos SATAg: " +getPosSATag(temp));
        				System.out.println("sign SATag: " +getForwardSATag(temp));
        				System.out.println("pos SATag end: " +getPosSATagEnd(temp));
        			
        			}
        		}
        	}
        	// only take 0 mismatches reads
        	// don't take secondary alignment reads
        	// don't take reads with no secondary alignment, all T-DNA reads are expected to have at least 2 alignments
        	// don't take any reads that are duplicates
        	// only take the T-DNA reads.
        	if (srec.getAttribute("SA") != null){
        		String SATag = (String) srec.getAttribute("SA");
        		String[] SAList = SATag.split(",|;");
        		int SALength = SAList.length;
        		String contigString = SATag.split(",")[0];
        		if((Translocation.getNMis0(srec)) && (srec.isSecondaryAlignment()==false) && (isDuplicate(srec)==false) && (srec.getFirstOfPairFlag()== sp.isFirstOfPairFlag())){
        				if(debug) {
            			//System.err.println("Still here");
            			if(srec.getReadName().contentEquals(testName)) {
                			System.err.println("Still here");
            			}
        			}
        			NM0Count++;
        			//chr has to be filled and unequal
        			//take the reverse complement if we are looking at reads going reverse
        			if(srec.getReadNegativeStrandFlag()==true) {
        				if(debug) {
                			//System.err.println("Still here !positiveStrand");
                			//System.err.println(srec.getReadString());
                		}
        				String cigar = srec.getCigarString();
        				//boolean currentNegativeStrandFlag = srec.getReadNegativeStrandFlag();
        				//System.out.println("BEFORE RC "+srec.getReadNegativeStrandFlag());
        				
        				srec.reverseComplement();
        				//the problem with reversing the strandflag, is that some things still depend on the strandedness even after the seq has been reversed.
        				//if(srec.getReadNegativeStrandFlag()==currentNegativeStrandFlag) {
        				//	srec.setReadNegativeStrandFlag(!currentNegativeStrandFlag);
        				//}
        				//System.out.println("AFTER RC "+srec.getReadNegativeStrandFlag());
        				//System.exit(0);
        				//bug so reverse it myself
        				Cigar tempCigar = srec.getCigar();
        				if(tempCigar.numCigarElements()>1 && cigar.equals(srec.getCigarString())) {
        					Cigar rev = new Cigar();
        					//reverse the cigar
        					for(int i = tempCigar.numCigarElements()-1;i>=0;i--) {
        						rev.add(tempCigar.getCigarElement(i));
        					}
        					srec.setCigar(rev);
        				}
        				//System.out.println("Made RC "+srec.getCigarString());
        			
        				if(srec.getReadName().contentEquals(testName)) {
        					System.out.println("AFTER RCing");
        					System.out.println("read name: " +srec.getReadName());
            				System.out.println("CIGAR: " +srec.getCigarString());
            				System.out.println("CIGAR length: " +srec.getCigarLength());
            				System.out.println("contig: " +srec.getContig());
            				System.out.println("contig SA tag 1: " +getContigSATag(srec));
            				System.out.println("add info: " +srec.toString());
            				System.out.println("sequence: " +srec.getReadString());
            				System.out.println("alignment start: " +srec.getAlignmentStart());
            				System.out.println("alignment end: " +srec.getAlignmentEnd());
            				System.out.println("Read negative strand flag: " +srec.getReadNegativeStrandFlag());
        					System.out.println("SATAg length: " +SALength);
        					System.out.println("pos SATAg: " +getPosSATag(srec));
        					System.out.println("sign SATag: " +getForwardSATag(srec));
        					System.out.println("pos SATag end: " +getPosSATagEnd(srec));
        					if (SALength > 6) {
        						System.out.println("pos SATAg 2: " +getPosSecondSATag(srec));
            					System.out.println("sign SATag 2: " +getForwardSecondSATag(srec));
            					System.out.println("contig SA tag 2: " +getContigSecondSATag(srec));
        					}
        					
        				}
        			}
        		
	        		if(debug) {
	        			if(srec.getReadName().contentEquals(testName)) {
	                		System.err.println("Still here2");
	            		}
	        		}
	        		if (SALength == 6) {
	        			if (!((srec.getContig().equals(sp.getChr())==true) && (sp.getChr().equals(getContigSATag(srec))==true))){
	        				if(			((srec.getReadNegativeStrandFlag()==false) && (srec.getContig().equals(sp.getChr())==true)  && (srec.getAlignmentStart()==primerStart)) //e.g. A00379:349:HM7WFDSXY:4:1218:1796:29528-3 & A00379:349:HM7WFDSXY:4:2376:6524:12477-4
	        						 || ((srec.getReadNegativeStrandFlag()==true)  && (srec.getContig().equals(sp.getChr())==true)	&& (srec.getAlignmentEnd()==primerEnd)) //A00379:349:HM7WFDSXY:4:1122:10800:11553-4 
	        						 
	        						 || ((srec.getReadNegativeStrandFlag()==true)  && (srec.getContig().equals(sp.getChr())==false) && (sp.getChr().equals(getContigSATag(srec))) && (getForwardSATag(srec)==false) && (getEndPosSATag(srec)==primerEnd)) //A00379:349:HM7WFDSXY:4:2232:31837:21746-2 & A00379:349:HM7WFDSXY:4:1157:11659:10755-2
	        						 || ((srec.getReadNegativeStrandFlag()==true)  && (srec.getContig().equals(sp.getChr())==false) && (sp.getChr().equals(getContigSATag(srec))) && (getForwardSATag(srec)==true)  && (getPosSATag(srec)==primerStart)) //A00379:349:HM7WFDSXY:4:1258:25102:13244-2 & A00379:349:HM7WFDSXY:4:1347:6406:9392-2
	        						 
	        						 || ((srec.getReadNegativeStrandFlag()==false) && (srec.getContig().equals(sp.getChr())==false) && (sp.getChr().equals(getContigSATag(srec))) && (getForwardSATag(srec)==true)  && (getPosSATag(srec)==primerStart)) //A00379:349:HM7WFDSXY:4:1507:17879:27038-1
	        						 || ((srec.getReadNegativeStrandFlag()==false) && (srec.getContig().equals(sp.getChr())==false) && (sp.getChr().equals(getContigSATag(srec))) && (getForwardSATag(srec)==false)  && (getPosSATagEnd(srec)==primerEnd))) //A00379:349:HM7WFDSXY:4:1250:15031:28902-1 (BL25_LZB1_RB_2_7137561)
	        					
	        					//maybe above there should be a filter to remove reads with cigarlength of 3, like in this case: A00379:349:HM7WFDSXY:4:1303:20781:4382-3 (BL30_LZB1_RB_3_22700768)
	        					//although checking whether it starts with the primer appears enough.
	        				{
	        					if(debug) {
	        						if(srec.getReadName().contentEquals(testName)) {
	        							System.err.println("Still here3");
	        						}
	        					}
	        					if(srec.getReadString().contains("N")) {
	        						Ncount++;
	        					}
	        					addTranslocation(srec, maxReadsPerTrans);		       					       			
	        				}
	        			}
	        		}//negative strand, alignment end is the end of the primer.
	        		if (SALength > 6) {
	        			if (!((srec.getContig().equals(sp.getChr())==true) && (sp.getChr().equals(getContigSATag(srec))==true) && (sp.getChr().equals(getContigSecondSATag(srec))==true))){
	        				if (	   ((srec.getReadNegativeStrandFlag()==false) && (srec.getContig().equals(sp.getChr())==true) && (srec.getAlignmentStart()==primerStart)) //e.g. A00379:349:HM7WFDSXY:4:1610:23384:24330-4
	        						|| ((srec.getReadNegativeStrandFlag()==true)  && (srec.getContig().equals(sp.getChr())==true) && (srec.getAlignmentEnd()==primerEnd)) // A00379:349:HM7WFDSXY:4:1226:10945:24862-3 & A00379:349:HM7WFDSXY:4:1223:26955:6433-3 & A00379:349:HM7WFDSXY:4:1309:17815:31908-3
	        						|| ((srec.getReadNegativeStrandFlag()==true)  && (srec.getContig().equals(sp.getChr())==true) && (srec.getCigarLength()==3) && (sp.getChr().equals(getContigSecondSATag(srec))==true) && (getForwardSecondSATag(srec)==true) && (getPosSATag2(srec)==primerStart)) //A00379:349:HM7WFDSXY:4:2642:1108:31344-3
	        						){
	        					if(debug) {
	        						if(srec.getReadName().contentEquals(testName)) {
	        							System.err.println("Still here4");
	        						}
	        					}
	        					if(srec.getReadString().contains("N")) {
	        						Ncount++;
	        					}
	        					addTranslocation(srec, maxReadsPerTrans);
	        				}
	        			}
	        		}
        		}
        		if(isDuplicate(srec)) {
       				duplicateFlag++;
       			}
        	}
        	if(count%1000000==0) {
        		System.out.println("Already processed "+count+" reads, NM0 reads: "+NM0Count+", Ncount: "+Ncount);
        		//break;
        	}
        	if(getTotalTranslocations()==maxTrans) {
        		System.err.println("Stopping because we have reached maxNr");
        		break;
        	}
        	if(count==maxReads) {
        		System.err.println("Reached "+maxReads+", stopping");
        		break;
        	}
        }
        System.out.println("Strand is forward: "+positiveStrand);
        System.out.println("Found "+count+" reads");
        System.out.println("Found "+NM0Count+" NM0Count");
        System.out.println("Found "+duplicateFlag+" duplicateFlag (should be >0)");
        r.close();
        System.out.println("Adding mates");
        addMates(sr);
        System.out.println("Added mates");
        System.out.println("Adding refGenomeParts");
        addRefGenomePart(rsf);
        System.out.println("Added refGenomeParts");
        
        System.out.println("hashAllRealPositions");
        this.hashAllRealPositions();
        System.out.println("hashAllRealPositions done");
        
        //}
        try {
			sr.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        if(debug) {
        	/*
        	String debugChr = "4";
        	int startDebug = 411102-5000;
    		int endDebug = startDebug+10000;
    		boolean forward = true;
    		Translocation tl = searchTranslocation(debugChr,startDebug, endDebug, forward);
    		if(tl!=null) {
    			System.out.println("DEBUG FROM HIER!");
    			System.out.println(tl.getReads());
    			//System.exit(0);
    			System.out.println(tl);
    			tl.printDebug();
    			System.out.println(tl.isForward());
    			System.exit(0);
    		}
    		else {
    			System.out.println("not found");
    		}
    		*/
		}
	}
	private int getTotalTranslocations() {
		int	count = 0;
		for(String key: trans.keySet()) {
			count+= trans.get(key).size();
		}
		return count;
	}
	/**
	 * @param srec
	 * @return the start of the first secondary alignment
	 */
	private static int getPosSATag(SAMRecord srec) {
		String SATag = (String) srec.getAttribute("SA");
		String intString = SATag.split(",|;")[1];
		return Integer.parseInt(intString);
	}
	private static int getPosSATag2(SAMRecord srec) {
		String SATag = (String) srec.getAttribute("SA");
		String intString = SATag.split(",|;")[7];
		return Integer.parseInt(intString);
	}
	private static int getEndPosSATag(SAMRecord srec) {
		String SATag = (String) srec.getAttribute("SA");
		int tempPos = Integer.parseInt(SATag.split(",|;")[1]);
		String SACigar = SATag.split(",|;")[3];
		int indexFirstM = SACigar.indexOf("M");
		int indexFirstS = SACigar.indexOf("S");
		String sASign = SATag.split(",|;")[2];
		int pos = -1;
		if ((indexFirstM > indexFirstS) && (sASign.equals("-")==true)) {
			int saLength = Integer.parseInt(SACigar.substring(indexFirstS+1, indexFirstM));
			pos = tempPos+saLength-1;
		}
		//add something in case the order of M and S is different or when the strand is +
		return pos;
	}
	private static String getContigSATag(SAMRecord srec) {
		String SATag = (String) srec.getAttribute("SA");
		return SATag.split(",|;")[0];
	}
	private static String getContigSecondSATag(SAMRecord srec) {
		String SATag = (String) srec.getAttribute("SA");
		return SATag.split(",|;")[6];
	}
	private static int getPosSecondSATag(SAMRecord srec) {
		String SATag = (String) srec.getAttribute("SA");
		String intString = SATag.split(",|;")[7];
		return Integer.parseInt(intString);
	}
	/**
	 * @param srec
	 * @return the end of the first secondary alignment.
	 */
	private static int getPosSATagEnd(SAMRecord srec) {
		String SATag = (String) srec.getAttribute("SA");
		String intString = SATag.split(",|;")[1];
		int pos = Integer.parseInt(intString);
		String cigarString = SATag.split(",|;")[3];
		int indexFirstM = cigarString.indexOf("M");
		int indexFirstS = cigarString.indexOf("S");
		//int indexFirstH = cigarString.indexOf("H");
		if (indexFirstS < indexFirstM) {
			int posSM = Integer.parseInt(cigarString.substring(indexFirstS+1, indexFirstM));
			int position = pos+posSM-1;
			return position;
		}
		if (indexFirstS > indexFirstM) {
			int posSM = Integer.parseInt(cigarString.substring(0, indexFirstM));
			int position = pos+posSM-1;
			return position;
		}
		int position =-1;
		return position;
	}

	/**
	 * @param srec
	 * @return the forward or reverse flags from the first entry in the SA field.
	 */
	private static boolean getForwardSATag(SAMRecord srec) {
		String SATag = (String) srec.getAttribute("SA");
		String signString = SATag.split(",|;")[2];
		if (signString.equals("+")){
				return true;
		}
		else {
			return false;
		}
	}
	private static boolean getForwardSecondSATag(SAMRecord srec) {
		String SATag = (String) srec.getAttribute("SA");
		String signString = SATag.split(",|;")[8];
		if (signString.equals("+")){
				return true;
		}
		else {
			return false;
		}
	}
	public ArrayList<Translocation> getTranslocations() {
		ArrayList<Translocation> al = new ArrayList<Translocation>();
		for(String key: trans.keySet()) {
			al.addAll(trans.get(key));
		}
		return al;
	}

	public void print(BufferedWriter bw, BufferedWriter bwOut, int minSupport) {
		System.out.println("print started");
		if(bwOut!=null) {
			printEventsComplete(bwOut, true);
		}
		else {
			printEventsComplete(bwOut, false);
		}
		System.out.println("printEventsComplete done");
		System.out.println("printContents started");
       	printContents(bw, minSupport);
       	System.out.println("printContents done");
	}

	private void printEventsComplete(BufferedWriter bwOut, boolean print) {
		for(String key: trans.keySet()){
			for(Translocation tl: trans.get(key)) {
				if(tl.getSams().size()>0) {
					try {
						//bwOut.write(tl.getIGVPos()+"\n");
						String s = tl.getReads();
						if(print) {
							bwOut.write(s+"\r\n");
						}
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				
			}
		}
		
	}
	public Translocation searchTranslocation(Translocation tl) {
		if(this.searchRealPositions.size()==0) {
			this.hashAllRealPositions();
		}
		String search = tl.getIGVPos()+tl.isForward();
		Translocation found = searchRealPositions.get(search);
		return found;
	}

	private Translocation searchTranslocation(String string, int startDebug, int endDebug, boolean forward) {
		ArrayList<Translocation> al = trans.get(string);
		if(al==null) {
			System.out.println("no translocation in chr "+string);
			return null;
		}
		for(Translocation tl: al) {
			int pos = tl.getPosition1();
			if(pos>startDebug && pos < endDebug && tl.isForward() == forward) {
				//System.out.println("FOUND!");
				return tl;
			}
		}
		return null;
	}
	/**
	 * Finds and prints the locations of the LB and RB nicks
	 */
	public void testLBRB() {
		ReferenceSequenceFile rsf = ReferenceSequenceFileFactory.getReferenceSequenceFile(sp.getRef());
	    //System.out.println(rsf.isIndexed());
		//System.out.println("["+options.getRefFile()+"]");
	    String chr = sp.getChr();
	    List<SAMSequenceRecord> list = rsf.getSequenceDictionary().getSequences();
	    for(SAMSequenceRecord s: list) {
	    	System.out.println(s.getSequenceName());
	    }
	    ReferenceSequence rs = rsf.getSequence(chr);
	    
	    //System.out.println(rs==null);
	    //should give the position of the first occurrence of LB in TDNA which is indexLB
	    String TDNA = rs.getBaseString().toUpperCase();
	    int indexLB = TDNA.indexOf(LB);
	    // here LB and RB are declared to be forward
	    boolean LBisForward = true;
	    boolean RBisForward = true;
	    //if LB never occurs, reverse complement. If LB does occur, shift position by 4 because of the nick, and declare LB reverse
	    if(indexLB == -1) {
	    	indexLB = TDNA.indexOf(Utils.reverseComplement(LB));
	    	if(indexLB>=0) {
	    		indexLB+=4;
	    	}
	    	LBisForward = false;
	    }
	    //if the LB can be found, it should be forward. The index should then be shifted by 22.
	    else {
	    	indexLB +=22;
	    }
	    //same for RB, but +3 or +23
	    int indexRB = TDNA.indexOf(RB);
	    if(indexRB == -1) {
	    	indexRB = TDNA.indexOf(Utils.reverseComplement(RB));
	    	if(indexRB>=0) {
	    		indexRB+=3;
	    	}
	    	RBisForward = false;
	    }
	    else {
	    	indexRB +=23;
	    }
	    //if neither borders are found, give an error, show which is wrong, and stop the program
	    if(indexLB == -1 || indexRB == -1) {
	    	System.err.println("LB or RB could not be found in sequence "+chr);
	    	System.err.println("LB: "+indexLB);
		    System.err.println("RB: "+indexRB);
		    System.exit(0);
	    }
	    System.out.println("Setting LB "+indexLB+" "+LBisForward);
	    System.out.println("Setting RB "+indexRB+" "+RBisForward);
	    sp.setTDNALBPos(indexLB, LBisForward);
	    sp.setTDNARBPos(indexRB, RBisForward);
	}
	public void hashAllRealPositions() {
		for(Translocation tl: this.getTranslocations()) {
			if(tl.getSams().size()>0) {
				String key = tl.getIGVPos()+tl.isForward();
				//that is not supposed to happen, but can happen if two events are really close
				if(searchRealPositions.containsKey(key)) {
					tl.setMultipleEvents();
					searchRealPositions.get(key).setMultipleEvents();
				}
				else {
					searchRealPositions.put(tl.getIGVPos()+tl.isForward(),tl);
				}
			}
		}
	}
}
