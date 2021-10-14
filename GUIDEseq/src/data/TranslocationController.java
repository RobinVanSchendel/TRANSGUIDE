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
	public static final String RB = "GTTTACCCGCCAATATATCCTGTCA"; //this is the default one
	private SamplePrimer sp;
	private HashMap<String, Translocation> searchRealPositions = new HashMap<String, Translocation>();
	private boolean debug = false;
	public static final String testName = "A00379:349:HM7WFDSXY:4:1250:15031:28902";
	public static final String testPosition = "1:25672880";
	ReferenceSequenceFile rsf = null;
	
	public TranslocationController(SamplePrimer sp) {
		this.sp = sp;
		ReferenceSequenceFile rsfOrig = ReferenceSequenceFileFactory.getReferenceSequenceFile(sp.getRef());
	    rsf = new BufferedReferenceSequenceFile(rsfOrig);
	    this.debug = sp.debug();
	}
	
	private Translocation getNearestTranslocation(SAMRecordWrap s) {
		int pos = s.getPosition2(sp);
		if (pos!=-1) {  //filter for faulty positions
			Translocation nearest = null;
			int minDis = Integer.MAX_VALUE;
			if(trans.get(s.getMateReferenceName()) == null) {
				return null;
			}
			if (sp.getChr().equals(s.getMateReferenceName())==true) {
				if(s.getReadName().contentEquals(testName)) { 
					System.err.println("Mate is primarily aligned to the plasmid, nearest translocation will be null");}
				return null;
			}
			for(Translocation tl: trans.get(s.getMateReferenceName())) {
				if(tl.isForward() != s.getMateNegativeStrandFlag()) {
					int tempPos = tl.getPosition();
					int distance = Math.abs(pos-tempPos);
					if(distance<minDis) {
						minDis = distance;
						nearest = tl; 
					}
					else if(minDis>=0 && minDis<MAXDIST && nearest!=null) {
						return nearest;
					}
				}
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
	public Translocation addTranslocation(SAMRecordWrap s, int maxReads) {
		Translocation nearest = getNearestTranslocation(s);
		if(nearest !=null) {
			if(nearest.getSams().size()<maxReads) {
				nearest.addSam(s);
			}
			return nearest;
		}
		else if (s.getPosition2(sp)!=-1){
			Translocation tl = new Translocation(s, sp, rsf);
			if(s.getReadName().contentEquals(testName)) { 
			System.out.println("addTranslocation - checkpoint1");}
			//sometimes the sam is not added due to filtering of secondary alignments
			if((tl.getNrSupportingReads()>0) && (sp.getChr().equals(s.getMateReferenceName())==false)) {
				ArrayList<Translocation> al = trans.get(s.getMateReferenceName());
				if(al==null) {
					al = new ArrayList<Translocation>();
					trans.put(s.getMateReferenceName(), al);
				}
				al.add(tl);
				if(s.getReadName().contentEquals(testName)) { 
					System.out.println("addTranslocation - checkpoint2");}
				return tl;
			}
			else {
				//System.out.println(tl.getContigMate()+" has 0 supporting reads");
			}
		}
		if(s.getReadName().contentEquals(testName)) { 
			System.err.println("addTranslocation - no translocation added");}
		return null;
	}
    

	public void printContents(BufferedWriter bw, long minSupport) {
		//sorting doesn't work at the moment
		//trans.sort(comparing(Translocation::getNrSupportingReads, Collections.reverseOrder()));
		try {
			int counter = 0;
			for(String key: trans.keySet()){
				int printNr = 0;
				for(Translocation tl: trans.get(key)) {
					if(tl.getSams().size()>0) {
						if(tl.getNrPartialAnchors()>=minSupport) {
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
				SAMRecordIterator sri = sr.query(tl.getContigMate(), tl.getPosition()-around, tl.getPosition()+around, false);
				while(sri.hasNext()) {
					SAMRecord srec = sri.next();
					if(srec.getContig()!=null) {
						if(!isDuplicate(srec) && tl.containsRecord(srec.getReadName())  && srec.getFirstOfPairFlag()!=sp.isFirstOfPairFlag() && !srec.getContig().equals(sp.getChr())) {
							SAMRecordWrap s = new SAMRecordWrap(srec);
							tl.addSam(s);
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
			}
		}
		System.out.println("Finished processing "+count+" mates ");
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
	public void addRefGenomePart(long minSupport) {
		for(String key: trans.keySet()){
			for(Translocation tl: trans.get(key)) {
				if(tl.getNrPartialAnchors()>=minSupport) {
						tl.addRefSequence();
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
		System.out.println(sp.getFile());
		SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(sp.getFile());
		SamReader sr2 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(sp.getFile());
		if(!sr.hasIndex()) {
			System.out.println("no index available for this bam file. Please run samtools index "+sp.getFile());
			System.exit(0);
		}
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
    	System.out.println("Primers "+primerStart+" "+primerEnd+" "+positiveStrand+" "+sp.getPrimer());
    	//System.out.println("positions "+start+" - "+end);
    	//System.out.println("["+options.getChr()+"]");
    	//System.out.println(start);
    	//System.out.println(end);
    	//System.out.println(sr.hasIndex());
    	SAMRecordIterator r = sr.iterator();
    	System.out.println("Primerstart:" +primerStart +", Primerend:" +primerEnd +", positivestrand:"+positiveStrand +", primerseq:"+sp.getPrimer());
	    
	    int count = 0;
	    int NM0Count = 0;
	    int duplicateFlag = 0;
	    int Ncount = 0;
	 
        while(r.hasNext()) {
        	count++;
        	SAMRecord srec = r.next();
        	if(debug) {
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
        			if (temp!= null) {
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
        			if (temp == null) {
        				System.out.println("mate not found");
        			}
        		}
        	}
        	// don't take secondary alignment reads
        	// don't take reads with no secondary alignment, all T-DNA reads are expected to have at least 2 alignments
        	// don't take any reads that are duplicates
        	// take the read containing the T-DNA primer
        	if (!srec.isSecondaryAlignment() &&
    			srec.getFirstOfPairFlag() == sp.isFirstOfPairFlag() &&
    			!isDuplicate(srec) &&
    			srec.hasAttribute("SA")){
	        		//System.out.println(srec.getAttribute("SA"));
	        		String SATag = (String) srec.getAttribute("SA");
	        		String[] SAList = SATag.split(",|;");
	        		int SALength = SAList.length;
    				if(debug) {
    					if(srec.getReadName().contentEquals(testName)) {
            			System.err.println("launch analysis - checkpoint 1");
    					}
    				}
        			NM0Count++;
        			//System.out.println(srec.getContig()+" "+srec.getReadNegativeStrandFlag()+" "+srec.getAlignmentStart()+" "+srec.getReadString());
        			
        			if(srec.getReadNegativeStrandFlag()==true) {
        				if(debug) {
                		}
        				//This has to be done via this method as the
        				//reverseComplement method from SAMRecord contains a bug
        				//that does not reverse the Cigar
            			//Translocation.reverseComplementSAMRecord(srec);
            			srec.reverseComplement();

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
	                		System.err.println("launch analysis - checkpoint 2");
	            		}
	        		}
	        		//Checks if the record is mapped at all to the T-DNA plasmid.
	        		//There may be instances where checking whether it starts with the primer isn't enough, when the match with the T-DNA is too short to be picked up by the mapper.
	        		if (
	        			((SALength == 6) && (!((srec.getContig().equals(sp.getChr())==true) && (sp.getChr().equals(getContigSATag(srec))==true)))) ||
	        			((SALength >6) && (!((srec.getContig().equals(sp.getChr())==true) && (sp.getChr().equals(getContigSATag(srec))==true) && (sp.getChr().equals(getContigSecondSATag(srec))==true))))
	        			)
	        		{
	        			if (srec.getReadString().startsWith(sp.getPrimer())) {
	        				if(debug) {
	    	        			if(srec.getReadName().contentEquals(testName)) {
	    	        				System.err.println("launch analysis - checkpoint 3");
	    	        			}
	    	        		}
	    	        		if(srec.getReadString().contains("N")) {
	    	        			Ncount++;
	    	        		}
	    	        		addTranslocation(new SAMRecordWrap(srec), maxReadsPerTrans);	
	        			}
	        		}		       	
        		}
        		if(isDuplicate(srec)) {
       				duplicateFlag++;
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
        System.out.println("Found "+this.getTotalTranslocations()+" translocations");
        System.out.println("Found "+count+" reads");
        System.out.println("Found "+NM0Count+" NM0Count");
        System.out.println("Found "+duplicateFlag+" duplicateFlag (should be >0)");
        r.close();        
        System.out.println("Adding mates");
        addMates(sr);
        System.out.println("Added mates");
        System.out.println("Splitting T-DNA reads and finding minimal junction");
        addTDNASplit();
        System.out.println("Adding refGenomeParts");
        addRefGenomePart(sp.getMinSupport());
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
		if(srec.hasAttribute("SA")) {
			String SATag = (String) srec.getAttribute("SA");
			String intString = SATag.split(",|;")[1];
			return Integer.parseInt(intString);
		}
		return -1;
	}

	private static String getContigSATag(SAMRecord srec) {
		if(srec.hasAttribute("SA")) {
			String SATag = (String) srec.getAttribute("SA");
			return SATag.split(",|;")[0];
		}
		return null;
	}
	private static String getContigSecondSATag(SAMRecord srec) {
		if(srec.hasAttribute("SA")) {
			String SATag = (String) srec.getAttribute("SA");
			return SATag.split(",|;")[6];
		}
		return null;
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
		if (indexFirstS < indexFirstM) {
			int posSM = Integer.parseInt(cigarString.substring(indexFirstS+1, indexFirstM));
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
				if (tl.getIGVPos() == testPosition) {
					System.out.println("Junction found and printing");
				}
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
	private void addTDNASplit() {
		for(String key: trans.keySet()){
			for(Translocation tl: trans.get(key)) {
				if(tl.getSams().size()>0) {
					tl.splitTDNAReads();
					tl.getMinimumJunction();
				}
				
			}
		}
		
	}

	public static int getSACigarLength(SAMRecord s) {
		String SATag = (String) s.getAttribute("SA");
		String cigarString = SATag.split(",|;")[3];
		int countS = countChar(cigarString, 'S');
		int countM = countChar(cigarString, 'M');
		//int countI = countChar(cigarString, 'I');
		//int countD = countChar(cigarString, 'D');
		int totalCount = countS+countM;
		return totalCount;
	}
	public static int getSASecondCigarLength(SAMRecord s) {
		String SATag = (String) s.getAttribute("SA");
		String cigarString = SATag.split(",|;")[9];
		int countS = countChar(cigarString, 'S');
		int countM = countChar(cigarString, 'M');
		//int countI = countChar(cigarString, 'I');
		//int countD = countChar(cigarString, 'D');
		int totalCount = countS+countM;
		return totalCount;
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

	public Translocation searchTranslocation(Translocation tl) {
		if(this.searchRealPositions.size()==0) {
			this.hashAllRealPositions();
		}
		String search = tl.getIGVPos()+tl.isForward();
		Translocation found = searchRealPositions.get(search);
		return found;
	}

	/**
	 * Finds and prints the locations of the LB and RB nicks
	 */
	public void testLBRB() {
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
