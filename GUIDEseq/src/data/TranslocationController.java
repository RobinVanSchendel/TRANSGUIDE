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
	//temporarily changed from 500 to 0.
	public static final int MAXDIST = 10;
	private int primerStart, primerEnd;
	public static final String LB = "GTTTACACCACAATATATCCTGCCA";
	public static final String RB = "GTTTACCCGCCAATATATCCTGTCA";
	private SamplePrimer sp;
	private HashMap<String, Translocation> searchRealPositions = new HashMap<String, Translocation>();
	boolean debug = true;
	
	public TranslocationController(SamplePrimer sp) {
		this.sp = sp;
	}
    
	//I think this whole block needs to be altered.
	private Translocation getNearestTranslocation(SAMRecord s) {
		if (Translocation.getPosition2(s, sp)!=-1) {
			//pos is ok
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
			//sometimes the sam is not added due to filtering of secondary alignments
			if(tl.getNrSupportingReads()>0) {
				ArrayList<Translocation> al = trans.get(s.getMateReferenceName());
				if(al==null) {
					al = new ArrayList<Translocation>();
					trans.put(s.getMateReferenceName(), al);
				}
				al.add(tl);
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

	public void addMates(SamReader sr) {
		int around = 500;
		int duplicates = 0;
		int count = 0;
		for(String key: trans.keySet()){
			for(Translocation tl: trans.get(key)) {
				count++;
				//System.out.println(tl.getContigMate());
				SAMRecordIterator sri = sr.query(tl.getContigMate(), tl.getPosition1()-around, tl.getPosition1()+around, false);
				int counter2 = 0;
				//System.out.println(tl.getContigMate()+":"+tl.getPosition1());
				String testName = "A00379:349:HM7WFDSXY:4:2223:3685:2206";
				while(sri.hasNext()) {
					SAMRecord srec = sri.next();
					//System.out.println(srec.getReadName());
					if(srec.getReadName().contentEquals(testName)) {
						//System.err.println("found it"+srec.getReadName()+" "+srec.getReadLength());
					}
					if(srec.getContig()!=null) {
						if(srec.getReadName().contentEquals(testName)) {
							//System.err.println("found it still here"+srec.getReadName()+" "+srec.getCigarString());
							//System.err.println("found it still here"+srec.getReadName()+" "+isDuplicate(srec));
							//System.err.println("found it still here"+srec.getReadName()+" "+tl.containsRecord(srec.getReadName()));
							
							
						}
						//if first read is in, why remove the second if it is duplicate?
						/*
						if(tl.containsRecord(srec.getReadName())) {
							System.out.println(srec.getDuplicateReadFlag()+" "+isDuplicate(srec));
							if(isDuplicate(srec)) {
								System.out.println(srec.getCigarString());
								tl.printMate(srec);
							}
						*/
						if(!isDuplicate(srec) && tl.containsRecord(srec.getReadName())) {
							//System.out.println("adding "+srec.getContig()+":" +srec.getAlignmentStart()+"-"+srec.getAlignmentEnd());
							//System.out.println("adding "+srec.getReadName() + " "+srec.getContig());
							if(srec.getReadName().contentEquals(testName)) {
								System.err.println("found it adding"+srec.getReadName()+" "+srec.getReadLength());
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
						//not sure if I want these
						/*
						else if(!isDuplicate(srec)) {
							if(tl.isForward()) {
								if(srec.getFirstOfPairFlag()!= sp.isFirstOfPairFlag() && !srec.getReadNegativeStrandFlag()) {
									tl.addSam(srec);
								}
								else if(srec.getFirstOfPairFlag()== sp.isFirstOfPairFlag() && srec.getReadNegativeStrandFlag()) {
									tl.addSam(srec);
								}
							}
							else {
								if(srec.getFirstOfPairFlag()!= sp.isFirstOfPairFlag() && srec.getReadNegativeStrandFlag()) {
									tl.addSam(srec);
								}
								else if(srec.getFirstOfPairFlag()== sp.isFirstOfPairFlag() && !srec.getReadNegativeStrandFlag()) {
									tl.addSam(srec);
								}
								
							}
						}
						*/
						if(isDuplicate(srec)) {
							duplicates++;
						}
					}
					//System.out.println(counter2);
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
    	System.out.println(primerStart+" "+primerEnd+" "+positiveStrand+" "+sp.getPrimer());
    	//System.out.println("positions "+start+" - "+end);
    	//System.out.println("["+options.getChr()+"]");
    	//System.out.println(start);
    	//System.out.println(end);
    	//System.out.println(sr.hasIndex());
	    SAMRecordIterator r = sr.query(sp.getChr(), start, end, false);
	    
	    int count = 0;
	    int NM0Count = 0;
	    int duplicateFlag = 0;
	    int Ncount = 0;
	    int NoCount = 0;
	 
	    String debugReadName = "A00379:349:HM7WFDSXY:4:2223:3685:2206";
        while(r.hasNext()) {
        	count++;
        	SAMRecord srec = r.next();
        	if(debug) {
        		//System.err.println("Still here");
        		if(srec.getReadName().contentEquals(debugReadName)) {
        			System.out.println(srec.getReadName());
        			System.out.println(srec.getCigarString());
        			System.out.println(srec.getContig());
        			System.out.println(srec.toString());
        			System.out.println(srec.getReadString());
        			System.out.println("MATE");
        			SAMRecord temp = sr2.queryMate(srec);
        			System.out.println(temp.getReadName());
        			System.out.println(temp.getCigarString());
        			System.out.println(temp.getContig());
        			System.out.println(temp.toString());
        			System.out.println(temp.getReadString());
        			System.out.println("EO MATE");
        		}
        	}
        	//only take 0 mismatches reads
        	//20200920 only take reads with max cigar length of 2
        	if(Translocation.getNMis0(srec) ) {//&& srec.getCigarLength()<=2) {
        		if(debug) {
            		//System.err.println("Still here");
            		if(srec.getReadName().contentEquals(debugReadName)) {
                		System.err.println("Still here");
            		}
        		}
        		NM0Count++;
        		//chr has to be filled and unequal
        		//take the reverse complement if we are looking at reads going reverse
        		if(!positiveStrand) {
        			if(debug) {
                		//System.err.println("Still here !positiveStrand");
                		//System.err.println(srec.getReadString());
                	}
        			String cigar = srec.getCigarString();
        			boolean currentNegativeStrandFlag = srec.getReadNegativeStrandFlag();
        			//System.out.println("BEFORE RC "+srec.getReadNegativeStrandFlag());
        			srec.reverseComplement();
        			if(srec.getReadNegativeStrandFlag()==currentNegativeStrandFlag) {
        				srec.setReadNegativeStrandFlag(!currentNegativeStrandFlag);
        			}
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
        		}
	        	if(srec.getContig() != null && !srec.getContig().equals(srec.getMateReferenceName())) {
	        		if(debug) {
	        			if(srec.getReadName().contentEquals(debugReadName)) {
	                		System.err.println("Still here2 ");
	            		}
                	}
	        		//read should start with the primer
		       		//if(srec.getReadString().startsWith(sp.getPrimerPart(20))) {
	        		//System.out.println(srec.getAlignmentStart() + " " + srec.getAlignmentEnd());
	        		//System.out.println(srec.getReadString());
	       			if((positiveStrand && srec.getAlignmentStart()==primerStart) ||
	       					(!positiveStrand && srec.getAlignmentEnd()==primerEnd)){
		       			if(debug) {
		       				if(srec.getReadName().contentEquals(debugReadName)) {
		                		System.err.println("Still here3");
		                		System.err.println("is duplicate "+srec.getDuplicateReadFlag());
		                		System.err.println("is sa "+srec.isSecondaryAlignment());
		                		System.err.println("is sa/su "+srec.isSecondaryOrSupplementary());
		            		}
	                	}
		       			//System.out.println("starts correct "+srec.getReadNegativeStrandFlag() );
		       			//System.out.println(positiveStrand);
		       			//System.out.println(srec.getDuplicateReadFlag());
		       			//System.out.println(srec.getFirstOfPairFlag() == sp.isFirstOfPairFlag());
		       			//no duplicates
		       			
		       			if(!isDuplicate(srec) && srec.getFirstOfPairFlag() == sp.isFirstOfPairFlag()) {// && srec.getReadNegativeStrandFlag() == positiveStrand) {
		       				//System.out.println("adding "+srec.getReadName()+" "+ isDuplicate(srec));
		       				if(debug) {
			       				if(srec.getReadName().contentEquals(debugReadName)) {
			                		System.err.println("Still here4");
			                		System.err.println(srec.isSecondaryAlignment());
			       				}
		       					
		       				}
		       				if(srec.getReadString().contains("N")) {
		       					Ncount++;
		       				}
		       				else {
		       					NoCount++;
		       				}
	       					addTranslocation(srec, maxReadsPerTrans);
		       			}
		       			if(isDuplicate(srec)) {
		       				duplicateFlag++;
		       			}
		       		}
		       		
	        	}
	        	/*
	        	else if (srec.getContig() != null && srec.getContig().equals(srec.getMateReferenceName())) {
	        		//System.out.println(positiveStrand);
	        		//System.out.println(srec.getAlignmentStart());
	        		if(srec.getReadString().startsWith(options.getPrimerPart(20))) {
	        			
	        			int dist = Math.abs(srec.getMateAlignmentStart()-primerStart);
		       			if(dist>500) {
	        			
			        		if(!srec.getDuplicateReadFlag() && srec.getFirstOfPairFlag() == sp.isFirstOfPairFlag() && !srec.getReadNegativeStrandFlag() == positiveStrand) {
			       				//System.out.println("adding");
			       				if(debug) {
			       					
			       				}
		       					Translocation tl = addTranslocation(srec);
		       					
		       					
		       					
		       					SAMRecord mate = sr2.queryMate(srec);

		       					//System.out.println(srec.toString());
		       					//System.out.println(srec.getReadString());
		       					if(mate.getCigarString().contentEquals("299M")) {
		       						System.out.println("####START");
		       						//System.out.println(srec.getMateAlignmentStart());
			       					//System.out.println(tl);
			       					System.out.println(mate.getReadString());
			       					System.out.println(mate.getCigarString());
			       					System.out.println(mate.getContig()+":"+mate.getAlignmentStart()+"-"+mate.getAlignmentEnd());
			       					System.out.println("####END");
		       					}
			       			}
		       			}
	        		}
	        	}
	        	*/
        	}
        	if(count%1000000==0) {
        		System.out.println("Already processed "+count+" reads, NM0 reads "+NM0Count+" Ncount "+Ncount+" Nocount: "+NoCount);
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
	    	//following 4 lines are obsolete?
		    //LBisForward = false;
	    	//RBisForward = true;
	    	//indexLB = 1;
		   // indexRB = TDNA.length();
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
