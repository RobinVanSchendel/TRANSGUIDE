package data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

import static java.util.Comparator.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class TranslocationController {
	private HashMap<String, ArrayList<Translocation>> trans = new HashMap<String, ArrayList<Translocation>>();
	private MyOptions options;
	public static final int MAXDIST = 500;
	private int primerStart, primerEnd;
	public static final String LB = "GTTTACACCACAATATATCCTGCCA";
	public static final String RB = "GTTTACCCGCCAATATATCCTGTCA";
	
	public TranslocationController(MyOptions options) {
		this.options = options;
	}

	public Translocation addTranslocation(SAMRecord s) {
		Translocation nearest = getNearestTranslocation(s);
		int pos = Translocation.getPosition(s, options);
		if(nearest !=null) {
			nearest.addSam(s);
			return nearest;
		}
		else {
			Translocation tl = new Translocation(s, options);
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

	private Translocation getNearestTranslocation(SAMRecord s) {
		int pos = Translocation.getPosition(s, options);
		Translocation nearest = null;
		int minDis = Integer.MAX_VALUE;
		if(trans.get(s.getMateReferenceName()) == null) {
			//System.out.println("don't yet have "+s.getMateReferenceName());
			return null;
		}
		for(Translocation tl: trans.get(s.getMateReferenceName())) {
			//System.out.println(tl.getContigMate()+":"+s.getMateReferenceName());
			//getMateNegativeStrandFlag == Forward
			//System.out.println("same Contig");
			//System.out.println(tl.isForward()+":"+!s.getMateNegativeStrandFlag());
			//if(s.getMateReferenceName().contentEquals("5")) {
				//System.out.println(s.getMateReferenceName()+" contains "+trans.get(s.getMateReferenceName()).size()+" search "+pos);
				//System.out.println(tl.isForward() != s.getMateNegativeStrandFlag());
				/*
				for(Translocation tl2: trans.get(s.getMateReferenceName())){
					System.out.println("\t"+tl2.getPosition());
					System.out.println(tl.isForward() != s.getMateNegativeStrandFlag());
					int tempPos = tl.getPosition();
					int distance = Math.abs(pos-tempPos);
					System.out.println("\tdistance:\t"+distance);
				}
				*/
			//}
			if(tl.isForward() != s.getMateNegativeStrandFlag()) {
				//System.out.println("same orientation");
				int tempPos = tl.getPosition();
				int distance = Math.abs(pos-tempPos);
				//System.out.println("distance: "+distance);
				//if(distance>20000000) {
				//	System.out.println(tl.getContigMate()+"\t"+tl.getPosition());
				//	System.out.println(pos);
				//}
				if(distance<minDis) {
					minDis = distance;
					nearest = tl; 
				}
				//small optimization
				if(minDis == 0) {
					//System.out.println("shortcut! "+minDis);
					//System.out.println("found it "+s.getMateReferenceName()+" "+ pos);
					return nearest;
				}
				//is that correct?
				else if(minDis>0 && minDis<MAXDIST && nearest!=null) {
					//System.out.println("shortcut "+minDis);
					//System.out.println("found it "+s.getMateReferenceName()+" "+ pos);
					return nearest;
				}
			}
		}
		//System.out.println(s.getMateReferenceName() + " "+pos);
		//System.out.println(minDis);
		//only do this when we have a nearest
		if (minDis<MAXDIST) {
			//System.out.println("found it "+s.getMateReferenceName()+" "+ pos);
			return nearest;
		}
		//System.out.println("Creating new "+s.getMateReferenceName()+" "+ pos);
		//System.out.println(minDis);
		return null;
	}

	public void printContents(long minSupport) {
		//sorting doesn't work at the moment
		//trans.sort(comparing(Translocation::getNrSupportingReads, Collections.reverseOrder()));
		System.out.println("File\t"+Translocation.getHeader());
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(options.getOutput()));
			writer.write("File\t"+Translocation.getHeader()+"\n");
			for(String key: trans.keySet()){
				for(Translocation tl: trans.get(key)) {
					if(tl.getTranslocationSequence().contentEquals("TTGAGCTTGGATCAGATTGTCGTTTCCCGCCTTCAGTTTAAACTATCAGTGTTTGAACAAATAACACATTGTGGTGTTTAATGAATCGTGGTGGGATATATTGGCTAGAGCAGCTTGAGCTTGAAATGGAAAGGAGTGAAGAGTAAAGAAG")) {
						//tl.printDebug();
						//System.out.println("hier!");
						//System.exit(0);
					}
					if(tl.getNrSupportingReads()>=minSupport) {
						String output = tl.toString();
						System.out.println(options.getBam().getName()+"\t"+output);
						writer.write(options.getBam()+"\t"+output+"\n");
					}
				}
			}
			writer.close();
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
				SAMRecordIterator sri = sr.query(tl.getContigMate(), tl.getPosition()-around, tl.getPosition()+around, false);
				int counter2 = 0;
				while(sri.hasNext()) {
					SAMRecord srec = sri.next();
					//System.out.println(srec.getReadName());
					if(!srec.getDuplicateReadFlag() && tl.containsRecord(srec.getReadName())) {
						//System.out.println("adding "+srec.getContig()+":" +srec.getAlignmentStart()+"-"+srec.getAlignmentEnd());
						//System.out.println("adding "+srec.getReadName());
						boolean added = tl.addSam(srec);
					}
					if(srec.getDuplicateReadFlag()) {
						duplicates++;
					}
					//System.out.println(counter2);
				}
				sri.close();
				if(count%1000==0) {
					System.out.println("Already processed "+count+" mates "+trans.get(key).size() + " from chr "+key);
				}
			}
		}
		System.out.println("Found "+duplicates+" duplicate mates (should be >0)");
	}

	public void addRefGenomePart(ReferenceSequenceFile rsf) {
		for(String key: trans.keySet()){
			for(Translocation tl: trans.get(key)) {
				if(tl.isOK()) {
					tl.addRefSequence(rsf);
				}
			}
		}
	}

	public void launchAnalysis() {
		SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(options.getBam());
		//SamReader sr2 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(options.getBam());
		if(!sr.hasIndex()) {
			System.out.println("no index available for this bam file. Please run samtools index "+options.getBam());
			System.exit(0);
		}
	    ReferenceSequenceFile rsf = ReferenceSequenceFileFactory.getReferenceSequenceFile(options.getRefFile());
	    //SAMRecordIterator r = sr.iterator();
	    String chr = options.getChr();
	    ReferenceSequence rs = rsf.getSequence(chr);
	    if(rs == null) {
	    	System.err.println("Chromosome "+chr+" is not part of the reference genome");
	    	System.exit(0);
	    }
	    else {
	    	//System.out.println("Chromosome "+chr+" was found");
	    }
	    String seq = rs.getBaseString().toUpperCase();
	    int index = seq.indexOf(options.getPrimer());
	    int start = -1;
	    int end = -1;
	    boolean positiveStrand = true;
	    if(index == -1) {
	    	index = seq.indexOf(Utils.reverseComplement(options.getPrimer()));
	    	positiveStrand = false;
	    }
	    else {
	    	//System.out.println("Primer found on + strand at position: "+index+" ("+options.getPrimer()+")");
	    }
	   
	    if(index == -1) {
	    	System.err.println("Primer could not be found");
	    	System.exit(0);
	    }
		start = index;
    	end = index+options.getPrimer().length();
    	primerStart = start;
    	primerEnd = end;
    	//System.out.println("positions "+start+" - "+end);
    	//System.out.println(options.getChr());
    	//System.out.println(start);
    	//System.out.println(end);
    	//System.out.println(sr.hasIndex());
	    SAMRecordIterator r = sr.query(options.getChr(), start, end, false);
	    
	    int count = 0;
	    int NM0Count = 0;
	    int duplicateFlag = 0;
	    boolean debug = false;
        while(r.hasNext()) {
        	count++;
        	SAMRecord srec = r.next();
        	if(debug) {
        		//System.err.println("Still here");
        	}
        	//only take 0 mismatches reads
        	//20200920 only take reads with max cigar length of 2
        	if(Translocation.getNMis0(srec) && srec.getCigarLength()<=2) {
        		if(debug) {
            		//System.err.println("Still here NM0");
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
        			srec.reverseComplement();
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
                		//System.err.println("Still here getting translocation");
                		//System.err.println(srec.getReadString());
                	}
	        		//read should start with the primer
		       		if(srec.getReadString().startsWith(options.getPrimerPart(20))) {
		       			if(debug) {
	                		//System.err.println("Still here primer part is ok");
	                	}
		       			//System.out.println("starts correct");
		       			//no duplicates
		       			
		       			if(!srec.getDuplicateReadFlag() && srec.getFirstOfPairFlag() == options.isFirstOfPairFlag() && !srec.getReadNegativeStrandFlag() == positiveStrand) {
		       				//System.out.println("adding");
		       				if(debug) {
		       					
		       				}
	       					addTranslocation(srec);
		       			}
		       			if(srec.getDuplicateReadFlag()) {
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
	        			
			        		if(!srec.getDuplicateReadFlag() && srec.getFirstOfPairFlag() == options.isFirstOfPairFlag() && !srec.getReadNegativeStrandFlag() == positiveStrand) {
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
        		System.out.println("Already processed "+count+" reads, NM0 reads "+NM0Count);
        		//break;
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
        printContents(options.getMinSupport());
        try {
			sr.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        if(debug) {
        	String debugChr = "5";
        	int startDebug = 24099423-5000;
    		int endDebug = startDebug+10000;
    		boolean forward = false;
    		Translocation tl = searchTranslocation(debugChr,startDebug, endDebug, forward);
    		if(tl!=null) {
    			System.out.println("DEBUG FROM HIER!");
    			System.out.println(tl);
    			tl.printDebug();
    			System.out.println(tl.isForward());
    			System.exit(0);
    		}
    		else {
    			System.out.println("not found");
    		}
		}
	}

	private Translocation searchTranslocation(String string, int startDebug, int endDebug, boolean forward) {
		ArrayList<Translocation> al = trans.get(string);
		if(al==null) {
			System.out.println("no translocation in chr "+string);
			return null;
		}
		for(Translocation tl: al) {
			int pos = tl.getPosition();
			if(pos>startDebug && pos < endDebug && tl.isForward() == forward) {
				//System.out.println("FOUND!");
				return tl;
			}
		}
		return null;
	}

	public void testLBRB() {
		ReferenceSequenceFile rsf = ReferenceSequenceFileFactory.getReferenceSequenceFile(options.getRefFile());
	    //SAMRecordIterator r = sr.iterator();
	    String chr = options.getChr();
	    ReferenceSequence rs = rsf.getSequence(chr);
	    String TDNA = rs.getBaseString().toUpperCase();
	    int indexLB = TDNA.indexOf(LB);
	    //maybe reverse complement?
	    if(indexLB == -1) {
	    	indexLB = TDNA.indexOf(Utils.reverseComplement(LB));
	    	if(indexLB>=0) {
	    		indexLB+=3;
	    	}
	    }
	    //nick position
	    else {
	    	indexLB +=22;
	    }
	    int indexRB = TDNA.indexOf(RB);
	    if(indexRB == -1) {
	    	indexRB = TDNA.indexOf(Utils.reverseComplement(RB));
	    	if(indexRB>=0) {
	    		indexRB+=2;
	    	}
	    }
	    else {
	    	indexRB +=23;
	    }
	    if(indexLB == -1 && indexRB == -1) {
	    	//System.err.println("LB or RB could not be found in sequence "+chr);
	    	//System.err.println("LB: "+indexLB);
		    //System.err.println("RB: "+indexRB);
		    //System.exit(0);
	    	indexLB = 1;
		    indexRB = TDNA.length();
	    }
	    options.setTDNALBPos(indexLB);
	    options.setTDNARBPos(indexRB);
	}
}
