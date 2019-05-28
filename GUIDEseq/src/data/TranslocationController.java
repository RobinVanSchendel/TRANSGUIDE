package data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

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
	private ArrayList<Translocation> trans = new ArrayList<Translocation>();
	private MyOptions options;
	public static final int MAXDIST = 500;
	
	public TranslocationController(MyOptions options) {
		this.options = options;
	}

	public void addTranslocation(SAMRecord s) {
		Translocation nearest = getNearestTranslocation(s);
		if(nearest !=null) {
			nearest.addSam(s);
		}
		else {
			Translocation tl = new Translocation(s, options);
			//sometimes the sam is not added due to filtering of secondary alignments
			if(tl.getNrSupportingReads()>0) {
				trans.add(tl);
			}
		}
	}

	private Translocation getNearestTranslocation(SAMRecord s) {
		int pos = Translocation.getPosition(s, options);
		Translocation nearest = null;
		int minDis = Integer.MAX_VALUE;
		for(Translocation tl: trans) {
			//System.out.println(tl.getContigMate()+":"+s.getMateReferenceName());
			if(tl.getContigMate().equals(s.getMateReferenceName())) {
				//getMateNegativeStrandFlag == Forward
				//System.out.println("same Contig");
				//System.out.println(tl.isForward()+":"+!s.getMateNegativeStrandFlag());
				if(tl.isForward() != s.getMateNegativeStrandFlag()) {
					//System.out.println("same orientation");
					int tempPos = tl.getPosition();
					int distance = Math.abs(pos-tempPos);
					//System.out.println("distance: "+distance);
					if(distance<minDis) {
						minDis = distance;
						nearest = tl; 
					}
				}
			}
		}
		//only do this when we have a nearest
		if (minDis<MAXDIST) {
			return nearest;
		}
		return null;
	}

	public void printContents(long minSupport) {
		trans.sort(comparing(Translocation::getNrSupportingReads, Collections.reverseOrder()));
		System.out.println(Translocation.getHeader());
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(options.getOutput()));
			writer.write("File\t"+Translocation.getHeader()+"\n");
			for(Translocation tl: trans) {
				if(tl.getNrSupportingReads()>=minSupport) {
					String output = tl.toString();
					System.out.println(output);
					writer.write(options.getBam()+"\t"+output+"\n");
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
		for(Translocation tl: trans) {
			SAMRecordIterator sri = sr.query(tl.getContigMate(), tl.getPosition()-around, tl.getPosition()+around, false);
			while(sri.hasNext()) {
				SAMRecord srec = sri.next();
				//System.out.println(srec.getReadName());
				if(!srec.getDuplicateReadFlag() && tl.containsRecord(srec.getReadName())) {
					boolean added = tl.addSam(srec);
				}
				if(srec.getDuplicateReadFlag()) {
					duplicates++;
				}
			}
			sri.close();
		}
		System.out.println("Found "+duplicates+" duplicate mates (should be >0)");
	}

	public void addRefGenomePart(ReferenceSequenceFile rsf) {
		for(Translocation tl: trans) {
			if(tl.isOK()) {
				tl.addRefSequence(rsf);
			}
		}
		
	}

	public void launchAnalysis() {
		SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(options.getBam());
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
	    	System.out.println("Chromosome "+chr+" was found");
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
	    if(index == -1) {
	    	System.err.println("Primer could not be found");
	    	System.exit(0);
	    }
		start = index;
    	end = index+options.getPrimer().length();
    	System.out.println("positions "+start+" - "+end);
	    SAMRecordIterator r = sr.query(options.getChr(), start, end, false);
	    int count = 0;
	    int NM0Count = 0;
	    int duplicateFlag = 0;
        while(r.hasNext()) {
        	count++;
        	SAMRecord srec = r.next();
        	//only take 0 mismatches reads
        	if(Translocation.getNMis0(srec)) {
        		NM0Count++;
        		//chr has to be filled and unequal
        		//take the reverse complement if we are looking at reads going reverse
        		if(!positiveStrand) {
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
        		}
	        	if(srec.getContig() != null && !srec.getContig().equals(srec.getMateReferenceName())) {
	        		//read should start with the primer
		       		if(srec.getReadString().startsWith(options.getPrimerPart(20))) {
		       			//no duplicates
		       			if(!srec.getDuplicateReadFlag() && srec.getFirstOfPairFlag() == options.isFirstOfPairFlag() && !srec.getReadNegativeStrandFlag() == positiveStrand) {
	       					addTranslocation(srec);
		       			}
		       			if(srec.getDuplicateReadFlag()) {
		       				duplicateFlag++;
		       			}
		       		}
		       		
	        	}
        	}
        }
        System.out.println("Strand is forward: "+positiveStrand);
        System.out.println("Found "+count+" reads");
        System.out.println("Found "+NM0Count+" NM0Count");
        System.out.println("Found "+duplicateFlag+" duplicateFlag (should be >0)");
        r.close();
        addMates(sr);
        addRefGenomePart(rsf);
        printContents(options.getMinSupport());
        try {
			sr.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
