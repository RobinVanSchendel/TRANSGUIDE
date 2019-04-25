package data;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import htsjdk.samtools.*;

public class SAMReader {
	public final static String str = "pCAMBIA3301_GUIDE";
	//HERE we used P7 primers in the TDNA, so the getFirstOfPairFlag is false
	public final static boolean getFirstOfPairFlag = false;
	public final static boolean forwardRB = true; 
	
	public static void main(String[] args) {
		File bamFile = new File("E:\\Project_GUIDESeq\\Rbplus.sorted.bam");
		
		final int minSupport = 5;
		
		
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);
        //SAMRecordIterator r = sr.iterator();
        SAMRecordIterator r = sr.query(str, 2276, 2278, false);
        ArrayList<SAMRecord> printable = new ArrayList<SAMRecord>();
        
        TranslocationController tc = new TranslocationController();
        int count = 0;
        while(r.hasNext()) {
        	SAMRecord srec = r.next();
        	if(Translocation.getNMis0(srec)) {
	        	if(srec.getContig() != null && !srec.getContig().equals(srec.getMateReferenceName())) {
		       		//System.out.println(srec.getContig()+":"+srec.getAlignmentStart()+"-"+srec.getAlignmentEnd());
		       		//System.out.println(srec.getReadString());
		       		//System.out.println(srec.getCigarString());
		       		if(!srec.getReadNegativeStrandFlag()) {
		       			//System.out.println("forward");
		       		}
		       		else {
		       			//System.out.println("reverse");
		       		}
		       		if(srec.getReadString().startsWith("CAAACTAGGATAAATTATCGCGCGCGGTGT")) {
		       			if(!srec.getDuplicateReadFlag() && srec.getFirstOfPairFlag() == getFirstOfPairFlag && !srec.getReadNegativeStrandFlag() == forwardRB) {
		       				//System.out.println(srec.getMateReferenceName()+"\t"+srec.getMateAlignmentStart()+"\t"+srec.getMateNegativeStrandFlag()+"\t"+srec.getCigarString()+"\t"+srec.getReadString()+"\t"+srec.isSecondaryAlignment()+"\t"+srec.getMappingQuality());
		       				printable.add(srec);
		       				//System.out.println("adding");
		       				//if(srec.getMateAlignmentStart()==4889469) {
		       					tc.addTranslocation(srec);
		       				//}
		       				count++;
		       			}
		       		}
	        	}
        	}
        }
        r.close();
        tc.addMates(sr);
        tc.printContents(minSupport);
        try {
			sr.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }

}
