package test;

import java.io.File;
import java.io.IOException;
import java.util.List;

import data.BufferedReferenceSequenceFile;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

public class ReferenceGenomeTester {

	public static void main(String[] args) {
		File testRef = new File("E:\\NGS\\Genome_Scan_104596\\Ref\\pUBC-YFP-mod0.fa");
		
		ReferenceSequenceFile rsf = ReferenceSequenceFileFactory.getReferenceSequenceFile(testRef);
	    BufferedReferenceSequenceFile bufferedRsf = new BufferedReferenceSequenceFile(rsf);
		
	    List<SAMSequenceRecord> chrs = rsf.getSequenceDictionary().getSequences();
	    boolean errors = false;
	    for(SAMSequenceRecord s: chrs) {
	    	//System.out.println(s.getSequenceName());
	    	for(int start = 100; start<s.getSequenceLength();start+=100000) {
			    int size = 50;
			    int end = start+size;
			    String seq = rsf.getSubsequenceAt(s.getSequenceName(), start, end).getBaseString();
			    String seq2 = bufferedRsf.getSubsequenceAt(s.getSequenceName(), start, end).getBaseString();
			    if(!seq.contentEquals(seq2)) {
			    	System.err.println(seq);
				    System.err.println(seq2);
				    System.err.println("====");
				    errors = true;
			    }
			    else {
			    	System.out.println("OK");
			    }
	    	}
	    }
	    try {
			rsf.close();
			bufferedRsf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    if(errors) {
	    	System.out.println("Something is not quite right in the BufferedReferenceSequenceFile");
	    }
	    else {
	    	System.out.println("All the same outcomes - Test Passed");
	    }
	    
	}

}
