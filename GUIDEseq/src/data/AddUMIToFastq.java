package data;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.zip.GZIPOutputStream;

import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.core.util.iter.StreamingIterator;
import org.jcvi.jillion.trace.fastq.FastqDataStore;
import org.jcvi.jillion.trace.fastq.FastqFileDataStoreBuilder;
import org.jcvi.jillion.trace.fastq.FastqQualityCodec;
import org.jcvi.jillion.trace.fastq.FastqRecord;
import org.jcvi.jillion.trace.fastq.FastqWriter;
import org.jcvi.jillion.trace.fastq.FastqWriterBuilder;

public class AddUMIToFastq {

	public static void main(String[] args) {
		String dir = "E:\\NGS\\Genome_Scan_104596\\Raw";
		ArrayList<File> R1s = getR1s(new File(dir));
		for(File f: R1s) {
			System.out.println(f.getName());
		
			File R1 = f;
			File R2 = new File(f.getAbsolutePath().replace("R1.fastq.gz", "R3.fastq.gz"));
			//umis are in R2
			File umi = new File(f.getAbsolutePath().replace("R1.fastq.gz", "R2.fastq.gz"));
			
			String outR1String = R1.getName().replace(".fastq.gz", ".umi.fastq");
			String outR2String = R2.getName().replace(".fastq.gz", ".umi.fastq");
			File outR1 = new File(R1.getParent()+File.separatorChar+outR1String);
			File outR2 = new File(R2.getParent()+File.separatorChar+outR2String);
			
			try {
				FastqWriter writerR1 = new FastqWriterBuilder(outR1).build();
				FastqWriter writerR2 = new FastqWriterBuilder(outR2).build();
			
			FastqDataStore datastoreR1 = new FastqFileDataStoreBuilder(R1)
					.qualityCodec(FastqQualityCodec.SANGER)
	                .hint(DataStoreProviderHint.ITERATION_ONLY)
	                .build();
			FastqDataStore datastoreR2 = new FastqFileDataStoreBuilder(R2)
					.qualityCodec(FastqQualityCodec.SANGER)
	                .hint(DataStoreProviderHint.ITERATION_ONLY)
	                .build();
			FastqDataStore umiStore = new FastqFileDataStoreBuilder(umi)
					.qualityCodec(FastqQualityCodec.SANGER)
	                .hint(DataStoreProviderHint.ITERATION_ONLY)
	                .build();
			
			StreamingIterator<FastqRecord> iterR1 = datastoreR1.iterator();
			StreamingIterator<FastqRecord> iterR2 = datastoreR2.iterator();
			StreamingIterator<FastqRecord> iterUmi = umiStore.iterator();
			
			int counter = 0;
			while(iterR1.hasNext()) {
				FastqRecord fqR1 = iterR1.next();
				FastqRecord fqR2 = iterR2.next();
				FastqRecord fqUmi = iterUmi.next();
				FastqRecord fqTestR1 = fqR1.toBuilder().comment("BC:"+fqUmi.getNucleotideSequence()).build();
				FastqRecord fqTestR2 = fqR2.toBuilder().comment("BC:"+fqUmi.getNucleotideSequence()).build();
				writerR1.write(fqTestR1);
				writerR2.write(fqTestR2);
				counter++;
				if(counter%100000==0) {
					System.out.println("Already processed "+counter+" reads");
				}
			}
			iterR1.close();
			iterR2.close();
			iterUmi.close();
			writerR1.close();
			writerR2.close();
			datastoreR1.close();
			datastoreR2.close();
			umiStore.close();
			compressGZIP(outR1);
			compressGZIP(outR2);
			
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}
	private static ArrayList<File> getR1s(File file) {
		ArrayList<File> ret = new ArrayList<File>();
		for(File f: file.listFiles()) {
			if(f.getName().endsWith("R1.fastq.gz")) {
				ret.add(f);
			}
		}
		return ret;
	}
	public static void compressGZIP(File input) throws IOException {
		File output = new File(input.getAbsolutePath()+".gz");
        try (GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(output))){
            try (FileInputStream in = new FileInputStream(input)){
                byte[] buffer = new byte[1024];
                int len;
                while((len=in.read(buffer)) != -1){
                    out.write(buffer, 0, len);
                }
            }
        }
    }

}
