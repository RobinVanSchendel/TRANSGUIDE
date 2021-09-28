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

public class ConsolidateReads {

	public static void main(String[] args) {
		String dir = "E:\\NGS\\Genome_Scan_104596\\Raw\\UMI";
		//int maxReads = 4000000;
		int maxReads = 4000000;
		if(args.length>1) {
			dir = args[1];
		}
		
		ArrayList<File> R1s = getR1s(new File(dir));
		for(File f: R1s) {
			//if(!f.getName().contains("HM7LNDSXY_104269-001-001-168_TAAGGCGA-GCGATCTA_L004_R1.umi.fastq.gz")) {
			//	continue;
			//}
			//if(getFileSizeMegaBytes(f)<100)
			System.out.println(f.getName());
			File R1 = f;
			File R2 = new File(f.getAbsolutePath().replace("R1.umi.fastq.gz", "R3.umi.fastq.gz"));
			
			String outR1String = R1.getName().replace(".fastq.gz", ".cons.fastq");
			String outR2String = R2.getName().replace(".fastq.gz", ".cons.fastq");
			File outR1 = new File(R1.getParent()+File.separatorChar+outR1String);
			File outR2 = new File(R2.getParent()+File.separatorChar+outR2String);
			
			File gz1 = new File(R1.getParent()+File.separatorChar+outR1String+".gz");
			//skip done files
			if(gz1.exists()) {
				System.out.println(gz1.getAbsolutePath()+" exists");
				continue;
			}
			
			//System.out.println(outR1.getAbsolutePath());
			//System.out.println(outR2.getAbsolutePath());
			//System.exit(0);
			
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
			
			StreamingIterator<FastqRecord> iterR1 = datastoreR1.iterator();
			StreamingIterator<FastqRecord> iterR2 = datastoreR2.iterator();
			
			int counter = 0;
			HashMap<String, Reads> reads = new HashMap<String, Reads>();
			int nrNonAdded = 0;
			
			while(iterR1.hasNext()) {
				FastqRecord fqR1 = iterR1.next();
				FastqRecord fqR2 = iterR2.next();
				String[] parts = fqR1.getId().split("BC:");
				String umi = parts[1]+fqR1.getNucleotideSequence().toString().substring(0, 6)+fqR2.getNucleotideSequence().toString().substring(0, 6);
				
				//System.out.println(umi);
				//System.exit(0);
				
				//FastqRecord fqTestR1 = fqR1.toBuilder().comment("BC:"+fqUmi.getNucleotideSequence()).build();
				//FastqRecord fqTestR2 = fqR2.toBuilder().comment("BC:"+fqUmi.getNucleotideSequence()).build();
				Reads r = reads.get(umi);
				if(r==null) {
					r = new Reads();
					reads.put(umi,r);
				}
				boolean hasAdded = r.addReads(fqR1, fqR2);
				if(!hasAdded) {
					nrNonAdded++;
				}
				//added
				else {
					counter++;
					if(counter%100000==0) {
						System.out.println("Already processed "+counter+" reads");
					}
				}
				//writerR1.write(fqTestR1);
				//writerR2.write(fqTestR2);
				if(counter>=maxReads) {
					break;
				}
			}
			System.out.println("We have "+counter+" reads");
			System.out.println("We have skipped "+nrNonAdded+" reads");
			System.out.println("We have "+reads.size()+" UMIs");
			//now consolidate each read and write
			for(String key: reads.keySet()) {
				FastqRecord fqR1 = reads.get(key).consolidateR1();
				FastqRecord fqR2 = reads.get(key).consolidateR2();
				writerR1.write(fqR1);
				writerR2.write(fqR2);
			}
			iterR1.close();
			iterR2.close();
			writerR1.close();
			writerR2.close();
			datastoreR1.close();
			datastoreR2.close();
			compressGZIP(outR1);
			//delete large file now
			outR1.delete();
			compressGZIP(outR2);
			//delete large file now
			outR2.delete();
			System.out.println(outR1.getAbsolutePath());
			System.out.println(outR2.getAbsolutePath());
			
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}
	private static ArrayList<File> getR1s(File file) {
		ArrayList<File> ret = new ArrayList<File>();
		if(file.isDirectory()) {
			for(File f: file.listFiles()) {
				if(f.getName().endsWith("R1.umi.fastq.gz")) {
					ret.add(f);
				}
			}
		}
		else {
			ret.add(file);
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
	private static double getFileSizeMegaBytes(File file) {
		return (double) file.length() / (1024 * 1024);
	}

}
