package data;

import java.io.File;
import java.util.ArrayList;
import java.util.Vector;

public class OutputCommands {

	public static void main(String[] args) {
		String dir = "E:\\NGS\\GUIDEseq_Exp6";
		//String dir = "E:\\NGS\\all_bams\\analysisTDNAonly151";
		//String dir = "E:\\NGS\\GUIDEseq\\sortedbams";
		//String refseq = "E:\\NGS\\GUIDEseq\\RefSeq\\Arabidopsis_thaliana.TAIR10.28_plus_pCAMBIA2201.fa";
		String refseq = "E:\\NGS\\GUIDEseq_Exp2\\RefSeq\\Arabidopsis_thaliana.TAIR10.28_pUBC.dna.genome.fa";
		//String refseq = "E:\\NGS\\GUIDEseq_Exp2\\RefSeq\\Arabidopsis_thaliana.TAIR10.28_pUBC_tdna.dna.genome.fa";
		ArrayList<File> files = searchSortedBam(new File(dir));
		PrimerController pc = new PrimerController(new File("Sample_Primer.txt"));
		//String chr= "pUBC-YFP-Dest";
		//String chr= "-c pUBC_tDNA_part";
		for(File f: files) {
			System.out.println(f.getName()+"\t"+pc.getPrimer(f));
		}
		for(File f: files) {
			//if(f.getName().contains("rmdup")) {
				//if(f.getName().contains("LZ")) {
					//String primer = getPrimer(f);
					String primer = pc.getPrimer(f);
					String chr = pc.getChr(f);
					//if(f.getName().startsWith("LB")) {
						Vector<String> v = new Vector<String>();
						v.add("-i");
						v.add("\""+f.getAbsolutePath()+"\"");
						v.add("-p");
						v.add(primer);
						v.add("-c");
						v.add(chr);
						v.add("-r");
						v.add(pc.getRefSeq(f).getAbsolutePath());
						v.add("-P7");
						String s = "-i \""+f.getAbsolutePath()+"\" "+primer+" "+chr+" -r "+refseq +" -P7";
						//System.out.println(s);
					//}
					//else {
						//System.out.println("-i \""+f.getAbsolutePath()+"\" "+rb);
					//}
				//}
						//if(f.getName().contains("LZ34-1-L_S5.sorted.bam")) {
							SAMReader.main(v.toArray(new String[v.size()]));
						//}
			//}
		}
	}

	private static String getPrimer(File f) {
		String part = f.getName().replace(".sorted.bam", "");
		part = part.replace("_rmdup", "");
		String primer = "";
		if(part.endsWith("LZ003")) {
			primer = "agctgataattcaattcgg";
		}
		else if(part.endsWith("LZ004")) {
			primer = "ctgataattcaattcggcgtta";
		}
		else if(part.endsWith("LZ007")) {
			primer = "ttgagcttggatcagatt";
		}
		else if(part.endsWith("LZ008")) {
			primer = "cttgagcttggatcaga";			
		}
		else if(part.startsWith("LZ30_1_S")) {
			primer = "GATAATTCAATTCGGCGTTA";
		}
		else if(part.startsWith("LZ30_3_2")) {
			primer = "AGCTGATAATTCAATTCGG";
		}
		else if(part.startsWith("LZ30_3_3")) {
			primer = "GATAATTCAATTCGGCGTTA";
		}
		else if(part.startsWith("LZ30_3_4")) {
			primer = "AGCTGATAATTCAATTCGG";
		}
		else if(part.startsWith("LZ30_3_5")) {
			primer = "TTGAGCTTGGATCAGATT";
		}
		else if(part.startsWith("LZ30_3_6")) {
			primer = "CTTGAGCTTGGATCAGA";
		}
		else if(part.contentEquals("LZ30_mix_3_S1")) {
			primer = "GATAATTCAATTCGGCGTTA";
		}
		else if(part.contentEquals("LZ30_mix_4_S2")) {
			primer = "CTTGAGCTTGGATCAGA";
		}
		else if(part.contains("-L_")) {
			primer = "GATAATTCAATTCGGCGTTA";
		}
		else if(part.contains("-R_")) {
			primer = "CTTGAGCTTGGATCAGA";
		}
		
		else {
			System.err.println("I do not know that ending "+part);
		}
		return primer;
	}

	private static ArrayList<File> searchSortedBam(File file) {
		ArrayList<File> files = new ArrayList<File>();
		if(file.isDirectory()) {
			for(File f: file.listFiles()) {
				if(f.getAbsolutePath().endsWith(".sorted.bam")) {
					files.add(f);
				}
			}
		}
		return files;
	}

}
