package data;

import java.io.File;
import java.util.ArrayList;

public class OutputCommands {

	public static void main(String[] args) {
		String dir = "E:\\NGS\\GUIDEseq\\sortedbams";
		String refseq = "E:\\NGS\\GUIDEseq\\RefSeq\\Arabidopsis_thaliana.TAIR10.28_plus_pCAMBIA2201.fa";
		ArrayList<File> files = searchSortedBam(new File(dir));
		String lb = "-p ttaattcagtacattaaaaacg -P7 -r "+refseq+" -c pCAMBIA3301_GUIDE";
		String rb = "-p caaactaggataaattatcg  -P7 -r "+refseq+" -c pCAMBIA3301_GUIDE";
		for(File f: files) {
			if(f.getName().contains("Total")) {
				if(f.getName().startsWith("LB")) {
					System.out.println("-i \""+f.getAbsolutePath()+"\" "+lb);
				}
				else {
					System.out.println("-i \""+f.getAbsolutePath()+"\" "+rb);
				}
			}
		}
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
