package data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class SimplifyTSVTable {

	public static void main(String[] args) {
		File fileIn = null;
		if(args.length==1) {
			String file = args[0];
			fileIn = new File(file);
		}
		else {
			String file = "C:\\Users\\rvanschendel\\git\\guideseq-t-dna\\GUIDEseq\\out.txt";
			fileIn = new File(file);
		}
		ArrayList<File> list = new ArrayList<File>();
		if(!fileIn.exists()) {
			System.err.println("File "+fileIn.getAbsolutePath()+" does noet exist or is a directory");
			System.exit(0);
		}
		if(fileIn.isDirectory()) {
			for(File tempFile: fileIn.listFiles()) {
				if(!tempFile.getName().contains("reduced") && !tempFile.getName().contains("stats")) {
					if(tempFile.getName().contains(".txt")) {
						list.add(tempFile);
						System.out.println(tempFile.getName());
					}
				}
			}
		}
		else if(fileIn.isFile()) {
			list.add(fileIn);
		}
		//System.exit(0);
		
		
		for(File f: list) {
			String outputFile = f.getAbsolutePath()+"_reduced.txt";
			if(f.getAbsolutePath().endsWith(".txt")) {
				outputFile = f.getAbsolutePath().replace(".txt", "_reducedCols.txt");
			}
			File outputF = new File(outputFile);
			long sizeBytes = f.length();
			long sizeKB = sizeBytes/1024;
			long sizeMB = sizeKB/1024;
			boolean keepUnique = true;
			String uniqueColumn = "IGVPos";
			try {
				Scanner s = new Scanner(f);
				int line = 0;
				String[] keepColumns = {"Sample","Genotype","DNA","DNANonMatching","DNAMatching","multipleEvents","tailSize","NrSupportingReads","NrAnchors","NrAnchorsIncludingPartial","NrSupportJunction","NrNotSupportJunction","NrSupportJunctionHQ","NrNotSupportJunctionGQ","countLBWeird","NrPrimaryAlignment","NrSecondaryAlignment","Chr","Position","RealPosition","realPositionCounter","DistanceToLBRB","LB/RB","IGVPos","isForward","CigarString","JunctionSequence","TDNASequence","GenomicSequence","getTDNASequenceNew","getGenomicSequenceNew"
						,"getGenomicSequenceIntNew","getGenomicSequenceIntMinMaxNew","getGenomicSequenceIntMedian"
						,"getGenomicSequenceIntMean","getGenomicSequenceIntSD","getMostRepeatedPreGenome","getMostRepeatedConsensus","Type","Homology","HomologyLength","Filler","FillerLength","FillerIsTemplated","getLargestMatchString","getSubS","getSubS2","getType","getLengthS","getPosS","getFirstHit","getFirstPos","TotalLength","RefSequence","error","warning","isOK","getTDNASequenceDiffReads","getTDNASequenceDiffReadsHighestContributor"
						};
				
				ArrayList<Integer> colKeepHM = null;
				FileWriter writer = new FileWriter(outputF);
				BufferedWriter bw = new BufferedWriter(writer);
				
				while(s.hasNextLine()) {
					String curLine = s.nextLine();
					if(line==0) {
						String[] cols = curLine.split("\t");
						int count = 0;
						//for(String col: cols) {
							//System.out.println(count+"\t"+col);
							colKeepHM = createKeepColumns(cols,keepColumns);
							count++;
						//}
						//System.out.println(String.join("\",\"", cols));
						//System.exit(0);
					}
					String[] cols = curLine.split("\t");
					StringBuffer output = new StringBuffer(); 
					for(int colNr=0;colNr<cols.length;colNr++) {
						if(colKeepHM.contains(colNr)) {
							if(output.length()>0) {
								output.append("\t");
							}
							output.append(cols[colNr]);
						}
					}
					output.append("\n");
					bw.write(output.toString());
					line++;
				}
				s.close();
				bw.close();
				long sizeBytesOut = outputF.length();
				long sizeKBytesOut = sizeBytesOut/1024;
				long sizeMBytesOut = sizeKBytesOut/1024;
				System.out.println(f.getAbsolutePath()+ " size: "+sizeMB+"Mb");
				System.out.println(outputF.getAbsolutePath()+" reduced to: "+sizeMBytesOut+"Mb");
				
				
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	private static ArrayList<Integer> createKeepColumns(String[] cols, String[] keepColumns) {
		ArrayList<Integer> colKeepHM = new ArrayList<Integer>();
		int colCounter = 0;
		for(String col: cols) {
			for(String kC: keepColumns) {
				if(col.contentEquals(kC)) {
					colKeepHM.add(colCounter);
				}
			}
			colCounter++;
		}
		return colKeepHM;
	}

}
