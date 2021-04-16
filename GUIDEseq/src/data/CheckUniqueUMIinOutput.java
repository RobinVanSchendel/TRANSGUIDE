package data;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;
import java.util.Vector;

public class CheckUniqueUMIinOutput {

	public static void main(String[] args) {
		File f = new File("out.txt");
		HashMap<String, IGVPos> hm = new HashMap<String, IGVPos> ();
		try {
			Scanner s = new Scanner(f);
			String header = s.nextLine();
			int umiColumn = getColumnNr(header,"umis");
			int igvColumn = getColumnNr(header,"IGVPos");
			int isForwardColumn = getColumnNr(header,"isForward");
			int sampleColumn = getColumnNr(header,"Sample");
			//System.out.println(umiColumn);
			//System.out.println(igvColumn);
			//System.out.println(header);
			int lineNr = 1;
			int hits = 0;
			while(s.hasNextLine()) {
				String line = s.nextLine();
				String sampleN = line.split("\t")[sampleColumn];
				//skip those for now
				if(sampleN.contains("tim")) {
					continue;
				}
				String umis = line.split("\t")[umiColumn];
				String igv = line.split("\t")[igvColumn]+line.split("\t")[isForwardColumn];
				if(!hm.containsKey(igv)) {
					IGVPos igvpos = new IGVPos(igv);
					hm.put(igv,igvpos);
				}
				IGVPos pos = hm.get(igv);
				Sample sample = new Sample(sampleN, umis);
				pos.add(sample);
				lineNr++;
			}
			for(String key: hm.keySet()) {
				IGVPos igvpos = hm.get(key);
				System.out.println(igvpos+"\t"+igvpos.hasMultipleSamples()+"\t"+igvpos.getNrSamples()+"\t"+igvpos.containsDuplicateUMI()+"\t"+igvpos.getnrUMIS()+"\t"+igvpos.getnrduplUMIS()+"\t"+igvpos.getTop2()+"\t"+igvpos.getSampleNames()+"\t"+igvpos.getMaxUMI());
			}
			//System.out.println(hits+" found");
			s.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private static void addUMI(HashMap<String, Vector<String>> hm, String igv, String umi) {
		Vector<String> v = new Vector<String>();
		if(hm.containsKey(igv)) {
			v = hm.get(igv);
		}
		String[] umis = umi.split("\\|");
		for(String umiString: umis) {
			v.add(umiString);
		}
		hm.put(igv,v);
	}

	private static String containsDuplicateUMI(Vector<String> vector, String umi) {
		String[] umis = umi.split("\\|");
		String retString = "";
		for(String umiString: umis) {
			if(vector.contains(umiString)) {
				if(retString.length()>0){
					retString+="|";
				}
				retString+=umiString;
			}
		}
		if(retString.length()==0) {
			return null;
		}
		return retString;
	}

	private static int getColumnNr(String header, String string) {
		String[] parts = header.split("\t");
		for(int i=0;i<parts.length;i++) {
			if(parts[i].contentEquals(string)) {
				return i;
			}
		}
		return -1;
	}

}
