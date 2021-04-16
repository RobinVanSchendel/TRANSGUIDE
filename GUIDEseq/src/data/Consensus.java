package data;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Map.Entry;
import java.util.stream.Collectors;

public class Consensus {
	private ArrayList<String> strings;
	public Consensus(ArrayList<String> list) {
		this.strings = list;
	}
	public Consensus() {
		strings = new ArrayList<String>();
	}
	public String getMostRepeatedString() {
		if(strings.size()==0) {
			return null;
		}
		String mostRepeatedWord 
	    = strings.stream()
	          .collect(Collectors.groupingBy(w -> w, Collectors.counting()))
	          .entrySet()
	          .stream()
	          .max(Comparator.comparing(Entry::getValue))
	          .get()
	          .getKey();
		return mostRepeatedWord;
	}
	public int getMostRepeatedStringNr() {
		String most = this.getMostRepeatedString();
		int count = 0;
		for(String s: strings) {
			if(s.contentEquals(most)) {
				count++;
			}
		}
		return count;
	}
	public double getMostRepeatedStringFraction() {
		return getMostRepeatedStringNr()/(double)strings.size();
	}
	public int getLongestString() {
		int maxSize = -1;
		for(String s: strings) {
			if(s.length()>maxSize) {
				maxSize = s.length();
			}
		}
		return maxSize;
	}
	public String getConsensusString() {
		StringBuffer str = new StringBuffer(); 
		for(int i = 0;i<getLongestString();i++) {
			Consensus c = new Consensus();
			for(String s: strings) {
				if(s.length()>i) {
					c.add(s.charAt(i));
				}
				//not sure if this is needed, but maybe it is
				else {
					//c.add(" ");
				}
			}
			if(c.getMostRepeatedStringFraction()==1) {
				str.append(c.getMostRepeatedString());
			}
			else {
				str.append(c.getMostRepeatedString().toLowerCase());
			}
		}
		return str.toString();
	}
	public String getConsensusStringReplaceByN(double minFreq) {
		StringBuffer str = new StringBuffer(); 
		for(int i = 0;i<getLongestString();i++) {
			Consensus c = new Consensus();
			for(String s: strings) {
				if(s.length()>i) {
					c.add(s.charAt(i));
				}
				//not sure if this is needed, but maybe it is
				else {
					//c.add(" ");
				}
			}
			if(c.getMostRepeatedStringFraction()>=minFreq) {
				str.append(c.getMostRepeatedString());
			}
			else {
				str.append('N');
			}
		}
		return str.toString();
	}
	private void add(char charAt) {
		strings.add(""+charAt);
	}
	public void add(String readPart) {
		if(readPart!=null) {
			strings.add(readPart);
		}
	}
	public int size() {
		return strings.size();
	}
	public String getStrings() {
		StringBuffer sb = new StringBuffer();
		for(String s: strings) {
			sb.append("[").append(s).append("]").append("\n");
		}
		return sb.toString();
	}

}
