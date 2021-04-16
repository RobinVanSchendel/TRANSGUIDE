package data;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Map.Entry;
import java.util.stream.Collectors;

public class ConsensusInt{
	private ArrayList<Integer> ints;
	public ConsensusInt(ArrayList<Integer> list) {
		this.ints = list;
	}
	public ConsensusInt() {
		ints = new ArrayList<Integer>();
	}
	public int getMostRepeatedInt() {
		if(ints.size()==0) {
			return Integer.MIN_VALUE;
		}
		int mostRepeatedWord 
	    = ints.stream()
	          .collect(Collectors.groupingBy(w -> w, Collectors.counting()))
	          .entrySet()
	          .stream()
	          .max(Comparator.comparing(Entry::getValue))
	          .get()
	          .getKey();
		return mostRepeatedWord;
	}
	public int getMostRepeatedIntNr() {
		int most = this.getMostRepeatedInt();
		int count = 0;
		for(int s: ints) {
			if(s == most) {
				count++;
			}
		}
		return count;
	}
	public double getMostRepeatedStringFraction() {
		return getMostRepeatedIntNr()/(double)ints.size();
	}
	public int getMin() {
		int current = Integer.MAX_VALUE;
		for(int i: ints) {
			if(i<current) {
				current = i;
			}
		}
		return current;
	}
	public int getMax() {
		int current = Integer.MIN_VALUE;
		for(int i: ints) {
			if(i>current) {
				current = i;
			}
		}
		return current;
	}
	public void add(int readPart) {
		ints.add(readPart);
	}
	public int size() {
		return ints.size();
	}
	public double getMedian() {
		if(ints.size()==0) {
			return 0;
		}
		Integer[] numArray = ints.toArray(new Integer[ints.size()]);
		Arrays.sort(numArray);
		double median;
		if (numArray.length % 2 == 0)
		    median = ((double)numArray[numArray.length/2] + (double)numArray[numArray.length/2 - 1])/2;
		else
		    median = (double) numArray[numArray.length/2];
		return median;
	}
	public double getMean() {
		double total = 0;
		for(int i: ints) {
			total+=i;
		}
		return total/ints.size();
	}
	public double getSD()
    {
        double standardDeviation = 0.0;
        double mean = getMean();

        for(int num: ints) {
            standardDeviation += Math.pow(num - mean, 2);
        }
        return Math.sqrt(standardDeviation/ints.size());
    }
}
