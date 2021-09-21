package data;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class Utils {
	public static int totalSAMRecordReadSize(SAMRecord srec) {
		if(!srec.getCigarString().contains("H")) {
			return srec.getReadLength();
		}
		else {
			int length = 0;
			//also include the 'H' element
			for(CigarElement ce: srec.getCigar()) {
				length+=ce.getLength();
			}
			return length;
		}
	}
	public static String longestCommonSubstring(String S1, String S2)
	{
		int Start = 0;
	    int Max = 0;
	    for (int i = 0; i < S1.length(); i++)
	    {
	        for (int j = 0; j < S2.length(); j++)
	        {
	            int x = 0;
	            while (Character.toUpperCase(S1.charAt(i + x)) == Character.toUpperCase(S2.charAt(j + x)))
	            {
	                x++;
	                if (((i + x) >= S1.length()) || ((j + x) >= S2.length())) break;
	            }
	            if (x > Max)
	            {
	                Max = x;
	                Start = i;
	            }
	         }
	    }
	    return S1.substring(Start, (Start + Max));
	}
	public static String reverseComplement(String dna){
		StringBuffer dnaRev = new StringBuffer(dna).reverse();
		StringBuffer revCom = new StringBuffer();
		for(char c: dnaRev.toString().toCharArray()){
			switch(c){
				case 'a':
					revCom.append('t');
					break;
				case 'A':
					revCom.append('T');
					break;
				case 't':
					revCom.append('a');
					break;
				case 'T':
					revCom.append('A');
					break;
				case 'c':
					revCom.append('g');
					break;
				case 'C':
					revCom.append('G');
					break;
				case 'g':
					revCom.append('c');
					break;
				case 'G':
					revCom.append('C');
					break;
				case 'N':
					revCom.append('N');
					break;
				case 'n':
					revCom.append('n');
					break;
				//for second hit only!
				case 'x':
					revCom.append('x');
					break;
				default:
					//System.out.println("Can't complement "+dna);
					//System.out.println("Can't complement "+c);
					//System.exit(0);
			}
		}
		return revCom.toString();
	}
	public static boolean cigarStringFollowsSMS(String cigarString) {
		Pattern p = Pattern.compile("\\d*[S]\\d*[M]\\d*[S]");
		Matcher m = p.matcher(cigarString);
		boolean b = m.matches();
		return b;
	}
	/**
	 * Creates a stringbuffer with cigarstrings. If the length of "s" is larger than "size", return the part of the string from the start position until "size".
	 * If the length is smaller than "size", add spaces to the end until that length has been reached. Then return the string plus those spaces. 
	 * @param cigarString
	 * @param size
	 * @return
	 */
	public static String getString(String str, int size) {
		StringBuffer s = new StringBuffer(str);
		if(s.length()>size) {
			return s.toString().substring(0, size);
		}
		while(s.length()<size) {
			s.append(" ");
		}
		return s.toString();
	}
}
