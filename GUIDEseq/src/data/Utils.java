package data;

public class Utils {
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
}
