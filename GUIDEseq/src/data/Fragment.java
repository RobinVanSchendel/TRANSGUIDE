package data;

public class Fragment {
private String frag;
	
	
	/**
	 * Creates a new Fragment based upon a String representing a sequence of nucleotides, 
	 * containing only the uppercase characters G, C, A and T.
	 * @param nucleotides
	 * @throws IllegalArgumentException if invalid characters are in the sequence of nucleotides 
	 */
	public Fragment(String nucleotides) throws IllegalArgumentException {
		
		boolean k = false;
		
		
		for(int i = 0; i < nucleotides.length(); i++){
			
			char lol = nucleotides.charAt(i);
			if(lol=='A'||lol=='G'||lol=='C'||lol=='T' || lol=='N'){
				k = true;
			}
			else{
				k = false;
			}
			
			if(k == false){
			throw new IllegalArgumentException("Dosent work: "+lol);
			
			}
			
			frag = nucleotides;
			
		}
		
		
		
	}
	
	/**
	 * Returns the length of this fragment.
	 * 
	 * @return the length of this fragment
	 */
	public int length() {
		
		return frag.length();
	}
	
	/**
	 * Returns a String representation of this fragment, exactly as was passed to the constructor.
	 * 
	 * @return a String representation of this fragment
	 */
	@Override
	public String toString() {
		return frag;
	}
	
	/**
	 * Return true if and only if this fragment contains the same sequence of nucleotides
	 * as another sequence.
	 */
	@Override
	public boolean equals(Object o) {
		
		boolean k = true;
		if(o instanceof Fragment == false){
			k = false;
		}
		
		else{
			Fragment oF = (Fragment) o;
			if(this.length()!=oF.length()) {
				return false;
			}
			else {
				for(int i = 0; i<frag.length(); i++){
					if(frag.charAt(i) == o.toString().charAt(i)){
						k = true;
					}
					else
						k = false;
				}
			}
			
		}
		
		return k;
	}
	
	/**
	 * Returns the number of nucleotides of overlap between the end of this fragment and
	 * the start of another fragment, f.
	 * 
	 * The largest overlap is found, for example, CAA and AAG have an maximum overlap of 2.
	 * 
	 * @param f the other fragment
	 * @return the number of nucleotides of overlap
	 */
	public int calculateOverlap(Fragment f) { 
		int overLap = 0;
	
		String underFrag = f.toString();
	    final int maxlen = Integer.min(frag.length(), underFrag.length());
	    
	    //shortcut
	    if(frag.contains(f.toString()) || f.toString().contains(frag)) {
	    	return Math.min(f.length(), length());
	    }
	    
	    for(int i = maxlen; i>0; i--){
	    	String part = f.toString().substring(0, i);
	    	//System.out.println("Trying\t"+part);
	    	if(this.frag.endsWith(part)){
	    		//System.out.println("hier!");
	    		if(i>overLap){
	    			overLap = i;
	    			//no need to search further
	    			break;
	    		}
	    	}
	    	
	    }
	    return overLap;
	    
	}
	
	
	/**
	 * Returns a new fragment based upon merging this fragment with another fragment f.
	 * 
	 * This fragment will be on the left, the other fragment will be on the right; the
	 * fragments will be overlapped as much as possible during the merge.
	 *  
	 * @param f the other fragment
	 * @return a new fragment based upon merging this fragment with another fragment 
	 */
	public Fragment mergedWith(Fragment f) {
		
		int val = this.calculateOverlap(f);
		Fragment g = new Fragment(this.frag+f.frag.substring(val));
		return g;
	}
}
