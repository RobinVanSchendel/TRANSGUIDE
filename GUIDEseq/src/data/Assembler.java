package data;

import java.util.ArrayList;
import java.util.List;

public class Assembler {
	public List<Fragment> Fragments = new ArrayList<Fragment>();
	/**
	 * Creates a new Assembler containing a list of fragments.
	 * 
	 * The list is copied into this assembler so that the original list
	 * will not be modified by the actions of this assembler.
	 * 
	 * @param fragments
	 */
	public Assembler(List<Fragment> fragments) {
		
		for(int i = 0; i<fragments.size(); i++){
			Fragments.add(fragments.get(i));
		}
		
	}
	public Assembler() {
		// TODO Auto-generated constructor stub
	}
	public void addFragment(Fragment f) {
		Fragments.add(f);
	}
	
	/**
	 * Returns the current list of fragments this assembler contains.
	 * @return the current list of fragments
	 */
	public List<Fragment> getFragments() {
		return Fragments;
	}
	public String toString() {
		String s = "";
		for(Fragment f: Fragments) {
			if(s.length()>0) {
				s+="\n";
			}
			s+=f.toString();
		}
		return s;
	}
	
	/**
	 * Attempts to perform a single assembly, returning true if an assembly was performed.
	 * 
	 * This method chooses the best assembly possible, that is, 
	 * it merges the two fragments with the largest overlap, breaking ties
	 * between merged fragments by choosing the shorter merged fragment.
	 * 
	 * Merges must have an overlap of at least 1.
	 * 
	 * After merging two fragments into a new fragment, the new fragment is inserted into
	 * the list of fragments in this assembler, and the two original fragments are removed
	 * from the list.
	 * 
	 * @return true iff an assembly was performed
	 */
	public boolean assembleOnce() {
	 
		//System.out.println("assembleOnce");
		int c = 0;
		int maxoverlap = -1;
		Fragment m1 = null, m2 = null, m = null;
		for(Fragment frag1 : Fragments){
				for(Fragment frag2 : Fragments){
					if(frag1 == null || frag2 == null){
						return false;
					}
					if(frag1!=frag2){
						c = frag1.calculateOverlap(frag2);
						//System.out.println("overlap "+c);
						//System.out.println("maxoverlap "+maxoverlap);
						if(c>maxoverlap){
							maxoverlap = c;
							m1 = frag1;
							m2 = frag2;
						}
					}
					
				}
		}
		
		if(maxoverlap>0){
			Fragments.remove(m1);
			Fragments.remove(m2);
			m = m1.mergedWith(m2);
			Fragments.add(m);
			return true;
		}
		
		else {
			return false;
		}
		
		
		

	}
	
	/**
	 * Repeatedly assembles fragments until no more assembly can occur.
	 */
	public void assembleAll() {
		
		while(assembleOnce() !=false){
			assembleOnce();
			
		}
			
		
	}
}