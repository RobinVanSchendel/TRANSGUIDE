package data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Vector;

public class IGVPos {
	private String igvpos;
	private HashMap<String, Sample> hm = new HashMap<String, Sample>();

	public IGVPos(String igv) {
		this.igvpos = igv;
	}

	public void add(Sample sample) {
		hm.put(sample.getName(),sample);
	}
	public boolean containsDuplicateUMI(){
		Vector<String> v = new Vector<String>();
		for(String key: hm.keySet()) {
			Sample s = hm.get(key);
			Vector<String> umis = s.getUMIs();
			for(String umi: umis) {
				if(v.contains(umi)) {
					return true;
				}
				v.add(umi);
			}
		}
		//so no
		return false;
	}
	public int getnrUMIS() {
		int nr = 0;
		for(String key: hm.keySet()) {
			Sample s = hm.get(key);
			Vector<String> umis = s.getUMIs();
			nr+=umis.size();
		}
		return nr;
	}
	public int getnrduplUMIS() {
		Vector<String> v = new Vector<String>();
		int nrDups = 0;
		for(String key: hm.keySet()) {
			Sample s = hm.get(key);
			Vector<String> umis = s.getUMIs();
			for(String umi: umis) {
				if(v.contains(umi)) {
					nrDups++;
				}
				else {
					v.add(umi);
				}
			}
		}
		//so no
		return nrDups;
	}
	
	public String toString() {
		return igvpos;
	}

	public boolean hasMultipleSamples() {
		return hm.size()>1;
	}

	public int getNrSamples() {
		return hm.size();
	}

	public String getTop2() {
		Vector<Integer> v = new Vector<Integer>();
		for(String key: hm.keySet()) {
			v.add(hm.get(key).getUMIs().size());
		}
		Collections.sort(v,Collections.reverseOrder());
		String ret = "";
		for(int i: v) {
			if(ret.length()>0) {
				ret+="[";
			}
			ret+=i;
		}
		return ret;
	}

	public String getSampleNames() {
		ArrayList<String> v = new ArrayList<String>();
		for(String key: hm.keySet()) {
			v.add(key+":"+hm.get(key).getUMIs().size());
		}
		return String.join("[", v);
	}

	public int getMaxUMI() {
		Vector<Integer> v = new Vector<Integer>();
		for(String key: hm.keySet()) {
			v.add(hm.get(key).getUMIs().size());
		}
		Collections.sort(v,Collections.reverseOrder());
		return v.get(0);
	}
	
}
