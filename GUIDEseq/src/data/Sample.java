package data;

import java.util.Vector;

public class Sample {
	private String name;
	private Vector<String> v = new Vector<String>();
	public Sample(String sample, String umis2) {
		this.name = sample;
		addUMIs(umis2);
	}
	private void addUMIs(String umi) {
		String[] umis = umi.split("\\|");
		for(String umiString: umis) {
			v.add(umiString);
		}
	}
	public String getName() {
		return name;
	}
	public Vector<String> getUMIs() {
		return v;
	}

}
