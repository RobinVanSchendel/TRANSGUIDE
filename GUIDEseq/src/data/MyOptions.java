package data;

import java.io.File;
import java.lang.reflect.Method;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;

public class MyOptions {
private CommandLine cmd;
private int RBpos = -1;
private int LBpos = -1;
private boolean LBisForward;
private boolean RBisForward;
private boolean isLB;
private File ref;
private String primer;
private File bam;
private String chr;
	
	public MyOptions(CommandLine cmd) {
		this.cmd = cmd;
	}

	public File getBam() {
		if(bam!=null) {
			return bam;
		}
		String fileName = cmd.getOptionValue("input");
		File f = new File(fileName);
		if(!f.exists()) {
			System.err.println("Bam file "+f.getAbsolutePath()+" does not exist");
			System.exit(0);
		}
		return f;
	}
	public boolean isLBOption() {
		if(cmd.hasOption("LB")) {
			this.isLB = true;
			return true;
		}
		if(!cmd.hasOption("RB")){
			System.err.println("No LB or RB specified");
			System.exit(0);
		}
		this.isLB = false;
		return false;
	}

	public File getRefFile() {
		if(this.ref!=null) {
			return ref;
		}
		String fileName = cmd.getOptionValue("reference");
		File f = new File(fileName);
		if(!f.exists()) {
			System.err.println("Reference file "+f.getAbsolutePath()+" does not exist");
			System.exit(0);
		}
		return f;
	}
	public String getChr() {
		if(this.chr!=null) {
			return chr;
		}
		return cmd.getOptionValue("chr");
	}
	public boolean isFirstOfPairFlag() {
		if(cmd.hasOption("P5")) {
			return true;
		}
		if(cmd.hasOption("P7")) {
			return false;
		}
		System.err.println("Either P5 or P7 has to be specified");
		System.exit(0);
		return false;
	}

	public String getPrimer() {
		if(primer!=null) {
			return primer.toUpperCase();
		}
		String primer =  cmd.getOptionValue("primer").toUpperCase();
		return primer.replace("TCAGACGTGTGCTCTTCCGATCT", "");
	}
	public String getPrimerPart(int number) {
		String primer =  cmd.getOptionValue("primer").toUpperCase();
		if(this.primer!=null) {
			primer = this.primer.toUpperCase();
		}
		if(primer.length()>number) {
			return primer.substring(0, number);
		}
		return primer;
	}

	public long getMinSupport() {
		if(cmd.hasOption("minsupport")) {
			try {
				return (Long) cmd.getParsedOptionValue("minsupport");
			} catch (ParseException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return 1;
	}
	public void setTDNARBPos(int pos, boolean rBisForward) {
		this.RBisForward = rBisForward;
		this.RBpos = pos;
	}
	/** hard-coded the RB position
	 * 
	 */
	public int getTDNARBPos() {
		return RBpos;
	}
	public void setTDNALBPos(int pos, boolean lBisForward) {
		this.LBisForward = lBisForward;
		this.LBpos = pos;
	}
	/** hard-coded the LB position
	 * 
	 */
	public int getTDNALBPos() {
		return LBpos;
	}
	public boolean getTDNALBisForward() {
		return this.LBisForward;
	}
	public boolean getTDNARBisForward() {
		return this.RBisForward;
	}
	public String printParameters() {
		if(cmd!= null) {
			Class<?> c = this.getClass();
	        Method[] methods = c.getDeclaredMethods();
	        StringBuffer sb = new StringBuffer();
	        for(Method m: methods) {
	       		if(!m.getName().equals("printParameters") && !m.getName().equals("getPrimerPart") 
	       				&& !m.getName().equals("getClass") && !m.getName().equals("hasOptions")
	       				&& !m.getName().startsWith("set")) {
	       			try {
	       				if(m.invoke(this) != null){
	       					if(sb.length()>0) {
	    	            		sb.append("\n");
	    	            	}
		       				String str = m.invoke(this).toString();
		       				String name = m.getName();
		       				name = name.replace("get","");
							sb.append(name+": "+str);//+":"+m.invoke(this));
	       				}
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} 
					//sb.append(m);//+":"+m.invoke(this));
	       		}
	        }
	        return sb.toString();
		}
		return null;
	}
	public String getOutput() {
		if(cmd.hasOption('o')) {
			return cmd.getOptionValue('o');
		}
		else {
			return getBam().getAbsolutePath()+"output.txt";
		}
	}

	public static int getMaxMismatches() {
		return 1;
	}

	public void setLB(boolean lb) {
		this.isLB = lb;
	}
	public boolean isLB() {
		return this.isLB;
	}

	public void setRef(File ref) {
		this.ref = ref;
	}

	public void setPrimer(String primer) {
		this.primer = primer;
	}

	public void setFile(File f) {
		this.bam = f;
		
	}

	public void setChr(String chr) {
		this.chr = chr;
	}

	public int getMaxTranslocations() {
		return -1;
	}

	public int getMaxReadsPerTrans() {
		return 1000;
	}

	public int getMaxReads() {
		return -1;
	}

	public String getSample() {
		return cmd.getOptionValue("n");
	}

	public String getDNAsample() {
		if(cmd.hasOption("dna")) {
			return cmd.getOptionValue("dna");
		}
		return "";
	}

	public String getGenotype() {
		if(cmd.hasOption("genotype")) {
			return cmd.getOptionValue("genotype");
		}
		return "";
	}

	public String getEcotype() {
		if(cmd.hasOption("ecotype")) {
			return cmd.getOptionValue("ecotype");
		}
		return "";
	}

	public String getP5() {
		if(cmd.hasOption("P5")) {
			return "P5";
		}
		else if(cmd.hasOption("P7")) {
			return "P7";
		}
		return "";
	}

	public boolean isUMI() {
		// TODO Auto-generated method stub
		return false;
	}
}
