package data;

import java.io.File;
import java.lang.reflect.Method;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;

public class MyOptions {
private CommandLine cmd;
private int RBpos = -1;
private int LBpos = -1;
	
	public MyOptions(CommandLine cmd) {
		this.cmd = cmd;
	}

	public File getBam() {
		String fileName = cmd.getOptionValue("input");
		File f = new File(fileName);
		if(!f.exists()) {
			System.err.println("Bam file "+f.getAbsolutePath()+" does not exist");
			System.exit(0);
		}
		return f;
	}

	public File getRefFile() {
		String fileName = cmd.getOptionValue("reference");
		File f = new File(fileName);
		if(!f.exists()) {
			System.err.println("Reference file "+f.getAbsolutePath()+" does not exist");
			System.exit(0);
		}
		return f;
	}
	public String getChr() {
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
		return cmd.getOptionValue("primer").toUpperCase();
	}
	public String getPrimerPart(int number) {
		String primer =  cmd.getOptionValue("primer").toUpperCase();
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
	public void setTDNARBPos(int pos) {
		this.RBpos = pos;
	}
	/** hard-coded the RB position
	 * 
	 */
	public int getTDNARBPos() {
		return RBpos;
	}
	public void setTDNALBPos(int pos) {
		this.LBpos = pos;
	}
	/** hard-coded the LB position
	 * 
	 */
	public int getTDNALBPos() {
		return LBpos;
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
		return 0;
	}
}
