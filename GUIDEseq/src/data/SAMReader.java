package data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class SAMReader {
	public static void main(String[] args) {
		for(String str: args) {
			//System.out.println(str);
		}
		Options optionsApache = createOptions();
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = null;
		try {
			cmd = parser.parse( optionsApache, args);
		} catch (ParseException e) {
			System.err.println(e);
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "GUIDEseqAnalyzer", optionsApache );
			System.exit(0);
		}
		MyOptions options = new MyOptions(cmd);
		String dir = "E:\\NGS\\GUIDEseq_Exp5\\300bpPlasmidMapped\\LZ30-1-L_S1.sorted.bam";
		boolean recursive = true;
		ArrayList<File> files = OutputCommands.searchSortedBam(new File(dir), recursive);
		PrimerController pc = new PrimerController(new File("Sample_Primer.txt"));
		
		boolean allFound = true;
		for(File f: files) {
			System.out.println(f.getName());
			SamplePrimer sp = pc.getSamplePrimer(f);
			if(sp==null) {
				allFound = false;
				System.err.println(f.getAbsolutePath()+ " does not have an entry");
			}
		}
		if(!allFound) {
			System.exit(0);
		}
		
		File out = new File("out.txt");
		int nr = 0;
		for(File f: files) {
			SamplePrimer sp = pc.getSamplePrimer(f);
			if(sp==null) {
				System.err.println(f.getAbsolutePath()+ " does not have an entry");
			}
			String primer = sp.getPrimer();
			String chr = sp.getChr();
			File ref = sp.getRef();
			options.setLB(sp.isLB());
			options.setRef(ref);
			options.setPrimer(primer);
			options.setFile(f);
			options.setChr(chr);
			
			BufferedWriter bw = null;
			if(nr==0) {
				try {
					bw = new BufferedWriter(new FileWriter(out,false));
					bw.write("File\t"+Translocation.getHeader()+"\n");
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			else {
				try {
					bw = new BufferedWriter(new FileWriter(out,true));
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			//System.out.println(options.printParameters());
	        TranslocationController tc = new TranslocationController(options);
	        tc.testLBRB();
	        tc.launchAnalysis(bw);
	        try {
				bw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	        nr++;
		}
		System.out.println("Output written to: "+out.getAbsolutePath());
    }

	private static Options createOptions() {
		Options o = new Options();
		//threads
		Option i = Option.builder( "i" )
				.longOpt("input")
                .hasArg()
                .argName("FILE")
                .desc("the input sorted and indexed bam file that is used to map the reads" )
                .required(true)
                .build();
		o.addOption(i);
		Option r   = Option.builder( "r" )
				.longOpt("reference")
                .hasArg()
                .argName("FILE")
                .desc("the reference file that is used to map the reads" )
                .required(true)
                .build();
		o.addOption(r);
		
		Option p = Option.builder("p")
				.longOpt("primer")
				.hasArg()
				.argName("STRING")
				.desc("The internal primer used in the GUIDEseq experiment")
				.required()
				.build();
		o.addOption(p);
		
		Option c = Option.builder("c")
				.longOpt("chr")
				.hasArg()
				.argName("STRING")
				.desc("The plasmid name that is integrated into the genome")
				.required()
				.build();
		o.addOption(c);
		
		Option p5   = Option.builder( "P5" )
				.longOpt("P5")
                .desc("Is the primer used as a P5 primer during sequencing")
                .optionalArg(true)
                .build();
				o.addOption(p5);
				
		Option p7   = Option.builder( "P7" )
				.longOpt("P7")
                .desc("Is the primer used as a P7 primer during sequencing")
                .optionalArg(true)
                .build();
		
				o.addOption(p7);
		
		Option out   = Option.builder( "o" )
				.longOpt("output")
                .hasArg()
                .argName("FILE")
                .desc("the file that the output is written to" )
                .build();
		o.addOption(out);
		
		Option m   = Option.builder( "m" )
				.longOpt("minsupport")
				.argName("NUMBER")
                .type(Number.class)
                .hasArg()
                .desc("the minimum number of support to be output (5 is standard)" )
                .build();
		o.addOption(m);
		
		
		return o;
	}

}
