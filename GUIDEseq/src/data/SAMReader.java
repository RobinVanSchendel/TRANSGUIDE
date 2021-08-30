package data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

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
		//String dir = "E:\\NGS\\Genome_Scan_104269\\umi_cons_bam\\";
		//String dir = "E:\\NGS_data\\bamfiles\\300bp_uncropped\\";
		String dir = "E:\\NGS_data\\bamfiles\\";
		boolean recursive = true;
		boolean combineFiles = false;
		//default == 1 (for now, because duplicate position filter needs to come before anchor filter)
		int minSupport = 1;
		//default 500
		int maxReadsPerTrans = 1000;
		//default -1 == all
		int maxTranslocation = -1;
		//default -1 == all
		int maxReads = -1;
		
		ArrayList<File> files = OutputCommands.searchSortedBam(new File(dir), recursive);
		PrimerController pc = new PrimerController(new File("Sample_Primer.txt"));
		//PrimerController pc = new PrimerController(new File("Sample_Primer_RBcutter.txt"));
		File out = new File("out.txt");
		
		BufferedWriter allEvents = null;
		BufferedWriter bw = null;
		try {
			allEvents = new BufferedWriter(new FileWriter("allEvents.txt"));
			bw = new BufferedWriter(new FileWriter(out));
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		boolean allFound = true;
		for(File f: files) {
			SamplePrimer sp = pc.getSamplePrimer(f);
			if(sp==null) {
				allFound = false;
				System.out.println(f.getAbsolutePath()+ " does not have an entry");
				//System.out.println(f.getAbsolutePath()+"\t"+f.getName());
			}
			else {
				System.out.println(f.getAbsolutePath()+"\t"+sp.toString());
			}
		}
		if(!allFound) {
			//System.exit(0);
		}
		

		
		HashMap<SamplePrimer, TranslocationController> hm = new HashMap<SamplePrimer, TranslocationController>();
		if(!combineFiles) {
			try {
				bw.write(SamplePrimer.getSampleStringHeader()+"\t"+Translocation.getHeader()+"\r\n");
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		
		for(File f: files) {
			SamplePrimer sp = pc.getSamplePrimer(f);
			if(sp==null) {
				continue;
			}
			
			sp.setFile(f);
			System.out.println(f.getName());
			System.out.println(sp.getDNAsample());
			//if(f.getName().contains("LZ45-1-L")) {
			//	System.out.println("contains");
			//if(sp.getRun()!=null && (sp.getRun().contentEquals("104269") || sp.getRun().contains("Exp6"))) {
				
				String primer = sp.getPrimer();
				String chr = sp.getChr();
				File ref = sp.getRef();
				options.setLB(sp.isLB());
				options.setRef(ref);
				options.setPrimer(primer);
				options.setFile(f);
				options.setChr(chr);
				
				//System.out.println(options.printParameters());
		        TranslocationController tc = new TranslocationController(sp);
		        hm.put(sp,tc);
		        tc.testLBRB();
		        //should be -1
		        tc.launchAnalysis(maxTranslocation,maxReadsPerTrans, maxReads);
		        if(!combineFiles) {
		        	hm.get(sp).print(bw, allEvents, minSupport);
					hm.put(sp,null);
		        }
			//System.out.println("hier! " +sp.getRun());
			//}
		}
		//reiterate all samples in hash
		if(combineFiles) {
			for(SamplePrimer ps: hm.keySet()) {
				String DNAsample = ps.getDNAsample();
				TranslocationController tc = hm.get(ps);
				System.out.println("Searching SamplePrimer "+ps.getSample());
				for(SamplePrimer ps2: hm.keySet()) {
					if(ps != ps2) {
						TranslocationController tc2 = hm.get(ps2);
						for(Translocation tl: tc.getTranslocations()) {
							Translocation tl2 = tc2.searchTranslocation(tl);
							if(tl2!=null) {
								if(DNAsample.contentEquals(ps2.getDNAsample())){
									tl.addFoundOtherSampleMatchingDNA(ps2.getSample());
								}
								else {
									tl.addFoundOtherSampleNonMatchingDNA(ps2.getSample());
								}
							}
						}
					}
				}
			}
			
			
			//print
			try {
				bw.write(SamplePrimer.getSampleStringHeader()+"\t"+Translocation.getHeader()+"\r\n");
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			for(SamplePrimer sp: hm.keySet()) {
				hm.get(sp).print(bw, allEvents, minSupport);
				hm.put(sp,null);
			}
		}
		System.out.println("Output written to: "+out.getAbsolutePath());
		try {
			allEvents.close();
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }

	private static Options createOptions() {
		Options o = new Options();
		//threads
		Option i = Option.builder( "i" )
				.longOpt("input")
                .hasArg()
                .argName("FILE")
                .desc("the input sorted and indexed bam file that is used to map the reads" )
                .build();
		o.addOption(i);
		Option r   = Option.builder( "r" )
				.longOpt("reference")
                .hasArg()
                .argName("FILE")
                .desc("the reference file that is used to map the reads" )
                .build();
		o.addOption(r);
		
		Option p = Option.builder("p")
				.longOpt("primer")
				.hasArg()
				.argName("STRING")
				.desc("The internal primer used in the GUIDEseq experiment")
				.build();
		o.addOption(p);
		
		Option c = Option.builder("c")
				.longOpt("chr")
				.hasArg()
				.argName("STRING")
				.desc("The plasmid name that is integrated into the genome")
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
