package data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
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

public class GUIDESeqRun {
	public static void main(String[] args) {
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
		
		try {
			bw.write(SamplePrimer.getSampleStringHeader()+"\t"+Translocation.getHeader()+"\r\n");
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
		SamplePrimer sp = new SamplePrimer(options);
		
		System.out.println(options.printParameters());
		long previousTime = System.currentTimeMillis();
        TranslocationController tc = new TranslocationController(sp);
        tc.testLBRB();
        tc.launchAnalysis(options.getMaxTranslocations(),options.getMaxReadsPerTrans(), options.getMaxReads());
        tc.print(bw, allEvents, (int)options.getMinSupport());
		System.out.println("Output written to: "+out.getAbsolutePath());
		Path path = Paths.get(out.getAbsolutePath());
		try {
			System.out.println("Output size is "+Files.size(path));
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		long currentTime = System.currentTimeMillis();
		double elapsedTime = (currentTime - previousTime) / 1000.0;
		System.out.println("Time in seconds : " + elapsedTime);
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
                .required()
                .argName("FILE")
                .desc("the input sorted and indexed bam file that is used to map the reads" )
                .build();
		o.addOption(i);
		Option r   = Option.builder( "r" )
				.longOpt("reference")
                .hasArg()
                .required()
                .argName("FILE")
                .desc("the reference file that is used to map the reads" )
                .build();
		o.addOption(r);
		
		Option lb   = Option.builder( "LB" )
				.longOpt("LB")
                .desc("is LB or RB" )
                .build();
		o.addOption(lb);
		
		Option rb   = Option.builder( "RB" )
                .desc("is LB or RB" )
                .build();
		o.addOption(rb);
		
		Option p = Option.builder("p")
				.longOpt("primer")
				.hasArg()
				.argName("STRING")
				.required()
				.desc("The internal primer used in the GUIDEseq experiment")
				.build();
		o.addOption(p);
		
		Option c = Option.builder("c")
				.longOpt("chr")
				.hasArg()
				.argName("STRING")
				.required()
				.desc("The plasmid name that is being tested")
				.build();
		o.addOption(c);
		
		Option p5   = Option.builder( "P5" )
                .desc("Is the primer used as a P5 primer during sequencing")
                .optionalArg(true)
                .build();
				o.addOption(p5);
				
		Option p7   = Option.builder( "P7" )
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
		
		Option n = Option.builder("n")
				.longOpt("name")
				.hasArg()
				.argName("STRING")
				.required()
				.desc("The sample name")
				.build();
		o.addOption(n);
		
		Option g = Option.builder("g")
				.longOpt("genotype")
				.hasArg()
				.argName("STRING")
				.required()
				.desc("The sample name")
				.build();
		o.addOption(g);
		
		Option d = Option.builder("d")
				.longOpt("dna")
				.hasArg()
				.argName("STRING")
				.required()
				.desc("The sample name")
				.build();
		o.addOption(d);
		
		Option e = Option.builder("e")
				.longOpt("ecotype")
				.hasArg()
				.argName("STRING")
				.desc("The ecotype")
				.build();
		o.addOption(e);
		
		Option run = Option.builder("u")
				.longOpt("run")
				.hasArg()
				.argName("STRING")
				.desc("The sequence run")
				.build();
		o.addOption(run);
		
		Option debug = Option.builder("debug")
				.desc("run in debug mode")
				.build();
		o.addOption(debug);
		
		
		return o;
	}

}
