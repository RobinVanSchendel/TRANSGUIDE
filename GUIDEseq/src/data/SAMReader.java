package data;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class SAMReader {
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
		
		System.out.println(options.printParameters());
        TranslocationController tc = new TranslocationController(options);
        tc.launchAnalysis();
    }

	private static Options createOptions() {
		Options o = new Options();
		//threads
		Option i   = Option.builder( "i" )
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
                .required(true)
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
