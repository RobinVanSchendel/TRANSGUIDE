package data;

import java.util.ArrayList;

import org.jcvi.jillion.core.qual.QualitySequence;
import org.jcvi.jillion.core.qual.QualitySequenceBuilder;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.jillion.trace.fastq.FastqRecord;

public class Reads {
	private ArrayList<FastqRecord> R1s;
	private ArrayList<FastqRecord> R2s;
	public Reads() {
		R1s = new ArrayList<FastqRecord>();
		R2s = new ArrayList<FastqRecord>();
	}
	/** Adding max 100 reads, otherwise we can not fit it in memory
	 * 
	 */
	public boolean addReads(FastqRecord r1, FastqRecord r2) {
		if(R1s.size()<=100) {
			R1s.add(r1);
			R2s.add(r2);
			return true;
		}
		return false;
	}
	public FastqRecord consolidateR2() {
		FastqRecord cons = FastqRecord(R2s);
		return cons;
	}
	public FastqRecord consolidateR1() {
		FastqRecord cons = FastqRecord(R1s);
		return cons;
	}
	private static FastqRecord FastqRecord(ArrayList<FastqRecord> al) {
		if(al.size()>1) {
			Consensus c = new Consensus();
			for(FastqRecord f: al) {
				c.add(f.getNucleotideSequence().toString());
			}
			//get maxQuality per read
			byte[] quals = new byte[(int) al.get(0).getLength()];
			//get the consensus read
			String read = c.getConsensusStringReplaceByN(0.8);
			
			for(int i=0;i<al.get(0).getLength();i++) {
				char currentChar = read.charAt(i);
				byte max = 0;
				for(FastqRecord f: al) {
					char readChar = f.getNucleotideSequence().toString().charAt(i);
					//only get the base quality for the nucleotide that matches the consensus
					//N will very likely result in quality of 0 (=! in phred code)
					if(currentChar == readChar) {
						byte value = f.getQualitySequence().get(i).getQualityScore();
						if(value>max) {
							max = value;
						}
					}
				}
				quals[i] = max;
			}
			
			NucleotideSequence ns = new NucleotideSequenceBuilder().append(read).build();
			QualitySequence qs = new QualitySequenceBuilder().append(quals).build();
			FastqRecord fqR1 = al.get(0).toBuilder()
					.qualities(qs)
					.basecalls(ns)
					.build();
			
			return fqR1;
		}
		else {
			return al.get(0);
		}
	}
}
