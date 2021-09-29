package data;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Log;

import java.io.IOException;
import java.util.HashMap;

/**
 * Buffers entire reference to enable efficient random lookup of sequences
 * @author Daniel Cameron
 *
 */
public class BufferedReferenceSequenceFile implements ReferenceSequenceFile {
	private static final Log log = Log.getInstance(BufferedReferenceSequenceFile.class);
	private final ReferenceSequenceFile underlying;
	/**
	 * Cached contigs
	 */
	private HashMap<String, ReferenceSequence> cache = new HashMap<String, ReferenceSequence>();
	private ReferenceSequence[] referenceIndexLookup;
	public BufferedReferenceSequenceFile(ReferenceSequenceFile underlying) {
		this.underlying = underlying;
		this.referenceIndexLookup = new ReferenceSequence[underlying.getSequenceDictionary().getSequences().size()];
	}
	public byte getBase(int referenceIndex, int position) {
		ReferenceSequence seq = referenceIndexLookup[referenceIndex];
		if (seq == null) {
			synchronized (referenceIndexLookup) {
				seq = addToCache(underlying.getSequenceDictionary().getSequence(referenceIndex).getSequenceName());
			}
		}
		return seq.getBases()[position - 1];
	}
	@Override
	public SAMSequenceDictionary getSequenceDictionary() {
		return underlying.getSequenceDictionary();
	}
	@Override
	public ReferenceSequence nextSequence() {
		return underlying.nextSequence();
	}
	@Override
	public void reset() {
		underlying.reset();
	}
	@Override
	public boolean isIndexed() {
		return underlying.isIndexed();
	}
	/**
	 * Updates the cache to include the new contig
	 * @param contig
	 */
	private synchronized ReferenceSequence addToCache(String contig) {
		ReferenceSequence seq = cache.get(contig);
		if (seq != null) {
			// already populated by another thread while we were waiting to enter
			// this synchronized block
			return seq;
		}
		seq = underlying.getSequence(contig);
		cache.put(contig,seq);
		referenceIndexLookup[underlying.getSequenceDictionary().getSequence(contig).getSequenceIndex()] = seq;
		return seq;
	}
	@Override
	public ReferenceSequence getSequence(String contig) {
		ReferenceSequence seq = cache.get(contig);
		if (seq == null) {
			seq = addToCache(contig);
		}
		return seq;
	}
	@Override
	public ReferenceSequence getSubsequenceAt(String contig, long start, long stop) {
        int length = (int)(stop - start + 1);
		ReferenceSequence fullContig = getSequence(contig);
		if (length > fullContig.length()) {
			throw new IllegalArgumentException("subsequence out of contig bounds");
		}
		if (start > stop + 1) {
			throw new IllegalArgumentException("start after stop");
		}
		byte[] target = new byte[length];
		System.arraycopy(fullContig.getBases(), (int) (start - 1), target, 0, target.length);
		return new ReferenceSequence(fullContig.getName(), fullContig.getContigIndex(), target);
	}
	@Override
	public void close() throws IOException {
		underlying.close();	
	}
}