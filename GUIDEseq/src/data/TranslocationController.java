package data;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import htsjdk.samtools.*;
import static java.util.Comparator.*;

public class TranslocationController {
	private ArrayList<Translocation> trans = new ArrayList<Translocation>();
	public static final int MAXDIST = 500;
	
	public void addTranslocation(SAMRecord s) {
		Translocation nearest = getNearestTranslocation(s);
		if(nearest !=null) {
			nearest.addSam(s);
		}
		else {
			Translocation tl = new Translocation(s);
			//sometimes the sam is not added due to filtering of secondary alignments
			if(tl.getSize()>0) {
				trans.add(tl);
			}
		}
	}

	private Translocation getNearestTranslocation(SAMRecord s) {
		int pos = Translocation.getPosition(s);
		Translocation nearest = null;
		int minDis = Integer.MAX_VALUE;
		for(Translocation tl: trans) {
			//System.out.println(tl.getContigMate()+":"+s.getMateReferenceName());
			if(tl.getContigMate().equals(s.getMateReferenceName())) {
				//getMateNegativeStrandFlag == Forward
				//System.out.println("same Contig");
				//System.out.println(tl.isForward()+":"+!s.getMateNegativeStrandFlag());
				if(tl.isForward() != s.getMateNegativeStrandFlag()) {
					//System.out.println("same orientation");
					int tempPos = tl.getPosition();
					int distance = Math.abs(pos-tempPos);
					//System.out.println("distance: "+distance);
					if(distance<minDis) {
						minDis = distance;
						nearest = tl; 
					}
				}
			}
		}
		//only do this when we have a nearest
		if (minDis<MAXDIST) {
			return nearest;
		}
		return null;
	}

	public void printContents(int minSupport) {
		trans.sort(comparing(Translocation::getSize, Collections.reverseOrder()));
		for(Translocation tl: trans) {
			if(tl.getSize()>minSupport) {
				System.out.println(tl.toString());
			}
		}
	}

	public void addMates(SamReader sr) {
		int around = 500;
		for(Translocation tl: trans) {
			SAMRecordIterator sri = sr.query(tl.getContigMate(), tl.getPosition()-around, tl.getPosition()+around, false);
			while(sri.hasNext()) {
				SAMRecord srec = sri.next();
				//System.out.println(srec.getReadName());
				if(tl.containsRecord(srec.getReadName())) {
					tl.addSam(srec);
				}
				
			}
			sri.close();
		}
	}
}
