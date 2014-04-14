package clans.model;

import java.util.*;

import clans.model.proteins.AminoAcidSequence;
import clans.model.proteins.MinimalAttractionValue;
import clans.model.proteins.MinimalHsp;

public class SelectedSubsetHandling {

	public static MinimalHsp[] get_blasthits(MinimalHsp[] blasthits,
			int[] selectednames) {
		// convert the old sequence numbering to the new and remove all
		// non-relevant blast hsp's
		Vector<MinimalHsp> tmpvec = new Vector<MinimalHsp>();
		int qnum, hnum;
		HashMap<Integer, Integer> tmphash = new HashMap<Integer, Integer>(
				(int) (selectednames.length / 0.8) + 1, 0.8f);

		for (int i = 0; i < selectednames.length; i++) {
			tmphash.put(selectednames[i], new Integer(i));
		}

		for (int i = 0; i < blasthits.length; i++) {
			qnum = blasthits[i].query;
			hnum = blasthits[i].hit;
			if (tmphash.containsKey(qnum) && tmphash.containsKey(hnum)) {
				tmpvec.addElement(new MinimalHsp(tmphash.get(qnum).intValue(),
						tmphash.get(hnum).intValue(), blasthits[i].val));
			}
		}

		MinimalHsp[] retarr = new MinimalHsp[tmpvec.size()];
		tmpvec.copyInto(retarr);
		return retarr;
	}

	public static MinimalAttractionValue[] get_myattvals(MinimalAttractionValue[] myattvals,
			int[] selectednames) {
		// convert the old sequence numbering to the new and remove all
		// non-relevant attvals
		Vector<MinimalAttractionValue> tmpvec = new Vector<MinimalAttractionValue>();
		String qnum, hnum;
		HashMap<String, Integer> tmphash = new HashMap<String, Integer>(
				(int) (selectednames.length / 0.8) + 1, 0.8f);
		for (int i = 0; i < selectednames.length; i++) {
			tmphash.put(String.valueOf(selectednames[i]), new Integer(i));
		}

		for (int i = 0; i < myattvals.length; i++) {
			qnum = String.valueOf(myattvals[i].query);
			hnum = String.valueOf(myattvals[i].hit);
			if (tmphash.containsKey(qnum) && tmphash.containsKey(hnum)) {
				tmpvec.addElement(new MinimalAttractionValue(((Integer) tmphash.get(qnum))
						.intValue(), ((Integer) tmphash.get(hnum)).intValue(),
						myattvals[i].att));
			}
		}
		MinimalAttractionValue[] retarr = new MinimalAttractionValue[tmpvec.size()];
		tmpvec.copyInto(retarr);
		return retarr;
	}

	public static float[][] get_mymovearr(float[][] mymovearr, int[] selectednames) {
		int elements = selectednames.length;
		float[][] retarr = new float[elements][3];
		for (int i = 0; i < elements; i++) {
			retarr[i] = mymovearr[selectednames[i]];
		}
		return retarr;
	}

	public static float[][] get_myposarr(float[][] myposarr, int[] selectednames) {
		int elements = selectednames.length;
		float[][] retarr = new float[elements][3];
		for (int i = 0; i < elements; i++) {
			retarr[i] = myposarr[selectednames[i]];
		}
		return retarr;
	}

	public static AminoAcidSequence[] get_sequences(AminoAcidSequence[] inaln,
			int[] selectednames) {
		int elements = selectednames.length;
		AminoAcidSequence[] retarr = new AminoAcidSequence[elements];
		for (int i = 0; i < elements; i++) {
			retarr[i] = inaln[selectednames[i]];
		}
		return retarr;
	}

	public static String[] get_namearr(String[] namearr, int[] selectednames) {
		int elements = selectednames.length;
		String[] retarr = new String[elements];
		for (int i = 0; i < elements; i++) {
			retarr[i] = namearr[selectednames[i]];
		}
		return retarr;
	}

	public static float[] get_weights(float[] weights, int[] selectednames) {
		if (weights == null) {
			return null;
		}
		int elements = selectednames.length;
		float[] retarr = new float[elements];
		for (int i = 0; i < elements; i++) {
			retarr[i] = weights[selectednames[i]];
		}// end for i
		return retarr;
	}
}
