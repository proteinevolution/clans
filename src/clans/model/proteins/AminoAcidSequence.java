package clans.model.proteins;

/**
 * Represents a single protein sequence with name and one-letter code sequence.
 */
public class AminoAcidSequence {

	public String name;
	public String seq;

	public AminoAcidSequence() {
		name = "";
		seq = "";
	}

	public AminoAcidSequence(String name, String seq) {
		this.name = name;
		this.seq = seq;
	}
	
	/**
	 * Returns the length of the sequence.
	 */
	public int length() {
		return seq.length();
	}
}
