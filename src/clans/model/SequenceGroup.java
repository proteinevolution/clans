package clans.model;

import java.awt.Color;

/**
 * 
 * @author tancred
 */
public class SequenceGroup {

	public String name = "default";
	public int[] sequences = new int[0];
	public float groupconf = -1;
	public String confvals = null;
	public float[] seqconf = null;
    public java.awt.Color color = Color.red;
    public int type = 0;
    public int size = 5;
    public int[][] polygon = null;
    public boolean hide = false;

    public SequenceGroup() {
    }

    public SequenceGroup(String name, int[] sequences, int size, int type, Color color) {
        this.name = name;
        this.sequences = sequences;
        this.size = size;
        this.type = type;
        this.color = color;
    }

    public void remove(int rmindex) {
        int[] tmp = sequences;
        int seqnum = sequences.length;
        sequences = new int[seqnum - 1];

        for (int i = seqnum; --i > rmindex;) {
            sequences[i - 1] = tmp[i];
        }

        for (int i = rmindex; --i >= 0;) {
            sequences[i] = tmp[i];
        }
    }
    
    /**
     * returns the number of entries in the group
     * @return number of entries in the group 
     */
    public int size() {
        return sequences.length;
    }

}
