package clans;

import java.awt.Color;

/**
 * 
 * @author tancred
 */
public class SequenceGroup {

    String name = "default";
    int[] sequences = new int[0];
    float groupconf = -1;
    String confvals = null;
    float[] seqconf = null;
    java.awt.Color color = Color.red;
    int type = 0;
    int size = 5;
    int[][] polygon = null;
    boolean hide = false;

    void remove(int rmindex) {
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

}
