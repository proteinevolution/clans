/*
 * seqgroup.java
 *
 * Created on March 31, 2004, 2:31 PM
 */
package clans;
/**
 *
 * @author  tancred
 */
public class SequenceGroup {
    
    /** Creates a new instance of seqgroup */
    public SequenceGroup() {
    }
    
    String name="default";
    int[] sequences=new int[0];
    float groupconf=-1;
    String confvals=null;
    float[] seqconf=null;
    java.awt.Color color=java.awt.Color.red;
    int type=0;
    int size=5;
    int[][] polygon=null;
    boolean hide=false;

    void remove(int rmindex){
        int[] tmp=sequences;
        int seqnum=java.lang.reflect.Array.getLength(sequences);
        sequences=new int[seqnum-1];
        for(int i=seqnum;--i>rmindex;){
            sequences[i-1]=tmp[i];
        }//end for i
        for(int i=rmindex;--i>=0;){
            sequences[i]=tmp[i];
        }//end for i
    }//end remove

}
