/*
 * zoom.java
 *
 * Created on August 20, 2003, 6:11 PM
 */
package clans;
import java.util.*;
/**
 *
 * @author  tancred
 */
public class zoomdata {
    
    /** Creates a new instance of zoom */
    public zoomdata() {
    }
    
    //--------------------------------------------------------------------------
    
    static MinimalHsp[] getblasthitsubset(MinimalHsp[] blasthits,int[]selectednames){
        //convert the old sequence numbering to the new and remove all non-relevant blast hsp's
        int elements=java.lang.reflect.Array.getLength(selectednames);
        Vector tmpvec=new Vector();
        String qnum, hnum;
        HashMap tmphash=new HashMap((int)(elements/0.8)+1,0.8f);
        for(int i=0;i<elements;i++){
            tmphash.put(String.valueOf(selectednames[i]),new Integer(i));
        }//end for i
        elements=java.lang.reflect.Array.getLength(blasthits);
        for(int i=0;i<elements;i++){
            qnum=String.valueOf(blasthits[i].query);
            hnum=String.valueOf(blasthits[i].hit);
            if(tmphash.containsKey(qnum)&&tmphash.containsKey(hnum)){
                tmpvec.addElement(new MinimalHsp(((Integer)tmphash.get(qnum)).intValue(),((Integer)tmphash.get(hnum)).intValue(),blasthits[i].val));
            }
        }//end for i
        MinimalHsp[] retarr=new MinimalHsp[tmpvec.size()];
        tmpvec.copyInto(retarr);
        return retarr;
    }//end getblasthitsubset
    
    //--------------------------------------------------------------------------
    
    static minattvals[] getmyattvalssubset(minattvals[] myattvals, int[] selectednames){
        //convert the old sequence numbering to the new and remove all non-relevant attvals
        int elements=java.lang.reflect.Array.getLength(selectednames);
        Vector tmpvec=new Vector();
        String qnum, hnum;
        HashMap tmphash=new HashMap((int)(elements/0.8)+1,0.8f);
        for(int i=0;i<elements;i++){
            tmphash.put(String.valueOf(selectednames[i]),new Integer(i));
        }//end for i
        elements=java.lang.reflect.Array.getLength(myattvals);
        for(int i=0;i<elements;i++){
            qnum=String.valueOf(myattvals[i].query);
            hnum=String.valueOf(myattvals[i].hit);
            if(tmphash.containsKey(qnum)&&tmphash.containsKey(hnum)){
                tmpvec.addElement(new minattvals(((Integer)tmphash.get(qnum)).intValue(),((Integer)tmphash.get(hnum)).intValue(),myattvals[i].att));
            }
        }//end for i
        minattvals[] retarr=new minattvals[tmpvec.size()];
        tmpvec.copyInto(retarr);
        return retarr;
    }//end getmyattvalssubset
    
    //--------------------------------------------------------------------------
    
    static float[][] getmymovearrsubset(float[][] mymovearr,int[] selectednames){
        int elements=java.lang.reflect.Array.getLength(selectednames);
        float[][] retarr=new float[elements][3];
        for(int i=0;i<elements;i++){
            retarr[i]=mymovearr[selectednames[i]];
        }//end for i
        return retarr;
    }//end getmymovearr
    
    //--------------------------------------------------------------------------
    
    static float[][] getmyposarrsubset(float[][] myposarr,int[] selectednames){
        int elements=java.lang.reflect.Array.getLength(selectednames);
        float[][] retarr=new float[elements][3];
        for(int i=0;i<elements;i++){
            retarr[i]=myposarr[selectednames[i]];
        }//end for i
        return retarr;
    }//end getmyposarrsubset
    
    //--------------------------------------------------------------------------
    
    static AminoAcidSequence[] getinalnsubset(AminoAcidSequence[] inaln,int[] selectednames){
        int elements=java.lang.reflect.Array.getLength(selectednames);
        AminoAcidSequence[] retarr=new AminoAcidSequence[elements];
        for(int i=0;i<elements;i++){
            retarr[i]=inaln[selectednames[i]];
        }//end for i
        return retarr;
    }//end getinalnsubset
    
    //--------------------------------------------------------------------------
    
    static String[] getnamearrsubset(String[] namearr,int[] selectednames){
        int elements=java.lang.reflect.Array.getLength(selectednames);
        String[] retarr=new String[elements];
        for(int i=0;i<elements;i++){
            retarr[i]=namearr[selectednames[i]];
        }//end for i
        return retarr;
    }//end getnamearrsubset
    
    //--------------------------------------------------------------------------
    
    static float[] getweightssubset(float[] weights,int[] selectednames){
        if(weights==null){
            return null;
        }
        int elements=java.lang.reflect.Array.getLength(selectednames);
        float[] retarr=new float[elements];
        for(int i=0;i<elements;i++){
            retarr[i]=weights[selectednames[i]];
        }//end for i
        return retarr;
    }//end getnamearrsubset
    
}//end class
