/*
 * saverunobject.java
 *
 * Created on September 19, 2003, 4:59 PM
 */
package clans;
/**
 *
 * @author  tancred
 */
public class saverunobject {
    
    /** Creates a new instance of saverunobject */
    public saverunobject() {
    }
    
    java.io.File file=null;
    minhsp[] blasthits=null;
    aaseq[] inaln=null;
    minattvals[] attvals=null;
    float[][] posarr=null;
    float[] weights=null;
    int ovalsize=4;
    int dotsize=2;
    int groupsize=4;
    double minpval=1;
    double cooling=1;
    double currcool=1;
    float attfactor=10;
    int attvalpow=1;
    float repfactor=10;
    int repvalpow=1;
    float dampening=0.2f;
    float zoom=1;
    double minattract=1;
    float maxmove=0.1f;
    double pval=-1;
    double[][] rotmtx={{1,0,0},{0,1,0},{0,0,1}};//default
    boolean cluster2d=false;
    boolean showinfo=true;
    boolean usescval=false;
    boolean complexatt=true;
    java.util.Vector seqgroupsvec=new java.util.Vector();
    String blastpath="blastall -p blastp";
    String formatdbpath="formatdb";
    java.util.ArrayList mapfiles=new java.util.ArrayList();
    java.util.ArrayList lookupfiles=new java.util.ArrayList();
    java.util.Vector affyfiles=null;
    java.awt.Color[] colorarr=null;
    float[] colorcutoffs=null;
    boolean usefoldchange=false;
    boolean avgfoldchange=false;
    String namesdmp_file=null;
    String nodesdmp_file=null;
    int rounds=0;
    
}//end class
