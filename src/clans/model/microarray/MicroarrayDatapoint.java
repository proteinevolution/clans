/*
 * datapoint.java
 *
 * Created on September 7, 2004, 10:13 AM
 */
package clans.model.microarray;
/**
 *
 * @author  tancred
 */
public class MicroarrayDatapoint {
    
    /** Creates a new instance of datapoint */
    public MicroarrayDatapoint() {
    }
    
    //contain the name of the datapoint, the values (replicates and avg+stdev) 
    //and the corresponding wild-type values (replicates and avg+stdev)
    public String name="";
    public String info="";
    public float value=0;
    public float stdev=0;
    public float presence=0;
    public float[] values=null;
    public float[] datpresences=null;
    public float[] wtvalues=null;
    public float[] wtpresences=null;
    public float wtval=0;
    public float wtstdev=0;
    public float wtpresence=0;
    public boolean use=true;
    
}
