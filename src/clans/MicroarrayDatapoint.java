/*
 * datapoint.java
 *
 * Created on September 7, 2004, 10:13 AM
 */
package clans;
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
    String name="";
    String info="";
    float value=0;
    float stdev=0;
    float presence=0;
    float[] values=null;
    float[] datpresences=null;
    float[] wtvalues=null;
    float[] wtpresences=null;
    float wtval=0;
    float wtstdev=0;
    float wtpresence=0;
    boolean use=true;
    
}
