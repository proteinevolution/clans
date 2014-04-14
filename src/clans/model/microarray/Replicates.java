/*
 * replicates.java
 *
 * Created on September 7, 2004, 10:00 AM
 */
package clans.model.microarray;
/**
 *
 * @author  tancred
 */
public class Replicates {
    
    /** Creates a new instance of replicates */
    public Replicates() {
    }
    public String abbreviation=null;
    public String name="no name";
    public String wtname="no wtname";
    public java.io.File[] replicate=null;
    public java.io.File[] wtreplicate=null;
    public int replicates=0;
    public int wtreplicates=0;
}
