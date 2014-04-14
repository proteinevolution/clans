/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package clans;

/**
 *
 * @author tancred
 */
public class HighScoringSegmentPair {
/** Creates a new instance of hsp */
    public HighScoringSegmentPair() {
        qname=new String("");
        qseq=new String("");
        hname=new String("");
        hseq=new String("");
        qstart=-1;
        qend=-1;
        hstart=-1;
        hend=-1;
        value=-1;
        
    }
    String qname;
    String qseq;
    String hname;
    String hseq;
    int qstart;
    int qend;
    int hstart;
    int hend;
    double value;
    
}
