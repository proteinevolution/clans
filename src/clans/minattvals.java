/*
 * minattvals.java
 *
 * Created on September 2, 2005, 2:18 PM
 */
package clans;
/**
 *
 * @author  tancred
 */
public class minattvals {
    
    /** Creates a new instance of minattvals */
    public minattvals() {
    }
    
    public minattvals(int qnum, int hnum, float att){
        this.query=qnum;
        this.hit=hnum;
        this.att=att;
    }
    
    int query=-1;
    int hit=-1;
    float att=0;
    
}
