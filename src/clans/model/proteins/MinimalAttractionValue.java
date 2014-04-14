/*
 * minattvals.java
 *
 * Created on September 2, 2005, 2:18 PM
 */
package clans.model.proteins;
/**
 *
 * @author  tancred
 */
public class MinimalAttractionValue {
    
    /** Creates a new instance of minattvals */
    public MinimalAttractionValue() {
    }
    
    public MinimalAttractionValue(int qnum, int hnum, float att){
        this.query=qnum;
        this.hit=hnum;
        this.att=att;
    }
    
    public int query=-1;
    public int hit=-1;
    public float att=0;
    
}
