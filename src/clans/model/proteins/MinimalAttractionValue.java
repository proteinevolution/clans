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
        this.query=-1;
        this.hit=-1;
    }

    public MinimalAttractionValue(int seqLeft, int seqRight) {
        // Attraction is supposed to be symmetrical
        // So just use a cannonical representation
        // This functionality has been moved into the constructor
        if(seqLeft < seqRight) {
            this.query = seqLeft;
            this.hit   = seqRight;
        } else {
            this.query = seqRight;
            this.hit   = seqLeft;
        }

        this.att=att;
    }

    public MinimalAttractionValue(int qnum, int hnum, float att) {
        this.query=qnum;
        this.hit=hnum;
        this.att=att;
    }
    
    public final int query;
    public final int hit;
    public float att=0;

	/**
	 * Only depends on query and hit
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + query;
		result = prime * result + hit;
		return result;
	}

	/**
	 * True if both objects have the same value for query and hit.
	 * The val field is ignored.
	 * This way a MinimalAttractionValue object can be used as key and value
	 * in hash table. This saves memory.
	 * @param obj
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;

		MinimalAttractionValue other = (MinimalAttractionValue) obj;

		return (query == other.query && hit == other.hit);
	}
}
