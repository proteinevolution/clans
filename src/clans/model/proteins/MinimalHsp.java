package clans.model.proteins;

public class MinimalHsp {

	public final int query;
	public final int hit;
	public double[] val = new double[0];

	/**
	 * 
	 */
	public MinimalHsp() {
		this.query = -1;
		this.hit = -1;
	}

	/**
	 * 
	 * @param i1
	 * @param i2
	 */
	public MinimalHsp(int i1, int i2) {
		this.query = i1;
		this.hit = i2;
	}

	/**
	 * 
	 * @param i1
	 * @param i2
	 * @param pv
	 */
	public MinimalHsp(int i1, int i2, double[] pv) {
		this.query = i1;
		this.hit = i2;
		this.val = pv;
	}

	/**
	 * 
	 * @param i1
	 * @param i2
	 * @param pv
	 */
	public MinimalHsp(int i1, int i2, double pv) {
		this.query = i1;
		this.hit = i2;
		this.val = new double[1];
		val[0] = pv;
	}

	/**
	 * 
	 * @param pv
	 */
	public void addpval(double pv) {
		int length = val.length;
		double[] tmp = new double[length + 1];
		for (int i = 0; i < length; i++) {
			tmp[i] = val[i];
		}// end for i
		tmp[length] = pv;
		val = tmp;
	}

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
	 * This way a MinimalHsp object can be used as key and value
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

		MinimalHsp other = (MinimalHsp) obj;

		return (query == other.query && hit == other.hit);
	}
}
