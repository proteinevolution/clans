package clans.model.proteins;

public class MinimalHsp {

	public int query = -1;
	public int hit = -1;
	public double[] val = new double[0];

	/**
	 * 
	 */
	public MinimalHsp() {
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
}
