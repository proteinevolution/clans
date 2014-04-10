package clans.misc;

public class MyMath {
	
	/**
	 * Changes the values within {@code matrix} to those of an identity matrix (diagonals = 1, rest = 0).
	 * 
	 * @param matrix
	 *            The matrix to be altered.
	 */
	public static void setTo3x3IdentityMatrix(double[][] matrix) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				matrix[i][j] = (i == j) ? 1 : 0;
			}
		}
	}
	
	public static double[][] create3x3IdentityMatrix() {
		double[][] matrix = new double[3][3];
		setTo3x3IdentityMatrix(matrix);
		return matrix;
	}
}