package a4fracture;

import no.uib.cipr.matrix.Vector;

public interface MatrixMult {

	/**  
	 * Computes product Av = thisMatrix * v 
	 * @param v
	 * @param Av result
	 */
    public void mult(Vector v, Vector Av);

}
