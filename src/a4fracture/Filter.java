package a4fracture;

import no.uib.cipr.matrix.Vector;

public interface Filter {

    /**
     * removes disallowed parts of v by projection
     * @param v
     */
    public void filter(Vector v);
    
}
