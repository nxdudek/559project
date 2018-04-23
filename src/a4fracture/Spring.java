package a4fracture;

import javax.vecmath.Vector2d;

/**
 * Spring class for 599 assignment 1
 * @author kry
 */
public class Spring {

    Particle p1 = null;
    Particle p2 = null;
    boolean unique;
        
    /**
     * All springs share the same stiffness coefficient
     */
    public static double k = 1;
    
    /**
     * All springs share the same camping coefficient
     */
    public static double b = 1;
    
    /**
     * Rest length
     */
    double l0 = 0;
    
    /** Empty constructor for doing something special for a mouse spring */ 
    public Spring() {
    	// do nothing
    }
    
    /**
     * Creates a spring connecting two particles.
     * The rest length should be set
     * @param p1
     * @param p2
     */
    public Spring( Particle p1, Particle p2 ) {
        this.p1 = p1;
        this.p2 = p2;
        computeRestLength();
        p1.springs.add(this);
        p2.springs.add(this);
    }
    
    /**
     * Computes the rest length of the connected particles
     */
    public void computeRestLength() {
        l0 = p1.p0.distance( p2.p0 );
    }
    
    /**
     * Sets the rest length of the connected particles with their current positions
     */
    public void setRestLength() {
        l0 = p1.p.distance( p2.p );
    }
    
    /**
     * Applies spring forces to the two particles
     */
    public void apply() {
        Vector2d force = new Vector2d();
        
        force.sub( p2.p, p1.p );
        double l = force.length();
        force.normalize();
        force.scale( (l-l0)*k );
        p1.addForce(force);
        force.scale(-1);
        p2.addForce(force);
        
        force.sub( p2.p, p1.p );
        force.normalize();
        Vector2d v = new Vector2d();
        v.sub(p2.v, p1.v);
        double rv = force.dot(v);
        force.scale( b * rv );
        p1.addForce(force);
        force.scale(-1);
        p2.addForce(force);            
    }
    
}
