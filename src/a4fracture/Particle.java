package a4fracture;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import no.uib.cipr.matrix.DenseMatrix;

import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;
import java.util.ArrayList;

/**
 * Particle class for 599 assignment 1
 * @author kry
 */
public class Particle {
    
    /** true means that the particle can not move */
	public boolean pinned = false;
        
    /** The mass of the particle */
	public double mass = 1;
    
    /** current position of the particle */
    public Point2d p = new Point2d();
    
    /** current velocity of the particle */
    public Vector2d v = new Vector2d();
    
    /** initial position of the particle */
    public Point2d p0 = new Point2d();
    
    /** initial velocity of the particle */
    public Vector2d v0 = new Vector2d();
    
    /** force acting on this particle */
    public Vector2d f = new Vector2d();
    
    /** position differential for implicit integration*/
    public Vector2d dx = new Vector2d();
    
    /** force differential for implicit integration */
    public Vector2d df = new Vector2d();
    
    /**
     * A list of springs to which this particle is attached
     * (used to recompute rest lengths when particles are moved)
     */
    public ArrayList<Spring> springs = new ArrayList<Spring>();
    
    /**
     * A list of adjacent triangles, necessary for processing fracture
     */
    ArrayList<FEMTriangle> tris = new ArrayList<FEMTriangle>();
    
    ArrayList<Vector2d> fplus = new ArrayList<Vector2d>();
    
    ArrayList<Vector2d> fminus = new ArrayList<Vector2d>();

    Vector2d fMinusSum = new Vector2d();

    Vector2d fPlusSum = new Vector2d();

    /** The separtaion tensor for deciding fracture points and directions */
    Matrix2d separationTensor = new Matrix2d();

    /** Eigenvalue decomposition of the separation tensor */
    Matrix2d evd;

    public int addTriangle( FEMTriangle tri ) {
    	int i = tris.size();
    	tris.add( tri );    	
    	//fplus.add( new DenseVector(2) );
    	//fminus.add( new DenseVector(2) );
    	return i;
    }

    public void resetTensors()
    {
        separationTensor.zero();
        Matrix2d posSum = new Matrix2d();
        Matrix2d negSum = new Matrix2d();
    }

    
    public void computeSeparationTensor() {
    	
    	// TODO: Objective 3: compute the separtaion tensor
    	// note thet fMiuns and fPlus should be set by the FEMTriangle
    	// you'll want to compute the sums and assmeble the matrix from
    	// outer products.  Note MTJ matrices have a rank1 method for adding
    	// an outer prodcut to a matrix. 

        separationTensor.zero();
        Matrix2d posSum = new Matrix2d();
        Matrix2d negSum = new Matrix2d();

        for (Vector2d v : fplus) {
            fPlusSum.add(v);
            // sums m(f) for all f in f+
            posSum.rank1(1, v);
        }

        for (Vector2d v : fminus) {
            fMinusSum.add(v);
            // sums m(f) for all f in f- and negates them
            negSum.rank1(-1, v);
        }

        separationTensor.rank1(-1, fPlusSum);
        separationTensor.rank1(1, fMinusSum);
        separationTensor.add(posSum);
        separationTensor.add(negSum);
        separationTensor.scale(.5);

        // reset everything to prepare for the next iteration
        fplus.clear();
        fminus.clear();
        fPlusSum.x = 0;
        fPlusSum.y = 0;
        fMinusSum.x = 0;
        fMinusSum.y = 0;

        // here is osme code to compte the EVD of the separation tensor
    	try {
			if ( !Double.isFinite(separationTensor.a) ||
				 !Double.isFinite(separationTensor.b) ||
			     !Double.isFinite(separationTensor.c) ||
			     !Double.isFinite(separationTensor.d) ) {
				evd = null;
				return;
			}
			evd = new Matrix2d(separationTensor);
			evd.evd();
		} catch ( Exception e ) {
            e.printStackTrace();
        }
    }
    
    public void drawSeparationTensor( GLAutoDrawable drawable, double s ) {
    	if ( evd == null ) return;
    	GL2 gl = drawable.getGL().getGL2();
    	gl.glBegin( GL.GL_LINES );
    	DenseMatrix evecs = new DenseMatrix(2, 2);
        evecs.set(0, 0, evd.v1.x);
        evecs.set(1, 0, evd.v1.y);
        evecs.set(0, 1, evd.v2.x);
        evecs.set(1, 1, evd.v2.y);
		double[] evalreal = new double[2]; // = evd.getRealEigenvalues();
        evalreal[0] = evd.ev1;
        evalreal[1] = evd.ev2;
		double s1 = s * evalreal[0];
		double s2 = s * evalreal[1];
		if ( evalreal[0] < 0 ) {
			gl.glColor3f(0.75f,0,0);
		} else {
			gl.glColor3f(0,0.75f,0);
		}
		gl.glVertex2d( p.x, p.y);
		gl.glVertex2d( p.x + evecs.get(0,0) * s1, p.y + evecs.get(1,0) * s1 );
		if ( evalreal[1] < 0 ) {
			gl.glColor3f(0.75f,0,0);
		} else {
			gl.glColor3f(0,0.75f,0);
		}
		gl.glVertex2d( p.x, p.y);
		gl.glVertex2d( p.x + evecs.get(0,1) * s2, p.y + evecs.get(1,1) *s2 );

    	gl.glEnd();
    }
    
    /**
     * Creates a partilce copy, useful for fracture events
     * @param
     */
    public Particle( Particle A ) {
    	p.set( A.p );
    	p0.set( A.p0 );
    	v.set( A.v );
    	v0.set( A.v0 );
    	// NOTE: Mass will need to be recomputed!
    }
    
    /**
     * Creates a particle with the given position and velocity
     * @param x
     * @param y
     * @param vx
     * @param vy
     */
    public Particle( double x, double y, double vx, double vy ) {
        p0.set(x,y);
        v0.set(vx,vy);
        reset();
    }
    
    /**
     * Resets the position of this particle
     */
    public void reset() {
        p.set(p0);
        v.set(v0);
        f.set(0,0);
    }
    
    /**
     * Adds the given force to this particle.
     * Note that you probably want to set the force to zero 
     * before accumulating forces. 
     * @param force
     */
    public void addForce( Vector2d force ) {
        f.add(force);
    }
    
    /**
     * Computes the distance of a point to this particle
     * @param x
     * @param y
     * @return the distance
     */
    public double distance( double x, double y ) {
        Point2d tmp = new Point2d( x, y );
        return tmp.distance(p);
    }
   
}
