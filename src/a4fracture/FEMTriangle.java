package a4fracture;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import no.uib.cipr.matrix.DenseMatrix;

import javax.vecmath.Vector2d;
import java.text.DecimalFormat;

/**
 * FEM triangle spring... a number of MTJ matrices are created for your convenience, 
 * but note that you may ultimately note that rolling your own 2D matrix class would
 * be much more efficient here.  Note that there are nice techniqeus for quickly 
 * computing SVD or EVD of a 2x2 matrix.
 * 
 * @author kry
 */
public class FEMTriangle {

	/**
	 * Particles in the triangle.  Note that order matters!
	 * (clockwise order give positive area)
	 */
	Particle A, B, C;
	
	/** index of this tri in each particle's triangle list */
	int Ai, Bi, Ci;
	
	/** Lam� parameter */
	static double mu;		
	
	/** Lam� parameter */
	static double lambda;   
	
	/** Area of the material space triangle */
	double area;
	
	/** used for computing lumped mass at particles */
	double density;
	
	/** deformed shape matrix */
	Matrix2d Ds = new Matrix2d();
	
	/** material space shape matrix */
	Matrix2d Dm = new Matrix2d();

	Matrix2d Dminv = new Matrix2d();
	Matrix2d DminvT = new Matrix2d();
	
	/** deformation gradient */
	Matrix2d F = new Matrix2d();
	
	/** Strain tensor */
	Matrix2d E = new Matrix2d();
	
	/** Piola stress tensor */
	Matrix2d P = new Matrix2d();
	
	/** Forces */
	Matrix2d H = new Matrix2d();
	
	/** rotation matrix for SVD co-rotated elasticity */
	Matrix2d R = new Matrix2d();

	Matrix2d stress = new Matrix2d();

	/** Eigenvalue decomposition of the stress tensor */
	Matrix2d evd;

	Matrix2d tensileStress = new Matrix2d();

	Matrix2d compressiveStress = new Matrix2d();

	public FEMTriangle( Particle A, Particle B, Particle C ) {
		this.A = A;
		this.B = B;
		this.C = C;
		Ai = A.addTriangle( this );
		Bi = B.addTriangle( this );
		Ci = C.addTriangle( this );

		// TODO: Objective 1:  Set up the material shape matrix and other quantities you need

		Dm.a = A.p.x - C.p.x;
		Dm.b = B.p.x - C.p.x;
		Dm.c = A.p.y - C.p.y;
		Dm.d = B.p.y - C.p.y;

		Dminv.inverse(Dm);
		DminvT.transpose(Dminv);

		double a = A.p.distance(B.p);
		double b = A.p.distance(C.p);
		double c = B.p.distance(C.p);
		double s = (a + b + c) / 2;
		area = Math.sqrt(s * (s - a) * (s - b) * (s - c));

	}
	
	/**
	 * Reset to clear MTJ matrices, incase they get bad values
	 */
	void reset() {
		Ds.zero();
		F.zero();
		E.zero();
		P.zero();
		H.zero();
	}

	public void addSprings(ParticleSystem system)
	{
		Spring s1 = new Spring(A, B);
		Spring s2 = new Spring(A, C);
		Spring s3 = new Spring(B, C);

		system.springs.add(s1);
		system.springs.add(s2);
		system.springs.add(s3);
	}
	
	/**
	 * Helper function in case you want to print out any MTJ matrices!
	 * @param M
	 */
	public void printMatrix( DenseMatrix M ) {
		final DecimalFormat df = new DecimalFormat("0.00");
		for ( int row = 0; row < M.numRows(); row ++ ) {
			for ( int col = 0; col < M.numColumns(); col++ ) {
				System.out.print( df.format(M.get(row, col)) + " " );
			}
			System.out.println();
		}
	}

	public boolean sharesTwo(FEMTriangle t)
	{
		int i = 0;
		if(A.equals(t.A) || A.equals(t.B) || A.equals(t.C))
			++i;
		if(B.equals(t.A) || B.equals(t.B) || B.equals(t.C))
			++i;
		if(C.equals(t.A) || C.equals(t.B) || C.equals(t.C))
			++i;

		if(i == 2)
			return true;
		else
			return false;
	}

	public Vector2d getCenter()
	{
		Vector2d c = new Vector2d();
		c.x = (A.p.x + B.p.x + C.p.x) / 3;
		c.y = (A.p.y + B.p.y + C.p.y) / 3;
		return c;
	}
	
	public void applyForce() {

		// TODO: Objective 1: compute the forces
		// add the contribution to the particle f vector
		Ds.a = A.p.x - C.p.x;
		Ds.b = B.p.x - C.p.x;
		Ds.c = A.p.y - C.p.y;
		Ds.d = B.p.y - C.p.y;

		F.mult(Ds, Dminv);

		Matrix2d U = new Matrix2d();
		Matrix2d V = new Matrix2d();
		Matrix2d Fhat = new Matrix2d();

		F.irving2004(U, Fhat, V);

		Matrix2d Vt = new Matrix2d();
		Vt.transpose(V);

		F.R = new Matrix2d();
		F.S = new Matrix2d();

		F.R.mult(U, Vt);
		F.S.mult(V, Fhat);
		F.S.mult(Vt);

		// 2μ(F−R) +λtr(RTF−I)R

		Matrix2d temp = new Matrix2d(F);
		Matrix2d temp2 = new Matrix2d();

		temp.sub(F.R);
		temp.scale(2 * mu);

		temp2.transAmult(F.R, F);
		temp2.a -= 1;
		temp2.d -= 1;

		P.a = F.R.a;
		P.b = F.R.b;
		P.c = F.R.c;
		P.d = F.R.d;

		P.scale(lambda * temp2.trace());
		P.add(temp);

		// final calculation for H
		H.mult(P, DminvT);
		H.scale(-area);

		A.f.x += H.a;
		A.f.y += H.c;

		B.f.x += H.b;
		B.f.y += H.d;

		C.f.x -= (H.a + H.b);
		C.f.y -= (H.c + H.d);

		// TODO: Objective 2: visualize eignevalues of the stress tensor
		// NOTE: this is partiall done for you here, with code 
		// below to draw... but you'll may want to look at the Particle
		// class for ideas on cleaning this up!

		stress.mult(F.S, F.S);
		stress.a -= 1;
		stress.d -= 1;
		stress.scale(.5);

		try {
			// this one frezes on NaN instad of dying
			if ( !Double.isFinite( stress.a) ||
				 !Double.isFinite( stress.b) ||
				 !Double.isFinite( stress.c) ||
				 !Double.isFinite( stress.d) ) {
				evd = null;
				return ;
			}
			evd = new Matrix2d(stress);
			evd.evd();
		} catch ( Exception e ) {
			e.printStackTrace();
		}
		
		// TODO: Objective 3: compute compressive and tensile stress and vertex contribs
		// see fplus and fminus lists in the particle

		compressiveStress.zero();
		compressiveStress.rank1(Math.max(0, evd.ev1), evd.v1);
		compressiveStress.rank1(Math.max(0, evd.ev2), evd.v2);

		tensileStress.zero();
		tensileStress.rank1(Math.min(0, evd.ev1), evd.v1);
		tensileStress.rank1(Math.min(0, evd.ev2), evd.v2);

		H.mult(compressiveStress, DminvT);
		H.scale(-area);
		A.fplus.add(new Vector2d(H.a, H.c));
		B.fplus.add(new Vector2d(H.b, H.d));
		C.fplus.add(new Vector2d(-H.a - H.b, -H.c - H.d));

		H.mult(tensileStress, DminvT);
		H.scale(-area);
		A.fminus.add(new Vector2d(H.a, H.c));
		B.fminus.add(new Vector2d(H.b, H.d));
		C.fminus.add(new Vector2d(-H.a - H.b, -H.c - H.d));
	}
	
	/**
	 * Use position differential set in particles, Particle.dx, to compute 
	 * and accumulate force differentials, Particle.df
	 */
	public void computeForceDifferentials( ) {
		
		// TODO: Optional: for raleigh beta damping and implicit integration
		
	}
	
	public void display( GLAutoDrawable drawable, double s) {
		GL2 gl = drawable.getGL().getGL2();
		Vector2d p = new Vector2d();
		p.add( A.p );
		p.add( B.p );
		p.add( C.p );
		p.scale(1.0/3.0);
		
		// TODO: Objective 2: draw the stress tensor eigenvectors
		// (see also code in Particle wrt separation tensor)

		double ev0 = 0;  // eigenvalue
		double ev0x = 0; // eigenvector x component
		double ev0y = 0; // eigenvector y component
		double ev1 = 0;  // eigenvalue
		double ev1x = 0; // eigenvector x component
		double ev1y = 0;

		if(evd.v1 != null && evd.v2 != null) {
			ev0 = evd.ev1;  // eigenvalue
			ev0x = evd.v1.x; // eigenvector x component
			ev0y = evd.v1.y; // eigenvector y component
			ev1 = evd.ev2;  // eigenvalue
			ev1x = evd.v2.x; // eigenvector x component
			ev1y = evd.v2.y; // eigenvector y component
		}

		double s0 = s * ev0;
		double s1 = s * ev1;
				
		gl.glBegin( GL.GL_LINES );
		if ( ev0 < 0 ) {
			gl.glColor3f(0.75f,0,0);
		} else {
			gl.glColor3f(0,0.75f,0);
		}
		gl.glVertex2d( p.x, p.y);
		gl.glVertex2d( p.x + ev0x * s0, p.y + ev0y * s0 );
		if ( ev1 < 0 ) {
			gl.glColor3f(0.75f,0,0);
		} else {
			gl.glColor3f(0,0.75f,0);
		}
		gl.glVertex2d( p.x, p.y);
		gl.glVertex2d( p.x + ev1x * s1, p.y + ev1y * s1 );
		gl.glEnd();
	}
	
}
