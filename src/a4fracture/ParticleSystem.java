package a4fracture;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.CollapsiblePanel;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.SceneGraphNode;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.vecmath.Vector2d;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Implementation of a simple particle system.
 * 
 * Note that the particle system implements Function, that is, it evaluates
 * its derivatives to return to the step method which is called by implementations
 * of the Integrator interface.
 * 
 * Note also that it is actually the updateParticles method in this class which 
 * should be calling the Integrator step method! 
 * 
 * @author kry
 */
public class ParticleSystem implements SceneGraphNode, Filter, MatrixMult {

    RobustCCD robustCCD = new RobustCCD();

    /** the particle list */
    public List<Particle> particles = new LinkedList<Particle>();
    
    /** the spring list for collisions (unused in this assignment) */
    public List<Spring> collidableEdges = new LinkedList<Spring>();
    
    /** the spring list, (empty unless you're doing something out of the scope of the assignment */
    public List<Spring> springs = new LinkedList<Spring>();
    
    /** leaf springs connect 3 particles with a hinge like spring */
    public List<FEMTriangle> femSprings = new LinkedList<FEMTriangle>();
      
    /** The name of the particle system. */
    public String name = "";
    
    Spring mouseSpring;
    boolean useMouseSpring = false;
    
    /**
     * Creates an empty particle system
     */
    public ParticleSystem() {
    	mouseSpring = new Spring();
    	mouseSpring.p2 = new Particle(0,0,0,0); // a special particle!
        // creates an empty system!
    }
    
    /**
     * Resets the positions of all particles to their initial states
     */
    public void resetParticles() {
        for ( Particle p : particles ) {
            p.reset();
        }
        for ( FEMTriangle tri : femSprings ) {
        	tri.reset(); // try to get rid of NaNs
        }
        time = 0;
    }
    
    /**
     * Deletes all particles, and as such removes all springs too.
     */
    public void clear() {        
        particles.clear();
        springs.clear();
        collidableEdges.clear();
        femSprings.clear();
        name = "";
    }

    public void addSpringsFromTriangles()
    {
        springs.clear();
        for(FEMTriangle t : femSprings)
        {
            t.addSprings(this);
        }
    }

    public void recomputeBorders()
    {
        for(Spring s : springs)
        {
            s.unique = true;
        }

        collidableEdges.clear();
        for(Spring s1 : springs)
        {
            for(Spring s2 : springs)
            {
                if(!s1.equals(s2))
                {
                    if(s1.p1.equals(s2.p1) && s1.p2.equals(s2.p2) || (s1.p1.equals(s2.p2) && s1.p2.equals(s2.p1))) {
                        s1.unique = false;
                        s2.unique = false;
                    }
                }
            }
        }

        for(Spring s : springs)
        {
            if(s.unique)
                collidableEdges.add(new Spring(s.p1, s.p2));
        }
    }

    public void newTriBorders(Particle p)
    {
        for(FEMTriangle f  : p.tris)
        {
            if(p.equals(f.A))
            {
                collidableEdges.add(new Spring(p, f.B));
                collidableEdges.add(new Spring(p, f.C));
            }
            else if(p.equals(f.B))
            {
                collidableEdges.add(new Spring(p, f.A));
                collidableEdges.add(new Spring(p, f.C));
            }
            else
            {
                collidableEdges.add(new Spring(p, f.A));
                collidableEdges.add(new Spring(p, f.B));
            }
        }
    }

    public double time = 0;
        
    /**
     * Advances time and updates the position of all particles
     * @param h 
     * @return true if update was successful
     */
    public void updateParticles( double h ) {
        Vector2d tmp = new Vector2d();

        // set the global spring properties
        Spring.k = k.getValue();
        Spring.b = b.getValue();
        double E = YoungModulus.getValue();
        double nu = PoissonRatio.getValue();
        FEMTriangle.mu = E/2/(1+nu);
        FEMTriangle.lambda = E*nu/(1+nu)/(1-2*nu);

        // collision detection
        robustCCD.restitution = restitution.getValue();
        robustCCD.H = H.getValue();

     /*   if ( repulsion.getValue() ) {
            robustCCD.applyRepulsion( h, this, H.getValue() );
        }*/

        if ( ! robustCCD.check( h, this ) ) {
            System.out.println("collision failed?");
        }


        // Update the velocity of the particles as per symplectic Euler
        if ( ! implicit.getValue() ) {

            computeForces();

	        for ( Particle p : particles ) {
	            if(p.mass == 0) p.mass = 1;

	            if ( p.pinned ) {            
	                p.f.set(0,0); // just to make sure!
	                p.v.set(0,0);
	            } else {
	                tmp.scale( h / p.mass, p.f );
	                p.v.add( tmp );            
	            }
	        }
        } else {
	        // TODO: Optional: implicit integration for stiff systems
        }

        // Finally, update the positions using the velocity 
        for ( Particle p : particles ) {
            if ( p.pinned ) continue;
            // symplectic Euler 
            tmp.scale( h, p.v );
            p.p.add( tmp );
            p.f.set(0,0);
        }

        time = time + h;
        return;
    }
    
    private void computeForces() {
        Vector2d tmp = new Vector2d();
    	double damping  = c.getValue();
        for ( Particle p : particles ) {
            p.f.set( 0, useg.getValue() ? g.getValue() * p.mass : 0 );
            tmp.scale( -damping, p.v );
            p.f.add( tmp );                        
        }
        if ( useMouseSpring ) {
        	Vector2d force = new Vector2d();            
            force.sub( mouseSpring.p2.p, mouseSpring.p1.p );            
            force.scale( Spring.k );
            mouseSpring.p1.addForce(force);
            force.scale(-1);
            mouseSpring.p2.addForce(force);            
        }
        for ( Spring s : springs ) { // shouldn't be any
            s.apply();
        }

        for ( FEMTriangle f : femSprings ) {
        	f.applyForce();
        }

        double alpha = RaleighAlpha.getValue();

        for ( Particle p : particles ) {
        	if ( useRAlpha.getValue() ) {
        		p.f.x += - alpha * p.mass * p.v.x;
        		p.f.y += - alpha * p.mass * p.v.y;
        	}
        }

        // TODO: optional: consider implementing Raleigh beta damping 

        // compute the separation tensors at each node
        for ( Particle p : particles ) {
        	if ( p.tris.size() == 0 ) continue;
        	p.computeSeparationTensor();
        }

        // TODO: Objective 4: process fracture events
        // search for the non-pinned vertex that has a
        // separation tensor eigenvalue that exceeds
        // the toughness threshold the most, and only
        // only process one fracture event per time step

        double maxEV = 0;
        // just initialze with any random particle. shouldnt matter since maxEV is still 0?
        Particle frac = new Particle(0, 0, 0, 0);

        // find largest eigenvalue of separation tensors of all particles and remember which particle
        for(Particle p : particles)
        {
            if(p.pinned) continue;

            // particle must be attached to at least 2 triangles to fracture
            if(p.tris.size() > 1) {
                if (p.evd != null) {
                    if (p.evd.ev1 > maxEV) {
                        maxEV = p.evd.ev1;
                        frac = p;
                    }
                    if (p.evd.ev2 > maxEV) {
                        maxEV = p.evd.ev2;
                        frac = p;
                    }
                }
                else
                    System.out.println("null evd for particle");
            }
        }

        // if it exceeds the material toughness
        if(maxEV > toughness.getFloatValue() && frac.tris.size() > 0)
        {
            List<FEMTriangle> top = new ArrayList<>();
            List<FEMTriangle> bot = new ArrayList<>();

            Vector2d fracLine;
            if(frac.evd.ev1 > frac.evd.ev2)
                fracLine = frac.evd.v1;
            else
                fracLine = frac.evd.v2;

            for(FEMTriangle tri : frac.tris)
            {
                Vector2d c = tri.getCenter();
                Vector2d lineToCenter = new Vector2d();

                lineToCenter.x = c.x - frac.p.x;
                lineToCenter.y = c.y - frac.p.y;

                if(fracLine.dot(lineToCenter) > 0)
                    top.add(tri);
                else
                    bot.add(tri);
            }

            Particle dup = null;
            if(bot.size() > 0 && top.size() > 0) {
                dup = createParticle(frac.p.x, frac.p.y, frac.v.x, frac.v.y);
                dup.resetTensors();
                for (FEMTriangle tri : bot) {
                    dup.addTriangle(tri);
                    frac.tris.remove(tri);

                    if(tri.A.equals(frac))
                        tri.A = dup;
                    else if(tri.B.equals(frac))
                        tri.B = dup;
                    else if(tri.C.equals(frac))
                        tri.C = dup;
                    else
                        System.out.println("bad");
                }
            }

            List<Spring> toRemove = new ArrayList<Spring>();
            for(Spring s : collidableEdges)
            {
                if(s.p1.equals(frac) || s.p2.equals(frac))
                    toRemove.add(s);
            }

            for(Spring s : toRemove)
            {
                collidableEdges.remove(s);
            }

            // inelegant code to remove hinge joints
            // can be commented out by adding a /* here, which goes until marked below

            if(dup != null)
            {
                // add all the particles in triangles that are affected by the most recent fracture
                List<Particle> adjacent = new ArrayList<>();
                for(FEMTriangle t : frac.tris)
                {
                    adjacent.add(t.A);
                    adjacent.add(t.B);
                    adjacent.add(t.C);
                }
                for(FEMTriangle t : dup.tris)
                {
                    adjacent.add(t.A);
                    adjacent.add(t.B);
                    adjacent.add(t.C);
                }

                // for each of these particles, check if its a hinge
                for(Particle p : adjacent) {
                    List<FEMTriangle> uncon = new ArrayList<>();
                    List<FEMTriangle> con = new ArrayList<>();
                    boolean modified = true;

                    // add all the adjacent triangles to this particle to an unconnected triangle list
                    // except for one, which is the 'connected' triangle to start with
                    uncon.addAll(p.tris);
                    con.add(uncon.get(0));
                    uncon.remove(0);

                    // if something was changed and there are still unconnected triangles, keep iterating over the list
                    while (modified && uncon.size() > 0) {
                        modified = false;
                        FEMTriangle toMove = null;
                        // check each of the unconnected triangles against each of hte connected ones
                        // i know this is n^2 but since there are like < 5 triangles per particle, it works ok
                        for (FEMTriangle t : uncon) {
                            for (FEMTriangle q : con) {
                                // sharesTwo just checks if any 2 particles in the two triangles are the same
                                // ie, if they share an edge
                                if (t.sharesTwo(q)) {
                                    // if they do, the unconnected one is connected to the root
                                    modified = true;
                                    toMove = t;
                                    break;
                                }
                            }
                            if (modified) break;
                        }
                        // if it is connected, move it to the other list and iterate again
                        if (toMove != null) {
                            uncon.remove(toMove);
                            con.add(toMove);
                        }
                    }

                    // if, after checking all the connections, there are unconnected particles
                    // then you need to do a second fracturing on the particle in question as well
                    if (uncon.size() > 0 && con.size() > 0) {
                        dup = createParticle(p.p.x, p.p.y, p.v.x, p.v.y);
                        dup.resetTensors();
                        for (FEMTriangle tri : uncon) {
                            dup.addTriangle(tri);
                            p.tris.remove(tri);

                            if (tri.A.equals(p))
                                tri.A = dup;
                            else if (tri.B.equals(p))
                                tri.B = dup;
                            else if (tri.C.equals(p))
                                tri.C = dup;
                            else
                                System.out.println("bad");
                        }
                    }
                }
            }
            // STAR SLASH BETWEEN THESE TWO COMMENTS

            // TO ALLOW HINGE JOINTS AGAIN

           // newTriBorders(dup);
           // newTriBorders(frac);

            addSpringsFromTriangles();
            recomputeBorders();
        }
    }

    DenseVector xdot;
    DenseVector xz;
    DenseVector deltax;
    DenseVector deltaxdot;
    DenseVector f;
    DenseVector rhs;
    ConjugateGradientMTJ cgMTJ;
    
    /**
     * Initializes variables for backward Euler integration
     */
    public void init() {
        int N = particles.size();
        if ( xdot == null || xdot.size() != 2*N ) {
        	xdot = new DenseVector( 2*N );
        	xz = new DenseVector( 2*N );
        	deltax = new DenseVector( 2*N );
        	deltaxdot = new DenseVector( 2*N );
        	f = new DenseVector( 2*N );
            rhs = new DenseVector( 2*N );
        	cgMTJ = new ConjugateGradientMTJ( 2*N );
            cgMTJ.setFilter(this);
        }
    }
    
    /**
     * Fills in the provided vector with the particle velocities.
     * Helper function for implicit integration
     * @param xd
     */
    private void getVelocity( DenseVector xd ) {
    	int j = 0;
        for ( Particle p : particles ) {
            if( p.pinned ) {
                xd.set( j, 0 );
                xd.set( j+1, 0 );
            } else {
                xd.set( j, p.v.x );
                xd.set( j+1, p.v.y );
            }
            j += 2;
        }       
    }

    /** Helper function for implicit integration */
    private void setVelocity( DenseVector xdot ) {
    	int j = 0;
        for ( Particle p : particles ) {
            if( p.pinned ) {
                p.v.set(0,0);
            } else {
                p.v.x = xdot.get(j);
                p.v.y = xdot.get(j+1);
            }
            j += 2;
        }
    }
    
    /** Helper function for implicit integration */
    private void getForce( DenseVector f ) {
    	int j = 0;
        for ( Particle p : particles ) {
        	f.set( j, p.f.x );
        	f.set( j+1, p.f.y );
        	j += 2;
        }
    }
    
    /**
     * Fills the provided vector with the current positions of the particles.
     * Helper function for implicit integration
     * @param x
     */
    private void getPosition( DenseVector x ) {
    	int j = 0;
    	for ( Particle p: particles ) {
            x.set( j, p.p.x );
            x.set( j+1, p.p.y );
            j += 2;
    	}	
    }
    
    /**
     * Creates a new particle and adds it to the system
     * @param x
     * @param y
     * @param vx
     * @param vy
     * @return the new particle
     */
    public Particle createParticle( double x, double y, double vx, double vy ) {
        Particle p = new Particle( x, y, vx, vy );
        particles.add( p );
        return p;
    }
    
    /**
     * Creates a new spring between two particles and adds it to the system.
     * @param p1
     * @param p2
     * @return the new spring
     */
    public Spring createSpring( Particle p1, Particle p2 ) {
        Spring s = new Spring( p1, p2 ); 
        springs.add( s );
        collidableEdges.add( s );
        return s;
    }
    
    @Override
    public void init(GLAutoDrawable drawable) {
        // do nothing
    }

    @Override
    public void display(GLAutoDrawable drawable) {
        GL2 gl = drawable.getGL().getGL2();

        if ( useMouseSpring ) {
        	gl.glColor4d( 1,0,0, 1 );
        	gl.glBegin( GL.GL_LINES );
            gl.glVertex2d( mouseSpring.p1.p.x, mouseSpring.p1.p.y );
            gl.glVertex2d( mouseSpring.p2.p.x, mouseSpring.p2.p.y );
            gl.glEnd();        
        }

        double evs = evScale.getValue();

        gl.glLineWidth( 1 );
        for ( FEMTriangle f : femSprings ) {
            gl.glColor4d( 0,0,0, 0.5 );
            gl.glBegin( GL.GL_LINE_LOOP );
            gl.glVertex2d( f.A.p.x, f.A.p.y );
            gl.glVertex2d( f.B.p.x, f.B.p.y );
            gl.glVertex2d( f.C.p.x, f.C.p.y );
            gl.glEnd();
            if ( drawStressTensor.getValue() ) {
                f.display(drawable, evs );
            }
        }
        
      //  gl.glColor4d( 0,0.5,0, 1 );
        gl.glLineWidth( 1 );
        gl.glBegin( GL.GL_LINES );
        /*
        for (Spring s : springs) {
            gl.glVertex2d( s.p1.p.x, s.p1.p.y );
            gl.glVertex2d( s.p2.p.x, s.p2.p.y );
        }        */
        gl.glColor4d( 1,0,0, 1 );
        for ( Spring e : collidableEdges ) {
        	gl.glVertex2d( e.p1.p.x, e.p1.p.y );
        	gl.glVertex2d( e.p2.p.x, e.p2.p.y );
        }
        gl.glEnd();



        double sts = stScale.getValue();
        
        if ( drawSeparationTensor.getValue() ) {
        	for ( Particle p : particles ) {
        		p.drawSeparationTensor(drawable, sts);
        	}
        }
        
        if ( drawParticles.getValue() ) {
            gl.glPointSize( pointSize.getFloatValue() );
            gl.glBegin( GL.GL_POINTS );
            for ( Particle p : particles ) {
                // transparency is used to get smooth edges on the particles
                double alpha = 1;//  0.75;
                if ( p.pinned ) {
                    gl.glColor4d( 1, 0, 0, alpha );
                } else {
                	gl.glColor4d( 0,0,0, alpha );
                }
                gl.glVertex2d( p.p.x, p.p.y );
            }
            gl.glEnd();
        }
        
    }    
    
    BooleanParameter implicit = new BooleanParameter( "Implicit integration", false );
    
    IntParameter newtonIterations = new IntParameter( "Newton root solve iterations", 1, 1, 20 );

    BooleanParameter drawParticles = new BooleanParameter( "draw Particles", true ) ;
    
    DoubleParameter pointSize = new DoubleParameter("point size", 5, 1, 25);
    
    DoubleParameter evScale = new DoubleParameter("ev scale", 0.1, 1e-3, 1e3 );
    DoubleParameter stScale = new DoubleParameter("st scale", 1e-7, 1e-10, 1 );
    BooleanParameter drawSeparationTensor = new BooleanParameter("draw separation tensor", false );
    BooleanParameter drawStressTensor = new BooleanParameter("draw stress tensor", false );
    
    BooleanParameter useg = new BooleanParameter( "use gravity", false );
    
    DoubleParameter g = new DoubleParameter( "gravity", 98, 0.01, 1000 );
    
    DoubleParameter k = new DoubleParameter( "spring stiffness", 100, 0.01, 100000 );
        
    DoubleParameter PoissonRatio  = new DoubleParameter( "Poisson Ratio", 0.3, -1, .5 );

    DoubleParameter YoungModulus = new DoubleParameter( "YoungModulus", 2000, 1, 1e10 );

    DoubleParameter toughness = new DoubleParameter("material toughness", 2, 1, 1e8 );
    
    DoubleParameter RaleighAlpha = new DoubleParameter("Raleigh alpha", 1e-2, 1e-3, 1e3 );
    DoubleParameter RaleighBeta = new DoubleParameter("Raleigh beta", 1e-2, 1e-3, 1e3 );
    BooleanParameter useRAlpha = new BooleanParameter("use Raleigh alpha", true );
    BooleanParameter useRBeta = new BooleanParameter("use Raleigh beta", false );
    
    DoubleParameter b = new DoubleParameter( "spring damping", 1, 0, 10 );
    
    DoubleParameter c = new DoubleParameter( "viscous damping", 1, 0, 10 );

    DoubleParameter restitution = new DoubleParameter( "restitution", .0001, 0, 1 );
    
    DoubleParameter H = new DoubleParameter( "min distance (H)", 2, 0.1, 10 );

    JTextArea comments = new JTextArea("<comments>");
        
    BooleanParameter showCommentsAndParameters = new BooleanParameter("show comments and parameters", true );
            
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        
        VerticalFlowPanel vfp0 = new VerticalFlowPanel();
        vfp0.setBorder( new TitledBorder("Viewing Parameters" ) );
        vfp0.add( drawParticles.getControls() );
        vfp0.add( pointSize.getSliderControls(false) );
        vfp0.add( comments );
        vfp0.add( showCommentsAndParameters.getControls() );
        CollapsiblePanel cp0 = new CollapsiblePanel( vfp0.getPanel() );
        cp0.collapse();
        vfp.add( cp0 );
        
        VerticalFlowPanel vfp1 = new VerticalFlowPanel();
        vfp1.setBorder( new TitledBorder("Simulation Parameters"));
        vfp1.add( implicit.getControls() );
        vfp1.add( newtonIterations.getControls() );
        vfp1.add( useg.getControls() );
        vfp1.add( g.getSliderControls(true) );
        vfp1.add( k.getSliderControls(true) );
        vfp1.add( b.getSliderControls(false) );
        vfp1.add( c.getSliderControls(false) );
        vfp1.add( YoungModulus.getSliderControls(true));
        vfp1.add( PoissonRatio.getSliderControls(false));
        vfp1.add( toughness.getSliderControls(true));
        vfp1.add( RaleighAlpha.getSliderControls(true));
        vfp1.add( RaleighBeta.getSliderControls(true));
        vfp1.add( useRAlpha.getControls() );
        vfp1.add( useRBeta.getControls() );
        vfp1.add( evScale.getSliderControls(true));
        vfp1.add( stScale.getSliderControls(true));
        vfp1.add( drawSeparationTensor.getControls() );
        vfp1.add( drawStressTensor.getControls() );
        
        
        vfp1.add( restitution.getSliderControls(false) );
        vfp1.add( H.getSliderControls(false) );
        CollapsiblePanel cp1 = new CollapsiblePanel( vfp1.getPanel() );
        cp1.collapse();
        vfp.add( cp1 );
        
        return vfp.getPanel();        
    }
    
    @Override
    public String toString() {
    	DecimalFormat df = new DecimalFormat("0.000");
        String s = "particles = " + particles.size() + "\n" + "time = " + df.format(time) + "\n";
        if ( showCommentsAndParameters.getValue() ) {
        	s +=
               comments.getText() + "\n" + 
               "mouse spring stiffness = " + k.getValue() + "\n" + 
               "Young's modulus = " + YoungModulus.getValue() + "\n" +
               "Poisson's ratio = " + PoissonRatio.getValue() + "\n" +
               "Raleigh alpha = " + RaleighAlpha.getValue() + " " + ((!useRAlpha.getValue())?"(unused)":"") + "\n" +
               "Raleigh beta = " + RaleighBeta.getValue() + " " + ((!useRBeta.getValue())?"(unused)":"")+ "\n" +
               "\n";
        }
        return s;
    }

	@Override
	public void mult(Vector v, Vector Av) {
		// for implicit integration
	}

	@Override
	public void filter(Vector v) {
		int i = 0;
		for ( Particle p : particles ) {
			if ( p.pinned ){
            	v.set( i, 0 );
            	v.set( i, 1 );
            }
			i += 2;
        }
		// could do other filters based on Witkin and Baraff
	}
    
}
