package a4fracture;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import com.jogamp.opengl.util.gl2.GLUT;
import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.swing.CollapsiblePanel;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.EasyViewer;
import mintools.viewer.Interactor;
import mintools.viewer.SceneGraphNode;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;

/**
 * Fracture and continuum elastic simulation application
 * @author kry
 */
public class FractureApp implements SceneGraphNode, Interactor {

    private EasyViewer ev;
    
    public ParticleSystem system;

    public TestSystems testSystems;

    private double grabThresh = 10;
    
    /**
     * Entry point for application
     * @param args
     */
    public static void main(String[] args) {
        new FractureApp();        
    }
        
    /**
     * Creates the application / scene instance
     */
    public FractureApp() {
        system = new ParticleSystem();
        testSystems = new TestSystems( system );
        ev = new EasyViewer( "A4 - FEM + Fracture", this, new Dimension(640,480), new Dimension(640,480) );
        // we add ourselves as an interactor to set up mouse and keyboard controls
        ev.addInteractor(this);
    }
     
    @Override
    public void init(GLAutoDrawable drawable) {
        GL gl = drawable.getGL();
        gl.glEnable( GL.GL_BLEND );
        gl.glBlendFunc( GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA );
        gl.glEnable( GL.GL_LINE_SMOOTH );
        gl.glEnable( GL2.GL_POINT_SMOOTH );
        system.init(drawable);
    }
        
    @Override
    public void display(GLAutoDrawable drawable) {

        GL2 gl = drawable.getGL().getGL2();        
        gl.glClearColor( 1, 1, 1, 1 );
        gl.glClear( GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT );
        
        // First advance the system (if it is running or wants to be stepped)
        if ( run.getValue() || stepRequest ) {
            for ( int i = 0; i < substeps.getValue(); i++ ) {
                system.updateParticles( stepsize.getValue() / substeps.getValue() );                
            }
        }
                
        // We're doing 2D drawing only, so we'll 
        // call this helper function to set up a 2D projection 
        // where the origin is at the top left corner and the 
        // units are pixels
        EasyViewer.beginOverlay(drawable);
        
        gl.glTranslated(screenTranslationx.getValue(), screenTranslationy.getValue(), 0);
        
        gl.glDisable( GL2.GL_LIGHTING );
        
        system.display( drawable );

        // Here we'll display some extra stuff for the interface        
        if ( mouseDown ) {
            if ( grabbed ) {
                gl.glPointSize( 15f );
                gl.glColor4d(0,1,0,0.95);
                gl.glBegin( GL.GL_POINTS );
                gl.glVertex2d( p1.p.x, p1.p.y );
                gl.glEnd();        
            }
        } else {
            if ( mouseInWindow ) {
                findCloseParticles( xcurrent, ycurrent );
                if ( p1 != null && d1 < grabThresh ) {
                    gl.glPointSize( 15f );
                    gl.glColor4d(0,1,0,0.95);
                    gl.glBegin( GL.GL_POINTS );
                    gl.glVertex2d( p1.p.x, p1.p.y );
                    gl.glEnd();        
                }
            }
        }

        // Finally we'll display a string with useful information
        // about the system and the current stepping        
        String text = system.toString();
        if ( displayStepInformation.getValue() ) {
                text +=
                    "\n" + "h = " + stepsize.getValue() + "\n" +
                      "substeps = " + (int)(double) substeps.getValue() ;
        }
        gl.glColor3f(0.5f,0.5f,0.5f);
        EasyViewer.printTextLines( drawable, text, 10, 15, 18, GLUT.BITMAP_9_BY_15 );
        //EasyViewer.printTextLines( drawable, text );
        EasyViewer.endOverlay(drawable);    

        // If we're recording, we'll save the step to an image file.
        // we'll also clear the step request here.
        if ( run.getValue() || stepRequest ) {
            stepRequest = false;        
            if ( record.getValue() ) {
                // write the frame
                File file = new File( "stills/" + dumpName + format.format(nextFrameNum) + ".png" );                                             
                nextFrameNum++;
                file = new File(file.getAbsolutePath().trim());
                ev.snapshot(drawable, file);
            }
        }
    }
    
    /** 
     * boolean to signal that the system was stepped and that a 
     * frame should be recorded if recording is enabled
     */
    private boolean stepRequest = false;
        
    /**
     * Base name of images to save
     */
    private String dumpName = "img";
    
    /**
     * Index for the frame number we're saving.
     */
    private int nextFrameNum = 0;
    
    /**
     * For formating the image file name when recording frames
     */
    private NumberFormat format = new DecimalFormat("00000");
    
    private BooleanParameter record = new BooleanParameter( "record each step to image file (press ENTER in canvas to toggle)", false );

    private BooleanParameter run = new BooleanParameter( "simulate (press SPACE in canvas to toggle)", false );
    
    private DoubleParameter stepsize = new DoubleParameter( "step size", 0.01, 1e-5, 1 );
    
    private DoubleParameter substeps = new DoubleParameter( "sub steps (integer)", 1, 1, 100);
        
    private BooleanParameter displayStepInformation = new BooleanParameter( "display stepping information overlay", true );
    
    private DoubleParameter screenTranslationx = new DoubleParameter( "screen x" , 0, -500,500 );
    private DoubleParameter screenTranslationy = new DoubleParameter( "screen x" , 0, -500,500 );
    
    
    @Override
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
    
        vfp.add( screenTranslationx.getSliderControls(false) );
        vfp.add( screenTranslationy.getSliderControls(false) );
        
        vfp.add( record.getControls() );        
        
        VerticalFlowPanel vfp0 = new VerticalFlowPanel();
        vfp0.setBorder( new TitledBorder("Numerical Integration Controls"));
        vfp0.add( run.getControls() );        
        vfp0.add( stepsize.getSliderControls(true) );
        vfp0.add( substeps.getControls() );
        vfp0.add( displayStepInformation.getControls() );
        CollapsiblePanel cp0 = new CollapsiblePanel( vfp0.getPanel() );
        cp0.collapse();
        vfp.add( cp0 );
        
        vfp.add( system.getControls() );
        
        vfp.add( testSystems.getControls() );
        
        return vfp.getPanel();
    }
    
    // Some member variables to help us keep track of close particles
    private Particle p1 = null;
    private Particle p2 = null;
    private double d1 = 0;
    private double d2 = 0;
    
    /**
     * Finds the two closest particles for showing potential spring connections
     * @param x 
     * @param y 
     */
    private void findCloseParticles( int x, int y ) {
        List<Particle> particles = system.particles;
        p1 = null;
        p2 = null;
        d1 = 0;
        d2 = 0;
        if ( particles.size() > 0 ) {
            for ( Particle p : particles ) {                
                double d = p.distance( x, y );
                if ( p1 == null || d < d1 ) {
                    p2 = p1; d2 = d1; p1 = p; d1 = d;
                } else if ( p2 == null || d < d2 ) {
                    p2 = p; d2 = d;
                }
            }      
        }
    }  
    
    private int xdown = 0;
    private int ydown = 0;
    private int xcurrent = 0;
    private int ycurrent = 0;
    private boolean mouseDown = false;
    private boolean mouseInWindow = false;
    private boolean grabbed = false;
    
    @Override
    public void attach(Component component) {
        component.addMouseMotionListener( new MouseMotionListener() {
            @Override
            public void mouseDragged(MouseEvent e) {
                xcurrent = e.getPoint().x;
                ycurrent = e.getPoint().y;
                if ( system.useMouseSpring ) {
                	system.mouseSpring.p2.p.set( xcurrent, ycurrent );
                	system.mouseSpring.p2.v.set( 0, 0 );
                } else if ( grabbed ) {
                    p1.p.set( xcurrent, ycurrent );
                    p1.v.set( 0, 0 ); 
                    if ( ! run.getValue() ) {
                        p1.p0.set( p1.p );
                        p1.v0.set( p1.v );
                        for ( Spring s : p1.springs ) {
                            s.computeRestLength();
                        }
                    }
                } else {
                    findCloseParticles(xcurrent, ycurrent);
                }
            }
            @Override
            public void mouseMoved(MouseEvent e) {
                xcurrent = e.getPoint().x;
                ycurrent = e.getPoint().y;
            }
        } );
        component.addMouseListener( new MouseListener() {
            @Override
            public void mouseClicked(MouseEvent e) {
                // do nothing
            }
            @Override
            public void mouseEntered(MouseEvent e) {
                mouseInWindow = true;
            }
            @Override
            public void mouseExited(MouseEvent e) {
                // clear the potential spring lines we're drawing
                mouseInWindow = false;
            }
            @Override
            public void mousePressed(MouseEvent e) {
                xdown = e.getPoint().x;
                ydown = e.getPoint().y;
                xcurrent = xdown;
                ycurrent = ydown;
                mouseDown = true;
                findCloseParticles(xcurrent, ycurrent);
                if ( p1 != null && d1 < grabThresh ) {
                    //p1.pinned = true;
                	if ( run.getValue() ) {
                		system.mouseSpring.p1 = p1;
                		system.mouseSpring.p2.p.set( xcurrent, ycurrent );
                		system.useMouseSpring = true;
                	} else {
                		grabbed = true;
                	}
                }
            }
            @Override
            public void mouseReleased(MouseEvent e) {
                mouseDown = false;
                if ( system.useMouseSpring ) {
                	system.useMouseSpring = false;
                } else if ( ! grabbed && ! run.getValue() ) {
                	// Let's not create new springs when we click in the window
//                    double x = e.getPoint().x;
//                    double y = e.getPoint().y;
//                    Particle p = system.createParticle( x, y, 0, 0 );
//                    if ( p1 != null && d1 < maxDist ) {
//                        system.createSpring( p, p1 );
//                    }
//                    if ( p2 != null && d2 < maxDist ) {
//                        system.createSpring( p, p2 );
//                    }  
                } else {
                    if ( p1 != null ) p1.pinned = ! p1.pinned;
                }
                grabbed = false;
            }
        } );
        component.addKeyListener( new KeyAdapter() {
            @Override
            public void keyPressed(KeyEvent e) {
                if ( e.getKeyCode() == KeyEvent.VK_SPACE ) {
                    run.setValue( ! run.getValue() ); 
                } else if ( e.getKeyCode() == KeyEvent.VK_S ) {                    
                    stepRequest = true;
                } else if ( e.getKeyCode() == KeyEvent.VK_R ) {
                    system.resetParticles();
                } else if ( e.getKeyCode() == KeyEvent.VK_C ) {                   
                    system.clear();
                    p1 = null;
                    p2 = null;
                } else if ( e.getKeyCode() == KeyEvent.VK_ESCAPE ) {
                    // quit the program
                    ev.stop();
                } else if ( e.getKeyCode() == KeyEvent.VK_ENTER ) {
                    // toggle recording of steps to png files
                    record.setValue( ! record.getValue() );
                } 
            }
        } );
    }
    
}
