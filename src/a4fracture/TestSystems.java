package a4fracture;

import mintools.parameters.BooleanParameter;
import mintools.swing.CollapsiblePanel;
import mintools.swing.FileSelect;
import mintools.swing.VerticalFlowPanel;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.vecmath.Matrix3d;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

/**
 * @author kry
 */
public class TestSystems {
        
    private BooleanParameter clearFirst = new BooleanParameter( "clear current system before creating new systems", true );
    
    private ParticleSystem system;
           
    /**
     * Creates a new test system 
     * @param system
     */
    public TestSystems( ParticleSystem system ) {
        this.system = system;
    }
    
    /**
     * Quick and dirty generic test generation button
     * @author kry
     */
    private class TestButton extends JButton implements ActionListener {
        private static final long serialVersionUID = 1L;
        private int testNumber;
        public TestButton( String name, int testNumber ) {
            super( name );
            this.testNumber = testNumber;
            addActionListener( this );
        }
        @Override
        public void actionPerformed(ActionEvent e) {
            createSystem(this.testNumber);
        }
    }
    
    /**
     * Gets the control panel for setting different systems.
     * @return the control panel
     */
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.setBorder( new TitledBorder("Particle System Test Systems"));

        vfp.add( clearFirst.getControls() );
        
        for ( int i = 0; i < tests.length; i++ ) {            
            vfp.add( new TestButton(tests[i], i) );             
        }

        CollapsiblePanel cp = new CollapsiblePanel( vfp.getPanel() );
        cp.collapse();
        return cp;   
    }
 
    public String[] tests = {            
            "single FEM triangle",
            "load FEM triangles from square.1.poly",
            "select FEM triangles to load from directory",
    };
    
    /**
     * Creates one of a number of simple test systems.
     *
     * Small systems are more useful for debugging!
     * 
     * @param which
     */
    public void createSystem( int which ) {
        if ( clearFirst.getValue() ) {
        	system.clear();
        }
        
        if ( which == 0 ) {
        	// looks like they need to be clockwise to be positive area... hopefully Triangle is consistent.
        	Particle p1 = new Particle( 100, 200, 0, 0 ); 
        	Particle p2 = new Particle( 120, 220, 0, 0 );
        	Particle p3 = new Particle( 140, 200, 0, 0 );        	
        	system.particles.add( p1 );
        	system.particles.add( p2 );
        	system.particles.add( p3 );
        	FEMTriangle tri = new FEMTriangle(p1, p3, p2);
        	tri.density = 1e-2;
        	system.femSprings.add( tri );
        	p1.mass = 1.0/3.0*tri.area*tri.density;
        	p2.mass = 1.0/3.0*tri.area*tri.density;
        	p3.mass = 1.0/3.0*tri.area*tri.density;

        	system.collidableEdges.add( new Spring(p1, p2) );
        	system.collidableEdges.add( new Spring(p2, p3) );
        	system.collidableEdges.add( new Spring(p3, p1) );

        } else if ( which == 1 ) {
        	// could actually do sommething that loads from the 
    		// output of Triangle, and can scale and place different meshed objects
    		// setting appropriate masses 
    		Matrix3d M = new Matrix3d();
    		M.setIdentity();
    		M.m02 = 300;
    		M.m12 = 300;
    		Triangle.loadTriangles(system, "dataFracture/square", M);    		
            system.name = tests[which];            
        } else if ( which == 2 ) {
        	File f = FileSelect.select("ele", "elements file", "load", "dataFracture", true);
        	String s = f.getName();
        	String fname = s.substring(0, s.length() - 6); // chop the ".1.ele" off
        	Matrix3d M = new Matrix3d();
    		M.setIdentity();
    		M.m02 = 300;
    		M.m12 = 300;
        	Triangle.loadTriangles(system, "dataFracture/"+fname, M);
        }

        system.addSpringsFromTriangles();
        system.recomputeBorders();
    }



}
