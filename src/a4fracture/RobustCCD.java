package a4fracture;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

import com.jogamp.opengl.math.VectorUtil;

/**
 * Implementation of a robust collision detection.
 * @author kry
 */
public class RobustCCD {

	double restitution;

	double H;
	int maxResolutionIterations = 200;
	double eps = .000001;

	/**
	 * Creates the new continuous collision detection and response object
	 */
	public RobustCCD() {
		// do nothing
	}

	/**
	 * Try to deal with contacts before they happen
	 * @param h
	 * @param system
	 */
	public void applyRepulsion( double h, ParticleSystem system , double H) {

		// TODO: OBJECTIVE 4 elastic repulsion forces

		// check for point-edge proximity between all particle-edge pairs
		// find the normal
		// take care to deal with segment end points carefully
		// compute an appropriate  impulse if the distance is less than H
		// make sure your impulse is in the correct direction!
		// don't apply the impulse if the relative velocity is separating fast enough (see Bridson et al. 2002)
		// distribute impulse to the three particles involved in the appropriate manner

		// again, for each particle and spring combination
		for (Particle p : system.particles)
		{
			for (Spring s : system.collidableEdges)
			{
				if(!p.equals(s.p1) && !p.equals(s.p2)) // make sure the particle isnt in the spring
				{
					// a bunch of code to get the distance between the particle and the spring
					double d = 0;
					double l2 = s.p1.distance(s.p2.p.x, s.p2.p.y);
					l2 = l2 * l2;
					
					if(l2 == 0)
						d = s.p1.distance(p.p.x, p.p.y);

					Vector2d v1 = new Vector2d(p.p.x, p.p.y);
					v1.x -= s.p1.p.x;
					v1.y -= s.p1.p.y;
					Vector2d v2 = new Vector2d(s.p2.p.x, s.p2.p.y);
					v2.x -= s.p1.p.x;
					v2.y -= s.p1.p.y;
					double t = Math.max(0, Math.min(1, v1.dot(v2) / l2));

					v2.scale(t);
					v2.add(s.p1.p);
					
					d = p.distance(v2.x, v2.y);
					// ultimately, d is the distance between the two
					// v1 is the closest point on the spring
					// t is parameterized distance of closest point along spring

					// n is normal vector
					Vector2d n = new Vector2d(s.p1.p.y - s.p2.p.y, s.p2.p.x - s.p1.p.x);
					n.normalize();

					// make sure the normal is pointing in the right direction. otherwise negate it.
					if(n.dot(new Vector2d(p.p.x - v1.x, p.p.y - v1.y)) < 0)
						n.negate();

					// set spring stiffness to use for repulsion spring.
					double k = 100;

					// if the particle is closer than the specified threshold distance
					if(d < H)
					{
						// and the velocity is appropriate (.1 * d * h taken from Bridson 2002)
						if(p.v.dot(n) > .1 * d * h)
						{
							// compute the spring force, using threshold as default spring length
							Vector2d l = new Vector2d();	
							l.sub(p.p, v2);

							Vector2d ldot = new Vector2d(0, 0);
							double left = k * (l.length() - H); // here
							double right = ldot.dot(l);
							right = .01 * right/l.length();
							double total = left + right;

							Vector2d FA = l;
							FA.scale(-total/l.length());
							FA.scale(h);

							Vector2d FB = new Vector2d(FA);
							FB.negate();

							// if the particles arent pinned, apply the forces.
							if(!p.pinned)
							{
								p.v.x += FA.x;
								p.v.y += FA.y;
							}

							if(!s.p1.pinned)
							{
								s.p1.v.x += t * FB.x;
								s.p1.v.y += t * FB.y;
							}

							if(!s.p2.pinned)
							{
								s.p2.v.x += (1 - t) * FB.x;
								s.p2.v.y += (1 - t) * FB.y;
							}
						}
					}
				}
			}
		}
	}

	/**
	 * Checks all collisions in interval t to t+h
	 * @param h
	 * @param system 
	 * @return true if all collisions resolved
	 */
	public boolean check( double h, ParticleSystem system ) {        

		// For each particle-edge pair, find the roots for when the three particles are
		// co-linear, and then pick the first root on (0,h] which corresponds to an 
		// actual collision.  Compute a collision response.  That is, compute an appropriate
		// collision normal, compute the impulse, and then apply the impulse to the associated
		// particles.  Be sure to deal with pinning constraints!  Repeat until all collisions
		// are resolved and it is safe to advacne time

		// You probably want to write other methods to help with the job, and
		// call them here.  Or alternatively you can write one large and ugly
		// monolithic function here.

		// TODO: OBJECTIVE 1 continuous collision detection    	
		// TODO: OBJECTIVE 2 compute collision impulses
		// TODO: OBJECTIVE 3 iterative to resolve all collisions

		// make sure to iterate with a maximum
		boolean collided = true;
		int i = maxResolutionIterations;
		while(collided && i > 0)
		{
			collided = false;
			i--;
			for (Particle p : system.particles) // iterate trhough all particles 
			{
				for (Spring s : system.collidableEdges) // check that particle against each spring
				{
					if(!p.equals(s.p1) && !p.equals(s.p2)) // but not a spring that the particle is part of
					{
						// find out if the spring-particle pair has an alpha
						// returns -1 if not
						double alpha = getCollisionAlpha(s.p1, s.p2, p, h);

						// if there is a valid one (in range 0 to h)
						if(alpha > 0 - eps)
						{
							// notify the loop that there has been a collision
							collided = true;
							// and find the resulting impulse
							double j = getj(s.p1, s.p2, p, h, alpha);

							// get the normal vector n
							Vector2d n = new Vector2d(s.p1.p.y - s.p2.p.y, s.p2.p.x - s.p1.p.x);
							n.normalize();

							// and apply the force to each unpinned particle
							if(!s.p1.pinned)
							{
								s.p1.v.x -= j * n.x * (1 - alpha) / s.p1.mass;
								s.p1.v.y -= j * n.y * (1 - alpha) / s.p1.mass;
							}
							if(!s.p2.pinned)
							{
								s.p2.v.x -= j * n.x * alpha / s.p2.mass;
								s.p2.v.y -= j * n.y * alpha / s.p2.mass;
							}    
							if(!p.pinned)
							{
								p.v.x += j * n.x / p.mass;
								p.v.y += j * n.y / p.mass;
							}
						}
					}
				}
			}
		}

		return true; // can return false if you give up
	}

	double getj(Particle A, Particle B, Particle C, double h, double alpha)
	{
		// omalpha stands for one minus alpha
		double omalpha = 1 - alpha;

		// velocities of the point of contact with the particle that hit the spring
		double Pvx = omalpha * A.v.x + alpha * B.v.x;
		double Pvy = omalpha * A.v.y + alpha * B.v.y;

		// relative velocity before collision
		Vector2d Vrn = new Vector2d(Pvx - C.v.x, Pvy - C.v.y);

		// normal vector to the spring
		Vector2d normal = new Vector2d(A.p.y - B.p.y, B.p.x - A.p.x);
		normal.normalize();

		// calculate j
		double j = Vrn.dot(normal);
		j = j * (1 + restitution);

		double t1 = (1 - omalpha) * (1 - omalpha);
		double t2 = alpha * alpha;
		double t3 = 1;

		// if particles are pinned, their contribution goes to 0 as mass approaches infinity
		if(A.pinned)
			t1 = 0;
		else
			t1 = t1/A.mass;

		if(B.pinned)
			t2 = 0;
		else
			t2 = t2/B.mass;

		if(C.pinned)
			t3 = 0;
		else
			t3 = t3/C.mass;

		j = j  / (t1 + t2 + t3);

		return j;
	}

	double getCollisionAlpha(Particle A, Particle B, Particle C, double h)
	{
		// basically solve for the whole giant ax + by + c = 0 equation, from
		// the 3x3 det(Ax-Cx Ay-Cy 0)  ... etc matrix

		// Ax By - Cx By - Ax Cy - Ay Bx + Ay Cx + Cy Bx
		double a = 0;
		double b = 0;
		double c = 0;

		a += A.v.x * B.v.y;
		a -= C.v.x * B.v.y;
		a -= A.v.x * C.v.y;
		a -= A.v.y * B.v.x;
		a += A.v.y * C.v.x;
		a += C.v.y * B.v.x;

		// Ax By - Cx By - Ax Cy - Ay Bx + Ay Cx + Cy Bx
		b += A.v.x * B.p.y;
		b += A.p.x * B.v.y;

		b -= C.v.x * B.p.y;
		b -= C.p.x * B.v.y;

		b -= A.v.x * C.p.y;
		b -= A.p.x * C.v.y;

		b -= A.v.y * B.p.x;
		b -= A.p.y * B.v.x;

		b += A.v.y * C.p.x;
		b += A.p.y * C.v.x;

		b += C.v.y * B.p.x;
		b += C.p.y * B.v.x;

		// Ax By - Cx By - Ax Cy - Ay Bx + Ay Cx + Cy Bx
		c += A.p.x * B.p.y;
		c -= C.p.x * B.p.y;
		c -= A.p.x * C.p.y;
		c -= A.p.y * B.p.x;
		c += A.p.y * C.p.x;
		c += C.p.y * B.p.x;

		// now that you have the components of the quadratic, you can solve it
		double root = b * b - 4 * a * c;
		double r1, r2;
		double alpha = -1;
		List<Double> collisions = new ArrayList<Double>();

		// case where a = 0 so its a linear equation
		if(a == 0)
		{
			if(b != 0)
				collisions.add(-c / b);
		}
		// case where you have 2 roots
		else if(root > 0)
		{
			r1 = Math.sqrt(root);
			r2 = Math.sqrt(root);

			r1 = -b + r1;
			r2 = -b - r2;

			r1 = r1 / (a * 2);
			r2 = r2 / (a * 2);

			collisions.add(r1);
			collisions.add(r2);
		}
		// case where you have 1 root
		else if(root == 0)
		{
			r1 = -b;
			r1 = r1 / (a * 2);

			collisions.add(r1);
		}
		else
		{
			// not a real collision
		}

		double d = 5;
		// find the collision with the smallest d value from the ones you just computed
		// d is time in this case, so in other words, the next collision you would have
		for(Double dc : collisions)
		{
			// only care if its in the next time step, ie between 0 and h
			if (dc > 0 - eps && dc <= h + eps)
			{
				if (dc < d)
					d = dc;
			}
		}

		if(d > 0 - eps && d <= h + eps)
		{
			// cA cB cC are the collision positions of A, B, and C, which are used to find the
			// resulting alpha value, parameterize the line
			Particle cA = new Particle(A.p.x + A.v.x * d, A.p.y + A.v.y * d, A.v.x, A.v.y);
			Particle cB = new Particle(B.p.x + B.v.x * d, B.p.y + B.v.y * d, B.v.x, B.v.y);
			Particle cC = new Particle(C.p.x + C.v.x * d, C.p.y + C.v.y * d, C.v.x, C.v.y);

			if (Math.abs(cA.p.y - cB.p.y) > Math.abs(cA.p.x - cB.p.x))
				alpha = (cA.p.y - cC.p.y) / (cA.p.y - cB.p.y);
			else
				alpha = (cA.p.x - cC.p.x) / (cA.p.x - cB.p.x);

			// make sure alpha is actually a point on the line. 
			if(alpha > 0 - eps && alpha < 1 + eps)
				return alpha;
		}		


		return -1;
	}

}
