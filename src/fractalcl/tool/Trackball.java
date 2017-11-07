package fractalcl.tool;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

import fr.meteo.synopsis.client.geometry.geometry3d.Camera3D;
import fr.meteo.synopsis.client.geometry.geometry3d.Geometry3D;
import fr.meteo.synopsis.client.geometry.geometry3d.Line3D;
import fr.meteo.synopsis.client.geometry.geometry3d.Point3D;
import fr.meteo.synopsis.client.geometry.geometry3d.Quaternion;
import fr.meteo.synopsis.client.geometry.geometry3d.Ray3D;
import fr.meteo.synopsis.client.geometry.geometry3d.Sphere3D;
import fr.meteo.synopsis.client.geometry.geometry3d.Vector3D;
import java.awt.geom.Point2D;

/**
 *
 * @author durands
 */
public class Trackball {
    
	// Some constants for convenience.
	private static final int X = 0;
	private static final int Y = 1;
	private static final int Z = 2;
	
	/**
	*  Storage for 3D point information passed between methods.
	**/
//	private final float[] op = new float[3], 
//                oq = new float[3], a = new float[3];
	
        
        // Trouve le point de la droite le plus proche de c
	public static Point3D closest(final Point3D c, final Line3D l, final Point3D result) {
            final Vector3D pc = new Vector3D(l.p0, c);
            return l.getParametricPosition(Vector3D.dot(pc, l.n), result);
        }
        
        // Trouve le point de la sphere le plus proche de c
	public static Point3D closest(final Point3D c, final Sphere3D sph, final Point3D result) {
            Ray3D ray2 = new Ray3D(sph.center, c);
            return ray2.getParametricPosition(sph.radius, result);
        }        

        public static Point3D intersect(final Sphere3D sphere, final Line3D ray, final Point3D result) {
            final double radius2 = sphere.radius*sphere.radius;
            final Vector3D dist = new Vector3D(sphere.center, ray.p0);
            final double B = Vector3D.dot(dist, ray.n);
            final double D = B*B - Vector3D.dot(dist, dist) + radius2;
            if (D < 0) {
                return null;
            }
            final double thc = Math.sqrt(D);
            if (B - thc > 0)
                return ray.getParametricPosition(B - thc, result);
            // Sinon on retourne le point le plus proche sur la sphere
            Point3D closestOnLine = closest(sphere.center, ray, result);
            return closest(closestOnLine, sphere, result);
        }
        
	/**
	*  Calculate a rotation matrix based on the axis and angle of rotation
	*  from the last 2 locations of the mouse relative to the Virtual
	*  Sphere cue circle.
	*
     * @param cam
	*  @param  pnt1       The 1st mouse location in the window.
	*  @param  pnt2       The 2nd mouse location in the window.
	*  @param  cueCenter  The center of the virtual sphere in the window.
	*  @param  cueRadius  The radius of the virtual sphere.
	*  @return A reference to the input rotMatrix is returned with the elements filled in.
	**/
	public static Camera3D makeRotationMtx(final Camera3D cam, 
                final Point2D pnt1, final Point2D pnt2, final Point3D cueCenter, final double cueRadius) {
            
            Ray3D ray0 = new Ray3D(new Point3D(), new Vector3D()),
                  ray1 = new Ray3D(new Point3D(), new Vector3D());
            
            cam.screenToRay(pnt1, ray0);
            cam.screenToRay(pnt2, ray1);

            Sphere3D sph = new Sphere3D(cueCenter, cueRadius);
            Point3D pt0 = new Point3D(), pt1 = new Point3D();
            
            intersect(sph, ray0, pt0);
            intersect(sph, ray1, pt1);
            
            Vector3D v1 = new Vector3D(sph.center, pt0, true),
                     v2 = new Vector3D(sph.center, pt1, true);
            
            Vector3D axis = Vector3D.cross(v1, v2, new Vector3D());
            double angle = Vector3D.dot(v1, v2);
            
            return cam.getRotated(new Quaternion(axis, angle), cueCenter);
	}
        
	public static Camera3D makeRotationMtx(final Camera3D cam, final Point3D pt0, final Point2D pnt2, final Point3D cueCenter) {

            double cueRadius = Geometry3D.distance(pt0, cueCenter);
            
            Ray3D ray1 = new Ray3D(new Point3D(), new Vector3D());
            cam.screenToRay(pnt2, ray1);
            Sphere3D sph = new Sphere3D(cueCenter, cueRadius);
            Point3D pt1 = new Point3D();
            
            intersect(sph, ray1, pt1);
            
            Vector3D v1 = new Vector3D(cueCenter, pt0, true),
                     v2 = new Vector3D(cueCenter, pt1, true);
            
            Vector3D axis = Vector3D.cross(v1, v2, new Vector3D());
            double angle = Vector3D.dot(v1, v2);
            
            return cam.getRotated(new Quaternion(axis, angle), cueCenter);
	}
	
	
	/**
	*  Project a 2D point on a circle to a 3D point on the +z hemisphere of a unit sphere.
	*  If the 2D point is outside the circle, it is first mapped to the nearest point on
	*  the circle before projection.
	*  Orthographic projection is used, though technically the field of view of the camera
	*  should be taken into account.  However, the discrepancy is neglegible.
	*
	*  @param  p         Window point to be projected onto the sphere.
	*  @param  cueCenter Location of center of virtual sphere in window.
	*  @param  cueRadius The radius of the virtual sphere.
	*  @param  v         Storage for the 3D projected point created by this method.
	**/
//	private static void pointOnUnitSphere(Point p, Point cueCenter, int cueRadius, float[] v) {
//		
//		/* Turn the mouse points into vectors relative to the center of the circle
//	 	*  and normalize them.  Note we need to flip the y value since the 3D coordinate
//	 	*  has positive y going up.  */
//		float vx = (p.x - cueCenter.x)/(float)cueRadius;
//		float vy = (cueCenter.y - p.y)/(float)cueRadius;
//		float lengthSqared = vx*vx + vy*vy;
//		
//		/* Project the point onto the sphere, assuming orthographic projection.
//	 	*  Points beyond the virtual sphere are normalized onto 
//	 	*  edge of the sphere (where z = 0).  */
//	 	float vz = 0;
//		if (lengthSqared < 1)
//			vz = (float)Math.sqrt(1.0 - lengthSqared);
//			
//		else {
//			float length = (float)Math.sqrt(lengthSqared);
//			vx /= length;
//			vy /= length;
//		}
//		
//		v[X] = vx;
//		v[Y] = vy;
//		v[Z] = vz;
//	}
	
	/**
	*  Computes a rotation matrix that would map (rotate) vectors op onto oq.
	*  The rotation is about an axis perpendicular to op and oq.
	*  Note this routine won't work if op or oq are zero vectors, or if they
	*  are parallel or antiparallel to each other.
	*
	*  <p>  Modification of Michael Pique's formula in 
	*       Graphics Gems Vol. 1.  Andrew Glassner, Ed.  Addison-Wesley.  </p>
	*
	*  @param  rotationMatrix  The 16 element rotation matrix to be filled in.
	*  @param  op              The 1st 3D vector.
	*  @param  oq              The 2nd 3D vector.
	**/
//	private void setRotationMatrix(float[] rotationMatrix, float[] op, float[] oq) {
//		
//		// Vector a is defined as a class variable to avoid wastefull memory allocations.
//
//		GLTools.crossProduct3D(op, oq, a);
//		float s = GLTools.length3D(a);
//		float c = GLTools.dotProduct3D(op, oq);
//		float t = 1 - c;
//
//		float ax = a[X];
//		float ay = a[Y];
//		float az = a[Z];
//		if (s > 0) {
//			ax /= s;
//			ay /= s;
//			az /= s;
//		}
//
//		float tax = t*ax;
//		float taxay = tax*ay, taxaz = tax*az;
//		float saz = s*az, say = s*ay;
//		rotationMatrix[0] = tax*ax + c;
//		rotationMatrix[1] = taxay + saz;
//		rotationMatrix[2] = taxaz - say;
//
//		float tay = t*ay;
//		float tayaz = tay*az;
//		float sax = s*ax;
//		rotationMatrix[4] = taxay - saz;
//		rotationMatrix[5] = tay*ay + c;
//		rotationMatrix[6] = tayaz + sax;
//
//		rotationMatrix[8] = taxaz + say;
//		rotationMatrix[9] = tayaz - sax;
//		rotationMatrix[10] = t*az*az + c;
//
//		rotationMatrix[3] = rotationMatrix[7] = rotationMatrix[11] = 
//		rotationMatrix[12] = rotationMatrix[13] = rotationMatrix[14] = 0;
//		rotationMatrix[15] = 1;
//	}
}