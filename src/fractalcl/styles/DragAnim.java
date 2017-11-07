/*
 * Permet de gerer l'animation d'un solide glissant sur une surface plane avec frottement
 */
package fractalcl.styles;

import fr.meteo.synopsis.client.geometry.geometry2d.Vec2;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Point2D;
import javax.swing.Timer;


/**
 *
 * @author sebastien.durand
 */
public final class DragAnim {

    // Constantes
    private static final long   
            DT_MAX = 100;               // Temps max pour considerer que deux evenements souris sont liés en ms
    private static final int    
            FRAME_DURATION = 10;        // Delai entre 2 frames de l'animation en Millisecondes
    public  static final double 
            MIN_VELOCITY = 25./1000.,   // 25 pixels par secondes
            MAX_VELOCITY = 5000./1000., // 5000 pixels par secondes
            R = .992; // .995;          // facteur de frottement (dt = 1 milliseconde)

    // Variables
    public interface DragAnimInterface {
        public void onMove(Point2D newPos, double dx, double dy);
        public void onEnd();
    }

    private final transient DragAnimInterface   dragInterface;
    private final transient Vec2            dragVelocity = new Vec2(0,0); // Vitesse en pixels par microseconde
    private transient long                      lastEventTime = 0; // Dernier temps lu en milliseconde
    private Point2D dragMem;
    private long dragWhen;
    
    private final Timer timerDragAnim = new Timer(FRAME_DURATION, new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
            if (dragVelocity.getLength() < MIN_VELOCITY) { // Plus assez de vitesse => Fin du mouvement
                dragInterface.onEnd();
                stop();
            } else {
                final double 
                        dt = System.nanoTime()/1000000. - lastEventTime,
                        dx = dragVelocity.x * dt,
                        dy = dragVelocity.y * dt;
                dragMem.setLocation(dragMem.getX()+dx, dragMem.getY()+dy);
                dragInterface.onMove(dragMem, dx, dy);
                dragVelocity.scale(Math.pow(R, dt));
                lastEventTime = (long)(System.nanoTime()/1000000.);
            }
        }
    });
                            
    public DragAnim(final DragAnimInterface parent) {
        this.dragInterface = parent;
    }
    
    public void start(final long t) {
        if ((t - lastEventTime > DT_MAX) || (dragVelocity.getLength() < MIN_VELOCITY)) { // Pas assez de vitesse => ignoré
            dragInterface.onEnd();
            stop();
        } else {
            lastEventTime = (long)(System.nanoTime()/1000000.);
            timerDragAnim.restart();
        }
    }
    
    public void animWithVelocity(final double vx, final double vy) {
        dragVelocity.set(vx,vy);
        lastEventTime = (long)(System.nanoTime()/1000000.);
        timerDragAnim.restart();
    }
    
    public void stop() {
        timerDragAnim.stop();      
        lastEventTime = 0;
        dragVelocity.set(0,0);
        dragMem = null;
        dragWhen = 0;
    }

    public void init(final Point2D p1, final long t1) {
        stop();
        dragMem = p1;
        dragWhen = t1;
    }
    
    public void updateVelocity(final Point2D p1, final long t1) {
        if (dragMem != null && dragWhen != 0) {
            updateVelocity(dragMem, p1, dragWhen, t1);
        }
        dragMem = p1;
        dragWhen = t1;
    }
    
    private void updateVelocity(final Point2D p0, final Point2D p1, final long t0, final long t1) {
        final double dt = t1 - t0;
        if ((dt < DT_MAX) && (dt > 0)){
            dragVelocity.set(p1, p0);
            dragVelocity.scale(1./dt); // Vitesse en pixel par milliseconde
            if (dragVelocity.getLength() > MAX_VELOCITY) { // Limitation de vitesse !
                dragVelocity.setLength(MAX_VELOCITY);    
            }
        } else {
            dragVelocity.set(0,0); // Vitesse nulle
        }
        lastEventTime = t1;
    }
}