/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.geom.Point2D;
import javax.swing.Timer;
import earth.tool.DragAnim;
import fr.meteo.synopsis.client.geometry.geometry3d.Camera3D;
import fr.meteo.synopsis.client.geometry.geometry3d.Point3D;
import fr.meteo.synopsis.client.geometry.geometry3d.Quaternion;
import fr.meteo.synopsis.client.geometry.geometry3d.Ray3D;
import fr.meteo.synopsis.client.geometry.geometry3d.Vector3D;
import fractalcl.tool.Trackball;
import javax.swing.JPanel;

/**
 *
 * @author sebastien.durand
 */
public class CameraController implements ComponentListener, MouseWheelListener, MouseListener, MouseMotionListener, KeyListener {
    protected static final String BUNDLE = "fr/meteo/synopsis/client/synmap/Bundle";

    /** Delai pour les timers en ms */
    private static final int DELAY = 500;
    
    /** The map */
    protected final CameraListener map;
            
    protected boolean bZooming = false,
                      bDragging = false,
                      bResizing = false,
                      bAllowDrag = true;

    // Window size    
    private int xSizeMem = 0,
                ySizeMem = 0;
    
    public Camera3D camera = new Camera3D(new Point3D(21,20.5,21), new Point3D(0,0,.2), new Vector3D(0,0,1), 3., 800, 600);
    
    // Dragging (avec animation)
    private Point2D dragMem = null;
    private long dragWhen = 0;
    Point3D rotStart, rotCenter;
    Camera3D cameraMem;
     
    private DragAnim dragAnim = new DragAnim(new DragAnim.DragAnimInterface() {
        @Override
        public void onMove(double dx, double dy) {
            map.onDragged(dx, dy);
        }
        @Override
        public void onEnd() {
            bDragging = false;
            map.onDragEnd();
        }
    });
    

    /** Timer de fin de détection du mouvement de souris pour mise à jour "réelle" */
    private final Timer timerWheel,
                        timerResize;

    /**
     * Constructor with parameter
     * @param map The map
     */
    public CameraController(final CameraListener map, final JPanel panel) {
        this.map = map;

        ActionListener actionWheel = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                onStopWheel();
            }
        };
        timerWheel = new Timer(DELAY, actionWheel);
        timerWheel.setRepeats(false);
        
        ActionListener actionResize = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                onStopResize();
            }
        };
        timerResize = new Timer(DELAY, actionResize);
        timerResize.setRepeats(false);
        
        panel.addMouseListener(this);
        panel.addMouseMotionListener(this);
        panel.addMouseWheelListener(this);
        panel.addComponentListener(this);
    }
    
    
// Dragging ///////////////////////////////////////////////   
    
    @Override
    public void mouseDragged(MouseEvent e) {
        double dx =0,dy=0;
        if (dragMem != null) {
            dx = (e.getX() - dragMem.getX())/(double)camera.getWidth();
            dy = (e.getY() - dragMem.getY())/(double)camera.getHeight();
        }
        
        if (bAllowDrag) {
            if (!bDragging) { // Cas particulier de la surcharge du mouse pressed par la classe fille
                mousePressed(e);
            } else {
                if (e.isShiftDown()) {
                    camera.moveLookDir(dy);
                    camera.moveSideDir(-dx);
                    
                } else if (e.isControlDown()) {
//                    Point3D rotC = rotCenter;
//                    if (rotC == null) {
//                        double dist = map.getZMapDistanceAt(0,0);
//                        Ray3D ray = camera.screenToRay(new Point2D.Double(0,0), new Ray3D());
//                        rotC = ray.getParametricPosition(dist, new Point3D());
//                    }
//                    if (cameraMem != null &&rotStart != null && rotCenter != null) {
//                        camera = Trackball.makeRotationMtx(cameraMem, rotStart, new Point2D.Double(e.getX(),e.getY()), rotC);
//                    } 
                    camera = camera.getRotated(new Quaternion(camera.getSideDir(), dy), camera.getPosition());
                    camera = camera.getRotated(new Quaternion(camera.getUpDir(), dx), camera.getPosition());
                } else {
                    camera.moveUpDir(dy);
                    camera.moveSideDir(-dx);
                }
                dragAnim.updateVelocity(dragMem, e.getPoint(), dragWhen, e.getWhen());
                map.onDragged(e.getX() - dragMem.getX(), e.getY() - dragMem.getY());
                dragMem = e.getPoint();
                dragWhen = e.getWhen();
            }
        }
    }
    
    public void setFocusPoint(Point2D ptScreen, double dist) {
        Ray3D ray = camera.screenToRay(ptScreen, new Ray3D());
        this.rotCenter = ray.getParametricPosition(dist, new Point3D());
    }

    @Override
    public void mousePressed(MouseEvent e) {
       // if (e.isShiftDown()) return;
        map.onMousePressed(e);
        
        e.getComponent().requestFocus();
        if (e.getButton() == MouseEvent.BUTTON1) {
            dragAnim.stop();
            bDragging = true;
            dragMem = e.getPoint();
            dragWhen = e.getWhen();
            
            cameraMem = new Camera3D(camera);
            double distance = map.getZMapDistanceAt(e.getX(),e.getY());
            Ray3D ray = camera.screenToRay(new Point2D.Double(e.getX(),e.getY()), new Ray3D());
            rotStart = ray.getParametricPosition(distance, new Point3D());
        }
    }
    
    @Override
    public void mouseReleased(MouseEvent e) {
        if (e.isShiftDown()) return;
        
        if (bDragging) {
            dragAnim.start(e.getWhen());
        }
    }
    
// Resize //////////////////////////////////////////////////
    
    public void onStartResize(int w, int h) {
        bResizing = true;
        camera.setScreenSize(w, h);
        onResize(w, h);
    }

    public void onResize(int w, int h) {
        camera.setScreenSize(w, h);
        map.onResize(w,h);
    }

    public void onStopResize() {
        bResizing = false;
        map.onStopResize();
    }
    
 // Wheel //////////////////////////////////////////////////

    private void onStartWheel(MouseWheelEvent e) {
        bZooming = true;
        onWheel(e);
    }    
    double zoom = 1.;
    private void onWheel(MouseWheelEvent e) {
        if (e.isShiftDown()) return;
      
        camera.moveLookDir(e.getWheelRotation());
        
        zoom *= (((e.getWheelRotation() < 0) ? 1.01 : 1./1.01));
        map.onWheel(e.getWheelRotation()/*zoom*/, e.getPoint());        
    }
          
    private void onStopWheel() {
        bZooming = false;
        map.onStopWheel(zoom);
    }

    
    /**
     * Call when the wheel move
     * @param e the event of the wheel movement
     */
    @Override
    public void mouseWheelMoved(MouseWheelEvent e) {
        if (e.isShiftDown()) return;
        
        if (!timerWheel.isRunning()) {
            onStartWheel(e);
        } else {
            onWheel(e);
        }
        timerWheel.restart();
    }
        
    @Override
    public void componentResized(ComponentEvent e) {
        int w = e.getComponent().getWidth();
        int h = e.getComponent().getHeight();
        
        if ((w > 0) && ((w != xSizeMem) || (h != ySizeMem))) {
            xSizeMem = w;
            ySizeMem = h;
            if (!timerResize.isRunning()) {
                onStartResize(w,h);
            } else {
                onResize(w,h);
            }
            timerResize.restart();
        }
    }

    
    @Override
    public void keyReleased(KeyEvent e) {
//        if (bDraggingKey) {
//            bDraggingKey = false;
//            projectionModificationEnd();
//        }
    }
    
    @Override
    public void keyTyped(KeyEvent ke) {}

    @Override
    public void keyPressed(KeyEvent ke) {
        int key = ke.getKeyCode();
        if (key == KeyEvent.VK_UP) {
            camera.moveUpDir(1.);
          //  projectionModificationBegin();
            bDragging = true;
            dragAnim.animWithVelocity(0, DragAnim.MAX_VELOCITY*.7);
        } else if (key == KeyEvent.VK_DOWN) {
            camera.moveUpDir(-1.);
          //  projectionModificationBegin();
            bDragging = true;
            dragAnim.animWithVelocity(0, -DragAnim.MAX_VELOCITY*.7);            
        } else if (key == KeyEvent.VK_LEFT) {
            camera.moveSideDir(1.);
          //  projectionModificationBegin();
            bDragging = true;
            dragAnim.animWithVelocity(DragAnim.MAX_VELOCITY*.7,0);            
        } else if (key == KeyEvent.VK_RIGHT) {
            camera.moveSideDir(-1.);
          //  projectionModificationBegin();
            dragAnim.animWithVelocity(-DragAnim.MAX_VELOCITY*.7,0);            
        }

    }

    @Override public void componentMoved(ComponentEvent e) {}
    @Override public void componentShown(ComponentEvent e) {}
    @Override public void componentHidden(ComponentEvent e) {}
    @Override public void mouseClicked(MouseEvent me) {
        this.map.onMouseClicked(me);
    }
    @Override public void mouseEntered(MouseEvent me) {}
    @Override public void mouseExited(MouseEvent me) {}
    @Override public void mouseMoved(MouseEvent me) {}



    
    
}
