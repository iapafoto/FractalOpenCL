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
import javax.swing.JPanel;

/**
 *
 * @author sebastien.durand
 */
public class PanelController implements ComponentListener, MouseWheelListener, MouseListener, MouseMotionListener, KeyListener {
    protected static final String BUNDLE = "fr/meteo/synopsis/client/synmap/Bundle";

    /** Delai pour les timers en ms */
    private static final int DELAY = 500;
    
    /** The map */
    protected final PanelControllerListener map;
            
    protected boolean bZooming = false,
                      bDragging = false,
                      bResizing = false,
                      bAllowDrag = true;

    // Window size    
    private int xSizeMem = 0,
                ySizeMem = 0;
    
    // Dragging (avec animation)
    private Point2D dragMem = null;
    private long dragWhen = 0;
    
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
    public PanelController(final PanelControllerListener map, final JPanel panel) {
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
        if (e.isShiftDown()) return;
        
        if (bAllowDrag) {
            if (!bDragging) { // Cas particulier de la surcharge du mouse pressed par la classe fille
                mousePressed(e);
            } else {
                dragAnim.updateVelocity(dragMem, e.getPoint(), dragWhen, e.getWhen());
                map.onDragged(e.getX() - dragMem.getX(), e.getY() - dragMem.getY());
                dragMem = e.getPoint();
                dragWhen = e.getWhen();
            }
        }
    }

    @Override
    public void mousePressed(MouseEvent e) {
        if (e.isShiftDown()) return;
                
        e.getComponent().requestFocus();
        if (e.getButton() == MouseEvent.BUTTON1) {
            dragAnim.stop();
            bDragging = true;
            dragMem = e.getPoint();
            dragWhen = e.getWhen();
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
        onResize(w, h);
    }

    public void onResize(int w, int h) {
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
          //  projectionModificationBegin();
            bDragging = true;
            dragAnim.animWithVelocity(0, DragAnim.MAX_VELOCITY*.7);
        } else if (key == KeyEvent.VK_DOWN) {
          //  projectionModificationBegin();
            bDragging = true;
            dragAnim.animWithVelocity(0, -DragAnim.MAX_VELOCITY*.7);            
        } else if (key == KeyEvent.VK_LEFT) {
          //  projectionModificationBegin();
            bDragging = true;
            dragAnim.animWithVelocity(DragAnim.MAX_VELOCITY*.7,0);            
        } else if (key == KeyEvent.VK_RIGHT) {
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
