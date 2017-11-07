/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.styles;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Window;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.geom.Path2D;
import javax.swing.JLabel;
import javax.swing.JWindow;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;

/**
 *
 * @author sebastien
 */
public class FlyingToolTip {
        
    public static interface FlyingToolTipListener {
        public void onPopupHide();
        public Point getLocationOnScreen();
        public void requestFocus();
        public Component getComponent();
        
    }
    
    /**
     * Composant gerant l'affichage du tooltip de couleur indiquant la valeur du slider selectionné
     */
   // public static class FlyingToolTip {
    final private FlyingToolTipListener source;
    final private JLabel label;
    
    private JWindow window; // La fenetre Popup
    Path2D path;

    // Detection des changements intervenants sur la fenetre parent => fermeture du popup
    private final ComponentListener componentListener = new ComponentListener() {
        @Override
        public void componentResized(ComponentEvent e) {
            dispose();
            source.onPopupHide();
        }
        @Override
        public void componentMoved(ComponentEvent e) {
            dispose();
            source.onPopupHide();
        }
        @Override
        public void componentShown(ComponentEvent e) {
        }
        @Override
        public void componentHidden(ComponentEvent e) {
            dispose();
            source.onPopupHide();
        }
    };

    public FlyingToolTip(final FlyingToolTipListener source) {
        this.source = source;
        label = new JLabel("test", JLabel.CENTER) {
            @Override
            public void paintComponent(Graphics g) {
                Graphics2D g2 = (Graphics2D) g.create();
                int w = getWidth(), h = getHeight();
                g.setColor(getBackground());
                g.fillRect(0,0,w,h);

                if (path != null) {                    
                    g2.setStroke(new BasicStroke(4f));
                    g2.setColor(getForeground());
                    g2.draw(path);
                }
                setOpaque(false);
                super.paintComponent(g2);
                setOpaque(true);
                g2.dispose();
            }
        };    
        label.setBorder(new EmptyBorder(2,4,10,4));
        label.setOpaque(true);
    }

    private void setText(final String txt) {
        if (txt != null) {
            label.setText("<html>"+txt+"</html>");
        } else {
            label.setText("---");
        }
        window.pack();
        refreshShape();
    }

    private void setPosition(final int x, final int y) {
        final Dimension dim = window.getSize();
        final Point pt = source.getLocationOnScreen();
        window.setLocation(pt.x + x-dim.width/2+1, pt.y-dim.height+3);
    }

    private void setColors(final Color background, final Color foreground) {
        label.setBackground(background);
        label.setForeground(foreground);
    }

    public void show(final int x, final String txt, final Color color) {
        if (window != null) {
            dispose();
        }            
        // On le recree pour tenir compte du mode multi ecran des previs
        final Window frame = SwingUtilities.getWindowAncestor(source.getComponent());
        window = new JWindow(frame);

        window.getContentPane().add(label);
        window.pack();

        frame.addComponentListener(componentListener);

        setColors(color, Color.black);
        setText(txt);
        setPosition(x, 0);
        window.setVisible(true);
        source.requestFocus();
    }

    public void move(final int x, final String txt, final Color color) {
        if (color != null) {
            setColors(color, Color.black);
        }
        setText(txt);
        setPosition(x, 0);
    }

    public void dispose() {            
        if (window != null) {
            final Window frame = SwingUtilities.getWindowAncestor(source.getComponent());
            frame.removeComponentListener(componentListener);
            window.setVisible(false);
            window.dispose();
            window = null;
        }
    }

    public void refreshShape() {
        int w = (int) window.getWidth(),     
            h = (int) window.getHeight();
        path = new Path2D.Double();
        path.moveTo(0,0);
        path.lineTo(w,0);
        path.lineTo(w,h-7);
        path.lineTo(w/2+7,h-7);
        path.lineTo(w/2,h);
        path.lineTo(w/2-7,h-7);
        path.lineTo(0,h-7);
        path.closePath();
        window.setShape(path);
    }
}

