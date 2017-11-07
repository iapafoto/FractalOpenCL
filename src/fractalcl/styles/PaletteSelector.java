/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * PaletteSelector.java
 *
 * Created on 26 oct. 2011, 16:40:03
 */
package fractalcl.styles;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.ComponentEvent;
import java.awt.event.MouseEvent;


/**
 *
 * @author sebastien.durand
 */
public class PaletteSelector extends SelectorPanel<Color> {

    int nbLine = 1;
    int nbEltByLine = 1;
    int selectedId = 0;
    int eltSize = 10;
    Color[] palette;

    public static Color[] defaultColors() {
        int nbcolors = 125;
        Color[] colors = new Color[nbcolors];
        double step = Math.pow(255, 3) / (nbcolors - 1);
        for (int i = 0; i < colors.length; i++) {
            double prog = step * i;
            double r = prog % 255;
            prog /= 255;
            double v = prog % 255;
            prog /= 255;
            double b = prog % 255;
            colors[i] = new Color((int) r, (int) v, (int) b);
        }
        return colors;
    }

    /** for preview in Netbean Designer **/
    public PaletteSelector() {
        super();
        Color[] colors = defaultColors();
        setPalette(colors, colors[55]);
        addComponentListener(this);
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 400, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 300, Short.MAX_VALUE)
        );
    }// </editor-fold>//GEN-END:initComponents
    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables

    @Override
    protected void paintComponent(Graphics g) {
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        g2.setColor(getBackground());
        g2.fillRect(0, 0, getWidth(), getHeight());

        int sx = getSizeX();
        int sy = getSizeY();

        int nbColor = palette.length;
        int row, col, k = 0;

        for (row = 0; row < nbLine; row++) {
            for (col = 0; col < nbEltByLine; col++) {
                if (k >= nbColor) {
                    break;
                }

                g2.setColor(palette[k]);
                g2.fillRect(col * sx + 1, row * sy + 1, sx - 1, sy - 1);
                g2.setColor(Color.black);
                g2.drawRect(col * sx, row * sy, sx, sy);
                k++;
            }
        }
        highlightColor(g2, selectedId, Color.white, 5);
        if (over) {
            highlightColor(g2, overId, Color.white, 3);
        }
    }

    private void highlightColor(Graphics2D g2, int colorId, Color color, int zoom) {
        if ((colorId < 0) || (colorId >= palette.length)) {
            return;
        }

        int sx = getSizeX();
        int sy = getSizeY();
        int row = colorId / nbEltByLine;
        int col = colorId % nbEltByLine;
        int dx = 0;
        int dy = 0;

        if (col == 0) {
            dx = zoom + 2;
        }
        if (row == 0) {
            dy = zoom + 2;
        }
        if (col == nbEltByLine - 1) {
            dx = -zoom - 2;
        }
        if (row == nbLine - 1) {
            dy = -zoom - 2;
        }

        g2.setColor(palette[colorId]);
        g2.fillRect(dx + col * sx - zoom, dy + row * sy - zoom, sx + 2 * zoom, sy + 2 * zoom);
        g2.setStroke(new BasicStroke(3));
        g2.setColor(color);
        g2.drawRect(dx + col * sx - zoom - 2, dy + row * sy - zoom - 2, sx + 2 * (zoom + 2), sy + 2 * (zoom + 2));
        g2.setColor(Color.black);
        g2.drawRect(dx + col * sx - zoom, dy + row * sy - zoom, sx + 2 * zoom, sy + 2 * zoom);
    }

    @Override
    public void mouseClicked(MouseEvent e) {
        if (this.isEnabled()) {
            selectedId = getIdOnPoint(e.getX(), e.getY());
            repaint();
            firePropertyChange(SELECTION_CHANGE_PROPERTY, null, palette[selectedId]);
        }
    }

    private int getSizeX() {
        return eltSize;
    }

    private int getSizeY() {
        return eltSize;
    }

    @Override
    protected int getIdOnPoint(int x, int y) {
        return Math.min(Math.max(0, (y / getSizeY()) * nbEltByLine + Math.min(nbEltByLine - 1, x / getSizeX())), palette.length - 1);
    }

    private void calculBestSize(int w, int h) {
        int max = (int) Math.ceil(Math.sqrt(palette.length));
        int n2, s = 0;

        eltSize = 0;

        for (int n1 = 1; n1 < max; n1++) {
            n2 = (int) Math.ceil(((double) palette.length) / ((double) n1));
            if (n2 < 1) {
                continue;
            }
            s = w / n1;
            if (n2 * s <= h) { // ca passe
                if (s > eltSize) {
                    nbEltByLine = n1;
                    nbLine = n2;
                    eltSize = s;
                }
            } else {
                s = h / n1;
                if ((n2 * s <= w) && (s > eltSize)) { // ca passe
                    nbEltByLine = n2;
                    nbLine = n1;
                    eltSize = s;
                }
            }
            s = w / n2;
            if (n1 * s <= h) { // ca passe
                if (s > eltSize) {
                    nbEltByLine = n2;
                    nbLine = n1;
                    eltSize = s;
                }
            } else {
                s = h / n2;
                if ((n1 * s <= w) && (s > eltSize)) { // ca passe
                    nbEltByLine = n1;
                    nbLine = n2;
                    eltSize = s;
                }
            }
        }
    }

    @Override
    public void componentResized(ComponentEvent e) {
        calculBestSize(getWidth(), getHeight());
        repaint();
    }

    @Override
    public void componentShown(ComponentEvent e) {
        calculBestSize(getWidth(), getHeight());
        repaint();
    }

    void setPalette(Color[] colors, Color selectedColor) {
        palette = colors;
        ColorTools.sort(palette); // on tri la palette
        calculBestSize(getWidth(), getHeight());
        setSelectedValue(selectedColor);
    }

    @Override
    public Color getSelectedValue() {
        return palette[selectedId];
    }

    @Override
    public void setSelectedValue(Color selectedColor) {
        selectedId = 0;
        if (selectedColor != null) {
            for (int id = 0; id < palette.length; id++) {
                if (selectedColor.equals(palette[id])) {
                    selectedId = id;
                    break;
                }
            }
        }
        repaint();
    }
}
