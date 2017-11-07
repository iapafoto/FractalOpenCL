/*
 * 
 * Meteo France 2012
 */
package fractalcl.styles;

import java.awt.BasicStroke;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.GeneralPath;
import javax.swing.Icon;

/**
 * Pour montrer les épaisseurs avec des jolis lignes courbes (repompe du composant fait par Seb)
 * @author fguerry
 */
public class ThicknessIcon implements Icon {

    private final float thickness;

    public ThicknessIcon(float thickness) {
        this.thickness = thickness;
    }

    public float getThickness() {
        return thickness;
    }

    @Override
    public void paintIcon(Component c, Graphics g, int x, int y) {
        int w = getIconWidth();
        int h = getIconHeight();

        GeneralPath path = new GeneralPath();
        path.moveTo(0, h / 2);
        double segment = w / 3d;
        path.curveTo(segment, -h * 0.3, 2 * segment, h * 1.3, w, h / 2);

        Graphics2D g2 = (Graphics2D) g.create();
        g2.translate(x, y);
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);
        g2.setColor(c.getForeground());
        g2.setStroke(new BasicStroke(thickness));
        g2.draw(path);
        g2.dispose();
    }

    @Override
    public int getIconWidth() {
        return 120;
    }

    @Override
    public int getIconHeight() {
        return 12;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final ThicknessIcon other = (ThicknessIcon) obj;
        if (Float.floatToIntBits(this.thickness) != Float.floatToIntBits(other.thickness)) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 67 * hash + Float.floatToIntBits(this.thickness);
        return hash;
    }
}
