package fractalcl.styles;

import java.awt.Color;
import java.awt.Graphics2D;
import javax.swing.JComponent;

/**
 *
 * @author abordier
 */
public abstract class AbstractSliderComponent_1 extends JComponent {

    /** la largeur d'un cursseur */
    protected static final int CURSOR_W = 6;
    /** La hauteur d'un cursseur */
    protected static final int CURSOR_H = CURSOR_W * 2;
    /** La hauteur du slider */
    public static final int SLIDER_H = CURSOR_W;
    /** La position X du slider */
    protected static final int SLIDER_X = CURSOR_W;
    /** La taille de l'arc de cerle */
    public static final int ARC_SIZE = 4;
    private double limitMin = 100;
    private double limitMax = 200;
    private boolean isLogarithmic = false;
    protected int borderLeftRight = CURSOR_W;
    
    public AbstractSliderComponent_1(final double limitMin, final double limitMax) {
        innerSetLimit(limitMin, limitMax);
        
    }

    private void innerSetLimit(final double limitMin, final double limitMax) {
        this.limitMin = Math.min(limitMin, limitMax);
        this.limitMax = Math.max(limitMin, limitMax);
    }

    public void setLimits(final double limitMin, final double limitMax) {
        innerSetLimit(limitMin, limitMax);
    }

    public void setLogarithmic(final boolean isLogarithmic) {
        this.isLogarithmic = isLogarithmic;
    }

    protected double validateValueInLimits(final double value) {
        return Math.min(Math.max(limitMin, value), limitMax);
    }

    protected int valueToPosition(final double value) {
        return borderLeftRight + (int) (sliderWidth() * valueAsCoef(value));
    }

    protected double valueAsCoef(final double value) {
        double coef;
        if (isLogarithmic) {
            coef = Math.log(1 + value - limitMin) / Math.log(1 + limitMax - limitMin);
        } else {
            coef = (value - limitMin) / (limitMax - limitMin);
        }
        return coef;
    }

    protected double positionToValue(final int pos) {
        int ws = sliderWidth();
        double result;
        if (isLogarithmic) {
            result = (limitMin - 1 + (Math.exp((Math.log(1 + limitMax - limitMin)) * ((double) (pos - borderLeftRight) / ws))));
        } else {
            result = (limitMin + ((double) ((pos - borderLeftRight) * (limitMax - limitMin)) / ws));
        }
        return result;
    }

    /** la largeur du slider */
    protected int sliderWidth() {
        return getWidth() - borderLeftRight * 2 - 2;
    }

    /** La position Y du slider */
    protected int sliderY() {
        return (getHeight() / 2) - (SLIDER_H / 2);
    }

    /** dessine un cursseur dont le centre en X est centerX */
    protected void drawCursor(final Graphics2D g2, final Color fillColor, final Color drawColor, final int centerX) {
        drawCursor(g2, fillColor, drawColor, centerX, getHeight() / 2);
    }

    /** dessine un cursseur dont le centre en X est centerX, et le centre en Y est centerY */
    public static void drawCursor(final Graphics2D g, final Color fillColor, final Color drawColor, final int centerX, final int centerY) {
        Graphics2D g2 = (Graphics2D) g.create();
        g2.setColor(fillColor);
        g2.fillRoundRect(centerX - CURSOR_W / 2, centerY - CURSOR_H / 2, CURSOR_W, CURSOR_H, ARC_SIZE, ARC_SIZE);

        g2.setColor(drawColor);
        g2.drawRoundRect(centerX - CURSOR_W / 2, centerY - CURSOR_H / 2, CURSOR_W, CURSOR_H, ARC_SIZE, ARC_SIZE);
        g2.dispose();
    }

    protected void fillSlider(final Graphics2D g2) {
        g2.fillRoundRect(borderLeftRight, sliderY(), sliderWidth(), SLIDER_H, ARC_SIZE, ARC_SIZE);
    }

    protected void drawSlider(final Graphics2D g2) {
        g2.drawRoundRect(borderLeftRight, sliderY(), sliderWidth(), SLIDER_H, ARC_SIZE, ARC_SIZE);
    }

    protected void fillRectInSlider(final Graphics2D g2, final int x, final int width) {
        g2.fillRect(x, sliderY(), width, SLIDER_H);
    }

    /**
     * Affiche une string sur le cursseur
     * @param g2 l'objet garphique
     * @param string la chaine à afficher
     * @param x la position du cursseur en X
     * @param onTop si <code>true</code> la chaine est affichée sur le cursseur, si <code>false</code> la chaine est affichée sous le cursseur
     */
    protected void printStringOnCursor(final Graphics2D g2, final String string, final int x, final boolean onTop) {
        int y;
        if (onTop) {
            y = sliderY() - (CURSOR_H / 2);
        } else {
            y = sliderY() + (CURSOR_H * 2);
        }
        y--;

        int stringLen = (int) g2.getFontMetrics().getStringBounds(string, g2).getWidth();
        int newX = Math.min(Math.max(CURSOR_W / 2, x - (stringLen / 2)), getWidth() - stringLen - (CURSOR_W / 2));

        g2.drawString(string, newX, y);
    }

    public double getMinimum() {
        return limitMin;
    }

    public double getMaximum() {
        return limitMax;
    }
}
