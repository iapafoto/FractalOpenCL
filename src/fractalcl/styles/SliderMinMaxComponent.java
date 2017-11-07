package fractalcl.styles;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import java.awt.Color;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.MouseEvent;

/**
 *
 * @author sebastien.durand
 */
public class SliderMinMaxComponent extends SliderComponent {

    public static final String MIN_VALUE = "MIN_VALUE";
    public static final String MAX_VALUE = "MAX_VALUE";
    private float valueMin = 130;
    private float valueMax = 160;
    private boolean sliderMinDragging = false;
    private boolean sliderMaxDragging = false;
    private boolean isMinFixed = false;  // Permet de fixer la valeur min (changement impossible)
    private boolean isMaxFixed = false;  // Permet de fixer la valeur max (changement impossible)

    public SliderMinMaxComponent() {
        this(100, 200, 130, 160, Color.cyan, Color.red);
    }

    public SliderMinMaxComponent(float limitMin, float limitMax, float valueMin, float valueMax) {
        this(limitMin, limitMax, valueMin, valueMax, Color.cyan, Color.red);
    }

    public SliderMinMaxComponent(float limitMin, float limitMax, float valueMin, float valueMax, Color colorMin, Color colorMax) {
        super(limitMin, limitMax, colorMin, colorMax);
        setValues(valueMin, valueMax);
    }

    
    @Override
    public void setLimits(float limitMin, float limitMax) {
        super.setLimits(limitMin, limitMax);
        // On fait en sorte que les valeurs restent valide

        this.valueMin = (float) validateValueInLimits(valueMin);
        this.valueMax = (float) validateValueInLimits(valueMax);
        repaint();
    }

    public void setValues(float valueMin, float valueMax) {
        if (valueMin > valueMax) {
            float mem = valueMin;
            valueMin = valueMax;
            valueMax = mem;
        }
        this.valueMin = (float) validateValueInLimits(valueMin);
        this.valueMax = (float) validateValueInLimits(valueMax);
        repaint();
    }

    @Override
    protected void paintComponent(Graphics g) {
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        int xMin = valueToPosition(valueMin);
        int xMax = valueToPosition(valueMax);

        GradientPaint gradient = new GradientPaint(0, 0, colorLimitMin, sliderWidth(), 0, colorLimitMax, true);

        // slider
        g2.setColor(Color.white);
        fillSlider(g2);
        g2.setPaint(gradient);
        fillRectInSlider(g2, xMin, xMax - xMin);
        g2.setColor(Color.black);
        drawSlider(g2);

        // cursor
        Color fillColor = ColorTools.gradiant(colorLimitMin, colorLimitMax, valueAsCoef(valueMin));
        drawCursor(g2, fillColor, Color.black, xMin);

        fillColor = ColorTools.gradiant(colorLimitMin, colorLimitMax, valueAsCoef(valueMax));
        drawCursor(g2, fillColor, Color.black, xMax);

        // Si le nombre est entier, on fait sauter la virgule
        final String minValueAsString;
        final String maxValueAsString;
        if (valueMin % 1f == 0) {
            minValueAsString = Integer.toString((int) valueMin);
        } else {
            minValueAsString = Float.toString(valueMin);
        }
        if (valueMax % 1f == 0) {
            maxValueAsString = Integer.toString((int) valueMax);
        } else {
            maxValueAsString = Float.toString(valueMax);
        }
        printStringOnCursor(g2, minValueAsString, xMin, true);
        printStringOnCursor(g2, maxValueAsString, xMax, false);
    }

    @Override
    public void mouseDragged(MouseEvent e) {
        if (!isMinFixed && sliderMinDragging) {
            setValueMin((float) positionToValue(e.getX()));
            repaint();
            firePropertyChange(VALUE_CHANGING, null, new Float[] {valueMin, valueMax});
        } else if (!isMaxFixed && sliderMaxDragging) {
            setValueMax((float) positionToValue(e.getX()));
            repaint();
            firePropertyChange(VALUE_CHANGING, null, new Float[] {valueMin, valueMax});
        }
    }

    public void setValueMin(float min) {
        valueMin = (float) validateValueInLimits((double) min);
        if (valueMin > valueMax) {
            valueMax = valueMin;
        }
    }

    public void setValueMax(float max) {
        valueMax = (float) validateValueInLimits((double) max);
        if (valueMax < valueMin) {
            valueMin = valueMax;
        }
    }

    @Override
    public void mousePressed(MouseEvent e) {
        if (Math.abs(e.getX() - valueToPosition(valueMin)) < Math.abs(e.getX() - valueToPosition(valueMax))) {
            if (!isMinFixed) {
                sliderMinDragging = true;
            }
        } else {
            if (!isMaxFixed) {
                sliderMaxDragging = true;
            }
        }
    }

    @Override
    public void mouseReleased(MouseEvent e) {
        if (!isMinFixed && sliderMinDragging) {
            this.firePropertyChange(MIN_VALUE, null, valueMin);
            this.firePropertyChange(VALUE_CHANGED, null, new Float[] {valueMin, valueMax});
            sliderMinDragging = false;
        }
        if (!isMaxFixed && sliderMaxDragging) {
            this.firePropertyChange(MAX_VALUE, null, valueMax);
            this.firePropertyChange(VALUE_CHANGED, null, new Float[] {valueMin, valueMax});
            sliderMaxDragging = false;
        }
    }
    
    public void setMinFixed(boolean minFixed) {
        this.isMinFixed = minFixed;
    }

    public void setMaxFixed(boolean maxFixed) {
        this.isMaxFixed = maxFixed;
    }

    public float getStep() {
        return mPrecision;
    }

    public float getValueMin() {
        return valueMin;
    }
    
    public float getValueMax() {
        return valueMax;
    }
}
