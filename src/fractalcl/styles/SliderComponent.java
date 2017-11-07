package fractalcl.styles;



import java.awt.Color;
import java.awt.Dimension;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

/**
 *
 * @author sebastien.durand
 */
public class SliderComponent extends AbstractSliderComponent implements MouseMotionListener, MouseListener {

    public static final String VALUE_CHANGED = "VALUE_CHANGED";
    public static final String VALUE_CHANGING = "VALUE_CHANGING";
    protected Color colorLimitMin = new Color(0, 255, 255);
    protected Color colorLimitMax = Color.red;
    protected boolean drawValue = false;
    private float value = 160;
    private boolean sliderDragging = false;
    /** Le nbre de décimales des valeurs manipulées */
    protected int mPrecision = -4;

    public SliderComponent() {
        this(100, 200, 160, new Color(200, 255, 210), new Color(200, 210, 255));
    }

    public SliderComponent(double limitMin, double limitMax, float value) {
        this(limitMin, limitMax, value, new Color(200, 255, 210), new Color(200, 210, 255));
    }

    public SliderComponent(double limitMin, double limitMax, float value, Color colorMin, Color colorMax) {
        this(limitMin, limitMax, colorMin, colorMax);
        setValue(value);
    }

    public SliderComponent(double limitMin, double limitMax, Color colorMin, Color colorMax) {
        super(limitMin, limitMax);
        setValue(value);
        colorLimitMin = colorMin;
        colorLimitMax = colorMax;

        Dimension dim = new Dimension(50, 48);
        this.setPreferredSize(dim);
        this.setMinimumSize(dim);

        this.addMouseListener(this);
        this.addMouseMotionListener(this);
    }

    public void setPrecision(int precision) {
        mPrecision = precision;
    }

    public void setLimits(float limitMin, float limitMax) {
        super.setLimits(limitMin, limitMax);

        // On fait en sorte que les valeurs restent valide
        this.value = (float) validateValueInLimits(value);
        repaint();
    }

    public float getValue() {
        return value;
    }

    @Override
    protected void paintComponent(Graphics g) {
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        int cursorX = valueToPosition(value);

        GradientPaint gradient = new GradientPaint(0, 0, colorLimitMin, sliderWidth(), 0, colorLimitMax, true);

        // slider
        g2.setPaint(gradient);
        fillSlider(g2);
        g2.setColor(Color.black);
        drawSlider(g2);

        // cursor
        Color fillColor = ColorTools.gradiant(colorLimitMin, colorLimitMax, valueAsCoef(value));
        drawCursor(g2, fillColor, Color.black, cursorX);

        if (drawValue) {
           // final String valueAsString = getValueAsString(); //Float.toString(Math.round(value * 10) / 10f);
            printStringOnCursor(g2, getValueAsString(), cursorX, true);
        }
    }

    @Override
    public void mouseDragged(MouseEvent e) {
        if (sliderDragging) {
            setValue((float) positionToValue(e.getX()));
            firePropertyChange(VALUE_CHANGING, null, value);
            repaint();
        }
    }

    @Override
    protected double positionToValue(final int pos) {
        double result = super.positionToValue(pos);
        if (mPrecision > 0) {
            final double factor = Math.pow(10, mPrecision);
            result = (Math.round(result * factor) / factor);
        } else if (mPrecision < 0) {
            final double factor = Math.pow(10, mPrecision);
            result = (Math.round(result / factor) * factor);
        } else { // int
            result = Math.round(result);
        }
        return result;
    }
    
    public void setPaintLabels(final boolean drawValue) {
        this.drawValue = drawValue;
    }

    public void setValue(float value) {
        this.value = (float) validateValueInLimits(value);
        repaint();
    }

    public String getValueAsString() {
        double result = getValue();
        
        if (mPrecision > 0) {
            final double factor = Math.pow(10, mPrecision);
            return Double.toString((Math.round(result * factor) / factor));
        } else if (mPrecision < 0) {
            final double factor = Math.pow(10, -mPrecision);
            return Integer.toString((int)(Math.round(result / factor) * factor));
        } else { // int
            return Integer.toString((int)Math.round(result));
        }
    }

    @Override
    public void mousePressed(MouseEvent e) {
        sliderDragging = true;
    }

    @Override
    public void mouseReleased(MouseEvent e) {
        if (sliderDragging) {
            this.firePropertyChange(VALUE_CHANGED, null, value);
            sliderDragging = false;
        }
    }

    @Override
    public void mouseMoved(MouseEvent e) {
    }

    @Override
    public void mouseClicked(MouseEvent e) {
    }

    @Override
    public void mouseEntered(MouseEvent e) {
    }

    @Override
    public void mouseExited(MouseEvent e) {
    }
}
