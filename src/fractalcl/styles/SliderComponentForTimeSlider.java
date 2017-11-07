package fractalcl.styles;


import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Rectangle2D;

/**
 *
 * @author sebastien.durand
 */
public class SliderComponentForTimeSlider extends AbstractSliderComponentForTimeSlider implements MouseMotionListener, MouseListener {

    public static final String VALUE_CHANGED = "VALUE_CHANGED";
    public static final String VALUE_CHANGING = "VALUE_CHANGING";
    protected Color colorLimitMin = new Color(0, 255, 255);
    protected Color colorLimitMax = Color.red;
    protected boolean drawValue = false;
    protected double value = 160;
    protected boolean sliderDragging = false;
    /** Le nbre de décimales des valeurs manipulées */
    protected int mPrecision = 0;
    /** permet de n'envoyer que les changements effectif d'etat lors du drag du slider */
    protected Object 
            draggingValueMem,
            draggedValueMem;

    private long slidingDelay = 0;

    
    public SliderComponentForTimeSlider() {
        this(100, 200, 160, new Color(200, 255, 210), new Color(200, 210, 255));
    }

    public SliderComponentForTimeSlider(final double limitMin, final double limitMax, final double value) {
        this(limitMin, limitMax, value, new Color(200, 255, 210), new Color(200, 210, 255));
    }

    public SliderComponentForTimeSlider(final double limitMin, final double limitMax, final double value, final Color colorMin, final Color colorMax) {
        this(100, limitMin, limitMax, colorMin, colorMax);
        setValue(value);
    }

    public SliderComponentForTimeSlider(final int slidingDelay, final double limitMin, final double limitMax, final Color colorMin, final Color colorMax) {
        super(limitMin, limitMax);
        setValue(value);
        colorLimitMin = colorMin;
        colorLimitMax = colorMax;

        Dimension dim = new Dimension(50, 48);
        this.setPreferredSize(dim);
        this.setMinimumSize(dim);
        this.addMouseListener(this);
        this.addMouseMotionListener(this);
        this.slidingDelay = slidingDelay;
                
 //       timerSliding = new Timer(slidingDelay, actionSliding);
 //       timerSliding.setRepeats(false);
    }

    public void setPaintLabels(final boolean drawValue) {
        this.drawValue = drawValue;
    }

    public void setPrecision(final int precision) {
        mPrecision = precision;
    }

    @Override
    public void setLimits(final double limitMin, final double limitMax) {
        super.setLimits(limitMin, limitMax);

        // On fait en sorte que les valeurs restent valide
        this.value = validateValueInLimits(value);
        repaint();
    }

    public double getValue() {
        return value;
    }

    @Override
    protected void paintComponent(Graphics g) {
        Graphics2D g2 = (Graphics2D) g.create();
     //   g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
     //   g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        drawSlider(g2);
        int cursorX = valueToPosition(value);
        drawCursor(g2, cursorX);

        if (drawValue) {
            printStringOnCursor(g2, valueAsString(value), cursorX, true);
        }
        g2.dispose();
    }

    @Override
    protected void drawSlider(final Graphics2D g2) {
        GradientPaint gradient = new GradientPaint(0, 0, colorLimitMin, sliderWidth(), 0, colorLimitMax, true);
        g2.setPaint(gradient);
        fillSlider(g2);
        g2.setColor(Color.black);
        super.drawSlider(g2);
    }
    
    protected void drawCursor(final Graphics2D g2, int cursorX) {
        Color fillColor = Color.red; //ColorTools.mix(colorLimitMin, colorLimitMax, valueAsCoef(value));
        drawCursor(g2, fillColor, Color.black, cursorX);
    }
    
    protected String valueAsString(double value) {
        if (Math.round(value) == value) {
            return Integer.toString((int)value);
        }
        return Double.toString(Math.round(value * 10) / 10f);
    }
    
    @Override
    protected double positionToValue(final int pos) {
        double result = super.positionToValue(pos);
        if (mPrecision > 0) {
            final double factor = Math.pow(10, mPrecision);
            result = (Math.round(result * factor) / factor);
        } else if (mPrecision < 0) {
            final double factor = Math.pow(10, -mPrecision);
            result = (Math.round(result / factor) * factor);
        } else { // int
            result = Math.round(result);
        }
        return result;
    }
        
    public void setValue(final double value) {
        this.value = validateValueInLimits(value);
        repaint();
    }

    @Override
    public void mouseDragged(MouseEvent e) {
        if (sliderDragging) {
            setValue(positionToValue(e.getX()));
            onDragging();
            repaint();
        }
    }

    @Override
    public void mousePressed(MouseEvent e) {
        sliderDragging = true;
        draggedValueMem = getValueForEvent();
        draggingValueMem = draggedValueMem;
    }

    @Override
    public void mouseReleased(MouseEvent e) {
        if (sliderDragging) {
            final Object dragValue = getValueForEvent();
       //     timerSliding.stop();
            if (draggedValueMem !=null && !draggedValueMem.equals(dragValue)) {
                this.firePropertyChange(VALUE_CHANGED, draggedValueMem, dragValue);
            }
            draggedValueMem = draggingValueMem = null;
            sliderDragging = false;
        }
    }

    long timeOfLastDraggingEvent = 0;

    /**
     * En cours de drag, on n'envoi pas les evenements trop raprochees dans le temps (slidingDelay)
     */
    public void onDragging() {
        final Object newValue = getValueForEvent(); 
        // Si la valeur a change on envoi l'evenement
        if (draggingValueMem == null || !draggingValueMem.equals(newValue)) {
            final long t = System.nanoTime() / 1000000;
            if (t - timeOfLastDraggingEvent > slidingDelay) {
                firePropertyChange(VALUE_CHANGING, draggingValueMem, newValue);
                draggingValueMem = newValue;
                timeOfLastDraggingEvent = t;
            }
        }
    }

    @Override
    public void mouseClicked(MouseEvent e) {
        double newValue = positionToValue(e.getX());
        if (newValue != value) {      
            final Object oldValue = getValueForEvent();
            setValue(newValue);
            this.firePropertyChange(VALUE_CHANGED, oldValue, getValueForEvent());
        }
    }

    @Override
    public void mouseMoved(MouseEvent e) {
    }

    @Override
    public String getToolTipText(MouseEvent e) {
        return valueAsString(this.positionToValue(e.getX()));
    }
    
    @Override
    public Point getToolTipLocation(MouseEvent event) {
        final String txt = getToolTipText(event);
        if (txt != null) {
            FontMetrics fm = this.getFontMetrics(getFont());
            Rectangle2D r = fm.getStringBounds(txt, getGraphics());
            return new Point(event.getX() - (int)r.getWidth()/2, -22);
        }
        return null;
    }
    
    @Override
    public void mouseEntered(MouseEvent e) {
    }

    @Override
    public void mouseExited(MouseEvent e) {
    }

    protected Object getValueForEvent() {
        return value;
    }



}
