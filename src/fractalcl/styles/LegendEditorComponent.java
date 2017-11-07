/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.styles;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.MouseEvent;

import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import javax.swing.JColorChooser;

/**
 *
 * @author sebastien.durand
 */
public class LegendEditorComponent extends AbstractSliderComponent implements MouseMotionListener, MouseListener {

    public static final String VALUES_CHANGED = "VALUES_CHANGED";    
    public static final String VALUES_CHANGING = "VALUES_CHANGING";
    private static final Color COLOR_DISABLED = new Color(180, 180, 180);
    private int sliderDraggingId = -1;
    private int removeId = -1;
    private int removeY = 0;
    private final List<Double> values = new ArrayList<Double>();
    private final List<Color> colors = new ArrayList<Color>();
    protected boolean isGradiant = false;
    protected boolean drawValue = false;
    protected float mStep = .1f;

    public LegendEditorComponent() {
        this(100, 200, new double[]{130, 160, 190}, new Color[]{Color.green, Color.blue});
    }

    public LegendEditorComponent(double limitMin, double limitMax, double[] values, Color[] colors) {
        super(limitMin, limitMax);
       
        //   setValues(values);

        Dimension dim = new Dimension(50, 48);
        this.setPreferredSize(dim);
        this.setMinimumSize(dim);

        this.addMouseListener(this);
        this.addMouseMotionListener(this);
    }

    @Override
    public void setLimits(double limitMin, double limitMax) {
        setLimits(limitMin, limitMax, true);
    }

    public void setLimits(double limitMin, double limitMax, boolean clipValues) {
        super.setLimits(limitMin, limitMax);

        if (clipValues) {
            // On fait en sorte que les valeurs restent valides
            for (int i = 0; i < values.size(); i++) {
                values.set(i, validateValueInLimits(values.get(i)));
            }
            repaint();
        }
    }

    public void setValues(Collection<Double> lst) {
        // On fait en sorte que les valeurs restent valides
        values.clear();
        for (Double val : lst) {
            values.add(validateValueInLimits(val));
        }
        // Et qu'elles soient dans le bon ordre
        Collections.sort(values);
        repaint();
    }

    @Override
    protected void paintComponent(Graphics g) {
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        int h = getHeight();
        int h2 = h / 2;

        List<Integer> tabX = valuesToPosition();

        // slider
        g2.setColor(Color.white);
        fillSlider(g2);
        for (int sliderId = 0; sliderId < tabX.size() - 1; sliderId++) {
            if (isGradiant && isEnabled()) {
                g2.setPaint(new GradientPaint(tabX.get(sliderId), 0, colors.get(sliderId), tabX.get(sliderId + 1), 0, colors.get(sliderId + 1)));
            } else {
                g2.setColor(getFillColor(sliderId));
            }
            fillRectInSlider(g2, tabX.get(sliderId), tabX.get(sliderId + 1) - tabX.get(sliderId));
        }
        g2.setColor(getDrawColor());
        drawSlider(g2);

        // cursors
        for (int sliderId = getFirstShownSliderId(); sliderId < getLastShownSliderId(); sliderId++) {
            int posY = (sliderId == removeId) ? Math.min(Math.max(CURSOR_W + 2, removeY), h - CURSOR_W - 2) : h2;
            drawCursor(g2, getFillColor(sliderId), getDrawColor(), tabX.get(sliderId), posY);
        }

        if (removeId != -1) {
            int posY = Math.min(Math.max(CURSOR_W + 2, removeY), h - CURSOR_W - 2);
            g2.setColor(Color.red);
            g2.setStroke(new BasicStroke(2f));
            g2.drawLine(tabX.get(removeId) - CURSOR_W / 2 - 1, posY - CURSOR_W / 2 - 1, tabX.get(removeId) + CURSOR_W / 2 + 1, posY + CURSOR_W / 2 + 1);
            g2.drawLine(tabX.get(removeId) - CURSOR_W / 2 - 1, posY + CURSOR_W / 2 + 1, tabX.get(removeId) + CURSOR_W / 2 + 1, posY - CURSOR_W / 2 - 1);
        }

        // Si le pas est compris entre 0 et 1, on affiche les décimales
        String valueAsString;
        for (int sliderId = 0, s = tabX.size(); sliderId < s; sliderId++) {
            if (sliderId == removeId) {
                continue;
            }

            if (sliderDraggingId != -1 && sliderId != sliderDraggingId) {
                g2.setColor(COLOR_DISABLED);
            } else {
                g2.setColor(getDrawColor());
            }
            if (mStep >= 0 && mStep < 1) {
                valueAsString = Float.toString((int) (values.get(sliderId) * 10) / 10f);
            } else {
                valueAsString = Integer.toString((int) (Math.floor(values.get(sliderId))));
            }
            printStringOnCursor(g2, valueAsString, tabX.get(sliderId), sliderId % 2 == 0);
        }
    }

    private Color getDrawColor() {
        return isEnabled() ? Color.black : COLOR_DISABLED;
    }

    private Color getFillColor(final int noColor) {
        if (this.isEnabled() && !colors.isEmpty()) {
            return noColor < colors.size() ? colors.get(noColor) : colors.get(colors.size() - 1);
        }
        return this.getBackground();
    }

    private int getFirstShownSliderId() {
        return values.isEmpty() ? 0 : 1;
    }

    private int getLastShownSliderId() {
        return values.size() - 1;
    }

    @Override
    public void mouseDragged(MouseEvent mouseEvent) {
        if (!isEnabled()) {
            return;
        }

        if (sliderDraggingId == -1) { // Start drag
            List<Integer> tabX = valuesToPosition();
            int dMin = Integer.MAX_VALUE;
            int bestId = -1;
            for (int sliderId = getFirstShownSliderId(); sliderId < getLastShownSliderId(); sliderId++) {
                if (Math.abs(mouseEvent.getX() - tabX.get(sliderId)) < dMin) {
                    dMin = Math.abs(mouseEvent.getX() - tabX.get(sliderId));
                    bestId = sliderId;
                }
            }
            if ((dMin > CURSOR_W * 2) && (bestId != -1)) {
                // On cree un nouveau point pres de  bestId
                if (mouseEvent.getX() > tabX.get(bestId)) {
                    bestId = createValueAt(bestId);
                } else {
                    bestId = createValueAt(bestId - 1);
                }

            }
            sliderDraggingId = bestId;
        }

        if (sliderDraggingId != -1) { // Drag
            if (sliderDraggingId != 0 && sliderDraggingId != values.size() - 1 && // imposible d'enlever la premiere et la derniere valeur
                    Math.abs(mouseEvent.getPoint().getY() - getHeight() / 2) > 12) {
                removeId = sliderDraggingId;
                removeY = (int) mouseEvent.getPoint().getY();
            } else {
                removeId = -1;
            }
            setValue(sliderDraggingId, positionToValue(mouseEvent.getX()));
            firePropertyChange(VALUES_CHANGING, null, values);
            repaint();
        }
    }

    public void setValue(int id, double value) {
        values.set(id, validateValueInLimits(value));
        // Le moyen le plus simple pour deplacer les valeurs precedentes et suivantes (l'id reste celui de la valeur n)
        Collections.sort(values);
    }

    public List<Double> getValues() {
        return values;
    }

    public Palette getColorMap() {
        return new Palette(colors, values);
//        // Colors
//        List<String> colorNames = new ArrayList<String>();
//        String name;
//        for (Color c : colors) {
//            name = null;
//            if (palette != null) {
//                name = palette.getName(c);
//            }
//            if (name == null) {
//                name = ColorTools.toHtmlColor(c);
//            }
//            colorNames.add(name);
//        }
//        return new SSColorMap("", values, colorNames);
    }

    public void setColorMap(Palette palette) {
//        // Values
//        List<Double> lstValue = palette.getLevels();
//
//        // Colors conversion de texte a Couleur
//        ImmutableList<String> lstColorName = palette.getColorNames();
//        List<Color> newColors = new ArrayList<Color>();
//        Color c;
//        for (String name : lstColorName) {
//            c = null;
//            if (name.startsWith("#")) { // Couleur HTML
//                c = ColorTools.fromHtmlColor(name);
//            } else if (palette != null) { // Couleur dans la palette
//                c = palette.getColor(name);
//            }
//            if (c == null) {
//                c = Color.red;
//            }
//            newColors.add(c);
//        }

        setValues(Arrays.asList(palette.getValues()), Arrays.asList(palette.getColors()));
    }

    @Override
    public void mousePressed(MouseEvent e) {
    }

    @Override
    public void mouseReleased(MouseEvent e) {
        if (!isEnabled()) {
            return;
        }

        if (sliderDraggingId != -1) {
            if (removeId != -1) {
                removeValue(removeId);
            }
            sliderDraggingId = -1;
            removeId = -1;
            repaint();
            onValueChanged();
        }
    }

    private void removeValue(int id) {
        values.remove(Math.min(id, values.size() - 1));
        colors.remove(Math.min(id, colors.size() - 1));
    }

    private int createValueAt(int id) {
        int idVal = Math.min(Math.max(1, id + 1), values.size() - 1);

        values.add(idVal, (values.get(idVal - 1) + values.get(idVal)) / 2f);

        if (colors.size() >= 2) {
            int idColor = Math.min(Math.max(1, id + 1), colors.size() - 1);
            colors.add(idColor, createColor(colors.get(idColor - 1), colors.get(idColor)));
        } else if (!colors.isEmpty()) {
            colors.add(colors.get(0));
        } else {
            colors.add(Color.green);
        }

        return idVal;
    }

    // Retourne une couleur 
    private Color createColor(final Color c1, final Color c2) {
        Color cMoy = ColorTools.gradiant(c1, c2, .5);
        return cMoy;
    }

    private void onValueChanged() {
        this.firePropertyChange(VALUES_CHANGED, null, values);
    }

    private List<Integer> valuesToPosition() {
        List<Integer> pos = new ArrayList<Integer>();
        for (Double v : values) {
            pos.add(valueToPosition(v));
        }
        return pos;
    }

    @Override
    protected double positionToValue(int pos) {
        double result = super.positionToValue(pos);
        if (mStep > 0) {
            result -= result % mStep;
        }
        return result;
    }

    @Override
    public void mouseMoved(MouseEvent e) {
    }

    @Override
    public void mouseClicked(MouseEvent e) {
        if (!isEnabled()) {
            return;
        }

        if (e.getButton() == MouseEvent.BUTTON1 && e.getClickCount() == 2) {
            List<Integer> tabX = valuesToPosition();
            for (int i = 0, s = tabX.size() - 1; i < s; i++) {
                if (e.getX() > tabX.get(i) && e.getX() <= tabX.get(i + 1)) {
                    editColor(i);
                    break;
                }
            }
        }
    }

    private void editColor(final int id) {
        Color newColor = null;
        newColor = JColorChooser.showDialog(this, "Choose color", colors.get(id));
        
        if (newColor != null) {
            colors.set(id, newColor);
            repaint();
            onValueChanged();
        }
    }

    @Override
    public void mouseEntered(MouseEvent e) {
    }

    @Override
    public void mouseExited(MouseEvent e) {
    }

    public void setStep(float step) {
        mStep = step;
    }

    public void setValues(Collection<Double> values, Collection<Color> colors) {
        this.colors.clear();
        this.colors.addAll(colors);
        setValues(values);
    }

    public void setGradiant(boolean gradiant) {
        this.isGradiant = gradiant;
    }
}
