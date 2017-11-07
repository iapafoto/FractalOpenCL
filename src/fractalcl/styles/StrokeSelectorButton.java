/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.styles;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.GeneralPath;

/**
 *
 * @author sebastien.durand
 */
public class StrokeSelectorButton extends SelectorButton<Stroke, StrokeSelector> {
    
    public StrokeSelectorButton() {  
        super(100,100);
        // Mode par defaut
        Stroke[] strokes = { new BasicStroke(1f),new BasicStroke(2f),new BasicStroke(3f),new BasicStroke(4f)};
        setAvailableStroke(strokes, strokes[0]);
        init();
    }   
    
    public StrokeSelectorButton(final Stroke[] tab, final Stroke selected) {
        super(100,100);
        setAvailableStroke(tab, selected);
        init();
    }
    
    public final void setAvailableStroke(final Stroke[] tab, final Stroke selected) {
        getSelectorPanel().setAvailableStroke(tab, selected);
    }
    
    public void selectStroke(final float lineWidth) {
        getSelectorPanel().selectedValueFromLineWidth(lineWidth);
    }   
    public void selectStroke(final float[] dashArray) {
        getSelectorPanel().selectedValueFromDashArray(dashArray);
    }   
    public void selectStroke(final Stroke selected) {
        getSelectorPanel().setSelectedValue(selected);
    }    
    public Stroke getStroke() {
        return getSelectorPanel().getSelectedValue();
    }
    
    @Override
    public void setForeground(final Color c) {
        if (selectorPanel != null) {
            selectorPanel.setForeground(c);
        }
        super.setForeground(c);
    }
    
    @Override
    protected void onSelectionChanged(final Stroke newValue) {
    }

    @Override
    protected void paintSelection(Graphics2D g2, int w, int h) {
        GeneralPath path = new GeneralPath();
        w = w-16;
        path.moveTo(8, h/2);
        path.curveTo(8+w/3, 6, 8+(2*w)/3, h, 8+w, h/2);

        g2.setColor(selectorPanel.getForeground());
        g2.setStroke(getSelectedValue());
        g2.draw(path);
    }

    @Override
    protected StrokeSelector createSelector() {
        return new StrokeSelector();
    }
}
