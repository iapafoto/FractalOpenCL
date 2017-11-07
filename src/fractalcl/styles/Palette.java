/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.styles;

import java.awt.Color;
import java.util.List;

/**
 *
 * @author Sebastien Durand
 */
public final class Palette {
    Color[] colors;
    public Double[] values;
    
    public Palette(List<Color> colors, List<Double> values) {
      this.colors = colors.toArray(new Color[colors.size()]);
      this.values = values.toArray(new Double[colors.size()]);
    }
    public Palette(Color[] colors, Double[] values) {
        this.colors = colors;
        this.values = values;
        /*
        double v0 = -Double.MAX_VALUE;
        Interval interval;
        for (int i=0; i<values.length; i++) {
            interval = new Interval(v0, values[i+1]);
            if (v0 != -Double.MAX_VALUE) {
                mapIntervalMax.put(v0, interval);
            }
            v0 = values[i+1];
            mapIntervalMin.put(v0, interval);
        }*/
    }
    
    public Color[] getColors() {
        return colors;
    }
    public Double[] getValues() {
        return values;
    }
    
    public Color getColorForInterval(double v0, double v1) {
        double min = Math.min(v0,v1);
        int index = getIndexForExactValue(min);
        if (index >= 0 && index < colors.length) {
            return colors[index];
        } 
        return null;
    }
    
    // TODO precalculer les vecteurs de couleurs
    public Color getColorForVal(Double val) {
        if (val == null) {
            return null;
        }
        for (int i=1, s=values.length; i<s; i++) {
            if (values[i] >= val) {
                return colors[i-1];
            }
        }
        return colors[colors.length-1];
    }

    
    public Integer getIndexOfPreviousOrEqual(final double val) {
        for (int i=1, s=values.length; i<s; i++) {
            if (values[i] >= val) {
                return i-1;
            }
        }
        return colors.length-1;
    }
    
    public int getIndexForExactValue(double value) {
        // TODO possibilitee d accelerer avec une map ou une recherche dchotomique
        for (int i=0, s=values.length; i<s; i++) {
            if (values[i] == value) {
                return i;
            }
        }
        return -1;
    }
    /*
    Double getOtherInterval(double value, double vOther) {
        int index = getIndexForExactValue(value);
        if (index != -1) {
            if (vOther < value) { // L'autre interval est vers devant
                if (index<values.length-1) {
                    return values[index+1];
                } else {
                    return Double.MAX_VALUE;
                }
            } else { // L'autre interval est vers derriere
                if (index>0) {
                    return values[index-1];
                } else {
                    return -Double.MAX_VALUE;
                }
            }
        }
        return null;
    }
*/
    public Double getNext(double value) {
        int index = getIndexForExactValue(value);
        if (index != -1) {
            if (index<values.length-1)
                return values[index+1];
            else return Double.MAX_VALUE;
        }
        return null;
    }
    
    public Double getPrevious(double value) {
        int index = getIndexForExactValue(value);
        if (index != -1) {
            if (index>0)
                return values[index-1];
            else return -Double.MAX_VALUE;
        }
        return null;
    }


}
