/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.styles;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 *
 * @author sebastien.durand
 */
public class ColorTools {
    // Couleur invisible (transparence à fond)
    public final static Color INVISIBLE = new Color(64, 64, 64, 0);

    
    public static boolean isInvisible(Color c) {
        return (c == null || c.getAlpha() == 0);
    }

    public static boolean isVisible(Color c) {
        return (c != null && c.getAlpha() != 0);
    }
    
    private static int _byte(int v) {
        return Math.min(Math.max(0, v), 255);
    }

    private static int _byte(double v) {
        return _byte((int) v);
    }

    
    /** Determination de la couleur sans transparence obtenue par l'application d'une couleur transparente sur une autre */
    public static Color applyColorOn(Color c0, Color c1) {
        if (c0 == null) return c1;
        if (c1 == null) return c0;
        float   k = ((float)(c0.getAlpha()) / 255f),
                k0 = 1f - k;
        return new Color(
                _byte((int)(k0*(float)c1.getRed()   + k*(float)c0.getRed())),
                _byte((int)(k0*(float)c1.getGreen() + k*(float)c0.getGreen())),
                _byte((int)(k0*(float)c1.getBlue()  + k*(float)c0.getBlue()))); 
    }
    
    public static Color brighter(Color c, int k) {
        return new Color(Math.min(255, c.getRed() + k), Math.min(255, c.getGreen() + k), Math.min(255, c.getBlue() + k));
    }

    public static Color darker(Color c, int k) {
        return new Color(Math.max(0, c.getRed() - k), Math.max(0, c.getGreen() - k), Math.max(0, c.getBlue() - k));
    }

    public static Color alpha(Color c, double a) {
        return new Color(c.getRed(), c.getGreen(), c.getBlue(), Math.min(255, Math.max(0, (int)(255.*a))));
    }
        
    public static Color alpha(Color c, int a) {
        return new Color(c.getRed(), c.getGreen(), c.getBlue(), Math.min(255, Math.max(0, a)));
    }

    public static Color brightnessContrast(Color c) {
        float[] hsb = new float[3];
        Color.RGBtoHSB(c.getRed(), c.getGreen(), c.getBlue(), hsb);
        return new Color(Color.HSBtoRGB(hsb[0], hsb[1], ((hsb[2] > .5f) ? 0.0f : 1.0f)));
    }

    /** 
     * Tris les couleur par Hue, Saturation puis Brightness
     * Attention : modifie la liste passée en paramétre 
     **/
    public static List<Color> sortColorByHSB(List<Color> lstColor) {
        Collections.sort(lstColor, new ComparatorColorHSB());
        return lstColor;
    }

    public static Color[] sortColorByHSB(Color[] colors) {
        Arrays.sort(colors, new ComparatorColorHSB());
        return colors;
    }

    private static class ComparatorColorHSB implements Comparator<Color> {

        @Override
        public int compare(Color c1, Color c2) {

            int r1 = c1.getRed(), g1 = c1.getGreen(), b1 = c1.getBlue();
            int r2 = c2.getRed(), g2 = c2.getGreen(), b2 = c2.getBlue();

            if ((r1 == g1) && (g1 == b1)) {     // Si c1 est gris
                if ((r2 == g2) && (g2 == b2)) { // si c2 est gris aussi
                    if (r1 > r2) {
                        return 1;
                    }
                    if (r1 < r2) {
                        return -1;
                    }
                    return 0;
                }
                return 1;
            }

            float[] hsb1 = new float[3];
            float[] hsb2 = new float[3];

            Color.RGBtoHSB(r1, g1, b1, hsb1);
            Color.RGBtoHSB(r2, g2, b2, hsb2);

            if (hsb1[0] > hsb2[0]) {
                return 1;
            }
            if (hsb1[0] < hsb2[0]) {
                return -1;
            }
            if (hsb1[2] > hsb2[2]) {
                return 1;
            }
            if (hsb1[2] < hsb2[2]) {
                return -1;
            }
            if (hsb1[1] > hsb2[1]) {
                return 1;
            }
            if (hsb1[1] < hsb2[1]) {
                return -1;
            }
            return 0;
        }
    }

    /**
     * Convertit une couleur au format héxadécimal
     * @param color la couleur à convertir
     * @return le code héxadécimal de la couleur au format transparency-red-green-blue : TTRRGGBB
     */
    public static String convertToHexa(Color color) {
        return Integer.toHexString(color.getRGB());
    }
    
    /**
     * Renvoie une couleur intermédiaire entre c0 et c1, positionnée à k1
     * entre les deux. k1 est compris entre 0 et 1. Si k1 vaut 0, la couleur
     * renvoiyée est c0, et si k1 vaut 1 la couleur renvoyée vaut c1.
     * @param c0
     * @param c1
     * @param k1
     * @return 
     */
    public static Color gradiant(final Color c0, final Color c1, double k1) {
        if (c0 == null) {
            return c1;
        }
        if (c1 == null) {
            return c0;
        }
        k1 = Math.min(Math.max(0., k1), 1.);
        double k0 = 1. - k1;
        return new Color(_byte(k0 * c0.getRed() + k1 * c1.getRed()),
                _byte(k0 * c0.getGreen() + k1 * c1.getGreen()),
                _byte(k0 * c0.getBlue() + k1 * c1.getBlue()),
                _byte(k0 * c0.getAlpha() + k1 * c1.getAlpha()));
    }

    public static Color gradiant(final Double[] tab, final Color[] c, final double v) {
        if (v <= tab[0]) {
            return c[0];
        }
        if (v > tab[tab.length - 1]) {
            return c[tab.length - 1];
        }
        int i = 0;
        for (; v > tab[i]; i++) {
        }
        return gradiant(c[i - 1], c[i], (v - tab[i - 1]) / (tab[i] - tab[i - 1]));
    }

    public static String toHtmlColor(Color c) {
        if (c != null) {
            return "#" + Integer.toHexString(c.getRGB())/*.substring(2)*/.toUpperCase();
        }
        return null;
    }
    
    public static String toHtmlColor6(Color c) {
        if (c != null) {
            return "#" + Integer.toHexString(c.getRGB()).substring(2).toUpperCase();
        }
        return null;
    }
    
    public static Color fromHtmlColor(String colorCode) {
        if (colorCode == null || colorCode.isEmpty()) {
            return null;
        }
        if (colorCode.startsWith("#")) {
            colorCode = colorCode.substring(1);
        }

        final long hexValue = Long.parseLong(colorCode.toLowerCase(), 16);
        return new Color((int) hexValue, colorCode.length() > 6);
    }

    public static Color[] fromHtmlColor(String[] htmlColors) {
        Color[] colors = new Color[htmlColors.length];
        for (int i = 0; i < htmlColors.length; i++) {
            colors[i] = fromHtmlColor(htmlColors[i]);
        }
        return colors;
    }

    public static String[] toHtmlColor(Color[] colors) {
        String[] htmlColors = new String[colors.length];
        for (int i = 0; i < colors.length; i++) {
            htmlColors[i] = toHtmlColor(colors[i]);
        }
        return htmlColors;
    }




    public static String[] convertToHexa(Color[] colors) {
        final String[] tab = new String[colors.length];
        for (int i = 0; i < colors.length; i++) {
            tab[i] = Integer.toHexString(colors[i].getRGB());
        }
        return tab;
    }

    // Retourne une couleur au hasard
    public static Color random() {
        return new Color((int) (Math.random() * 255.), (int) (Math.random() * 255.), (int) (Math.random() * 255.));
    }

    /**
     * Renvoie un générateur de couleurs arc en ciel, qui cycle après le nombre indiqué de couleurs...
     * @param cycle
     * @return 
     */
    public static ColorGenerator createGenerator(int cycle) {
        return new ColorGenerator(cycle, 1);
    }

    /**
     * Renvoie un générateur de couleurs, qui contient au moins le nombre indiqué de couleurs,
     * et qui les renvoie dans un ordre permettant de les différencier correctement les
     * unes après les autres
     * @param cycle
     * @return 
     */
    public static ColorGenerator createScatteredGenerator(int cycle) {
        int dst = (int) Math.sqrt(cycle);
        if (dst < 2) {
            dst = 2;
        }
        while (cycle % dst != 0) {
            cycle++;
        }
        int increment = cycle / dst - 1;
        if (increment == 0) {
            increment = 1;
        }
        return new ColorGenerator(cycle, increment);
    }

    /**
     * Générateur de couleur, qui permet de créer des couleurs type arc en ciel parmi un nombre spécifié d'avance
     */
    public static class ColorGenerator {

        private int number;
        private int increment;
        private int current;

        public ColorGenerator(int cycle, int increment) {
            this.number = cycle;
            this.increment = increment;
        }

        public Color getNext() {
            Color col = Color.getHSBColor(current / (float) number, 0.5f, 0.7f);
            current += increment;
            current = current % number;
            return col;
        }
    }


    /** Modifie la palette pour la rendre plus jolie */
    public static void sort(Color[] palette) {
        if (palette == null || palette.length < 2) {
            return;
        }

        List<Color> colorsToAdd = new ArrayList<Color>();
        Collections.addAll(colorsToAdd, palette);

        // On garde la premiere couleur
        colorsToAdd.remove(palette[0]);

        // On place la couleur la plus eloignee de toutes les autres en premier
        // pour eviter de casser un beau degradé en 2 ensuite
        int idLonely = findLonelyColor(palette);
        Color colorLonely = new Color(palette[idLonely].getRGB());
        palette[idLonely] = palette[0]; // On conserve la premiere couleur
        palette[0] = colorLonely;

        // On ajoute chaque fois la couleur la plus proche restante pour avoir une bonne continuitee
        for (int i = 1; i < palette.length; i++) {
            palette[i] = findClosestColor(palette[i - 1], colorsToAdd);
            colorsToAdd.remove(palette[i]);
        }
    }

    public static int findLonelyColor(final Color[] palette) {
        if (palette == null || palette.length == 0) {
            return 0;
        }
        double dist, distMax = 0;
        int bestId = 0;
        for (int i = 0; i < palette.length; i++) {
            dist = findClosestColorDistance(palette[i], palette);
            if (dist > distMax) {
                distMax = dist;
                bestId = i;
            }
        }
        return bestId;
    }

    public static double findClosestColorDistance(final Color color, final Color[] palette) {
        double d, dMin = Double.MAX_VALUE;
        for (Color c : palette) {
            if (c != color) {
                d = calculateColorDistance(color, c);
                if (d < dMin) {
                    dMin = d;
                }
            }
        }
        return dMin;
    }

    public static Color findClosestColor(final Color color, final Color[] palette) {
        double d, dMin = Double.MAX_VALUE;
        Color best = null;
        for (Color c : palette) {
            d = calculateColorDistance(color, c);
            if (d < dMin) {
                best = c;
                dMin = d;
            }
        }
        return best;
    }

    public static Color findClosestColor(final Color color, final Collection<Color> palette) {
        return findClosestColor(color, palette, Double.MAX_VALUE);
    }

    public static Color findClosestColor(final Color color, final Collection<Color> palette, final double limit) {
        double d, dMin = limit;
        Color best = null;
        for (Color c : palette) {
            d = calculateColorDistance(color, c);
            if (d < dMin) {
                best = c;
                dMin = d;
            }
        }
        return best;
    }

    /** Calcul la distance entre 2 couleurs (cf: http://www.compuphase.com/cmetric.htm)  */
    public static double calculateColorDistance(final Color c1, final Color c2) {
        long rm= ((long) c1.getRed() + (long) c2.getRed()) / 2;
        long r = (long) c1.getRed()   - (long) c2.getRed();
        long g = (long) c1.getGreen() - (long) c2.getGreen();
        long b = (long) c1.getBlue()  - (long) c2.getBlue();
        return Math.sqrt((((512 + rm) * r * r) >> 8) + 4 * g * g + (((767 - rm) * b * b) >> 8));
    }


}
