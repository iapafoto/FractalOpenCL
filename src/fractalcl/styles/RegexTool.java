/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.styles;

import java.util.regex.Pattern;

/**
 *
 * @author visme
 */
public class RegexTool {

    private static final Pattern integerPattern = Pattern.compile("[-]?\\d+");
    private static final Pattern doublePattern = Pattern.compile("[-]?\\d+(\\.?\\d+)?");

    private RegexTool() {
    }

    /**
     * Renvoie vrai si la chaîne passée en paramètres correspond à un entier
     * @param integer la chaîne à tester
     * @return vrai si la chaîne passée en paramètres correspond à un entier
     */
    public static boolean isInteger(String integer) {
        return integerPattern.matcher(integer).matches();
    }

    /**
     * Renvoie vrai si la chaîne passée en paramètres correspond à un double
     * @param integer la chaîne à tester
     * @return vrai si la chaîne passée en paramètres correspond à un double
     */
    public static boolean isDouble(String dbl) {
        return doublePattern.matcher(dbl).matches();
    }
}
