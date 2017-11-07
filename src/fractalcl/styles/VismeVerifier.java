/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.styles;

import java.awt.Color;
import javax.swing.BorderFactory;
import javax.swing.InputVerifier;
import javax.swing.JComponent;

/**
 *
 * @author visme
 */
public abstract class VismeVerifier extends InputVerifier {

    @Override
    public boolean verify(JComponent input) {
        boolean isValid = isValid(input);
        if (!isValid) {
            input.setBorder(BorderFactory.createLineBorder(Color.red));
        } else {
            input.setBorder(BorderFactory.createLineBorder(Color.BLACK));
        }
        return isValid;
    }

    protected abstract boolean isValid(JComponent input);
}
