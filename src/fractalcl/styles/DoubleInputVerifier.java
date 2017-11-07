package fractalcl.styles;

import javax.swing.JComponent;
import javax.swing.text.JTextComponent;


/**
 *
 * @author visme
 */
public class DoubleInputVerifier extends VismeVerifier {

    @Override
    public boolean isValid(JComponent input) {
        if (input instanceof JTextComponent) {
            String text = ((JTextComponent) input).getText();
            return RegexTool.isDouble(text);
        }
        return false;
    }
}
