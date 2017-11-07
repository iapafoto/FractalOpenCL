package fractalcl.styles;

import javax.swing.InputVerifier;
import javax.swing.JComponent;
import javax.swing.text.JTextComponent;


/**
 *
 * @author visme
 */
public class LongInputVerifier extends InputVerifier {

    @Override
    public boolean verify(JComponent input) {
        if (input instanceof JTextComponent) {
            String text = ((JTextComponent) input).getText();
            return RegexTool.isInteger(text);
        }
        return false;
    }
}
