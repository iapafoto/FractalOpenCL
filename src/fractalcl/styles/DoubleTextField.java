package fractalcl.styles;

import java.awt.event.KeyEvent;
import javax.swing.JTextField;

/**
 *
 * @author visme
 */
public class DoubleTextField extends JTextField {

    private int[] mSpecialKeys = {KeyEvent.VK_ENTER, KeyEvent.VK_BACK_SPACE, KeyEvent.VK_DELETE, KeyEvent.VK_LEFT, KeyEvent.VK_RIGHT, KeyEvent.VK_UP, KeyEvent.VK_DOWN, KeyEvent.VK_TAB, KeyEvent.VK_END, KeyEvent.VK_BEGIN, KeyEvent.VK_HOME};
    private int max_length;

    public DoubleTextField() {
        this(0, -1);
    }

    public DoubleTextField(int max_length) {
        this(max_length, -1);
    }

    /**
     * @param max_length
     * @param id  1=TemperatureTextFiled, 2= PressureTextFiled
     */
    public DoubleTextField(int max_length, int id) {
        super();
        getDocument().putProperty("id", id);
        this.max_length = max_length;
    }

    @Override
    public void processKeyEvent(KeyEvent ev) {

        char c = ev.getKeyChar();
        // Dans tous les cas, on accepte les touches de contrôle
        if (isSpecialKey(ev.getKeyCode())) {
            super.processKeyEvent(ev);
            if (ev.getKeyCode() == KeyEvent.VK_ENTER) {
                firePropertyChange("value", null, getText());
            }
        } else {
            // si on est au max_length, on arrête
            if (max_length > 0 && ((getSelectedText() == null || getSelectedText().isEmpty()) && getDocument().getLength() >= max_length)) {
                ev.consume();
                return;
            }
            // Un '-' ne peut être qu'au début
            if (c == '-' && getCaretPosition() == 0 && !getText().contains("-")) {
                super.processKeyEvent(ev);
            } else {
                // Un '.' ne peut pas être au début et on ne peut pas en avoir plusieurs
                if (c == '.' && (getCaretPosition() != 0 && !getText().contains("."))) {
                    super.processKeyEvent(ev);
                } else {
                    // On accepte tous les chiffres
                    if (Character.isDigit(c)) {
                        super.processKeyEvent(ev);
                    } else {
                        // Et le reste, poubelle !
                        ev.consume();
                    }
                }
            }
        }
    }

    private boolean isSpecialKey(int keycode) {
        for (final int acceptedKey : mSpecialKeys) {
            if (keycode == acceptedKey) {
                return true;
            }
        }
        return false;
    }

    public void setText(float flt) {
        if (flt % 1 == 0) {
            setText(Long.toString(Math.round(flt)));
        } else {
            setText(Float.toString(flt));
        }
    }

    public void setText(double dbl) {
        if (dbl % 1 == 0) {
            setText(Long.toString(Math.round(dbl)));
        } else {
            setText(Double.toString(dbl));
        }
    }

    public void setText(int i) {
        setText(Integer.toString(i));
    }
}