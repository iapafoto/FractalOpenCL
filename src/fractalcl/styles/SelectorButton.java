package fractalcl.styles;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowFocusListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import javax.swing.BorderFactory;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JToggleButton;
import javax.swing.SwingUtilities;

/**
 *
 * @author sebastien.durand
 */
public abstract class SelectorButton<Type, Selector extends SelectorPanel<Type>> extends javax.swing.JToggleButton {

    private JFrame popupWindow;
    protected Selector selectorPanel;
    protected int popupWidth = 160;
    protected int popupHeight = 100;
    private int x[] = new int[3];
    private int y[] = new int[3];

    public SelectorButton() {
        selectorPanel = createSelector();
    }

    public SelectorButton(int popupWidth, int popupHeight) {
        this();
        this.popupWidth = popupWidth;
        this.popupHeight = popupHeight;
    }

    protected Selector getSelectorPanel() {
        return selectorPanel;
    }

    /**
     * Initialiser.
     */
    protected void init() {
        final JToggleButton _this = this;

        setText("");
        setFocusable(false);

        addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                // set popup window visibility
                if (!popupWindow.isVisible()) {
                    // set location relative to button
                    Point location = getLocation();
                    SwingUtilities.convertPointToScreen(location, getParent());
                    location.translate(0, getHeight()
                            + (getBorder() == null ? 0
                            : getBorder().getBorderInsets(_this).bottom));
                    popupWindow.setBounds(location.x, location.y, popupWidth, popupHeight);

                    // show the popup if not visible
                    popupWindow.setVisible(true);
                    popupWindow.requestFocus();
                } else {
                    // hide it otherwise
                    popupWindow.setVisible(false);
                }
            }
        });

        // use frame
        popupWindow = new JFrame();
        popupWindow.setUndecorated(true);
        // add some components to window
        popupWindow.getContentPane().setLayout(new BorderLayout());
        ((JComponent) popupWindow.getContentPane()).setBorder(BorderFactory.createEtchedBorder());
        popupWindow.getContentPane().add(selectorPanel, BorderLayout.CENTER);
        popupWindow.pack();

        popupWindow.addWindowFocusListener(new WindowFocusListener() {

            @Override
            public void windowGainedFocus(WindowEvent e) {
            }

            @Override
            public void windowLostFocus(WindowEvent e) {
                SwingUtilities.invokeLater(new Runnable() {

                    @Override
                    public void run() {
                        if (popupWindow.isVisible()) {
                            popupWindow.setVisible(false);
                        }
                    }
                });
            }
        });
        
        selectorPanel.addListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                onSelectionChanged((Type) evt.getNewValue());
                popupWindow.setVisible(false);
                repaint();
            }
        });
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        int w = getWidth(), h = getHeight(), h2 = h / 2;
        x[0] = w - 6;
        x[1] = w - 9;
        x[2] = w - 12;
        y[0] = y[2] = h2 - 1;
        y[1] = h2 + 2;
        g2.setColor(Color.black);
        g2.fillPolygon(x, y, 3);

        g2.setColor(Color.gray);
        g2.drawLine(w - 18, 3, w - 18, h - 4);

        paintSelection(g2, w - 18, h);
    }

    public Type getSelectedValue() {
        return selectorPanel.getSelectedValue();
    }

    public void setSelectedValue(Type value) {
        selectorPanel.setSelectedValue(value);
        repaint();
    }

    public void addListener(PropertyChangeListener listener) {
        selectorPanel.addListener(listener);
    }

    public void removeListeners() {
        selectorPanel.removeListeners();
        selectorPanel = null;
    }

    protected abstract void onSelectionChanged(Type newValue);

    protected abstract void paintSelection(Graphics2D g2, int w, int h);

    protected abstract Selector createSelector();
}
