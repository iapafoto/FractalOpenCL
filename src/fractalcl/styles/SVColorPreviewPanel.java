/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.styles;

import java.awt.Color;
import java.awt.Graphics;
import javax.swing.JPanel;

/**
 *
 * @author sebastien.durand
 */
public class SVColorPreviewPanel extends JPanel {
    public SVColorPreviewPanel() {
        super();
    }
    
    @Override
    public void paintComponent(Graphics g) {
        int w =getWidth();
        int h = getHeight();
        for (int x=0;x<w; x+=8) {
            for (int y=0;y<h; y+=8) {
                g.setColor((((x>>3)+(y>>3))%2==0)?Color.white:Color.gray);
                g.fillRect(x,y,8,8);
            }
        }
        g.setColor(getBackground());
        g.fillRect(0,0,w,h);
    }
}
