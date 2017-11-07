/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.editor;

import java.awt.Point;
import javax.swing.JTextPane;
import javax.swing.text.Element;

/**
 *
 * @author durands
 */
public class JEditorTextPane extends JTextPane {
    

    @Override
    public boolean getScrollableTracksViewportWidth() {
        return getUI().getPreferredSize(this).width  <= getParent().getSize().width;
    }
    
    public int lineIdFromPos(Point pt) {
        final int rowStartOffset = viewToModel( pt );
        final Element root = getDocument().getDefaultRootElement();
        return root.getElementIndex(rowStartOffset) +1;
    }
    
}

