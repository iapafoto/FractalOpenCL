/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.editor;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import javax.swing.JPanel;
import javax.swing.JTextPane;
import javax.swing.text.BadLocationException;
import javax.swing.text.Utilities;

/**
 *
 * @author durands
 */
public class RightEditorPane extends JPanel implements MouseMotionListener, MouseListener {
    Map<Integer, List<CompileInfo>> lstInfo = null;
    int totalNbLines;
    JTextPane textPane;
    
    public RightEditorPane() {
        super();
        this.addMouseListener(this);
        this.addMouseMotionListener(this);
    }
            
    public void setCompileInfo(Map<Integer, List<CompileInfo>> lstInfo, JTextPane textPane) {
        this.textPane = textPane;
        this.lstInfo = lstInfo;
                                        
        this.totalNbLines = textPane.getDocument().getDefaultRootElement().getElementCount();
            repaint();
        }
    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        
        if (lstInfo != null && !lstInfo.isEmpty()) {
            int y,w = getWidth();
            double lineSize = (double)(getHeight()-48) / (double)totalNbLines;
            Color c;
            for (Entry<Integer, List<CompileInfo>> e : lstInfo.entrySet()) {
                int worst = e.getValue().get(0).kind; /*10;
                for (CompileInfo info : e.getValue()) {
                    if (info.kind < worst) worst = info.kind;
                }*/
                c = (worst == 0) ? Color.red : (worst == 1) ? Color.orange : Color.blue; 
                g.setColor(c.brighter().brighter());
                y = 16+(int)(lineSize*e.getKey());
                g.fillRect(0,y-1,w,3);
                g.setColor(this.getBackground());
                g.drawRect(0,y-1,w,3);
            }
        }                
    }
    
    @Override
    public void mouseMoved(MouseEvent e) {
        List<CompileInfo> lst = getInfosAt(e.getY());
        String txt = "";
        if (lst != null) {
            txt += "<html>";
            for(int i=0; i<lst.size(); i++) {
                if (i!=0) txt += "<br/>";
                txt += lst.get(i).txt;
            }
            txt += "</html>";
        }
        setToolTipText(txt);
    }
    
    public static int getPos(int line, int column, JTextPane editor) {
//       return editor.getDocument().getDefaultRootElement().getElement(line).getStartOffset() + column;
        try {
            int off=0;
            int rn = 0;
            while( rn<line-1) {
                off=Utilities.getRowEnd(editor, off)+1;
                rn++;
            }
            return off + column-1;
        } catch (BadLocationException e) {
        }
        return 0;
    }
    
    @Override
    public void mouseClicked(MouseEvent e) {
        List<CompileInfo> lst = getInfosAt(e.getY());
        if (lst != null) {
            CompileInfo first = lst.get(0);
            textPane.setCaretPosition(getPos(first.line, first.ch, textPane));
            textPane.requestFocus();
        }
    }
    
    private List<CompileInfo> getInfosAt(int y) {
        if (lstInfo == null) return null;
        
        double lineSize = (double)(getHeight()-48) / (double)totalNbLines;
        int yInfo;
        int dMin = 10;
        int best = -1;
        for (Integer index : lstInfo.keySet()) {
            yInfo = 16+(int)(lineSize*index);
            if (Math.abs(yInfo-y)<dMin) {
                dMin = Math.abs(yInfo-y);
                best = index;
            }
        }
        List<CompileInfo> lst = this.lstInfo.get(best);
        return lst;
    }

    @Override
    public void mouseDragged(MouseEvent e) {
    }

    @Override
    public void mousePressed(MouseEvent e) {
    }

    @Override
    public void mouseReleased(MouseEvent e) {
    }

    @Override
    public void mouseEntered(MouseEvent e) {
    }

    @Override
    public void mouseExited(MouseEvent e) {
    }

}
