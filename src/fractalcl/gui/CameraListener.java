/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.gui;
    
import java.awt.Point;
import java.awt.event.MouseEvent;


/**
 *
 * @author sebastien.durand
 */
public interface CameraListener {

    public void onStopWheel(double zoom);
    /* Deplacement rapide de la carte (scroll, zoom, etc) */
  //  public ViewPort onViewPortModificationBegin(Object master);
  //  public void onViewPortModification   (Object master, ViewPort viewPort);
  //  public void onViewPortModificationEnd(Object master, ViewPort viewPort);

    /* Deplacement instantané de la carte vers un autre état */
  //  public void setViewPort(ViewPort viewPort);
  //  public ViewPort getViewPort();
    
  //  public void setCursor(Cursor defaultCursor);
  //  public void setInfo(String INFO_START);
  //  public void repaint();
    public void onMousePressed(MouseEvent me);
    
    public void onMouseClicked(MouseEvent me);
    
    public void onWheel(double zoom, Point pt);

    public void onStopResize();

    public void onDragged(double d, double d0);

    public void onDragEnd();

    public void onResize(int w, int h);

    public double getZMapDistanceAt(int x, int y);
}


