/***************************************************************************
* ShowGPL.java                                                             *
*--------------------------------------------------------------------------*
*                                                                          *
* created on   : 20 mar 2005                                               *
* copyright    : (C) 2005, Salvatore La Bua      <SLB(at)Shogoki.it>       *
* copyright    : (C) 2005, Calogero Crapanzano   <calosan(at)libero.it>    *
* copyright    : (C) 2005, Pietro Amato          <lovedstone(at)libero.it> *
*                                                                          * 
***************************************************************************/

/***************************************************************************
*                                                                          *
*   This program is free software; you can redistribute it and/or modify   *
*   it under the terms of the GNU General Public License as published by   *
*   the Free Software Foundation; either version 2 of the License, or      *
*   (at your option) any later version.                                    *
*                                                                          *
***************************************************************************/


package fractalcl.loader.dicom;


import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Font;
// import java.awt.event.ActionEvent;
// import java.awt.event.ActionListener;
// import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

// import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;


/**
 * The <code>ShowGPL</code> class is called by the <code>DicomReaderGUI</code>
 *  class to show the GNU General Public License to the end user.
 * 
 * @author Salvatore La Bua    <i>< slabua (at) gmail.com ></i>
 * 
 * @version 1.3.2-1
 */
public class ShowGPL extends JFrame {
    
    private final static String version = DicomReader.getAppVersion();
    private final static String appName = "DicomReader " + version + " - Dicom file decoder";
    
    //private static BufferedReader license = null;
    
    private JTextArea resultArea;
    
    
    /**
     * The class constructor puts all graphics application's components.
     */
    public ShowGPL() {
        super("DicomReader is under the GNU General Public License");
        
        // Main container
        Container c = getContentPane();
        c.setLayout(new BorderLayout(10, 10));
        
        // resultPanel resultArea
        JPanel resultPanel = new JPanel();
        resultPanel.setLayout(new BorderLayout(10, 10));
        
        resultArea = new JTextArea();
        resultArea.setLineWrap(true);
        resultArea.setWrapStyleWord(true);
        resultArea.setOpaque(false);
        resultArea.setEditable(false);
        
        //resultArea.setFont(new Font("Lucida Console", Font.PLAIN, 12));
        //resultArea.setFont(new JLabel().getFont());
        
        JScrollPane scrollPanel = new JScrollPane(resultArea);
        resultPanel.add(scrollPanel);
        
        c.add(resultPanel, BorderLayout.CENTER);
        
        // bottomPanel
        JPanel bottomPanel = new JPanel();
        JLabel bottomLabel = new JLabel(appName);
        bottomLabel.setEnabled(false);
        bottomPanel.add(bottomLabel);
        
        c.add(bottomPanel, BorderLayout.SOUTH);
        
        
    } // End of ShowGPL constructor
    
    /**
     * The <code>showResult</code> method opens the GPL file and shows it to the
     *  user.
     */
    public void showResult() {
        String line = null;
        String text = "";
	    
        try {
	        
	        BufferedReader license = new BufferedReader(
	                new InputStreamReader(getClass().getResourceAsStream("data/resources/gpl/COPYING")));
	        
	        while ((line = license.readLine()) != null)
	            text += line + "\n";
	        
	        license.close();
	        
	        resultArea.setText(text);
	        
	        resultArea.setCaretPosition(0);
	        
        } catch (IOException ioe) {
		    System.err.println(ioe);
		    JOptionPane.showMessageDialog(ShowGPL.this, ioe, "Exception", JOptionPane.ERROR_MESSAGE);
		} // End try-catch block
	    
    } // End of showResult method
    
    /**
     * The <code>display</code> method sets the window's dimensions and shows it
     *  in the center of the screen.<br />
     * It executes the <code>showResult</code> method to show the GPL license file.
     */
    public static void display() {
      final ShowGPL window = new ShowGPL();
      
      window.setSize(580, 550);
      window.setResizable(false);
      window.setLocationRelativeTo(null);
      window.setVisible(true);
      window.addWindowListener(
              new WindowAdapter() {
                  public void windowClosing(WindowEvent e) {
                      System.err.println("Closing...");
                      window.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
                  }
              }
      );
      
      window.showResult();
      
    } // End of display method
    
} // End of ShowGPL class