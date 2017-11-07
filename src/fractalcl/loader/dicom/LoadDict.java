/***************************************************************************
* LoadDict.java                                                            *
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


import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Hashtable;
import java.util.StringTokenizer;

import javax.swing.JOptionPane;


/**
 * The <code>LoadDict</code> class opens the Dicom tags dictionary.<br />
 *  It will be needed by the <code>DicomReader</code> class operations on Dicom
 *  tags.
 * 
 * @author Salvatore La Bua    <i>< slabua (at) gmail.com ></i>
 * 
 * @version 1.3.2-1
 * 
 * @see DicomReader
 */
public class LoadDict {
    
    // Hash table initialization
    private Hashtable table = new Hashtable();
    
    
    /**
     * <code>LoadDict</code> class constructor.
     */
    public LoadDict() {
        
        try {
            
            BufferedReader dictFile = new BufferedReader(
//                    new InputStreamReader(getClass().getResourceAsStream("/resources/dict/Dicom.dic")));
                    new InputStreamReader(getClass().getResourceAsStream("Dicom.dic")));
            
            String line = null;
            StringTokenizer tkn;
            
            while ((line = dictFile.readLine()) != null) {
                if (!line.startsWith("#")) {
                    
                    String tag = null;
                    DicomValues dicomvalue = new DicomValues();
                    tkn = new StringTokenizer(line);
                    tag = tkn.nextToken();
                    
                    dicomvalue.vr = tkn.nextToken();
                    dicomvalue.name = tkn.nextToken();
                    dicomvalue.vm = tkn.nextToken();
                    dicomvalue.version = tkn.nextToken();
                    
                    table.put(tag, dicomvalue);
                    
                } // End of if block
            } // End of while block
            
            dictFile.close();
            
        } catch (IOException ioe) {
    	    System.err.println(ioe);
    	    JOptionPane.showMessageDialog(null, ioe, "Exception", JOptionPane.ERROR_MESSAGE);
    	    System.exit(-1);
    	} // End of try-catch block
        
    } // End of LoadDict(dictFileName) constructor
    
    /**
     * The method returns the tag's <i>Name</i> field value.
     * 
     * @param tag <code>String</code> Current tag's <i>ID</i>.
     * 
     * @return <code>String</code> Tag's <i>Name</i> field value.
     */
    public String getName(String tag) {
        return ((DicomValues)table.get(tag)).name;
    } // End of getName(tag) method
    
    /**
     * The method returns the tag's <i>VR</i> field value.
     * 
     * @param tag <code>String</code> Current tag's <i>ID</i>.
     * 
     * @return <code>String</code> Tag's <i>VR</i> field value.
     */
    public String getVR(String tag) {
        return ((DicomValues)table.get(tag)).vr;
    } // End of getVR(tag) method
    
    /**
     * The method returns the tag's <i>VM</i> field value.
     * 
     * @param tag <code>String</code> Current tag's <i>ID</i>.
     * 
     * @return <code>String</code> Tag's <i>VM</i> field value.
     */
    public String getVM(String tag) {
        return ((DicomValues)table.get(tag)).vm;
    } // End of getVM(tag) method
    
    /**
     * The method returns the tag's <i>Version</i> field value.
     * 
     * @param tag <code>String</code> Current tag's <i>ID</i>.
     * 
     * @return <code>String</code> Tag's <i>Version</i> field value.
     */
    public String getVersion(String tag) {
        return ((DicomValues)table.get(tag)).version;
    } // End of getVersion(tag) method
    
    /**
     * The method returns <code>true</code> if current tag exists in the Dicom
     *  dictionary.
     * 
     * @param tag <code>String</code> Current tag's <i>ID</i> to check.
     */
    public boolean isContained(String tag) {
        return table.containsKey(tag);
    } // End of isContain(tag) method
    
    /**
     * The private class <code>Dicomvalue</code> is an utility class where current
     *  tag's field values are stored.
     */
    private class DicomValues {
        String name = null;
        String vr = null;
        String vm = null;
        String version = null;
    } // End of Dicomvalue private class
    
} // End of LoadDict class