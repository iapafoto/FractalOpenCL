/***************************************************************************
* DicomReader.java                                                         *
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


import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;

import javax.swing.JOptionPane;


/**
 * The <code>DicomReader</code> class decodes Dicom files writing headers file
 *  and images files.<br />
 * The image(s) files is simply composed of ascii pixel values ordered by rows
 *  and columns as in a trivial image.
 * 
 * @author Salvatore La Bua    <i>< slabua (at) gmail.com ></i>
 * 
 * @version 1.3.2-1
 * 
 * @see LoadDict
 */
public class DicomReader {
    
    // Application version
    private final static String version = "1.3.2-1";
    
    // Debug flag
    private static boolean DEBUG = true;
    
    // Files
    private static File dicomF = null;
    private static File asciiF = null;
    
    // Utility Flags
    private boolean isDicomFile = true;
    //private boolean asciiSaved = false;
    //private boolean imageSaved = false;
    //private boolean headersSaved = false;
    
    // Streams
    private FileInputStream dicomFile = null;
  //  private BufferedWriter headersFile = null;
    
    // Temporary Buffers
    private final byte[] buff1 = new byte[1];
    private final byte[] buff2 = new byte[2];
    private final byte[] buff4 = new byte[4];
    private final byte[] buff8 = new byte[8];
    private int tempBuff = 0;

    // Dicom tag's attributes
    private String tagID = null;
    private String tagName = null;
    private String tagVR = null;
   // private final String content = null;

    // Utility variable used to handle the "SQ ReferencedImageSequence" tag
    private int seqN = 0;

    // Image class object that will contain image's parameters
    private final Image image = new Image();
	
	
	/**
	 * Four parametrs <code>DicomReader</code> class constructor.<br />
	 * The constructor simply initializes variables' values to file names passed
	 *  as parameters to the class.
	 * 
	 * @param dfn <code>File</code> Dicom file name to decode.
	 * @param hfn <code>File</code> Headers file name.
	 * @param afn <code>File</code> Ascii image file name.
	 * @param ifn <code>File</code> PGM-P2 image file name.
	 */
	public DicomReader(File dfn, File afn) {
	    dicomF = dfn;
	    asciiF = afn;
	} // End of DicomReader(dfn, hfn, afn, ifn) constructor
	
	/**
	 * The <code>decode</code> method is the main method of the class.<br />
	 * It will use other private methods to do the more simplest functions such
	 *  as checks, informations reading/writing...
	 */
    public void decode() {
        // Dicom dictionary
        LoadDict dict = new LoadDict();
        
        try {
            // Streams opening
//            headersFile = new BufferedWriter(new FileWriter(headersF));
            
            // Dicom file type checking (with or without DICM string header)
            dicomFile = new FileInputStream(dicomF);
            dicomFile = this.checkDICM(dicomFile);
            
            int bytesToRead;
            // Starting tags interpretation from Dicom file
            while (true) {
                // Reading next tag
                tagID = this.nextTagID(dicomFile);
                // If tagID doesn't exist in the Dicom dictionary, it will be
                //  ignored with its values
                while (!dict.isContained(tagID)) {
               //     System.out.println("Unknown tagID: " + tagID);
                   //if (tagID.equals("(2d00,2d00)")) {
                   //     bytesToRead = 0;                        
                   // } else {
                        bytesToRead = this.readBytes(dicomFile, tagID, null);
                   // }
              //      System.out.println("tagID: " + tagID + " - Unknown - read: " + bytesToRead + " bytes");
                    byte[] value = new byte[bytesToRead];
                    if (dicomFile.read(value) <0) {
                        return;
//                        tagID = null;
//                        break;
                    }
                    tagID = this.nextTagID(dicomFile);
                } // End of while block
              //  if (tagID == null) continue;
                
                
                // Getting tag's name and type
                tagName = dict.getName(tagID);
                tagVR = dict.getVR(tagID);
                
           //     System.out.println("tagID: " + tagID + " - " + tagName);
                if (tagID.equals("(0008,0081)")) {
                    int debug=1;
                }
                // "ox PixelData" tag check (last tag)
                if (tagName.equals("PixelData")) {
                    // Skip 4 bytes
                    dicomFile.read(buff4);
                    // Interrupts the while block when reading last tag
                    break;
                } // End of if block
                
                // "OB FileMetaInformationVersion" tag check
                if (tagID.equals("(0002,0001)")) {
                    // Skip 4 bytes
                    dicomFile.read(buff4);
                }                
                // Reading the correct number of bytes by the current tag's type
                bytesToRead = this.readBytes(dicomFile, tagID, tagVR);
              //  System.out.println("bytesToRead = " + bytesToRead);
                // "SQ ReferencedImageSequence" tag check
                if (tagID.equals("(0008,1140)")) {
                    seqN = ((bytesToRead / 86) * 2);
                    // Skip 8 bytes
                    dicomFile.read(buff8);
                    tagID = this.nextTagID(dicomFile);
                    tagName = dict.getName(tagID);
                    tagVR = dict.getVR(tagID);
                    bytesToRead = this.readBytes(dicomFile, tagID, tagVR);
                } // End of if block

                // "SQ ReferencedImageSequence" sub-tags handling
                if (seqN != 0)
                    bytesToRead = bytesToRead - 1;

                if (bytesToRead <= 0) {
                    int test=1;
                }
                // Reading tag's value
                byte[] value = new byte[bytesToRead];
                dicomFile.read(value);
                String sValue = "";

                // "SQ ReferencedImageSequence" sub-tags handling
                if (seqN != 0) {
                    dicomFile.read(buff1);
                }
                
                // Handling tags' values by the current tag's type
                if (value.length > 0 && (tagVR.equals("US") || tagVR.equals("UL") || tagVR.equals("xs"))) {
                    int nValue = (((0xff & value[1]) << 8) | (0xff & value[0]));
                    sValue = Integer.toString(nValue);	
                }
                // "SQ ReferencedImageSequence" tag check
                else if (tagID.equals("(0008,1140)")) {
                    seqN = bytesToRead / 86;
                    // Skip 8 bytes
                    dicomFile.read(buff8);

                }
                // "OB FileMetaInformationVersion" tag check
                else if (tagID.equals("(0002,0001)")) {
                    int nValue = ((0xff & value[1]) | ((0xff & value[0]) << 8));
                    sValue = Integer.toString(nValue);	
                } else {
                    sValue = new String(value);
                } 
                // "SQ ReferencedImageSequence" sub-tags handling
                if (seqN != 0) {
                    if ((seqN % 2) == 0) {
                        seqN--;
                    } else {
                        if (seqN != 1) {
                            // Skip 8 bytes
                            dicomFile.read(buff8);
                        } // End of if block
                        seqN--;
                    } // End of if block
                } // End of if block
                // Storing tags that are image-related
                this.storeImageParam(tagID, sValue); 
                
             
            } // End of while block
        } catch (IOException ioe) {
            System.err.println(ioe);
            JOptionPane.showMessageDialog(null, ioe, "Exception", JOptionPane.ERROR_MESSAGE);
            System.exit(-1);
        } // End of try-catch block
    } // End of decode() method
    

    
    /**
     * The <code>chechDICM</code> method checks if the chosen file has the "DICM"
     *  header and acts to handle the file recognized.
     * 
     * @param df <code>FileInputStream</code> File to check.
     * 
     * @return dicomFile <code>FileInputStream</code> File eventually updated.
     */
    private FileInputStream checkDICM(FileInputStream df) {
        try {
            // Reading first four bytes of the Dicom file to decode.
            df.read(buff4);
//            tempBuff = ((0xff&buff4[3])<<32)|((0xff&buff4[2])<<16)|((0xff&buff4[1])<<8)|(0xff&buff4[0]);
            tempBuff = ((0xff&buff4[3])<<24)|((0xff&buff4[2])<<16)|((0xff&buff4[1])<<8)|(0xff&buff4[0]);
            df.close();
            df = new FileInputStream(dicomF);
            
            // If tempBuff contains a zero data value, we can have a Dicom file
            //  with "DICM" header (in this case first 128 bytes of the file are
            //  all zeros), either we can have a file that isn't a Dicom file.
            if (tempBuff == 0) {
                int counter = 0;
                df.read(buff1);
                counter++;
                
                while ((tempBuff = buff1[0]) != 'M') {
                    if (counter == 132) {
                        isDicomFile = false;
                        String errorString = "This file not seems to be a Dicom file.\n\n DICM string not found at 132nd byte.";
                    	System.err.println(errorString);
                    	JOptionPane.showMessageDialog(null, errorString, "Wrong File", JOptionPane.ERROR_MESSAGE);
                    	System.exit(-1);
                    	
                    } // End of if block
                    
                    df.read(buff1);
                    counter++;
                } // End of while block
                
            // If first four bytes aren't zeros, we can have a Dicom "DICM"
            //  headerless file (in this case we assume first four bytes identify
            //  a valid Dicom tag),  either we can have a file that isn't a Dicom
            //  file.
            } else {
                LoadDict dict = new LoadDict();
                
                String tagID = this.nextTagID(df);
                df.close();
                df = new FileInputStream(dicomF);
                
                if (!dict.isContained(tagID)) {
                    isDicomFile = false;
                    String errorString = "This file not seems to be a Dicom file." + "\n\n" +
                    					 "DICM stringless Dicom file:" + "\n" +
                    					 "First tag is not a valid Dicom tag.";
                	System.err.println(errorString);
                	JOptionPane.showMessageDialog(null, errorString, "Wrong File", JOptionPane.ERROR_MESSAGE);
                	System.exit(-1);
                } // End of if block
                
            } // End of if-else block
            
        } catch (IOException ioe) {
		    System.err.println(ioe);
		    JOptionPane.showMessageDialog(null, ioe, "Exception", JOptionPane.ERROR_MESSAGE);
		    System.exit(-1);
		} // End of try-catch block
        
        return df;
        
    } // End of checkDICM(df) method
    
    /**
     * The <code>nextTagID</code> method contructs the identifying string of the
     *  Dicom tag that will be searched into the Dicom dictionary to obtain its
     *  description. 
     * 
     * @param df <code>FileInputStream</code> File to check.
     *  
     * @return tagID <code>String</code> Current tag representative string.
     */
    private String nextTagID(FileInputStream df) {
        FileInputStream dicomFile = df;
        
        String group = null;
        String number = null;
        String tagID = null;
        
        try {
            
            // Extracting group from first two bytes of tagID
            dicomFile.read(buff2);
            tempBuff = (((0xff & buff2[1]) << 8) |
                    	 (0xff & buff2[0]));
            
            group  = Integer.toString((tempBuff & 0x0000f000) >> 12, 16);
            group += Integer.toString((tempBuff & 0x00000f00) >>  8, 16);
            group += Integer.toString((tempBuff & 0x000000f0) >>  4, 16);
            group += Integer.toString((tempBuff & 0x0000000f),       16);
            
            // Extracting number from last two bytes of tagID
            dicomFile.read(buff2);
            tempBuff = (((0xff & buff2[1]) << 8) |
                    	 (0xff & buff2[0]));
            
            number  = Integer.toString((tempBuff & 0x0000f000) >> 12, 16);
            number += Integer.toString((tempBuff & 0x00000f00) >>  8, 16);
            number += Integer.toString((tempBuff & 0x000000f0) >>  4, 16);
            number += Integer.toString((tempBuff & 0x0000000f),       16);
            
            // Constructing the tagID string
            tagID = ("(" + group + "," + number + ")");
            
        } catch (IOException ioe) {
		    System.err.println(ioe);
		    JOptionPane.showMessageDialog(null, ioe, "Exception", JOptionPane.ERROR_MESSAGE);
		    System.exit(-1);
		} // End of try-catch block
        
        return tagID;
        
    } // End of nextTagID(df) method
    
    /**
     * The <code>readBytes</code> method determines how many bytes must be read
     *  from time to time during the Dicom file's headers analysis.
     * 
     * @param df <code>FileInputStream</code> File to check.
     * @param tID <code>String</code> Current tag representative string.
     * 
     * @return readBytes <code>int</code> Number of bytes to read.
     */
    private int readBytes(FileInputStream df, String tID, String tVR) {
        FileInputStream dicomFile = df;
        String tagID = tID;
        String tagVR = tVR;
        int bytesToRead = 0;
        int bytesTemp = 0;
        
        String group = tagID.substring(1, 5);
        String number = tagID.substring(6, 10);
        
        try {
            
            dicomFile.read(buff4);
            if (tagID.equals("(0009,10e3 )")) {
                int debugv=1;
            }
            if (tagID.equals("(0011,1001)")) {
		        bytesToRead = (((0xff & buff4[3]) << 32) |
                                        ((0xff & buff4[2]) << 16) |
                                        ((0xff & buff4[1]) <<  8) |
                                         (0xff & buff4[0])         );
		    } else {
		    	bytesToRead = (((0xff & buff4[3]) << 8) | (0xff & buff4[2]) );
	           	
		    	if ((bytesToRead == 65535) /*FF FF FF FF*/ ) {
		    	    bytesTemp = (((0xff & buff4[1]) << 8) | (0xff & buff4[0]) );
		    	    
		    	    if (bytesTemp == bytesToRead) {
		    	        bytesToRead = 0;
		    	    } // End of if block
		    	    
		    	} else if (bytesToRead == 0) {
		        	bytesTemp = (((0xff & buff4[1]) << 8) | (0xff & buff4[0]) );
		        	if (bytesTemp == 16975) { // OB
                                    dicomFile.read(buff4);
		        	    bytesToRead = buff4[0];
                                } else if (bytesTemp == 20819 /*SQ*/ ) {
		        	    dicomFile.read(buff4);
		        	    bytesToRead = 0;
                                }
		        	else if ((bytesTemp == 16708) /*DA*/ || (bytesTemp == 20300) /*LO*/ ||
		        	           (bytesTemp == 19796) /*TM*/ || (bytesTemp == 18515) /*SH*/ ||
		        	           (bytesTemp == 21315) /*CS*/ || (bytesTemp == 21321) /*IS*/ ||
		        	           // Altri tag
		        	           (bytesTemp == 18773) /*UI*/ || (bytesTemp == 20048) /*PN*/ ||
                                        
		    //      (bytesToRead == 33792) ||  // SDU	  
                                        (bytesTemp == 21587) ||
		        	           (bytesTemp == 21580) /*LT*/ || (bytesTemp == 21316) /*DS*/   ) {
		        	    bytesToRead = 0;
	           		} else {
                                    bytesToRead = bytesTemp;
                                }		        	
		    	} // End of else-if block
                       // else if (bytesToRead == 21331 ) {
		       // 	    dicomFile.read(buff4);
		       // 	    bytesToRead = 0;
		       // 	}
		    } // End of if-else block
	        
        } catch (IOException ioe) {
		    System.err.println(ioe);
		    JOptionPane.showMessageDialog(null, ioe, "Exception", JOptionPane.ERROR_MESSAGE);
		    System.exit(-1);
		} // End of try-catch block
        
        return bytesToRead;
        
    } // End of readBytes(df, tID) method
    
    /**
     * The <code>lastTagHandling</code> method handles the last header's tag
     *  that requires a different treatment than others.
     * 
     * @return sValue <code>String</code> Last tag's numeric value.
     */
    private String lastTagHandling() {
        tempBuff = (((0xff & buff4[3]) << 32) |
                    ((0xff & buff4[2]) << 16) |
                    ((0xff & buff4[1]) <<  8) |
                     (0xff & buff4[0])         );
        
        return Integer.toString(tempBuff);        
    } // End of lastTagHandling method
    
    /**
     * The <code>storeImageParam</code> method sets the <code>Image</code>
     *  object's fields.<br />
     * It will contain some informations to store the image, such as dimensions,
     *  color depth and so on.
     * 
     * @param tID <code>String</code> Current tag representative string.
     * @param sV <code>String</code> Value to store into the object's attribute.
     * 
     * @see Image
     */
    private void storeImageParam(String tID, String sV) {
        if (tID.equals("(0028,0010)")) {
    	    image.rows = Integer.parseInt(sV);
    	} else if (tID.equals("(0028,0011)")) {
    	    image.cols = Integer.parseInt(sV);
    	} else if (tID.equals("(0028,0100)")) {
    	    image.colorDepth = Integer.parseInt(sV);  // BitsAllocated, not BitsStored!
    	} else if (tID.equals("(0028,1052)")) {
    	    image.rescaleIntercept = Integer.parseInt(sV.trim());  
    	} else if (tID.equals("(0028,1053)")) {
    	    image.rescaleSlope = Integer.parseInt(sV.trim());  
    	} else if (tID.equals("(0054,0081)")) { // || tagID.equals("(0028,0008)")) || tagID.equals("(0054,0053)")) {
            image.slices = Integer.parseInt(sV); // NumberOfFrames (or NumberOfSlices or NumberOfFramesInRotation)
        } else if (tagName.equals("PixelData")) {
            image.totBytes = Integer.parseInt(sV);
        } 
    } 

    public void storeData() {
        String imgFileName = asciiF.getAbsolutePath();        
        final int totSlices = image.slices;
        
        if (totSlices == 1) {
            this.saveImage2(imgFileName);            
        } else {
            int startExt = imgFileName.lastIndexOf(".");
            String ext = null;
            if (startExt != -1) {
                ext = imgFileName.substring(startExt, imgFileName.length());
            }
            for (int currentSlice = 0; currentSlice < totSlices; currentSlice++) {                
                if (startExt != -1) {
                    imgFileName = imgFileName.substring(0, startExt + 1) + (currentSlice + 1) + ext;
                } else {
                    imgFileName = imgFileName + "." + (currentSlice + 1);
                } // End of if-else block
                
                this.saveImage2(imgFileName);
            } // End of for block    
        } // End of if-else block
    } // End of storeData method

    private double[][] getValuesArray() {        
        int rows = image.rows;
        int cols = image.cols;
        int bpp = image.colorDepth;
        // pour passer en echelle de Hounsfield
        int mv = image.rescaleIntercept; 
        int k = image.rescaleSlope;
        
        final double[][] buff = new double[cols][]; 
        final byte[] buffrows = new byte[(bpp == 16?2:1)*rows*cols];
        try {
            dicomFile.read(buffrows);
            for (int x = 0; x < cols; x++) {                
                buff[x] = new double[rows]; 
                for (int y = 0; y < rows; y++) {
                    // Color depth check
                    if (bpp == 16) {                         
                        tempBuff = (((0xff & buffrows[x*2+y*cols*2+1]) << 8) | (0xff & buffrows[x*2+y*cols*2]));
                    } else if (bpp == 8) {
                	//dicomFile.read(buff1);
                        tempBuff = (0xff & buffrows[x+y*cols]);
                    } // End of if-else block
                    if (tempBuff!= 63536 && tempBuff!=64511 && tempBuff!=65535) {
                        
                        buff[x][y] = k*tempBuff + mv; // // On passe en echelle de Hounsfield
                    } else {
                        buff[x][y] = -1000;
                    }
                	
                } // End of for block
            } // End of for block
        } catch (IOException ioe) {
            return buff;
        } // End of try-catch block
        return buff;
    } // end of saveImage() method

    
    private void saveImage2(String asciiFN) {        
        int rows = image.rows;
        int cols = image.cols;
        int bpp = image.colorDepth;
        
        BufferedImage img = new BufferedImage(cols, rows, BufferedImage.TYPE_INT_RGB);
        final int[] buff = ((DataBufferInt) img.getRaster().getDataBuffer()).getData(); 
        
        try {
            int min = Integer.MAX_VALUE, max = Integer.MIN_VALUE;
            for (int j = 0; j < cols; j++) {                
                for (int i = 0; i < rows; i++) {
                    // Color depth check
                    if (bpp == 16) { 
                        dicomFile.read(buff2);
                        tempBuff = (((0xff & buff2[1]) <<  8) | (0xff & buff2[0]));
                    } else if (bpp == 8) {
                	dicomFile.read(buff1);
                        tempBuff = (0xff & buff1[0]);
                    } // End of if-else block
                    if (tempBuff!= 63536) {
                       if (tempBuff>max) max = tempBuff;
                        if (tempBuff<min) min = tempBuff;
                        int c = (254*tempBuff)/2358;
                        buff[j*rows+i] = (c<<16) + (c<<8) + (c);
                    }
                	
                } // End of for block
            } // End of for block
            System.out.print("limits : min="+min+" max="+max);
        } catch (IOException ioe) {
            System.err.println(ioe);
            JOptionPane.showMessageDialog(null, ioe, "Exception", JOptionPane.ERROR_MESSAGE);
            System.exit(-1);
        } // End of try-catch block
         try {
            File outputfile = new File(asciiFN+".png");
            ImageIO.write(img, "png", outputfile);
        } catch (IOException e) {

        }
    } // end of saveImage() method

    
    /**
     * The method returns the <code>isDicomFile</code> variable value.<br />
     * <code>True</code> is returned if the file to convert is recognized as a
     *  valid Dicom file.
     * 
     * @return isDicomFile <code>boolean</code> <code>isDicomFile</code> variable value.
     */
    public boolean getIsDicomFile() {
        return isDicomFile;
    } // End of getIsDicomFile method
    

    /**
     * The method returns the <code>version</code> variable value.<br />
     * 
     * @return version <code>String</code> Application's version number.
     */
    public static String getAppVersion() {
        return version;
    } // End of getAppVersion method
    
    /**
     * Four parameters static method of the <code>decode</code> method.
     * 
     * @param dfn <code>File</code> Dicom file name to decode.
	 * @param hfn <code>File</code> Headers file name.
	 * @param afn <code>File</code> Ascii image file name.
	 * @param ifn <code>File</code> PGM-P2 image file name.
	 */
    public static double[][] decodeToArray(File dfn) {
        DicomReader instance = new DicomReader(dfn, null);
        instance.decode();
        return instance.getValuesArray();
    }
    
    public static void orderDir(File folder) {
        int k=folder.listFiles().length;
        double[][][] buff = new double[folder.listFiles().length][][];
        Map<String, File> reorder = new HashMap<String, File>();
        for (final File fileEntry : folder.listFiles()) {
            if (!fileEntry.isDirectory()) {
                    long time = fileEntry.lastModified();
                    reorder.put(time+fileEntry.getName(), fileEntry);
            }
        }
        List<String> names = new ArrayList<String>(reorder.keySet());
        Collections.sort(names);
        int n=1000;
        String path1;
        String path2;
        for (String key : names) {
            try {            
                path1 = "c:/scanner/"+n+".dicom";
                path2 = "c:/scanner/png/"+n+".png";
                Files.copy((reorder.get(key)).toPath(), new File(path1).toPath());
                decodeToImg(new File(path1), new File(path2));
                n++;
            } catch (IOException ex) {
                Logger.getLogger(DicomReader.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    
    Collection<File> orderFileByNumber(File[] listFiles) {
        Map<Integer, File> map = new LinkedHashMap();  
        for (final File fileEntry : listFiles) {
            if (!fileEntry.isDirectory()) {
                String name = fileEntry.getName();
                String numberOnly= name.replaceAll("[^0-9]", "");
                Integer id = Integer.parseInt(numberOnly);
                map.put(id, fileEntry);
            }
        }
        return map.values();
    }
    
    public static double[][][] decodeDirToArray(File folder) {
        int k=folder.listFiles().length;
        double[][][] buff = new double[folder.listFiles().length][][];
        for (final File fileEntry : folder.listFiles()) {
            if (!fileEntry.isDirectory()) {
                DicomReader instance = new DicomReader(fileEntry, new File(fileEntry.getAbsolutePath() + ".png"));
                instance.decode();
                buff[--k] = instance.getValuesArray();
            //    instance.saveImage2(fileEntry.getAbsolutePath() + ".png");
            }
        }
/*
        double[][][] buff = new double[200][][];
        for (int i=0;i<buff.length;i++) {
            DicomReader instance = new DicomReader(new File("C:/IMAGES/IM000"+(100+i*4)), null);
            instance.decode();
            buff[i] = instance.getValuesArray();
        }         
*/
        return buff;
    }
    
    public static void decodeToImg(File dfn, File output) {
        DicomReader instance = new DicomReader(dfn, output);
        instance.decode();
        instance.storeData();
    } 
    /**
     * The <code>main</code> methods of the class simply calls the <code>decode</code>
     *  static method.
     * 
     * @param dicomF <code>File</code> Dicom file name to decode.
	 * @param headersF <code>File</code> Headers file name.
	 * @param asciiF <code>File</code> Ascii image file name.
	 * @param imageF <code>File</code> PGM-P2 image file name.
	 */
    public static void main(String[] args) {
        // Variables initialization
        int argc = 0;
		int nArgs = args.length;
		String thisArg = null;
		
		File dicomF = null;
		File asciiF = null;
		
		// Arguments check and initialization
		while (argc < nArgs) {
			thisArg = args[argc++].trim();
			if (thisArg != null) {
				if (dicomF == null) {
					dicomF = new File(thisArg);
				} else if (asciiF == null) {
					asciiF = new File(thisArg);
				}  // End of if-else block
			} // End of if block
		} // End of while block
        
		// First three parameters presence control
		if (dicomF != null && asciiF != null) {
                    decodeToImg(dicomF, asciiF);
		} else {
		    String message = "Missed parameters:" + "\n" +
		    				 " File dfn, File hfn, File afn, File ifn";
		    System.err.println(message);
		    JOptionPane.showMessageDialog(null, message, "Missed parameters",
		            								JOptionPane.ERROR_MESSAGE);
		    System.exit(-1);
		} // End of if-else block
    } // End of main method
    
    /**
     * The <code>Image</code> private class is an utility class.<br />
     * Its object will contain needed values to the correct image storing process.
     */
    private class Image {
        int rows = 0;
        int cols = 0;
        int colorDepth = 0;
        int slices = 1;
        int totBytes = 0;
        // pour passer en echelle de Hounsfield
        int rescaleIntercept = 0; 
        int rescaleSlope = 1;
    } // End of Image private class
    
} // End of DicomReader class 