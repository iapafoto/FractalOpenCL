/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.gui;

import com.nativelibs4java.util.IOUtils;
import fr.meteo.synopsis.client.geometry.geometry2d.Geometry;
import fr.meteo.synopsis.client.geometry.geometry3d.Camera3D;
import fr.meteo.synopsis.client.geometry.geometry3d.Point3D;
import fr.meteo.synopsis.client.geometry.geometry3d.Vector3D;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.beans.PropertyChangeEvent;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import fractalcl.api.OpenCLWithJOCL;
import fractalcl.editor.CompileInfo;
import fractalcl.editor.JCompilableCodeEditor;
import fractalcl.styles.SliderComponent;
import fractalcl.styles.SliderMinMaxComponent;
import java.awt.event.MouseEvent;
import java.awt.geom.Point2D;
import java.awt.image.DataBufferFloat;
import javax.imageio.ImageIO;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.JTextPane;
import javax.swing.Timer;
import javax.swing.filechooser.FileFilter;
import javax.swing.text.BadLocationException;
import javax.swing.text.Utilities;
import org.joda.time.DateTime;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

/**
 *
 * @author durands
 */
public abstract class CLEditor extends javax.swing.JPanel implements KeyListener, JCompilableCodeEditor.Compilator {
    
    private final String savePath = "C:\\Users\\durands\\Desktop\\Fractales\\Save";
    
    // TODO ajouter eclairage curseur via 
    // https://www.shadertoy.com/view/ltjGDd
    //  float sphLight( vec3 P, vec3 N, vec4 L) {
    //  vec3 oc = L.xyz - P;
    //  return max(0., L.w*dot(N,oc)/dot(oc,oc));
    //}
    public final static long SECONDE = 1000, MINUTE = 60 * SECONDE, HOUR = 60 * MINUTE, DAY = 24 * HOUR;
    protected String lastFile = null;

   // protected Point cursor = new Point(0,0);
  //  protected BufferedImage lowQualityPicture = null;

    //protected double zoom = 1;
    public ImagePanel screen = new ImagePanel();

    
    protected float[] zBuffer;
    //
   // public float cam_fov = 3.f;
   // public float cam_dof = 1.f;
    
  //  protected BufferedImage imageTempHQ = null;

    // Image HQ incrementale
    protected Timer timerHQ;
    protected float[] sumR, sumG, sumB;
    protected int nbCumulFrameHQ = 0;
        
    public final static int subPictureHQHeight = 16;
    
    protected int NB_HQ_FRAMES = 256*subPictureHQHeight;


    
// OPENCL Stuff    
//    public OpenCLWithJavaCL clManager = new OpenCLWithJavaCL();
    public OpenCLWithJOCL clManager = new OpenCLWithJOCL();
    
    
    public void onLoadCode() {
        //Handle open button action.
        JFileChooser fc = new JFileChooser();
        fc.setCurrentDirectory(new File(savePath));
        fc.setFileFilter(new FileFilter() {

            @Override
            public boolean accept(File f) {
                if (f.isDirectory()) {
                    return true;
                }
                String extension = "";

                int i = f.getName().lastIndexOf('.');
                if (i > 0) {
                    extension = f.getName().substring(i + 1);
                }
                if (extension != null) {
                    if (extension.equals("png")) {// || extension.equals("conf") || extension.equals("png")) {
                        return true;
                    }
                }
                return false;
            }

            @Override
            public String getDescription() {
                return "OpenCL file"; //To change body of generated methods, choose Tools | Templates.
            }
        });

        ImagePreviewPanel preview = new ImagePreviewPanel();
        fc.setAccessory(preview);
        fc.addPropertyChangeListener(preview);
       // fc.setFileView(new ThumbnailFileView());
        int returnVal = fc.showOpenDialog(this);
        
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            try {
                File file = fc.getSelectedFile();
                lastFile = file.getPath();
                //This is where a real application would open the file.
                if (file.getName().endsWith("cl")) {
                    this.codePane.setText(IOUtils.readText(file));
                } else {
                    if (file.getName().endsWith("png")) {
                        String path = file.getPath();
                        path = path.replace(".png", ".conf");
                        file = new File(path);
                    }
                    JSONParser jsonParser = new JSONParser();
                    JSONObject json = (JSONObject)jsonParser.parse(IOUtils.readText(file));
                    this.restoreConfigurationFromJson(json);
                }
                // log.append("Opening: " + file.getName() + "." + newline);
            } catch (IOException | ParseException ex) {
                Logger.getLogger(CLEditor.class.getName()).log(Level.SEVERE, null, ex);
            }
        } else {
            //  log.append("Open command cancelled by user." + newline);
        }
    }

    public void onSaveCode() {

          
        JFileChooser fc = new JFileChooser();

        if (lastFile != null) {
            fc.setCurrentDirectory(new File(lastFile));
        }

        int returnVal = fc.showSaveDialog(this);

        if (returnVal == JFileChooser.APPROVE_OPTION) {
            saveToFile(fc.getSelectedFile(), codePane.getText());
        }
    }

    final void saveToFile(File file, String txt) {
          BufferedWriter writer = null;
            try {
                writer = new BufferedWriter(new FileWriter(file));
                writer.write(txt);
            } catch (IOException e) {
            } finally {
                try {
                    if (writer != null) {
                        writer.close();
                    }
                } catch (IOException e) {
                }
            }
    }

    @Override
    public void keyTyped(KeyEvent e) {
        if (e.isControlDown() && e.getKeyCode() == KeyEvent.VK_S) {
            if (lastFile != null) {
                BufferedWriter writer = null;
                try {
                    writer = new BufferedWriter(new FileWriter(lastFile));
                    writer.write(codePane.getText());
                } catch (IOException ee) {
                } finally {
                    try {
                        if (writer != null) {
                            writer.close();
                        }
                    } catch (IOException ee) {
                    }
                }
            }
        }
    }

    @Override
    public void keyPressed(KeyEvent e) {
    }

    @Override
    public void keyReleased(KeyEvent e) {
    }

    SliderComponent[] sliders;
    SliderMinMaxComponent[] slidersMinMax;
    
   // ListenableFuture<Boolean> futurePicture = null;
    /**
     *
     * @param clFileName
     */
    public CLEditor(final String clFileName) {

        clManager.initOpenCL();
        initComponents();

        codePane.setCompilator(this);

      //  screenContainer.add(BorderLayout.CENTER, screen);

        try {
            // TODO read it on disk
            String src = IOUtils.readText(CLEditor.class.getResource("../opencl/" + clFileName));
            codePane.setText(src);
        } catch (IOException ex) {
            Logger.getLogger(CLEditor.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        sliders = new SliderComponent[] {
            slider1, slider2, slider3, slider4
        }; 
        for (SliderComponent slider : sliders) {
            slider.setLimits(-10., 10);
        }
        
        slidersMinMax = new SliderMinMaxComponent[] {
            sliderMinMax1, sliderMinMax2, sliderMinMax3, sliderMinMax4
        }; 
        
        for (SliderComponent slider : sliders) {
            slider.addPropertyChangeListener(SliderComponent.VALUE_CHANGING, (PropertyChangeEvent evt) -> {
                if (clManager.isKernelOk()) {
                    updateKernelArgsAfterCtrlChanged();
                    regenerateFast();
                }
            });
        }
         
        for (SliderMinMaxComponent slider : slidersMinMax) {
            slider.addPropertyChangeListener(SliderComponent.VALUE_CHANGING, (PropertyChangeEvent evt) -> {
                if (clManager.isKernelOk()) {
                    updateKernelArgsAfterCtrlChanged();
                    regenerateFast();
                }
            });
            slider.addPropertyChangeListener(SliderComponent.VALUE_CHANGED, (PropertyChangeEvent evt) -> {
                if (clManager.isKernelOk()) {
                    updateKernelArgsAfterCtrlChanged();
                    regenerateHDWithReinit();
                }
            });
        }
       
        timerHQ = new Timer(30, (ActionEvent e) -> {
            timerHQ.stop();
            
            final BufferedImage imageHQ = clManager.createFloatBufferedImage(screen.img.getWidth(), subPictureHQHeight, 4); 
            //new BufferedImage(screen.img.getWidth(), subPictureHQHeight, BufferedImage.TYPE_INT_ARGB);
            
            int nbHQSubPicture = screen.img.getHeight()/subPictureHQHeight;                
            int subPictureId = nbCumulFrameHQ%nbHQSubPicture; 
        //    clManager.setArg(PARAM_SUB_PICTURE, subPictureId);

            clManager.regenerate(true, imageHQ, null/*zBuffer*/, CLEditor.this);
        //    clManager.setArg(PARAM_SUB_PICTURE, 0);
            
            sumPictureHQ(imageHQ, subPictureId, nbHQSubPicture);
                
//            if (futurePicture != null && !futurePicture.isDone() && !futurePicture.isCancelled()) {
//                futurePicture.cancel(true);
//            }
//            futurePicture = calculPicture();
//            Futures.addCallback(futurePicture, new FutureCallback<Boolean>() {
//                @Override
//                public void onSuccess(Boolean v) {
                    nbCumulFrameHQ++;
                    screen.repaint();
                    if (nbCumulFrameHQ < NB_HQ_FRAMES) {
                        timerHQ.restart(); // Les 50 ms s'ecoulent apres la generation de l'image
                    }
//                }
//                @Override
//                public void onFailure(Throwable thrwbl) {
//                    //throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//                }
//            }, SwingExecutor.getInstance());
        });
        timerHQ.setRepeats(false);
    }
    
    protected void sumPictureHQ(BufferedImage input, int subPictureId, int nbSubPicture) {
        
      //  int subPictureId = nbCumulFrameHQ%nbHQSubPicture; 
        final float[] c1 = ((DataBufferFloat) input.getRaster().getDataBuffer()).getData();

        final int w = (int) input.getWidth(),
                  h = (int) input.getHeight();
        // TODO: tout ca pourrait etre fait en Kernel
        // Accumulation
        int offSub = w * h * subPictureId;
        if (nbCumulFrameHQ < nbSubPicture) {
            for (int off = 0; off < w * h; off++) {
                sumR[offSub+off] = c1[off*4+0];
                sumG[offSub+off] = c1[off*4+1];
                sumB[offSub+off] = c1[off*4+2];
            }
        } else {
            for (int off = 0; off < w * h; off++) {
                sumR[offSub+off] += c1[off*4+0];
                sumG[offSub+off] += c1[off*4+1];
                sumB[offSub+off] += c1[off*4+2];
            }
        }

        // Preparation de l'image somme
        final int[] c2 = ((DataBufferInt) screen.img.getRaster().getDataBuffer()).getData();
        float k = nbCumulFrameHQ/nbSubPicture + 1;
        
        for (int off = offSub; off <offSub + w * h; off++) {
            c2[off] = (0xff << 24)
                    + ((((int)(100.*sumR[off] / k)) & 0xff) << 16)
                    + ((((int)(100.*sumG[off] / k)) & 0xff) << 8)
                    + ((((int)(100.*sumB[off] / k)) & 0xff));
        }
    }
    /*
    protected void sumPictureHQ(BufferedImage input, int subPictureId, int nbSubPicture) {
        
      //  int subPictureId = nbCumulFrameHQ%nbHQSubPicture; 
        final int[] c1 = ((DataBufferInt) input.getRaster().getDataBuffer()).getData();

        final int w = (int) input.getWidth(),
                  h = (int) input.getHeight();
        // TODO: tout ca pourrait etre fait en Kernel
        // Accumulation
        int offSub = w * h * subPictureId;
        if (nbCumulFrameHQ < nbSubPicture) {
            for (int off = 0; off < w * h; off++) {
                sumR[offSub+off] = ((c1[off] >> 16) & 0xff);
                sumG[offSub+off] = ((c1[off] >> 8) & 0xff);
                sumB[offSub+off] = ((c1[off]) & 0xff);
            }
        } else {
            for (int off = 0; off < w * h; off++) {
                sumR[offSub+off] += ((c1[off] >> 16) & 0xff);
                sumG[offSub+off] += ((c1[off] >> 8) & 0xff);
                sumB[offSub+off] += ((c1[off]) & 0xff);
            }
        }

        // Preparation de l'image somme
        final int[] c2 = ((DataBufferInt) screen.img.getRaster().getDataBuffer()).getData();
        int k = nbCumulFrameHQ/nbSubPicture + 1;
        
        for (int off = offSub; off <offSub + w * h; off++) {
            c2[off] = (0xff << 24)
                    + (((sumR[off] / k) & 0xff) << 16)
                    + (((sumG[off] / k) & 0xff) << 8)
                    + (((sumB[off] / k) & 0xff));
        }
    }
    */
    BufferedImage imageLQ = null;

    public void regenerate() {
        int w = screen.img.getWidth();
        int h = screen.img.getHeight();
        //  regenerateFast();
        if (zBuffer == null || zBuffer.length != w*h) {      
            zBuffer = new float[w*h];
        }
     //   if (w*h < 800*600) { // blindage pour pas exploser si image trop grande
     //       clManager.regenerate(false, screen.img, zBuffer, this);
     //   }
        screen.repaint();
    }

    public void regenerateHD() {
        regenerate();
        timerHQ.restart();
        //clManager.regenerate(true, screen.img, zBuffer, this);
        //sumPictureHQ();
        //screen.repaint();
    }
    

    
    protected void onSavePicture() {
        final String basePath = savePath + "\\fractal_" + System.currentTimeMillis();
        
        System.out.println("Save picture to: " + basePath);
        
        saveToFile(new File(basePath + ".cl"), codePane.getText());
        saveToFile(new File(basePath + ".conf"), getConfigurationAsJson());
                
        try {
            File outputfile = new File(basePath + ".png");
            ImageIO.write(screen.img, "png", outputfile);
        } catch (IOException ex) {
            Logger.getLogger(CLEditor.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public String getConfigurationAsJson() {
        final Camera3D cam = getCamera();
        JSONObject 
                json = new JSONObject(),
                jsonCam = new JSONObject(),
                jsonSlider = new JSONObject();

        jsonCam.put("width", cam.getWidth());
        jsonCam.put("height", cam.getHeight());
        jsonCam.put("focal", cam.getFocal());
        jsonCam.put("focdist", cam.getFocusDistance());
        jsonCam.put("aperture", cam.getAperture());
        jsonCam.put("pos",  toJsonArray(new double[] { cam.getPosition().x, cam.getPosition().y, cam.getPosition().z }));
        jsonCam.put("look", toJsonArray(new double[] { cam.getLookDir().x,  cam.getLookDir().y,  cam.getLookDir().z }));
        jsonCam.put("up",   toJsonArray(new double[] { cam.getUpDir().x,    cam.getUpDir().y,    cam.getUpDir().z }));
        json.put("camera", jsonCam);
        
        jsonSlider.put("sliders", toJsonArray(getSliderValues()));
       // jsonSlider.put("slider2", getSliderValue(1));
        jsonSlider.put("sliderMins", toJsonArray(getSliderMins()));        
        jsonSlider.put("sliderMaxs", toJsonArray(getSliderMaxs()));
        json.put("slider", jsonSlider);
        
        json.put("code", codePane.getText());
        
        return json.toString();
    }
    
    public void restoreConfigurationFromJson(JSONObject json) {
        JSONObject jsonCam = (JSONObject)json.get("camera");
        
        final Camera3D cam = new Camera3D(
                jsonToPoint3D((JSONArray)jsonCam.get("pos")),
                jsonToVector3D((JSONArray)jsonCam.get("look")),
                jsonToVector3D((JSONArray)jsonCam.get("up")),
                jsonToDouble(jsonCam.get("focal")),                
                jsonToDouble(jsonCam.get("focdist")),
                jsonToDouble(jsonCam.get("aperture")),
                jsonToInt(jsonCam.get("width")),
                jsonToInt(jsonCam.get("height"))     
        );
        
        final JSONObject jsonSlider = (JSONObject)json.get("slider");
        
        if (jsonSlider.containsKey("slider1")) { // retrocompatibilite
            sliders[0].setValue(jsonToFloat(jsonSlider.get("slider1")));
            sliders[1].setValue(jsonToFloat(jsonSlider.get("slider2")));
        } else {
            float[] values = jsonToFloatArray((JSONArray)jsonSlider.get("sliders"));
            for (int i=0; i< values.length; i++) 
                sliders[i].setValue(values[i]);
        }
        
        float[] mins = jsonToFloatArray((JSONArray)jsonSlider.get("sliderMins"));
        float[] maxs = jsonToFloatArray((JSONArray)jsonSlider.get("sliderMaxs"));
        
        for (int i=0; i<mins.length; i++) {
            slidersMinMax[i].setValues(mins[i], maxs[i]);
        }
        
        screen.controller.camera = cam;
        cam.setScreenSize(screen.getWidth(),screen.getHeight() );
 //       screen.setSize(cam.getWidth(), cam.getHeight());
        codePane.setText(json.get("code").toString());
    }
    
    double jsonToDouble(Object json) {
        return Double.valueOf(json.toString());  
    }    
    float jsonToFloat(Object json) {
        return Float.valueOf(json.toString());  
    }
    int jsonToInt(Object json) {
        return Integer.valueOf(json.toString());  
    }
    
    Point3D jsonToPoint3D(JSONArray json) {
        double[] arr = jsonToDoubleArray(json);
        return new Point3D(arr[0], arr[1], arr[2]);
    }
    
    Vector3D jsonToVector3D(JSONArray json) {
        double[] arr = jsonToDoubleArray(json);
        return new Vector3D(arr[0], arr[1], arr[2]);
    }
    
    double[] jsonToDoubleArray(JSONArray json) {
        double[] arr = new double[json.size()];
        for (int i=0; i<json.size(); i++) 
            arr[i] = Double.valueOf(json.get(i).toString());        
        return arr;
    }
    
    float[] jsonToFloatArray(JSONArray json) {
        float[] arr = new float[json.size()];
        for (int i=0; i<json.size(); i++) 
            arr[i] = Float.valueOf(json.get(i).toString());        
        return arr;
    }
        
    JSONArray toJsonArray(double[] vals) {
        JSONArray arr = new JSONArray();
        for (double v : vals) arr.add(Double.toString(v));        
        return arr;
    }
    
    JSONArray toJsonArray(float[] vals) {
        JSONArray arr = new JSONArray();
        for (float v : vals) arr.add(Float.toString(v));       
        return arr;
    }
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jTabbedPane1 = new javax.swing.JTabbedPane();
        jPanel1 = new javax.swing.JPanel();
        jSplitPane1 = new javax.swing.JSplitPane();
        jScrollPane2 = new javax.swing.JScrollPane();
        errorPane = new javax.swing.JTextPane();
        jPanel3 = new javax.swing.JPanel();
        codePane = new fractalcl.editor.JCompilableCodeEditor();
        toolbarPanel = new javax.swing.JPanel();
        jToggleButton1 = new javax.swing.JToggleButton();
        jToggleButton2 = new javax.swing.JToggleButton();
        jToggleButton3 = new javax.swing.JToggleButton();
        jToggleButton4 = new javax.swing.JToggleButton();
        jToggleButton5 = new javax.swing.JToggleButton();
        jPanel2 = new javax.swing.JPanel();
        slider1 = new fractalcl.styles.SliderComponent();
        slider2 = new fractalcl.styles.SliderComponent();
        slider3 = new fractalcl.styles.SliderComponent();
        slider4 = new fractalcl.styles.SliderComponent();
        sliderMinMax1 = new fractalcl.styles.SliderMinMaxComponent();
        sliderMinMax2 = new fractalcl.styles.SliderMinMaxComponent();
        sliderMinMax3 = new fractalcl.styles.SliderMinMaxComponent();
        sliderMinMax4 = new fractalcl.styles.SliderMinMaxComponent();

        jTabbedPane1.setBackground(new java.awt.Color(0, 0, 0));

        jSplitPane1.setBackground(new java.awt.Color(0, 0, 0));
        jSplitPane1.setBorder(null);
        jSplitPane1.setDividerSize(8);
        jSplitPane1.setOrientation(javax.swing.JSplitPane.VERTICAL_SPLIT);
        jSplitPane1.setResizeWeight(1.0);

        jScrollPane2.setBorder(null);

        errorPane.setBorder(null);
        errorPane.setForeground(new java.awt.Color(255, 102, 0));
        jScrollPane2.setViewportView(errorPane);

        jSplitPane1.setBottomComponent(jScrollPane2);

        jPanel3.setLayout(new java.awt.BorderLayout());

        codePane.setBackground(new java.awt.Color(0, 0, 0));
        codePane.setForeground(new java.awt.Color(255, 255, 255));
        jPanel3.add(codePane, java.awt.BorderLayout.CENTER);

        toolbarPanel.setBackground(new java.awt.Color(0, 0, 0));
        toolbarPanel.setDoubleBuffered(false);
        toolbarPanel.setMinimumSize(new java.awt.Dimension(100, 54));
        toolbarPanel.setPreferredSize(new java.awt.Dimension(296, 40));

        jToggleButton1.setText("Compile");
        jToggleButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jToggleButton1ActionPerformed(evt);
            }
        });

        jToggleButton2.setText("Load");
        jToggleButton2.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jToggleButton2ActionPerformed(evt);
            }
        });

        jToggleButton3.setText("Save");
        jToggleButton3.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jToggleButton3ActionPerformed(evt);
            }
        });

        jToggleButton4.setText("SavePng");
        jToggleButton4.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jToggleButton4ActionPerformed(evt);
            }
        });

        jToggleButton5.setText("Stop");
        jToggleButton5.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jToggleButton5ActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout toolbarPanelLayout = new javax.swing.GroupLayout(toolbarPanel);
        toolbarPanel.setLayout(toolbarPanelLayout);
        toolbarPanelLayout.setHorizontalGroup(
            toolbarPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, toolbarPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jToggleButton1, javax.swing.GroupLayout.PREFERRED_SIZE, 69, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jToggleButton5)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jToggleButton4)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 276, Short.MAX_VALUE)
                .addComponent(jToggleButton2)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jToggleButton3)
                .addContainerGap())
        );
        toolbarPanelLayout.setVerticalGroup(
            toolbarPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(toolbarPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(toolbarPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jToggleButton1)
                    .addComponent(jToggleButton2)
                    .addComponent(jToggleButton3)
                    .addComponent(jToggleButton4)
                    .addComponent(jToggleButton5))
                .addContainerGap(20, Short.MAX_VALUE))
        );

        jPanel3.add(toolbarPanel, java.awt.BorderLayout.SOUTH);

        jSplitPane1.setLeftComponent(jPanel3);

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jSplitPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 625, Short.MAX_VALUE)
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jSplitPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 509, Short.MAX_VALUE)
        );

        jTabbedPane1.addTab("tab1", jPanel1);

        jPanel2.setLayout(new javax.swing.BoxLayout(jPanel2, javax.swing.BoxLayout.PAGE_AXIS));
        jPanel2.add(slider1);
        jPanel2.add(slider2);
        jPanel2.add(slider3);
        jPanel2.add(slider4);
        jPanel2.add(sliderMinMax1);
        jPanel2.add(sliderMinMax2);
        jPanel2.add(sliderMinMax3);
        jPanel2.add(sliderMinMax4);

        jTabbedPane1.addTab("tab3", jPanel2);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jTabbedPane1, javax.swing.GroupLayout.Alignment.TRAILING)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jTabbedPane1)
        );
    }// </editor-fold>//GEN-END:initComponents

    private void jToggleButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jToggleButton1ActionPerformed
        if (onLoadProgram()) {
            onScreenResize();
            if (clManager.isKernelOk()) {
                updateKernelArgsAfterCompil();
                nbCumulFrameHQ = 0;
                regenerateHDWithReinit();
            }
        }
    }//GEN-LAST:event_jToggleButton1ActionPerformed

    private void jToggleButton2ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jToggleButton2ActionPerformed
        onLoadCode();
    }//GEN-LAST:event_jToggleButton2ActionPerformed

    private void jToggleButton3ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jToggleButton3ActionPerformed
        onSaveCode();
    }//GEN-LAST:event_jToggleButton3ActionPerformed

    private void jToggleButton4ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jToggleButton4ActionPerformed
        onSavePicture();
    }//GEN-LAST:event_jToggleButton4ActionPerformed

    private void jToggleButton5ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jToggleButton5ActionPerformed
        this.timerHQ.stop();
    }//GEN-LAST:event_jToggleButton5ActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private fractalcl.editor.JCompilableCodeEditor codePane;
    private javax.swing.JTextPane errorPane;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JSplitPane jSplitPane1;
    private javax.swing.JTabbedPane jTabbedPane1;
    private javax.swing.JToggleButton jToggleButton1;
    private javax.swing.JToggleButton jToggleButton2;
    private javax.swing.JToggleButton jToggleButton3;
    private javax.swing.JToggleButton jToggleButton4;
    private javax.swing.JToggleButton jToggleButton5;
    private fractalcl.styles.SliderComponent slider1;
    private fractalcl.styles.SliderComponent slider2;
    private fractalcl.styles.SliderComponent slider3;
    private fractalcl.styles.SliderComponent slider4;
    protected fractalcl.styles.SliderMinMaxComponent sliderMinMax1;
    protected fractalcl.styles.SliderMinMaxComponent sliderMinMax2;
    protected fractalcl.styles.SliderMinMaxComponent sliderMinMax3;
    protected fractalcl.styles.SliderMinMaxComponent sliderMinMax4;
    private javax.swing.JPanel toolbarPanel;
    // End of variables declaration//GEN-END:variables

    static class IntRef {

        IntRef() {
        }

        IntRef(int val) {
            v = val;
        }
        int v = 0;
    }

    static String getBetween(String txt, String s1, String s2, int[] pos, int posEnd) {
        final int[] pos2 = {pos[0]};

        final String s = getBetween(txt, s1, s2, pos2);
        if (pos2[0] > posEnd) {
            return null;
        }
        pos[0] = pos2[0];
        return s;
    }

    static String getBetweenFromLast(String txt, String s1, String s2, int[] pos) {
        int pos2 = txt.indexOf(s2, pos[0]);
        if (pos2 >= 0) {
            // pos1 += s1.length();
            int pos1 = txt.lastIndexOf(s1, pos2);
            if (pos1 >= 0) {
                pos[0] = pos2 + s2.length();
                return txt.substring(pos1 + s1.length(), pos2).trim();
            }
        }
        return null;
    }

    static String getBetween(String txt, String s1, String s2, int[] pos) {
        int pos1 = txt.indexOf(s1, pos[0]);
        if (pos1 >= 0) {
            pos1 += s1.length();
            int pos2 = txt.indexOf(s2, pos1);
            if (pos2 >= 0) {
                pos[0] = pos2 + s2.length();
                return txt.substring(pos1, pos2).trim();
            }
        }
        return null;
    }

    public void regenerateHDWithReinit() {
        timerHQ.stop();
//        if (futurePicture != null && !futurePicture.isDone() && !futurePicture.isCancelled()) {
//            futurePicture.cancel(true);
//        }
        int size = (int) (screen.img.getWidth() * screen.img.getHeight());

        // Reinitialisation des accumulateurs
        if (sumR == null || sumR.length != size) {
            sumR = new float[size];
            sumG = new float[size];
            sumB = new float[size];
        }
        nbCumulFrameHQ = 0;

     //   regenerateHD();
     
        timerHQ.restart();
    }

    protected void regenerateFast() {
//        if (futurePicture != null && !futurePicture.isDone() && !futurePicture.isCancelled()) {
//            futurePicture.cancel(true);
//        }
        timerHQ.stop();
        regenerate();
    }

    protected abstract void updateKernelArgsAfterCompil();

    protected abstract void updateKernelArgsAfterResize();

    protected abstract void updateKernelArgsAfterCtrlChanged();

    protected void updateKernelArgsAfterTimeChanged(DateTime t) {
    }

    protected void updateKernelArgsAfterBBoxChanging(double[] bbox) {
    }

    protected void updateKernelArgsAfterMapBBoxChanged(double[] bbox) {
    }

    protected void postProcessPicture(final BufferedImage img) {
    }

    public void updateKernelArgsRenderingMode(boolean antialiasing) {
    }

    protected void updateKernelArgsAfterBBoxChanged(double[] d) {
        updateKernelArgsAfterBBoxChanging(d);
        regenerateHDWithReinit();
    }

    public static int getRow(int pos, JTextPane editor) {
        int rn = (pos == 0) ? 1 : 0;
        try {
            int offs = pos;
            while (offs > 0) {
                offs = Utilities.getRowStart(editor, offs) - 1;
                rn++;
            }
        } catch (BadLocationException e) {
        }
        return rn;
    }

    public static int getColumn(int pos, JTextPane editor) {
        try {
            return pos - Utilities.getRowStart(editor, pos) + 1;
        } catch (BadLocationException e) {
        }
        return -1;
    }

    public float getSliderValue(int id) {
        return (float) sliders[id].getValue();
    }
    
    public float[] getSliderValues() {
        float[] values = new float[sliders.length];
        for (int i=0; i<sliders.length; i++)
            values[i] = getSliderValue(i);
        return values;
    }
    
    public float[] getSliderMins01() {
        return new float[]{
            (float) Geometry.invMix(sliderMinMax1.getMinimum(), sliderMinMax1.getMaximum(), sliderMinMax1.getValueMin()),
            (float) Geometry.invMix(sliderMinMax2.getMinimum(), sliderMinMax2.getMaximum(), sliderMinMax2.getValueMin()),
            (float) Geometry.invMix(sliderMinMax3.getMinimum(), sliderMinMax3.getMaximum(), sliderMinMax3.getValueMin()),
            (float) Geometry.invMix(sliderMinMax4.getMinimum(), sliderMinMax4.getMaximum(), sliderMinMax4.getValueMin())
        };
    }

    public float[] getSliderMaxs01() {
        return new float[]{
            (float) Geometry.invMix(sliderMinMax1.getMinimum(), sliderMinMax1.getMaximum(), sliderMinMax1.getValueMax()),
            (float) Geometry.invMix(sliderMinMax2.getMinimum(), sliderMinMax2.getMaximum(), sliderMinMax2.getValueMax()),
            (float) Geometry.invMix(sliderMinMax3.getMinimum(), sliderMinMax3.getMaximum(), sliderMinMax3.getValueMax()),
            (float) Geometry.invMix(sliderMinMax4.getMinimum(), sliderMinMax4.getMaximum(), sliderMinMax4.getValueMax())
        };
    }
    
    public float[] getSliderMins() {
        return new float[]{
            (float) sliderMinMax1.getValueMin(),
            (float) sliderMinMax2.getValueMin(),
            (float) sliderMinMax3.getValueMin(),
            (float) sliderMinMax4.getValueMin()};
    }
        
    public float[] getSliderMaxs() {
        return new float[]{
            (float) sliderMinMax1.getValueMax(),
            (float) sliderMinMax2.getValueMax(),
            (float) sliderMinMax3.getValueMax(),
            (float) sliderMinMax4.getValueMax()};

    }

    public static  float[] float4(Vector3D v) {
        return float4(v.x, v.y, v.z);
    }
    
    public static  float[] float4(Point3D p) {
        return float4(p.x, p.y, p.z);
    }
    
    public static float[] float4(final double v1, final double v2, final double v3) {
        return new float[]{(float) v1, (float) v2, (float) v3, 0f};
    }

    public static float[] float3(final double v1, final double v2, final double v3) {
        return new float[]{(float) v1, (float) v2, (float) v3};
    }

    public static float[] float4(final double v1, final double v2, final double v3, final double v4) {
        return new float[]{(float) v1, (float) v2, (float) v3, (float) v4};
    }
    
    public Camera3D getCamera() {
        return screen.controller.camera;
    }
    
    public class ImagePanel extends JPanel implements CameraListener {

        public BufferedImage img = null;
        public CameraController controller;

        ImagePanel() {
            super();
            controller = new CameraController(this, this) {};
        }

        @Override
        public void paintComponent(Graphics g) {

            g.setColor(Color.black);
            g.fillRect(0, 0, getWidth(), getHeight());
            if (img != null) {
                g.drawImage(img, (getWidth() - img.getWidth()) / 2, (getHeight() - img.getHeight()) / 2, img.getWidth(), img.getHeight(), null);
            }
        }

        @Override
        public Dimension getPreferredSize() {
            return new Dimension(256, 256);
        }

        @Override
        public void onStopWheel(double zoomLevel) {
            if (clManager.isKernelOk()) {
                updateKernelArgsAfterCtrlChanged();
                regenerateHDWithReinit();
            }
        }

        @Override
        public void onWheel(double zoomLevel, Point pt) {
            if (clManager.isKernelOk()) {
                updateKernelArgsAfterCtrlChanged();
                regenerateFast();
            }
        }

        @Override
        public void onResize(int w, int h) {
        }

        @Override
        public void onStopResize() {
            onScreenResize();
            if (clManager.isKernelOk()) {
                updateKernelArgsAfterResize();
                regenerateHDWithReinit();
            }
        }
        
        @Override
        public void onDragged(double dx, double dy) {
            if (clManager.isKernelOk()) {
                updateKernelArgsAfterCtrlChanged();
                regenerateFast();
            }
        }

        @Override
        public void onDragEnd() {
            if (clManager.isKernelOk()) {
                updateKernelArgsAfterCtrlChanged();
                regenerateHDWithReinit();
            }
        }
        @Override
        public void onMousePressed(MouseEvent me) {
            timerHQ.stop();
        }
        
        @Override
        public void onMouseClicked(MouseEvent me) {
            if (me.getButton() == MouseEvent.BUTTON3/*getClickCount() == 2*/ && zBuffer != null) {
                int dx = (getWidth() - img.getWidth()) / 2, 
                    dy = (getHeight() - img.getHeight()) / 2;

                float dist = zBuffer[me.getX()-dx+img.getWidth()*(me.getY()-dy)];
                if (dist >0) {
                    CLEditor.this.getCamera().setFocusDistance(dist);
                    controller.setFocusPoint(new Point2D.Double(me.getX(),me.getY()), dist);

                    updateKernelArgsAfterCtrlChanged();
                    regenerateHDWithReinit();
                }
            }
        }

        @Override
        public double getZMapDistanceAt(int x, int y) {
            if (zBuffer != null && x+getWidth()*y < zBuffer.length)
                return zBuffer[x+getWidth()*y];
            return 0.;
        }
    }

    public void initSliderMinMax(int id, float min, float max, float vmin, float vmax) {
        initSliderMinMax(id, min, max, vmin, vmax, false);
    }

    public void initSlider(int id, float min, float max, float v) {
        sliders[id].setLimits(min, max);
        sliders[id].setValue(v);
        sliders[id].setLogarithmic(false);
    }
    
    public void initSliderMinMax(int id, float min, float max, float vmin, float vmax, boolean isLogarithmic) {
        slidersMinMax[id].setLimits(min, max);
        slidersMinMax[id].setValues(vmin, vmax);
        slidersMinMax[id].setLogarithmic(isLogarithmic);
    }

    public boolean onLoadProgram() {
        String codeSrc = codePane.getText();
        String[] error = {""};
        if (clManager.createProgram(codeSrc, error)) {
            clManager.createKernels();
            this.updateKernelArgsAfterCompil();
            return true;
        } else {
            final Map<Integer, List<CompileInfo>> errors = displayError(codeSrc, error[0]);
        }
        return false;
    }

    @Override
    public Map<Integer, List<CompileInfo>> buildProgram(final String codeSrc) {
        errorPane.setText("");
        String[] error = {""};
        if (clManager.createProgram(codeSrc, error)) {
            errorPane.setForeground(Color.green);
            errorPane.setText("OK");
            toolbarPanel.setBackground(Color.black);
            return null;
        } else {
            return displayError(codeSrc, error[0]);
        }
    }
    
    private Map<Integer, List<CompileInfo>> displayError(String src, String errorMsg) {                 
        String[] lines = errorMsg.split("\n");
        errorPane.setForeground(Color.orange);
        errorPane.setText(errorMsg);
        Map<Integer, List<CompileInfo>> lstInfo = new LinkedHashMap<>();
        List<CompileInfo> lst;
        int kind = -1;
        String[] infoKind = {"error", "warning", "note"};
        String includeFile;
        for (String line : lines) {
            int pos = -1;
            for (int i = 0; i < infoKind.length; i++) {
                pos = line.lastIndexOf(infoKind[i] + ":");
                if (pos > 0) {
                    kind = i;
                    break;
                }
            }
            if (kind >= 0) {
                int p2 = line.lastIndexOf(":", pos);
                int p1 = line.lastIndexOf(":", p2 - 1);
                int p0 = line.lastIndexOf(":", p1 - 1);

                try {
                    // TODO parse line like:     C:\Users\durands\.javacl\includes\includes1481186168803523593\javaclsimple/opencl/Interpol.cl:43:62: note: passing argument to parameter 'buf' here"
                    int lineId = Integer.parseInt(line.substring(p0 + 1, p1));
                    int chr = Integer.parseInt(line.substring(p1 + 1, p2));

                    if (line.charAt(0) >= '0' && line.charAt(0) <= '9') {
                        includeFile = null;
                    } else {
                        includeFile = line.substring(0, p0);
                        int p4 = 0;
                        do {
                            p4 = src.indexOf("#include", p4);
                            if (p4 >= 0) {
                                String file = getBetween(src, "\"", "\"", new int[]{p4});
                                if (file != null && !file.isEmpty() && includeFile.endsWith(file)) {
                                    lineId = codePane.getLineIdOfCaret(p4);
                                    line = line.substring(p0 + 1);
                                    break;
                                }
                            }
                        } while (p4 >= 0);
                    }

                    lst = lstInfo.get(lineId);
                    if (lst == null) {
                        lst = new ArrayList<>();
                        lstInfo.put(lineId, lst);
                    }

                    lst.add(new CompileInfo(kind, lineId, chr, line, includeFile));
                } catch (Exception e2) {

                }
            }
        }

        if (!lstInfo.isEmpty()) {
            toolbarPanel.setBackground(Color.red);
        }

        return lstInfo;
    }
    
    int PARAM_ZBUFF_OUT = 6;
    
    public void onScreenResize() {
        final int width = screen.getWidth(), //(screen.getWidth() / 64) * 64,
                  height = screen.getHeight();//(screen.getHeight() / 64) * 64;

        if (screen.img == null || screen.img.getWidth() != width || screen.img.getHeight() != height) {
            screen.img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        }
        
        if (zBuffer == null || zBuffer.length != width*height) {
            zBuffer = new float[width*height]; 
        }
    //    if (imageTempHQ == null || imageTempHQ.getWidth() != width || imageTempHQ.getHeight() != height) {
    //        imageTempHQ = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
    //    }
        clManager.onScreenResize(screen.img, zBuffer);//, imageTempHQ);
      //  regenerateFast();
    }
    
    
    
//    
//    ListenableFuture<Boolean> calculPicture() {
//        final ListenableFuture<Boolean> res = ApplicationExecutor.getInstance().getExecutor().submit(new Callable<Boolean>() {
//            @Override
//            public Boolean call() throws Exception {
//                clManager.regenerate(true, screen.img, zBuffer, CLEditor.this);
//                sumPictureHQ();
//                return true;
//            }
//        });
//        return res;
//    }
//    
    
    
    
}
