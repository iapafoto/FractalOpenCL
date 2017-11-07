package fractalcl;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
import fractalcl.gui.CLEditor;
import fr.meteo.synopsis.client.geometry.geometry2d.Geometry;
import fr.meteo.synopsis.client.geometry.geometry3d.Camera3D;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferedImage;
import fractalcl.loader.DataLoader;
import java.util.Random;
import javax.swing.JFrame;
import org.joda.time.DateTime;

/**
 *
 * @author durands
 */
public class Fractal extends CLEditor {
    
    private final static int 
            PARAM_CAM_FROM = 0,
            PARAM_CAM_LOOK_DIR = 1,
            PARAM_CAM_UP_DIR = 2,
            PARAM_SLIDER = 3,
            PARAM_SLIDER2 = 4,
            PARAM_IMG_OUT = 5,
            PARAM_SUB_PIXEL = 6,
            PARAM_PALETTE = 7,
            PARAM_ZBUFF_OUT = 8,
            PARAM_CAM_CONF = 9;
    
    private final static BufferedImage buffPalette = DataLoader.getResourceImg("demPalette.png");
    
    public Fractal() {
        super("Fractal.cl");

        clManager.createInput(PARAM_PALETTE, buffPalette);

        initSlider(0, -5, 5, 0.70968f);
        initSlider(1, -5, 5, 0.92784f);
    
        initSliderMinMax(0, -5.f, 5, -0.63248f, 0.63248f);
        initSliderMinMax(1, -5.f, 5, -0.78632f, 0.78632f);
        initSliderMinMax(2, -5.f, 5, -0.875f,   0.875f);
        initSliderMinMax(3, -5.f, 5, .5f, .9f);

        //slidersMinMax1.setVisible(true);
        //sliderMinMax2.setVisible(true);
        //sliderMinMax3.setVisible(true);
        //sliderMinMax4.setVisible(true);

     //   final int width = (screen.getWidth() / 64) * 64,
     //             height = (screen.getHeight() / 64) * 64;
        zBuffer = new float[2000*1000];
    }

    int imgId = 0;


    @Override
    protected void updateKernelArgsAfterCompil() {
        updateKernelArgsAfterCtrlChanged();
      
        clManager.setImageOutOnArg(PARAM_IMG_OUT, false);
        clManager.setImageInOnArg(PARAM_PALETTE);
        clManager.setZBufferOnArg(PARAM_ZBUFF_OUT);
    }

    @Override
    protected void updateKernelArgsAfterResize() {
        
        clManager.setImageOutOnArg(PARAM_IMG_OUT, false);
        clManager.setZBufferOnArg(PARAM_ZBUFF_OUT);
//        clManager.setArg(PARAM_IMG_OUT, clManager.imageOut);
    }

    @Override
    protected void updateKernelArgsAfterBBoxChanged(final double[] bbox) {
        if (!clManager.isKernelOk()) {
            return;
        }
    }

    @Override
    protected void updateKernelArgsAfterTimeChanged(final DateTime t) {
        if (t != null) {
            Long time = t.getMillis();
           // clManager.setArg(PARAM_SUNPOS, sunPos);
        }
      //  clManager.setImageInOnArg(PARAM_BUFFER_2D);
    }

    @Override
    protected void updateKernelArgsAfterBBoxChanging(double[] bbox) {
        updateKernelArgsAfterCtrlChanged();
    }

    @Override
    protected void updateKernelArgsAfterCtrlChanged() {
        float[] sliderMin = getSliderMins();
        float[] sliderMax = getSliderMaxs();
        
        sliderMin[3] = getSliderValue(0);
        sliderMax[3] = getSliderValue(1);

        final Camera3D camera = getCamera();
        
        clManager.setArg(PARAM_CAM_FROM, float4(camera.getPosition()));
        clManager.setArg(PARAM_CAM_LOOK_DIR, float4(camera.getLookDir()));
        clManager.setArg(PARAM_CAM_UP_DIR, float4(camera.getSideDir()));
        clManager.setArg(PARAM_CAM_CONF, new float[] {(float)camera.getFocal(), (float)camera.getFocusDistance(), (float)camera.getAperture(), 0.f});
        clManager.setArg(PARAM_SUB_PIXEL, new float[]{.5f, .5f, 0.f, 0.f});
        clManager.setArg(PARAM_SLIDER, sliderMin);
        clManager.setArg(PARAM_SLIDER2, sliderMax);
    }
    
    float map(float a0, float b0, float a1, float b1, float v) {
        return Geometry.mix(a1, b1, (v - a0) / (b0 - a0));
    }
    
    float randomGenerator(long seed) {
        Random generator = new Random(seed);
        return generator.nextFloat();
    }
    float randSubX=0.f, randSubY=0.f, randsubDofX=0.f, randsubDofY=0.f;
    
    @Override
    public void updateKernelArgsRenderingMode(boolean antialiasing) {
        int nbHQSubPicture = screen.img.getHeight()/subPictureHQHeight;                
        int subPictureId = nbCumulFrameHQ%nbHQSubPicture; 
        long seed = nbCumulFrameHQ / nbHQSubPicture;
        
        if (subPictureId == 0) {
            randSubX = (float)Math.random();
            randSubY = (float)Math.random();
            do {
                randsubDofX = 2.f*((float)Math.random()-1.f);
                randsubDofY = 2.f*((float)Math.random()-1.f);
            } while(randsubDofX*randsubDofX + randsubDofY*randsubDofY > 1.f);
        }
        float 
                subX = (antialiasing) ? randSubX : .5f,
                subY = (antialiasing) ? randSubY : .5f;
        float subDofX = 0.f, subDofY = 0.f;
        
        if (antialiasing) {
            subDofX = randsubDofX;
            subDofY = randsubDofY;
            subY += subPictureId*subPictureHQHeight;
        }
        
        clManager.setArg(PARAM_SUB_PIXEL, new float[]{subX, subY, subDofX, subDofY});
   //     clManager.setArg(PARAM_QUALITY, antialiasing ? 1 : 0);
    }

    public static void main(String[] args) throws Exception {

        GraphicsEnvironment env = GraphicsEnvironment.getLocalGraphicsEnvironment();
        GraphicsDevice[] devices = env.getScreenDevices();
        GraphicsDevice device = devices[0];

        JFrame frame = new JFrame("Fractal 3D Editor", device.getDefaultConfiguration());
        JFrame framePicture = new JFrame("Fractal3D", device.getDefaultConfiguration());

        // Init IHM ----------------------------------------   
        Fractal ch = new Fractal();
        ch.setPreferredSize(new Dimension(800, 600));
        //  frame.setLayout(new BorderLayout());
        frame.add(ch, BorderLayout.CENTER);
        // -------------------------------------------------  
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        framePicture.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
        boolean isFullScreen = false; //device.isFullScreenSupported();
        
        frame.setUndecorated(isFullScreen);
        frame.setResizable(!isFullScreen);
        if (isFullScreen) {
            // Full-screen mode
            device.setFullScreenWindow(frame);
            // device.setDisplayMode(new DisplayMode(800,600,32,60));
            //   frame.createBufferStrategy(2);
            frame.validate();
        } else {
            // Windowed mode
            frame.pack();
            frame.setVisible(true);    
        }
        
        framePicture.add(ch.screen, BorderLayout.CENTER);
        framePicture.setResizable(true);
        framePicture.pack();
        framePicture.setVisible(true);
        
        ch.onLoadProgram();
    }


    
    
}
