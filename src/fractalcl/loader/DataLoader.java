package fractalcl.loader;

import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.swing.ImageIcon;

//
//
///**
// *
// * @author durands
// */
public final class DataLoader {
    
    public static ImageIcon getResourceIcon(String name) {
        return new ImageIcon(DataLoader.class.getClassLoader().getResource(name));
    }

    public static BufferedImage getResourceImg(final String name) {
        return toBufferedImage(getResourceIcon(name));
    }
    public static BufferedImage getImgARGB(final String name) {
        try {
            BufferedImage img = ImageIO.read(new File(name));
            if (img != null) {
                BufferedImage img2 = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_INT_ARGB);
                Graphics2D g2 = (Graphics2D) img2.getGraphics();
                g2.drawImage(img, null, 0, 0);
                g2.dispose();
                return img2;
            }
        } catch (IOException ex) {
            Logger.getLogger(DataLoader.class.getName()).log(Level.SEVERE, null, ex);
        }
        return null;
    }

    public static BufferedImage toBufferedImage(final ImageIcon icon) {
        if (icon == null) return null;
        final int w = icon.getIconWidth(), h = icon.getIconHeight();
        final BufferedImage bi = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB); 
        final Graphics g = bi.createGraphics();
        icon.paintIcon(null, g, 0,0);
        g.dispose();
        return bi;
    }
}

//    
//    public static float[] getDem(final ImageIcon icon, final float delta, final float scale) {
//        if (icon == null) return null;
//        final BufferedImage bi = toBufferedImage(icon);       
//        return getDem(bi, delta, scale, null);
//    }
//    
//    
//    public static float[] getDemBathy(final ImageIcon icon, final float delta, final float scale) {
//        if (icon == null) return null;
//        final BufferedImage bi = toBufferedImage(icon);       
//        return getDemBathy2(bi, delta, scale, null);
//    }
//  
//    
//    public static float[] getDemBathy2(final BufferedImage bi, final float delta, final float scale, final float[] outMinMax) {
//        if (bi.getRaster().getDataBuffer() instanceof DataBufferInt) {
//            final int[] buff = ((DataBufferInt)bi.getRaster().getDataBuffer()).getData();   
//            int size = buff.length;
//            final float[] demBuff = new float[size];
//            Color c;
//            if (outMinMax != null && outMinMax.length>1) {
//                outMinMax[0] = Float.MAX_VALUE;
//                outMinMax[1] = -Float.MAX_VALUE;
//            }
//            for (int off=0; off<size; off++) {
//                c = new Color(buff[off]);
//                demBuff[off] = delta + ((float)(c.getRed()-.6f*c.getBlue())*scale);
//                if (outMinMax != null && outMinMax.length>1) {
//                    if (demBuff[off] < outMinMax[0]) outMinMax[0] = demBuff[off];
//                    if (demBuff[off] > outMinMax[1]) 
//                        outMinMax[1] = demBuff[off];
//                }
//            }
//            return demBuff;   
//        }
//
//        return null;
//    }
//    
//    
//    public static float[] getDem(final BufferedImage bi, final float delta, final float scale, final float[] outMinMax) {
//        if (bi.getRaster().getDataBuffer() instanceof DataBufferInt) {
//            final int[] buff = ((DataBufferInt)bi.getRaster().getDataBuffer()).getData();   
//            int size = buff.length;
//            final float[] demBuff = new float[size];
//            Color c;
//            if (outMinMax != null && outMinMax.length>1) {
//                outMinMax[0] = Float.MAX_VALUE;
//                outMinMax[1] = -Float.MAX_VALUE;
//            }
//            for (int off=0; off<size; off++) {
//                c = new Color(buff[off]);
//                demBuff[off] = delta + ((float)(c.getRed()+c.getGreen()+c.getBlue())*scale/3f);
//                if (outMinMax != null && outMinMax.length>1) {
//                    if (demBuff[off] < outMinMax[0]) outMinMax[0] = demBuff[off];
//                    if (demBuff[off] > outMinMax[1]) 
//                        outMinMax[1] = demBuff[off];
//                }
//            }
//            return demBuff;   
//        } else if (bi.getRaster().getDataBuffer() instanceof DataBufferByte) {
//            final byte[] buff = ((DataBufferByte)bi.getRaster().getDataBuffer()).getData();   
//           
//            Color c;
//            if (outMinMax != null && outMinMax.length>1) {
//                outMinMax[0] = Float.MAX_VALUE;
//                outMinMax[1] = -Float.MAX_VALUE;
//            }
//            if (bi.getColorModel().getPixelSize() == 32) {
//                int size = buff.length/4;
//                final float[] demBuff = new float[size];
//                for (int off=0; off<size; off++) {
//                  //  c = new Color(buff[off], buff[off+1], buff[off+2], buff[off+3]);
//                    if ((buff[(off<<2)]&0xff) == 0)
//                        demBuff[off] = 0;
//                    else
//                        demBuff[off] = delta + ((float)((int)(buff[(off<<2)+1]&0xff)+(int)(buff[(off<<2)+2]&0xff)+(int)(buff[(off<<2)+3]&0xff))*scale/3f);
//                    if (outMinMax != null && outMinMax.length>1) {
//                        if (demBuff[off] < outMinMax[0]) 
//                            outMinMax[0] = demBuff[off];
//                        if (demBuff[off] > outMinMax[1]) 
//                            outMinMax[1] = demBuff[off];
//                    }
//                }
//                return demBuff;   
//            } else if (bi.getColorModel().getPixelSize() == 8) {
//                int size = buff.length;
//                final float[] demBuff = new float[size];
//                for (int off=0; off<size; off++) {
//                    demBuff[off] = delta + ((float)((int)(buff[off]&0xff))*scale);
//                    if (outMinMax != null && outMinMax.length>1) {
//                        if (demBuff[off] < outMinMax[0]) 
//                            outMinMax[0] = demBuff[off];
//                        if (demBuff[off] > outMinMax[1]) 
//                            outMinMax[1] = demBuff[off];
//                    } 
//                }
//                return demBuff;   
//            } 
//            
//        }
//        return null;
//    }
//
//    
//    
//    public static float[] getDicom(final String file, final int[] outSize, final float[] outMinMax) {
//        final double[][][] grib = DicomReader.decodeDirToArray(new File(file));
//        return to1D(grib, outSize,outMinMax);
//    }
//
//    public static float[] getDicom(final String file, 
//            final int xmin, final int ymin, final int zmin, final int xmax, final int ymax, final int zmax, 
//            final int[] outSize, final float[] outMinMax) {
//        final double[][][] grib = DicomReader.decodeDirToArray(new File(file));
//        return to1D(grib, xmin,ymin,zmin, xmax,ymax,zmax, outSize,outMinMax);
//    }
//
//    static long 
//            MINUTE = 60*1000, HOUR = 60*MINUTE, DAY = 24*HOUR;
//    
//    public static Map<Long,float[]> getGrib3DAnimFake(final int nb, final String file, final String param, final int[] outSize, final float[] outMinMax) {
//        final double[][][] grib = TextGribLoader.loadGrib3D(file, param);
//        final Map<Long,float[]> anim = new LinkedHashMap<>();
//        long t = (System.currentTimeMillis()/HOUR)*HOUR;
//        float[] next, last = to1D(grib, outSize,outMinMax);
//        anim.put(t, last);
//        for (int i=1;i<nb;i++) {
//            next = new float[last.length];
//            for (int off=0;off<last.length; off++) {
//                next[off] = last[off] + (float)(3.*Math.cos(i*6.28/(double)nb));
//            }
//            t += 15*MINUTE;
//            anim.put(t, next);
//            last = next;
//        }
//        return anim;
//    }
//   
//    public static Map<Long,float[]> getGrib3DAnim(final File folder, final String param, final int[] outSize, final float[] outMinMax, final float[] bbox) {
//        final Map<Long,float[]> anim = new LinkedHashMap<>();
//        String path;
//        int p0, p1;
//       // final int[] outSizeFrame = new int[3];
//        outMinMax[0] = Float.MAX_VALUE;
//        outMinMax[1] = -Float.MAX_VALUE;
//        final float[] outMinMaxFrame = new float[2];
//        for (final File fileEntry : folder.listFiles()) {
//            if (!fileEntry.isDirectory()) {
//                path = fileEntry.getPath();
//                System.out.println("Loading grib file : " + path);
//                int pos = path.lastIndexOf(File.separatorChar);
//                p0 = param == null ? path.indexOf("_", pos) : path.indexOf("_"+param+"_", pos);
//                if (p0>=0) {
//                    try {
//                        p0 += param == null ? 1 :(2+param.length());
//                        p1 = path.lastIndexOf(".");
//                        long time = Long.parseLong(path.substring(p0, p1));
//                        float[] grib3D = getGrib3D(path, outSize, outMinMaxFrame, bbox);
//                        if (outMinMaxFrame[0]<outMinMax[0]) outMinMax[0] = outMinMaxFrame[0];
//                        if (outMinMaxFrame[1]>outMinMax[1]) outMinMax[1] = outMinMaxFrame[1];
//                        anim.put(time, grib3D);
//                        
//                    } catch (IOException ex) {
//                        Logger.getLogger(DataLoader.class.getName()).log(Level.SEVERE, null, ex);
//                    }
//                }
//            }
//        }
//        return anim;
//    }
//    
//    public static float[] getGrib3DCSV(final String file, final String param, final int[] outSize, final float[] outMinMax) {
//        final double[][][] grib = TextGribLoader.loadGrib3D(file, param);
//        return to1D(grib, outSize,outMinMax);
//    }
//    
//    public static float[] getGrib3D(final String file, final int[] outSize, final float[] outMinMax, final float[] bbox) throws IOException {
//        final Grib3D grib = GribReader.decodeGrib3D(file);
//        if (bbox != null && bbox.length>=4) {
//            double[] bboxd = grib.getBBox();
//            for (int i=0;i<4;i++) bbox[i] = (float)bboxd[i];
//        }
//        return to1D(grib.getData(), outSize, outMinMax);
//    }
//
//    public static float[] getGrib2D(final String file, final int[] outSize, final float[] outMinMax, final float[] bbox) {
//        try {
//            final List<Grib2D> gribs = GribReader.decodeGrib(file);
//            if (gribs != null && !gribs.isEmpty()) {
//                Grib2D grib = gribs.get(0);
//                return getGrib2D(grib, outSize, outMinMax, bbox);
//            }
//        } catch (IOException ex) {
//            Logger.getLogger(DataLoader.class.getName()).log(Level.SEVERE, null, ex);
//        }
//        return null;
//    }
//
//    public static float[] getGrib2D(Grib2D grib, final int[] outSize, final float[] outMinMax, final float[] bbox) {
//        if (bbox != null && bbox.length>=4) {
//            double[] bboxd = grib.getBBox();
//            bbox[0] = (float)bboxd[1];
//            bbox[1] = (float)bboxd[0];
//            bbox[2] = (float)bboxd[3];
//            bbox[3] = (float)bboxd[2];
//        }
//        return to1D(grib.getData(), outSize, outMinMax);
//    }
//    
//    public static float[][] getGrib2D(final String file, final float[] bbox) {
//        try {
//            final List<Grib2D> gribs = GribReader.decodeGrib(file);
//            if (gribs != null && !gribs.isEmpty()) {
//                Grib2D grib = gribs.get(0);
//                if (bbox != null && bbox.length>=4) {
//                    double[] bboxd = grib.getBBox();
//                    bbox[0] = (float)bboxd[1];
//                    bbox[1] = (float)bboxd[0];
//                    bbox[2] = (float)bboxd[3];
//                    bbox[3] = (float)bboxd[2];
//                }
//                return grib.getData();
//            }
//        } catch (IOException ex) {
//            Logger.getLogger(DataLoader.class.getName()).log(Level.SEVERE, null, ex);
//        }
//        return null;
//    }
// 
//    public static float[] to1D(final float[][] grib, final int[] outSize, final float[] outMinMax) {
//        return to1D(grib, 0,0, grib[0].length, grib.length, outSize, outMinMax);
//    }
//    
//    public static float[] to1DYX(final float[][] grib, final int[] outSize, final float[] outMinMax) {
//        return to1DYX(grib, 0,0, grib[0].length, grib.length, outSize, outMinMax);
//    }
//    
//    public static float[] to1D(final float[][][] grib, final int[] outSize, final float[] outMinMax) {
//        return to1D(grib, 0,0,0, grib[0][0].length, grib[0].length, grib.length, outSize, outMinMax);
//    }
//    
//    public static float[] to1D(final double[][][] grib, final int[] outSize, final float[] outMinMax) {
//        return to1D(grib, 0,0,0, grib[0][0].length, grib[0].length, grib.length, outSize, outMinMax);
//    }
//    
//    public static float[] to1D(final double[][][] grib, final int xmin, final int ymin, int zmin, final int xmax, final int ymax, int zmax, 
//            final int[] outSize, final float[] outMinMax) {
//        
//        zmax = Math.min(grib.length,zmax);
//        
//        // Find min and max
//        double min = Double.MAX_VALUE, max = -Double.MAX_VALUE;
//        final int sz = zmax-zmin, sy = ymax-ymin, sx = xmax-xmin;
//
//        int off = 0;
//        double v;
//        final float[] data = new float[sx*sy*sz];
//        for(int z=zmin; z<zmax; z++) {
//            if (grib[z].length == 0) break;
//            for(int y=ymin; y<ymax; y++) {
//                for(int x=xmin; x<xmax; x++) {
//                    v = grib[z][y][x];
//                    data[off++] = (float)v;
//                    if (v<min) min = v;
//                    if (v>max) max = v;   
//                }
//            }
//        }
//        if (outMinMax != null && outMinMax.length>=2) {
//            outMinMax[0] = (float)min;
//            outMinMax[1] = (float)max;
//        }
//        if (outSize != null && outSize.length>=3) {
//            outSize[0] = sx;
//            outSize[1] = sy;
//            outSize[2] = sz;
//        }
//        return data;
//    }
//    
//    public static float[] to1D(final float[][][] grib, 
//            final int xmin, final int ymin, final int zmin, 
//            final int xmax, final int ymax, final int zmax, 
//            final int[] outSize, final float[] outMinMax) {
//        // Find min and max
//        float min = Float.MAX_VALUE, max = -Float.MAX_VALUE;
//        final int sz = zmax-zmin, sy = ymax-ymin, sx = xmax-ymin;
//        int off = 0;
//        Float v;
//        final float[] data = new float[sx*sy*sz];
//        for(int z=zmin; z<zmax; z++) {
//            for(int y=ymin; y<ymax; y++) {
//                for(int x=xmin; x<xmax; x++) {
//                    v = grib[z][y][x];
//                    data[off++] = v;
//                    if (v<min) min = v;
//                    if (v>max) max = v;   
//                }
//            }
//        }
//        if (outMinMax != null && outMinMax.length>=2) {
//            outMinMax[0] = min;
//            outMinMax[1] = max;
//        }
//        if (outSize != null && outSize.length>=3) {
//            outSize[0] = sx;
//            outSize[1] = sy;
//            outSize[2] = sz;
//        }
//        return data;
//    }
//
//    public static float[] to1D(final float[][] grib, 
//            final int xmin, final int ymin, 
//            final int xmax, final int ymax, 
//            final int[] outSize, final float[] outMinMax) {
//        // Find min and max
//        float min = Float.MAX_VALUE, max = -Float.MAX_VALUE;
//        final int sy = ymax-ymin, sx = xmax-ymin;
//        int off = 0;
//        Float v;
//        final float[] data = new float[sx*sy];
//        for(int y=ymin; y<ymax; y++) {
//            for(int x=xmin; x<xmax; x++) {
//                v = grib[y][x];
//                data[off++] = v;
//                if (v<min) min = v;
//                if (v>max) max = v;   
//            }
//        }
//        
//        if (outMinMax != null && outMinMax.length>=2) {
//            outMinMax[0] = min;
//            outMinMax[1] = max;
//        }
//        if (outSize != null && outSize.length>=2) {
//            outSize[0] = sx;
//            outSize[1] = sy;
//        }
//        return data;
//    }
//    
//    public static float[] to1DYX(final float[][] grib, final int xmin, final int ymin, final int xmax, final int ymax, final int[] outSize, final float[] outMinMax) {
//        // Find min and max
//        float min = Float.MAX_VALUE, max = -Float.MAX_VALUE;
//        final int sy = ymax-ymin, sx = xmax-ymin;
//        int off = 0;
//        Float v;
//        final float[] data = new float[sx*sy];
//        for(int x=xmax-1; x>=xmin; x--) {
//            for(int y=ymin; y<ymax; y++) {
//                v = grib[y][x];
//                data[off++] = v;
//                if (v<min) min = v;
//                if (v>max) max = v;   
//            }
//        }
//        
//        if (outMinMax != null && outMinMax.length>=2) {
//            outMinMax[0] = min;
//            outMinMax[1] = max;
//        }
//        if (outSize != null && outSize.length>=2) {
//            outSize[0] = sy;
//            outSize[1] = sx;
//        }
//        return data;
//    }
//    
//    public static float[] decale(float[] buf, float add, float mult, float[] out) {
//        for (int i=0; i<buf.length; i++) {
//            out[i] = (buf[i]+add)*mult; 
//        }
//        return buf;
//    }
//    
//    public static float[] mapValues(float[] buf, float min, float max, float[] out) {
//        float[] minMax = new float[2];
//        getMinMax(buf, minMax);
//        return decale(buf, min-minMax[0], max/(minMax[1]-minMax[0]), out);
//    }
//    
//    public static void getMinMax(float[] demBuff, float[] outMinMax) {
//        if (outMinMax != null && outMinMax.length>=2) {
//            float v, min = Float.MAX_VALUE, max = -Float.MAX_VALUE;
//            for(int off=0; off<demBuff.length; off++) {
//                v = demBuff[off];
//                if (v<min) min = v;
//                if (v>max) max = v;   
//            }
//            outMinMax[0] = min;
//            outMinMax[1] = max;
//        }
//    }
//   
//    public static Map<Long, float[]> getDemAnim(File folder, Object object, int[] demSize, float[] outMinMax, float[] dataBBox) {
//        final Map<Long,float[]> anim = new LinkedHashMap<>();
//        String path;
//        int p0, p1;
//
//        
//        if (outMinMax != null && outMinMax.length>1) {
//            outMinMax[0] = Float.MAX_VALUE;
//            outMinMax[1] = -Float.MAX_VALUE;
//        }
//        long time = System.currentTimeMillis();
//       // time = (time / DAY) * DAY;
//        
//        final float[] outMinMaxFrame = new float[2];
//        BufferedImage img = null;
//        for (final File fileEntry : folder.listFiles()) {
//            if (!fileEntry.isDirectory()) {
//                path = fileEntry.getPath();
//                
//                p0 = path.lastIndexOf("_");
//                p1 = path.lastIndexOf(".");
//                if (p0>=0) {
//                    try {
//                        p0 += 1;
//                        time = Long.parseLong(path.substring(p0,p1));
//
//                        img = ImageIO.read(fileEntry);
//                        float[] dem = getDem(img, 0, 1, outMinMaxFrame);      
//                        
//                        if (outMinMaxFrame[0]<outMinMax[0]) outMinMax[0] = outMinMaxFrame[0];
//                        if (outMinMaxFrame[1]>outMinMax[1]) outMinMax[1] = outMinMaxFrame[1];
//                        
//                        anim.put(time, dem);
//                        
//                    } catch (IOException ex) {
//                        Logger.getLogger(DataLoader.class.getName()).log(Level.SEVERE, null, ex);
//                    }
//                }
//                
////                p0 = path.indexOf("_");
////                if (p0>=0) {
////                    try {
////                        p0 += 1;
////                        long day = Long.parseLong(path.substring(p0,p0+2));
////                        long hour = Long.parseLong(path.substring(p0+3,p0+5));
////                        img = ImageIO.read(fileEntry);
////                        float[] dem = getDem(img, 0, 1, outMinMaxFrame);      
////                        
////                        if (outMinMaxFrame[0]<outMinMax[0]) outMinMax[0] = outMinMaxFrame[0];
////                        if (outMinMaxFrame[1]>outMinMax[1]) outMinMax[1] = outMinMaxFrame[1];
////                        
////                        anim.put(time+(day-24)*DAY+hour*HOUR, dem);
////                        
////                    } catch (IOException ex) {
////                        Logger.getLogger(DataLoader.class.getName()).log(Level.SEVERE, null, ex);
////                    }
////                }
//            }
//        }
//        if (img != null && demSize != null && demSize.length > 1) {
//            demSize[0] = img.getWidth();
//            demSize[1] = img.getHeight();
//        }
//        return anim;
//    }
//
// 
//    
//    
//    
//    void saveGrib3DToFile(String path, float[] data, final int[] outSize, final float[] outMinMax, final float[] bbox) {
//        try{
//            final RandomAccessFile rFile = new RandomAccessFile(path, "rw");
//            rFile.writeInt(outSize, 0, 3);      // nlon, nlat, nalti 
//            rFile.writeFloat(outMinMax, 3, 2);  // min max
//            rFile.writeFloat(bbox, 5, 4);       // bbox
//            rFile.writeFloat(data, 9, data.length);   // data   
//            rFile.close();
//        }
//        catch (IOException ex) {
//            System.err.println(ex.getMessage());
//        }
//    }
//    
//    //READ
//    float[] loadGrib3DToFile(String path, final int[] outSize, final float[] outMinMax, final float[] bbox) {
//        
//        try{
//            RandomAccessFile rFile = new RandomAccessFile("/sdcard/demo.data", "rw");
//            int count = (int) (rFile.length() / 4);
//            float[] data = new float[count-9];
//            rFile.readInt(outSize, 0, 3);      // nlon, nlat, nalti 
//            rFile.readFloat(outMinMax, 3, 2);  // min max
//            rFile.readFloat(bbox, 5, 4);       // bbox
//            rFile.readFloat(data, 9, count);   // data   
//            rFile.close();
//            return data;
//        } catch (IOException ex) {
//            System.err.println(ex.getMessage());
//        }
//        return null;
//    }
//    
//    public interface Grib3DListener {
//        public void onNewGribReceived(final SimpleGrib newSimpleGrib);
//        public void onGroundGribReceived(Grib2D simpleGrib);
//    }
//    
//    public static class SimpleGrib {
//        public long date;
//        public float hmin;
//        public float hmax;
//        public final int[] outSize = new int[4];
//        public final float[] outMinMax = new float[2];
//        public final float[] bbox = new float[4];
//        public float[] data;
//    }
//
//    public static void getGribGroundAsync(final String layerUrl, final Grib3DListener gribListener) {
//        
//        // recuperation du grib d'altitude
//        final ListenableFuture<Grib2D> futureGrib2DGround = Futures.transform(WebLoader.futureGetLayerGribGround(layerUrl), new AsyncFunction<byte[], Grib2D>() {
//            @Override
//            public ListenableFuture<Grib2D> apply(byte[] grib) throws Exception {
//                final List<Grib2D> gribs2D = GribReader.decodeGrib(grib);
//                if (gribs2D != null && !gribs2D.isEmpty()) {
//                    return Futures.immediateFuture(gribs2D.get(0));
//                }
//                return null;
//            }
//        });
//        
//        Futures.addCallback(futureGrib2DGround, new FutureCallback<Grib2D>() {
//            @Override
//            public void onSuccess(Grib2D simpleGrib) {
//                gribListener.onGroundGribReceived(simpleGrib);
//            }
//            @Override
//            public void onFailure(Throwable thrwbl) {  
//                int test = 0;
//            }
//        }, SwingExecutor.getInstance());  // Appel du listener dans l'EDT
//    }
//    
//    public static void getGrib3DStepByStepOnFile2(final String pathDataDir, final Grib3DListener gribListener) {
//        final File folder = new File(pathDataDir);
//
//        ListeningExecutorService executor = ApplicationExecutor.getInstance().getExecutor();
//        executor.submit(() -> {
//            int p0, p1, pos;
//            for (final File fileEntry : folder.listFiles()) {
//                if (!fileEntry.isDirectory()) {
//                    final String path = fileEntry.getPath();
//                    System.out.println("loading file : " + path);
//                  //  pos = path.lastIndexOf(File.separatorChar);
//                    p0 = path.lastIndexOf("_");
//                    if (p0>=0) {
//                        p1 = path.lastIndexOf(".");
//                        long time = Long.parseLong(path.substring(p0+1, p1));
//                        try {
//                            Grib3D grib = GribReader.decodeGrib3D(path);
//                            SimpleGrib simpleGrib = grib3DToSimpleGrib(time, grib);
//                            SwingUtilities.invokeLater(() -> gribListener.onNewGribReceived(simpleGrib));
//                            System.out.println("loading file : " + path + "OK");
//                        } catch (IOException ex) {
//                            Logger.getLogger(DataLoader.class.getName()).log(Level.SEVERE, null, ex);
//                            System.out.println("loading file : " + path + "NOK");
//                        }
//                    } else {
//                        System.out.println("loading file : " + path + "NOK");
//                    }
//                }
//            }
//        });
//    }
//  
//    public static void getGrib3DStepByStepOnFile(final String pathDataDir, final Grib3DListener gribListener) {
//        File folder = new File(pathDataDir);
//        Map<Long, ListenableFuture<Grib3D>> mapFuturesGrib3D = new HashMap();
//        int p0, p1, pos;
//        ListeningExecutorService executor = ApplicationExecutor.getInstance().getExecutor();
//        
//        for (final File fileEntry : folder.listFiles()) {
//            if (!fileEntry.isDirectory()) {
//                final String path = fileEntry.getPath();
//                pos = path.lastIndexOf(File.separatorChar);
//                p0 = path.indexOf("_", pos);
//                if (p0>=0) {
//                    p1 = path.lastIndexOf(".");
//                    long time = Long.parseLong(path.substring(p0+1, p1));
//                    mapFuturesGrib3D.put(time, executor.submit(() -> GribReader.decodeGrib3D(path)));
//                }
//            }
//        }
//        
//        mapTimeGribToGribListener(mapFuturesGrib3D, gribListener);
//    }
//    
//    public static void getGrib3DStepByStepOnUrl(final String layerUrl, final Grib3DListener gribListener) {
//        final ListenableFuture<Map<Long, ListenableFuture<Grib3D>>> mapFutures = WebLoader.futureGetGrib3D(layerUrl);
//        Futures.addCallback(mapFutures, new FutureCallback<Map<Long, ListenableFuture<Grib3D>>>() {
//            @Override
//            public void onSuccess(Map<Long, ListenableFuture<Grib3D>> v) {
//                mapTimeGribToGribListener(v, gribListener);
//            }
//            @Override
//            public void onFailure(Throwable thrwbl) {     
//                int test = 0;
//            }
//        });
//    }
//
//    public static void mapTimeGribToGribListener(final Map<Long, ListenableFuture<Grib3D>> mapFutures, final Grib3DListener gribListener) {
//        for (Entry<Long, ListenableFuture<Grib3D>> e : mapFutures.entrySet()) {
//            grib3DToSimpleGrib(e.getKey(), e.getValue(), gribListener);  
//        }
//    }
//     public static SimpleGrib grib3DToSimpleGrib(Long date, Grib3D grib3D) {
//        final SimpleGrib simple = new SimpleGrib();
//        simple.date = date;
//        final float[] bbox = new float[4];
//        if (bbox != null && bbox.length>=4) {
//            double[] bboxd = grib3D.getBBox();
//            for (int i=0;i<4;i++) simple.bbox[i] = (float)bboxd[i];
//        }
//        simple.data = to1D(grib3D.getData(), simple.outSize, simple.outMinMax);
//        simple.hmin = grib3D.getDx();
//        simple.hmax = simple.hmin + grib3D.getDx()*grib3D.getData().length;
//        return simple;
//    }
//
//    public static void grib3DToSimpleGrib(final Long date, final ListenableFuture<Grib3D> futureGrib3D, final Grib3DListener gribListener) {
//        
//        ListenableFuture<SimpleGrib> futureSimpleGrib = Futures.transform(futureGrib3D, new AsyncFunction<Grib3D,SimpleGrib>() {
//            @Override
//            public ListenableFuture<SimpleGrib> apply(Grib3D grib3D) throws Exception {
//                SimpleGrib simple = grib3DToSimpleGrib(date, grib3D);
//                return Futures.immediateFuture(simple);
//            }
//        });
//
//        Futures.addCallback(futureSimpleGrib, new FutureCallback<SimpleGrib>() {
//            @Override
//            public void onSuccess(SimpleGrib simpleGrib) {
//                gribListener.onNewGribReceived(simpleGrib);
//            }
//            @Override
//            public void onFailure(Throwable thrwbl) {
//                int test = 0;
//            }
//        }, SwingExecutor.getInstance());  // Appel du listener dans l'EDT
//    }
//    
//}
//
//
//    
