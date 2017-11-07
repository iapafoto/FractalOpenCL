/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.api;

import com.nativelibs4java.opencl.CLBuffer;
import com.nativelibs4java.opencl.CLBuildException;
import com.nativelibs4java.opencl.CLContext;
import com.nativelibs4java.opencl.CLDevice;
import com.nativelibs4java.opencl.CLEvent;
import com.nativelibs4java.opencl.CLImage2D;
import com.nativelibs4java.opencl.CLKernel;
import com.nativelibs4java.opencl.CLMem;
import com.nativelibs4java.opencl.CLPlatform;
import com.nativelibs4java.opencl.CLProgram;
import com.nativelibs4java.opencl.CLQueue;
import com.nativelibs4java.opencl.JavaCL;
import java.awt.image.BufferedImage;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import fractalcl.gui.CLEditor;
import org.bridj.Pointer;
import org.jocl.CL;
import org.jocl.Sizeof;
import org.jocl.cl_mem;

/**
 *
 * @author durands
 */
public class OpenCLWithJavaCL {

    private CLContext clContext = null;
    private CLQueue clQueue = null;
    private CLKernel clKernel = null;
    private CLImage2D imageOut;
    private CLProgram program;

    Map<Integer, CLImage2D> inputsImg = new HashMap();
    Map<Integer, CLBuffer<Float>> inputsFloatBuffer = new HashMap();
    Map<Integer, Map<Long, CLBuffer<Float>>> inputsTimedFloatBuffer = new HashMap<>();
    
    public void setImageInOnArg(int id) {
        CLImage2D clImg = inputsImg.get(id);
        if (clImg != null) {
            clKernel.setArg(id, clImg);
        }
    }

    public void createInput(int id, BufferedImage bi) {
        CLImage2D clImg = inputsImg.get(id);
        if (clImg != null) {
            clImg.release();
        }
        if (bi != null) {
            clImg = clContext.createImage2D(CLMem.Usage.Input, bi, false);
            inputsImg.put(id, clImg);
        }
    }

    public void setFloatBufferOnArg(int id) {
        CLBuffer<Float> buff = inputsFloatBuffer.get(id);
        if (buff != null) {
            clKernel.setArg(id, buff);
        }
    }
    
    public Set<Long> getTimeOnArg(int id) {
        Map<Long, CLBuffer<Float>> map = inputsTimedFloatBuffer.get(id);
        if (map != null) {
            return map.keySet();
        }
        return new HashSet();
    }
    
    public void setTimedFloatBufferOnArg(int id, Long time) {
        Map<Long, CLBuffer<Float>> map = inputsTimedFloatBuffer.get(id);
        if (map != null && !map.isEmpty()) {
            CLBuffer<Float> buff = null;
            if (time == null) {
                buff = map.values().iterator().next();
            }else {
                buff = map.get(time);
            }
            if (buff != null) {
                clKernel.setArg(id, buff);
            }
        }
    }
    
    public void createTimedInput(int id, long time, final float[] array) {
        Map<Long, CLBuffer<Float>> map = inputsTimedFloatBuffer.get(id);
        if (map == null) {
            map = new HashMap();
            inputsTimedFloatBuffer.put(id, map);
        }
        CLBuffer<Float> buff = map.get(time);
        if (buff != null) {
            buff.release();
        }
        if (array != null) {
            buff = clContext.createBuffer(CLMem.Usage.Input, Pointer.pointerToFloats(array), true);
            map.put(time, buff);
        }
    }
    
    public void createInput(int id, final float[] array) {
        CLBuffer<Float> buff = inputsFloatBuffer.get(id);
        if (buff != null) {
            buff.release();
        }
        if (array != null) {
            buff = clContext.createBuffer(CLMem.Usage.Input, Pointer.pointerToFloats(array), true);
            inputsFloatBuffer.put(id, buff);
        }
    }

    public boolean isKernelOk() {
        return clKernel != null;
    }

    public void onScreenResize(BufferedImage img) {
        if (imageOut == null || imageOut.getWidth() != img.getWidth() || imageOut.getHeight() != img.getHeight()) {
            if (imageOut != null) {
                imageOut.release();
            }
            if (clContext != null) {
                imageOut = clContext.createImage2D(CLMem.Usage.Output, img, true);
            }
        }
    }

    public boolean createProgram(String codeSrc, String[] error) {
        try {
            program = clContext.createProgram(codeSrc);
            if (program != null) {
                program.build();
                return true;
            }
        } catch (CLBuildException e) {
            error[0] = e.getMessage();
        }
        return false;
    }

    public boolean createKernels() {
        CLKernel[] kernels = program.createKernels();
        if (kernels != null && kernels.length > 0) {
            clKernel = kernels[0];
            return true;
        }
        return false;
    }

    public void initOpenCL() {
        CLDevice bestGPUDevice = null, bestCPUDevice = null;
        for (CLPlatform p : JavaCL.listPlatforms()) {
            System.out.println("---- PLATEFORM: " + p.getName() + "---------");

            for (CLDevice device : p.listGPUDevices(false)) {
                System.out.println("GPU DEVICE: " + device.getName() + " " + (device.isAvailable() ? "available" : "not available"));
                if (device.isAvailable()) {
                    System.out.println("   version:      " + device.getVersion());
                    System.out.println("   OpenCL:       " + device.getOpenCLCVersion());
                    System.out.println("   Nb Comp Unit: " + device.getMaxComputeUnits());
                    System.out.println("   Clock Freq:   " + device.getMaxClockFrequency());
                    System.out.println("   MaxParameterSize:   " + device.getMaxParameterSize());
                    System.out.println("   MaxConstantArgs:    " + device.getMaxConstantArgs());
                    System.out.println("   MaxSamplerss:       " + device.getMaxSamplers());
                    System.out.println("   NativeVectorWidthFloat:   " + device.getNativeVectorWidthFloat());
                    System.out.println("   PreferredVectorWidthFloat:" + device.getPreferredVectorWidthFloat());
                    System.out.println("   LocalMemSize:             " + device.getLocalMemSize());

                    if (bestGPUDevice == null) {
                        bestGPUDevice = device;
                    }
                }

            }
            for (CLDevice device : p.listCPUDevices(false)) {
                System.out.println("CPU DEVICE: " + device.getName() + " " + (device.isAvailable() ? "available" : "not available"));
                if (device.isAvailable()) {
                    System.out.println("   version:      " + device.getVersion());
                    System.out.println("   OpenCL:       " + device.getOpenCLCVersion());
                    System.out.println("   Nb Comp Unit: " + device.getMaxComputeUnits());
                    System.out.println("   Clock Freq:   " + device.getMaxClockFrequency());
                    System.out.println("   MaxParameterSize:   " + device.getMaxParameterSize());
                    System.out.println("   MaxConstantArgs:    " + device.getMaxConstantArgs());
                    System.out.println("   MaxSamplerss:       " + device.getMaxSamplers());
                    System.out.println("   NativeVectorWidthFloat:   " + device.getNativeVectorWidthFloat());
                    System.out.println("   PreferredVectorWidthFloat:" + device.getPreferredVectorWidthFloat());
                    System.out.println("   LocalMemSize:             " + device.getLocalMemSize());
                    if (bestCPUDevice == null) {
                        bestCPUDevice = device;
                    }
                }
            }
            // bestCPUDevice = p.getBestDevice();
        }
        
        //    if (bestCPUDevice != null) {
        //    clContext = JavaCL.createContext(null, bestCPUDevice);
        //    clQueue = clContext.createDefaultQueue();
        //    System.out.println("choose => CPU DEVICE: " + bestCPUDevice.getName() + " (" + bestCPUDevice.getVersion()+")");
        //    } else
          if (bestGPUDevice != null) {
            clContext = JavaCL.createContext(null, bestGPUDevice);
            clQueue = clContext.createDefaultQueue();
            System.out.println("choose => GPU DEVICE: " + bestGPUDevice.getName() + " (" + bestGPUDevice.getVersion() + ")");
        }
        clKernel = null;
    }

    public void regenerate(boolean antialiasing, BufferedImage imgout, CLEditor editor) {
        if (isKernelOk() && !clQueue.isOutOfOrder()) {
            try {
                clQueue.finish();
                long t0 = System.nanoTime();
                editor.updateKernelArgsRenderingMode(antialiasing);

                final int   w = (int) imageOut.getWidth(),
                            h = (int) imageOut.getHeight();

                CLEvent evt = clKernel.enqueueNDRange(clQueue, new int[]{w, h}/*, new int[]{32,32}*/);
                if (evt != null) {
                    clQueue.enqueueWaitForEvents(evt);
                    imageOut.read(clQueue, imgout, false, evt);
                
                    long dt = System.nanoTime() - t0;
                    System.out.println("Image (" + w + "x" + h + ") generated in " + (dt / 1000000) + "ms");
                }

            } catch (Exception e) {
                int test;
            }

        }
    }

    public CLBuffer<Float> floatArrayToArgInput(final float[] array) {

        return clContext.createBuffer(CLMem.Usage.Input, Pointer.pointerToFloats(array), true);
    }

    public CLBuffer<Float> floatArrayToArgOutput(final float[] array) {
        return clContext.createBuffer(CLMem.Usage.Output, Pointer.pointerToFloats(array), true);
    }

    public void setArg(int i, int[] arg) {
        clKernel.setArg(i, arg);
    }
    public void setArg(int id, float[] f) {
        clKernel.setArg(id, f);
    }
    public void setArg(int i, float arg) {
        clKernel.setArg(i, arg);
    }
    public void setArg(int i, int arg) {
        clKernel.setArg(i, arg);
    }
    public void setArg(int i, CLBuffer<Float> arg) {
        clKernel.setArg(i, arg);
    }

    public void setImageOutOnArg(int id, boolean antialiasing) {
        if (clKernel != null && imageOut != null) {
            clKernel.setArg(id, imageOut);
        }
//        if (clKernel != null && imageOutHD != null && imageOut != null) {
//            clKernel.setArg(i, antialiasing ? imageOutHD : imageOut);
//        }
    }


}
