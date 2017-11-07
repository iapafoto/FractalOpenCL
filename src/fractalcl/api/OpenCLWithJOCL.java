/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fractalcl.api;

import java.awt.Transparency;
import java.awt.color.ColorSpace;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.ComponentColorModel;
import java.awt.image.DataBuffer;
import java.awt.image.DataBufferFloat;
import java.awt.image.DataBufferInt;
import java.awt.image.PixelInterleavedSampleModel;
import java.awt.image.Raster;
import java.awt.image.SampleModel;
import java.awt.image.WritableRaster;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import fractalcl.gui.CLEditor;


import org.jocl.*;
import static org.jocl.CL.CL_CONTEXT_PLATFORM;
import static org.jocl.CL.CL_DEVICE_NAME;
import static org.jocl.CL.CL_DEVICE_TYPE_ALL;
import static org.jocl.CL.CL_DEVICE_VENDOR;
import static org.jocl.CL.CL_DRIVER_VERSION;
import static org.jocl.CL.CL_MEM_COPY_HOST_PTR;
import static org.jocl.CL.CL_MEM_READ_ONLY;
import static org.jocl.CL.CL_MEM_USE_HOST_PTR;
import static org.jocl.CL.CL_MEM_WRITE_ONLY;
import static org.jocl.CL.CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE;
import static org.jocl.CL.CL_QUEUE_PROFILING_ENABLE;
import static org.jocl.CL.CL_RGBA;
import static org.jocl.CL.CL_UNSIGNED_INT8;
import static org.jocl.CL.clBuildProgram;
import static org.jocl.CL.clCreateBuffer;
import static org.jocl.CL.clCreateCommandQueue;
import static org.jocl.CL.clCreateContext;
import static org.jocl.CL.clCreateKernel;
import static org.jocl.CL.clCreateProgramWithSource;
import static org.jocl.CL.clGetDeviceIDs;
import static org.jocl.CL.clGetDeviceInfo;
import static org.jocl.CL.clGetPlatformIDs;




/**
 *
 * @author durands
 */
public class OpenCLWithJOCL {

    private cl_context context;
    private cl_command_queue commandQueue;
    private cl_program program;
    private cl_kernel kernel;
    private cl_mem imageOut;
    private cl_mem zBufferOut;
    
    Map<Integer, cl_mem> inputsImg = new HashMap();
    Map<Integer, cl_mem> inputsOutputsFloatBuffer = new HashMap();
    Map<Integer, Map<Long, cl_mem>> inputsTimedClMem = new HashMap();
    
    public void setImageInOnArg(final int id) {
        cl_mem clImg = inputsImg.get(id);
        if (clImg != null) {
            CL.clSetKernelArg(kernel, id, Sizeof.cl_mem, Pointer.to(clImg));
        }
    }

    public void createInput(final int id, final BufferedImage bi) {
        cl_mem clImg = inputsImg.get(id);
        if (clImg != null) {
            CL.clReleaseMemObject(clImg);
        }
        if (bi != null) {
            clImg = initImageMemFloat(bi);
            inputsImg.put(id, clImg);
        } else {
            inputsImg.remove(id);
        }
    }
    
    public void createInputAs3DTexture(int id, float[] demBuff, int[] demSize) {
        cl_mem clImg = inputsImg.get(id);
        if (clImg != null) {
            CL.clReleaseMemObject(clImg);
        }
        if (demBuff != null) {
            clImg = initImage3DMemFloat(demBuff, demSize);
            inputsImg.put(id, clImg);
        } else {
            inputsImg.remove(id);
        }
    }    
    
    public void createInputAs2DTexture(int id, float[] demBuff, int[] demSize) {
        cl_mem clImg = inputsImg.get(id);
        if (clImg != null) {
            CL.clReleaseMemObject(clImg);
        }
        if (demBuff != null) {
            clImg = initImage2DMemFloat(demBuff, demSize);
            inputsImg.put(id, clImg);
        } else {
            inputsImg.remove(id);
        }
    } 
    
     public void createTimedInputAs3DTexture(final int id, long time, final float[] array, final int[] size) {
        Map<Long, cl_mem> map = inputsTimedClMem.get(id);
        if (map == null) {
            map = new HashMap();
            inputsTimedClMem.put(id, map);
        }
        cl_mem buff = map.get(time);
        if (buff != null) {
            CL.clReleaseMemObject(buff);
        }
        if (array != null) {
            buff = initImage3DMemFloat(array, size);
            map.put(time, buff);
        }
    }
    
    public Set<Long> getTimeOnArg(final int id) {
        Map<Long, cl_mem> map = inputsTimedClMem.get(id);
        if (map != null) {
            return map.keySet();
        }
        return new HashSet();
    }
    
    public void setTimedFloatBufferOnArg(final int id, final Long time) {
        Map<Long, cl_mem> map = inputsTimedClMem.get(id);
        if (map != null && !map.isEmpty()) {
            cl_mem buff = null;
            if (time == null) {
                buff = map.values().iterator().next();
            }else {
                buff = map.get(time);
            }
            if (buff != null) {
                CL.clSetKernelArg(kernel, id, Sizeof.cl_mem, Pointer.to(buff));
            }
        }
    }
    
    public void createTimedInput(final int id, long time, final float[] array) {
        Map<Long, cl_mem> map = inputsTimedClMem.get(id);
        if (map == null) {
            map = new HashMap();
            inputsTimedClMem.put(id, map);
        }
        cl_mem buff = map.get(time);
        if (buff != null) {
            CL.clReleaseMemObject(buff);
        }
        if (array != null) {
            buff = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Sizeof.cl_float * array.length, Pointer.to(array), null);
            map.put(time, buff);
        }
    }
    
    public void setFloatBufferOnArg(final int id) {
        cl_mem buff = inputsOutputsFloatBuffer.get(id);
        if (buff != null) {
            CL.clSetKernelArg(kernel, id, Sizeof.cl_mem, Pointer.to(buff));
        }
    }
    
    public void createInput(final int id, final float[] array) {
        cl_mem buff = inputsOutputsFloatBuffer.get(id);
        if (buff != null) {
            CL.clReleaseMemObject(buff);
        }
        if (array != null) {
            buff = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Sizeof.cl_float * array.length, Pointer.to(array), null);
            inputsOutputsFloatBuffer.put(id, buff);
        } else {
            inputsOutputsFloatBuffer.remove(id);
        }
    }
    public void createOutput(final int id, final float[] array) {
        cl_mem buff = inputsOutputsFloatBuffer.get(id);
        if (buff != null) {
            CL.clReleaseMemObject(buff);
        }
        if (array != null) {
            buff = clCreateBuffer(context, CL_MEM_WRITE_ONLY, Sizeof.cl_float * array.length, Pointer.to(array), null);
            inputsOutputsFloatBuffer.put(id, buff);
        } else {
            inputsOutputsFloatBuffer.remove(id);
        }
    }
    
    public boolean isKernelOk() {
        return kernel != null;
    }
    int imageOutWidth = 0, imageOutHeight = 0;
   // int imageOutHDWidth = 0, imageOutHDHeight = 0;

    public void onScreenResize(final BufferedImage img, float[] zBufferArray) {//, BufferedImage imageTempHQ) {
        if (imageOut == null || imageOutWidth != img.getWidth() || imageOutHeight != img.getHeight()) {
            if (imageOut != null) {
                CL.clReleaseMemObject(imageOut);
            }
            if (context != null) {
                imageOutWidth = img.getWidth();
                imageOutHeight = img.getHeight();
                imageOut = initImageOut(img.getWidth(), img.getHeight());

                if (zBufferOut != null) CL.clReleaseMemObject(zBufferOut);
                zBufferOut = clCreateBuffer(context, CL_MEM_WRITE_ONLY, Sizeof.cl_float * zBufferArray.length, Pointer.to(zBufferArray), null);            
            }
        }
    }

    public boolean createProgram(final String programSource, final String[] error) {
       try {
            // Create the program
            System.out.println("Creating program...");
            program = clCreateProgramWithSource(context, 1, new String[]{ programSource }, null, null);

            // Build the program
            System.out.println("Building program...");
            clBuildProgram(program, 0, null, null, null, null);
        } catch (CLException e) {
            error[0] = e.getMessage();
            return false;
        }
        return true;
    }

    public boolean createKernels() {
          // Create the kernel
        System.out.println("Creating kernel...");
        kernel = clCreateKernel(program, "render", null);  

        return true;
    }

    public void initOpenCL() {
        
        final int platformIndex = 1;
        final long deviceType = CL_DEVICE_TYPE_ALL;
        final int deviceIndex = 0;

        // Enable exceptions and subsequently omit error checks in this sample
        CL.setExceptionsEnabled(true);

        // Obtain the number of platforms
        int numPlatformsArray[] = new int[1];
        clGetPlatformIDs(0, null, numPlatformsArray);
        int numPlatforms = numPlatformsArray[0];

        // Obtain a platform ID
        cl_platform_id platforms[] = new cl_platform_id[numPlatforms];
        clGetPlatformIDs(platforms.length, platforms, null);
        cl_platform_id platform = platforms[platformIndex];

        // Initialize the context properties
        cl_context_properties contextProperties = new cl_context_properties();
        contextProperties.addProperty(CL_CONTEXT_PLATFORM, platform);
        
        // Obtain the number of devices for the platform
        int numDevicesArray[] = new int[1];
        clGetDeviceIDs(platform, deviceType, 0, null, numDevicesArray);
        int numDevices = numDevicesArray[0];
        
        // Obtain a device ID 
        cl_device_id devices[] = new cl_device_id[numDevices];
        clGetDeviceIDs(platform, deviceType, numDevices, devices, null);
        cl_device_id device = devices[deviceIndex];

        // CL_DEVICE_NAME
        String deviceName = getString(device, CL_DEVICE_NAME);
        System.out.println("--- Info for device "+deviceName+": ---");
        System.out.printf("CL_DEVICE_NAME: \t\t\t%s\n", deviceName);

        // CL_DEVICE_VENDOR
        String deviceVendor = getString(device, CL_DEVICE_VENDOR);
        System.out.printf("CL_DEVICE_VENDOR: \t\t\t%s\n", deviceVendor);

        // CL_DRIVER_VERSION
        String driverVersion = getString(device, CL_DRIVER_VERSION);
        System.out.printf("CL_DRIVER_VERSION: \t\t\t%s\n", driverVersion);
        
        // Create a context for the selected device
        context = clCreateContext(contextProperties, 1, new cl_device_id[]{device}, null, null, null);
                   
        // Create a command-queue
        System.out.println("Creating command queue...");
        long properties = 0;
        properties |= CL_QUEUE_PROFILING_ENABLE;
        properties |= CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE;
        commandQueue = clCreateCommandQueue(context, device, properties, null);
    }

    public void regenerate(final boolean antialiasing, final BufferedImage imgout, final float[] zBuffer, final CLEditor editor) {
        if (isKernelOk()) {
            final int w = imgout.getWidth(),//imageOutWidth,
                      h = imgout.getHeight(); //imageOutHeight;

            long t0 = System.nanoTime();
            CL.clFinish(commandQueue);
            
            editor.updateKernelArgsRenderingMode(antialiasing);
           
            CL.clEnqueueNDRangeKernel(commandQueue, kernel, 2, null, new long[]{w, h}, null, 0, null, null);

            // Read the pixel data into the output image
//          DataBufferInt dataBufferDst = (DataBufferInt)imgout.getRaster().getDataBuffer();
            DataBufferFloat dataBufferDst = (DataBufferFloat)imgout.getRaster().getDataBuffer();
            float dataDst[] = dataBufferDst.getData();
     
            CL.clEnqueueReadImage(
                commandQueue, imageOut, true, new long[3], new long[]{w,h,1},
                w * 4* Sizeof.cl_float /*Sizeof.cl_uint*/, 0, Pointer.to(dataDst), 0, null, null);
            
            if (!antialiasing && zBuffer != null) {
                CL.clEnqueueReadBuffer(commandQueue, zBufferOut, CL.CL_TRUE, 0,
                    w*h * Sizeof.cl_float, Pointer.to(zBuffer), 0, null, null); 
            }
            
            long dt = System.nanoTime() - t0;
            System.out.println("Image (" + w + "x" + h + ") generated in " + (dt / 1000000) + "ms");
        }
    }
    public void setArg(int id, int[] arg) {
        CL.clSetKernelArg(kernel, id, arg.length * Sizeof.cl_int, Pointer.to(arg));
    }
    public void setArg(int id, float[] arg) {
        CL.clSetKernelArg(kernel, id, arg.length * Sizeof.cl_float, Pointer.to(arg));        
    }
    public void setArg(int id, float arg) {
        CL.clSetKernelArg(kernel, id, Sizeof.cl_float, Pointer.to(new float[] {arg}));
    }
    public void setArg(int id, int arg) {
        CL.clSetKernelArg(kernel, id, Sizeof.cl_int, Pointer.to(new int[] {arg}));
    }
    
    public void setImageOutOnArg(int id, boolean antialiasing) {
        if (kernel != null && imageOut != null) {
            CL.clSetKernelArg(kernel, id, Sizeof.cl_mem, Pointer.to(imageOut));
        }
    }
    public void setZBufferOnArg(int id) {
        if (kernel != null && zBufferOut != null) {
            CL.clSetKernelArg(kernel, id, Sizeof.cl_mem, Pointer.to(zBufferOut));
        }
    }    
    /*
    void readImageOut(final cl_mem imgMem, final BufferedImage img) {
        // Read the pixel data into the output image
        DataBufferInt dataBufferDst = (DataBufferInt)img.getRaster().getDataBuffer();
        int dataDst[] = dataBufferDst.getData();
        
        CL.clEnqueueReadImage(
            commandQueue, imgMem, true, new long[3], new long[]{img.getWidth(), img.getHeight(), 1},
            img.getWidth() * Sizeof.cl_uint, 0, Pointer.to(dataDst), 0, null, null);
    }
    */
    private cl_mem initImageMem(BufferedImage img) {
        cl_image_format imageFormat = new cl_image_format();
        imageFormat.image_channel_order = CL_RGBA;
        imageFormat.image_channel_data_type = CL_UNSIGNED_INT8;
        
        // Create the memory object for the input image
        DataBufferInt dataBufferSrc = (DataBufferInt)img.getRaster().getDataBuffer();
        int dataSrc[] = dataBufferSrc.getData();

        cl_mem inputImageMem = CL.clCreateImage2D(
            context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
            new cl_image_format[]{imageFormat}, img.getWidth(), img.getHeight(),
            img.getWidth() * Sizeof.cl_uint, Pointer.to(dataSrc), null);
        
        return  inputImageMem;
    }
    
    private cl_mem initImage2DMemFloat(float[] data, int[] size) {
        cl_image_format imageFormat = new cl_image_format();
        imageFormat.image_channel_order = CL.CL_INTENSITY;
        imageFormat.image_channel_data_type = CL.CL_FLOAT;
        
        // Create the memory object for the input image
        int[] errorCode = {0};
        cl_mem inputImageMem = CL.clCreateImage2D(
            context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
            new cl_image_format[]{imageFormat}, size[0], size[1],
            size[0] * Sizeof.cl_float, Pointer.to(data), errorCode);
        
        return inputImageMem;
    }
    
    private cl_mem initImage3DMemFloat(float[] data, int[] size) {
        cl_image_format imageFormat = new cl_image_format();
        imageFormat.image_channel_order = CL.CL_INTENSITY;
        imageFormat.image_channel_data_type = CL.CL_FLOAT;
        
        // Create the memory object for the input image
        int[] errorCode = {0};
        cl_mem inputImageMem = CL.clCreateImage3D(
            context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
            new cl_image_format[]{imageFormat}, size[0], size[1], size[2],
            size[0] * Sizeof.cl_float, size[0] * size[1] * Sizeof.cl_float, Pointer.to(data), errorCode);
        
        return  inputImageMem;
    }
    
    private cl_mem initImageMemFloat(BufferedImage img) {
        cl_image_format imageFormat = new cl_image_format();
        imageFormat.image_channel_order = CL_RGBA;
        imageFormat.image_channel_data_type = CL.CL_FLOAT;
        
        // Create the memory object for the input image
        DataBufferInt dataBufferSrc = (DataBufferInt)img.getRaster().getDataBuffer();
        int dataSrc[] = dataBufferSrc.getData();

        float[] src = new float[dataSrc.length*4];
        int r, g, b, a, off4;
        for (int off=0; off<dataSrc.length; off++) {
            a = (dataSrc[off]>> 24) & 0xFF;
            r = (dataSrc[off]>> 16) & 0xFF;
            g = (dataSrc[off]>> 8) & 0xFF;
            b = (dataSrc[off]) & 0xFF;
            off4 = off*4;
            src[off4] =   (((float)a)/256.f);
            src[off4++] = (((float)r)/256.f);
            src[off4++] = (((float)g)/256.f);
            src[off4++] = (((float)b)/256.f);
        }
        
        cl_mem inputImageMem = CL.clCreateImage2D(
            context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
            new cl_image_format[]{imageFormat}, img.getWidth(), img.getHeight(),
            img.getWidth() * Sizeof.cl_float * 4, Pointer.to(src), null);
        
        return  inputImageMem;
    }
  
    private cl_mem initImageOut(int w, int h) {
        cl_image_format imageFormat = new cl_image_format();
        imageFormat.image_channel_order = CL_RGBA;
        imageFormat.image_channel_data_type = CL.CL_FLOAT; //CL_UNSIGNED_INT8; 
        return CL.clCreateImage2D(context, CL_MEM_WRITE_ONLY,  new cl_image_format[]{imageFormat}, w, h, 0, null, null);
    }
    
    public BufferedImage createFloatBufferedImage(int w, int h, int bands) {
        // Define dimensions and layout of the image
        //int bands = 4; // 4 bands for ARGB, 3 for RGB etc
        int[] bandOffsets = {0, 1, 2, 3}; // length == bands, 0 == R, 1 == G, 2 == B and 3 == A

        // Create a TYPE_FLOAT sample model (specifying how the pixels are stored)
        SampleModel sampleModel = new PixelInterleavedSampleModel(DataBuffer.TYPE_FLOAT, w, h, bands, w  * bands, bandOffsets);
        // ...and data buffer (where the pixels are stored)
        DataBuffer buffer = new DataBufferFloat(w * h * bands);

        // Wrap it in a writable raster
        WritableRaster raster = Raster.createWritableRaster(sampleModel, buffer, null);

        // Create a color model compatible with this sample model/raster (TYPE_FLOAT)
        // Note that the number of bands must equal the number of color components in the 
        // color space (3 for RGB) + 1 extra band if the color model contains alpha 
        ColorSpace colorSpace = ColorSpace.getInstance(ColorSpace.CS_sRGB);
        ColorModel colorModel = new ComponentColorModel(colorSpace, true, false, Transparency.TRANSLUCENT, DataBuffer.TYPE_FLOAT);

        // And finally create an image with this raster
        return new BufferedImage(colorModel, raster, colorModel.isAlphaPremultiplied(), null);
    }


    /**
     * Returns the value of the device info parameter with the given name
     *
     * @param device The device
     * @param paramName The parameter name
     * @return The value
     */
    private static String getString(cl_device_id device, int paramName)
    {
        // Obtain the length of the string that will be queried
        long size[] = new long[1];
        clGetDeviceInfo(device, paramName, 0, null, size);

        // Create a buffer of the appropriate size and fill it with the info
        byte buffer[] = new byte[(int)size[0]];
        clGetDeviceInfo(device, paramName, buffer.length, Pointer.to(buffer), null);

        // Create a string from the buffer (excluding the trailing \0 byte)
        return new String(buffer, 0, buffer.length-1);
    }


}
