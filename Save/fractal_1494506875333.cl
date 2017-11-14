#define TEXTURE_IN_FLOAT
//#define TEXTURE_OUT_FLOAT

#define DE_RAY pseudoKleinianSlider
#define DE_COLOR pseudoKleinianSliderColor

// #define DE kleinian
// #define DE deSpace
// #define DE apollonian

#define WITH_SHADOWS
//#define WITH_SUN
#define WITH_AO
//#define ONLY_AO
//#define WITH_VIGNETING
//#define WITH_DEPTH_OF_FIELD
//#define WITH_ISOLINE

#ifdef ONLY_AO
    #define BACK_COLOR 0.1f*(float3)(.8f, .8f, .8f) 
#else 
 //#define BACK_COLOR 0.1f*(float3)(.8f, .8f, .8f) 
    #define BACK_COLOR 3.f*(float3)(.1f, .16f, .34f) 
//#define BACK_COLOR (float3)(1.f)
//    #define BACK_COLOR (float3)(.8f, .65f, .54f) 
#endif

#define FOCUSBLUR .0002f

#define PRECISION_FACTOR 1e-4
#define MIN_DIST_RAYMARCHING .01f
#define MAX_DIST_RAYMARCHING 11.f

#define MIN_DIST_SHADOW 10.f*PRECISION_FACTOR
#define MAX_DIST_SHADOW 1.f
#define PRECISION_FACTOR_SHADOW 3e-4

#define MIN_DIST_AO 4.f*PRECISION_FACTOR
#define MAX_DIST_AO .05f
#define PRECISION_FACTOR_AO 3e-4

#define NB_ITERATION 7

#define NB_RAY_ITER 256

#define LIGHT_VEC normalize((float3)(.2,.7, 1.6) )

const sampler_t sampler_linear = CLK_NORMALIZED_COORDS_FALSE |  CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;



typedef struct Context {
    const float4 mins;
    const float4 maxs;
} Context;


inline float3 reflect(const float3 i, const float3 n){
  return i - 2.f * n * dot(n,i);
}

// -------------------------------------------------------------------

inline float hash1(const float seed) {
    float fractptr;
    return fract(sin(seed)*43758.5453123f, &fractptr);
}
inline float2 hash2(const float seed) {
    float2 fractptr;
    return fract(sin((float2)(seed*43758.5453123f,(seed+.1)*22578.1459123f)), &fractptr);
}
inline float3 hash3(const float seed) {
    float3 fractptr;
    return fract(sin((float3)(seed,seed+.1f,seed+.2f))*(float3)(43758.5453123f,22578.1459123f,19642.3490423f), &fractptr);
}

inline float3 randomCosWeightedHemiSphereVector(const float3 n, const float seed) {
    const float r1 = 2.f*M_PI_F*hash1(seed), r2 = hash1(seed+.1f);
    const float3 u = normalize(cross(n, (float3)(0.f,1.f,1.f))), v = cross(u,n);
    return normalize(sqrt(r2)*(cos(r1)*u + sin(r1)*v) + sqrt(1.f-r2)*n);
}

inline float3 randomSphereDirection(const float seed) {
    float2 r = 2.f*M_PI_F*hash2(seed);
    return (float3)(sin(r.x)*(float2)(sin(r.y),cos(r.y)),cos(r.x));
}

inline float3 randomHemisphereDirection(const float3 n, const float seed) {
    float2 r = 2.f*M_PI_F*hash2(seed);
    float3 dr =  (float3)(sin(r.x)*(float2)(sin(r.y),cos(r.y)),cos(r.x));
    float k = dot(dr,n);
    return k == 0.f ? n : normalize(k*dr);
}


int box(const float3 ro, const float3 rd, const float3 sz, float* tN, float* tF, float3* n) {
    const float3 m = 1.f/rd,  k = fabs(m)*sz,  a = -m*ro-k*.5f, b = a+k;
    *n = -sign(rd) * step(a.yzx,a.xyz) * step(a.zxy,a.xyz);
    *tN = max(max(a.x,a.y),a.z);
    *tF = min(min(b.x,b.y),b.z);
    return *tN<*tF && *tF>0.f ? 1.f : 0.f;
}

inline float map(const float a0, const float b0, const float a1, const float b1, const float v) {
    return mix(a1,b1,(v-a0)/(b0-a0));
}

// -----------------------------------------------------

float screen(const float3 p) {
	const float d1=length(p.yz-(float2)(.25f,0.f))-.5f,	
                    d2=length(p.yz-(float2)(.25f,2.f))-.5f;	
	return min(max(d1,fabs(p.x-.3f)-.01f),max(d2,fabs(p.x+2.3f)-.01f));
}


//----------------------------------------------------------------------------------------

//mat2 rmx;
//int rotater=-1f;
//__constant float3 CSize2 = (float3)(0.63248f,0.78632f,0.875f);

float sdBox( float3 p, float3 b ) {
  float3 d = fabs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

float sdSponge(float3 z){
    z /= .2;
    //folding
    for(int n=0;n <4;n++) {
       z = fabs(z);
       z.xy = (z.x<z.y) ? z.yx : z.xy;
       z.xz = (z.x<z.z) ? z.zx : z.xz;
       z.zy = (z.y<z.z) ? z.yz : z.zy;	 
       z = z*3.-2.;
       z.z += (z.z<-1.) ? 2. : 0.;
    }
    //distance to cube
    z = fabs(z) - (float3)(1.);
    float dis = min(max(z.x, max(z.y,z.z)),0.) + length(max(z,0.)); 
    //scale cube size to iterations
    return dis *.2f* pow(3.f, -3.f); 
}

float pseudoKleinianSlider(float3 p, const Context* ctx) {//knighty's pseudo kleinian
    float k,r2, scale=1.f, orb = 1.f;
    for(int i=0;i<NB_ITERATION;i++) {
       // if(i==rotater)p.xy=p.xy*rmx;
        p = 2.f*clamp(p, ctx->mins.xyz, ctx->maxs.xyz)-p;
        r2 = dot(p,p);
  //float d = sdBox(p, (float3)(.1f*scale));
        k = max(ctx->mins.w/(dot(p,p)),1.f);
        p *= k;
        scale *= k;
    }
    float rxy=length(p.xy);
	float d2 = ((rxy*p.z) / length(p))/scale;
    //float d0 = max(-rxy+ctx->maxs.w, /*fabs*/(rxy*p.z) / length(p))/scale;


p /= scale;
p *= 20.f;
	float d1 = sdSponge(p)/20.f;

   return .1f*max(d2, d1);
   //return .2f*sdBox(p, (float3)(.005f*scale))/scale;
  //  float rxy=length(p.xy);
  //  return .25f*max(d - ctx->maxs.w, d)/scale;
}


float3 pseudoKleinianSliderColor(float3 p, const Context* ctx) {//knighty's pseudo kleinian
    float k,r2, scale = 1.f, orb = 1.f;
    for(int i=0;i<7;i++) {
        p = 2.f*clamp(p, ctx->mins.xyz, ctx->maxs.xyz)-p;
        r2 = dot(p,p);
        orb = min(orb, r2);
        k = max(ctx->mins.w/r2,1.f);
        p *= k;
        scale *= k;
    }
 //float rxy = sdBox(p, (float3)(.1f*scale));
    return (float3)(0.f, sqrt(orb), orb);
}




inline float rayIntersect(const float3 ro, const float3 rd, const Context* ctx, const float precision, const float mind, const float maxd) {
    float h, t = mind;
    for(int i=0; i<NB_RAY_ITER; i++ ) {
        h = DE_RAY(ro+rd*t, ctx);
        if (h<precision*t || t>maxd) 
            return t;
        t += h;
    }
    return -1.f;
}

float3 trace(const float3 ro, const float3 rd, const Context* ctx ) {
    const float d = rayIntersect(ro, rd, ctx, PRECISION_FACTOR, MIN_DIST_RAYMARCHING, MAX_DIST_RAYMARCHING);
    if (d>0.f) {
        return (float3)(d, DE_COLOR(ro+rd*d, ctx).yz);
    }
    return (float3)(-1.f, 1.f, 0.f);
}

#ifdef WITH_AO

float calcAO(const float3 ro, const float3 n, const Context* ctx) {
    const float seed = hash1(ro.x*(ro.y*32.56)+ro.z*147.2 + ro.y);
    const float3 rd = randomHemisphereDirection(n, seed);    
    const float d = rayIntersect(ro, rd, ctx, PRECISION_FACTOR_AO, MIN_DIST_AO, MAX_DIST_AO);
    if (d>0.f) {
        return 1.f-clamp((MAX_DIST_AO-d)/(MAX_DIST_AO),0.f,1.f);
    }
    return 1.f;
}

#endif

#ifdef WITH_SHADOWS

float shadow(const float3 ro, const float3 rd,const Context* ctx ) {
    const float 
        seed = hash1(ro.x*(ro.y*32.56)+ro.z*147.2 + ro.y),
        d = rayIntersect(ro, rd, ctx, PRECISION_FACTOR_SHADOW, MIN_DIST_SHADOW, MAX_DIST_SHADOW);
    if (d>0.f) {
        return smoothstep(0.f, MAX_DIST_SHADOW, d);
    }
    return 1.f;
}

#endif



#ifdef WITH_AO

float3 forwardSF(float i, float n) {
    float r;
    const float PI  = 3.141592653589793238f;
    const float PHI = 1.618033988749894848f;
    float phi = 2.0*PI*fract(i/PHI,&r);
    float zi = 1.f - (2.f*i+1.f)/n;
    float sinTheta = sqrt( 1.f - zi*zi);
    return (float3)( cos(phi)*sinTheta, sin(phi)*sinTheta, zi);
}

float calcAO2( float3 pos, float3 nor, const Context* ctx ) {
    float ao = 0.f;
    float3 w;
    for(int i=0; i<16; i++ ) {
        w = forwardSF( (float)i, 16.f );
	w *= sign( dot(w,nor) );
        float h = ((float)i)/15.f;
        ao += clamp( DE_RAY(pos + nor*0.01f + w*h*0.15f, ctx)*2.f, 0.f, 1.f );
    }
    ao /= 16.f;
    return clamp( ao*16.f, 0.f, 1.f );
}


float calcAO3(float3 pos, float3 nor,const Context* ctx) {
    float totao = 0.0;
    for( int aoi=0; aoi<16; aoi++ ) {
        float3 aopos = -1.0+2.0*hash3(pos.x*pos.y+pos.z+(float)(aoi)*213.47);
        aopos *= sign( dot(aopos,nor) );
        aopos = pos + nor*0.01f + aopos*0.04f;
        float dd = clamp( DE_RAY( aopos, ctx )*8.0, 0.0, 1. );
        totao += dd;
    }
    totao /= 16.0;
	
    return clamp( totao*totao*100.0, 0.0, 1.0 );
}

float calcAO4( const float3 pos, const float3 nor, const Context* ctx ) {
    float3 aopos;
    float hr, dd, 
          occ = 0.f,
          sca = 1.f;
    for( int i=0; i<5; i++ ) {
        hr = MIN_DIST_AO + MAX_DIST_AO*(float)(i)/4.f;
        aopos =  nor * hr + pos;
        dd = DE_RAY(aopos, ctx);
        occ += -(dd-hr)*sca;
        sca *= 0.95f;
    }
    return clamp( 1.f - 3.f*occ, 0.f, 1.f );    
}

#endif


float3 getRay(float3 ro, float3 look, float2 uv){
    float3 f = normalize(look - ro);
    float3 r = normalize((float3)(f.z,0.f,-f.x));
    float3 u = cross (f,r);
    return normalize(f + uv.x * r + uv.y * u);
}

float3 calcNormal( float3 pos, float t,const Context* ctx ){
    const float precis = PRECISION_FACTOR * t * 0.57f;
    const float3 e = (float3)(precis, -precis, 0.f);

    return normalize(e.xyy*DE_RAY(pos + e.xyy, ctx) + 
		e.yyx*DE_RAY(pos + e.yyx, ctx) + 
		e.yxy*DE_RAY(pos + e.yxy, ctx) + 
                     e.xxx*DE_RAY(pos + e.xxx, ctx) );
}

/*
float3 calcNormal(const float3 p, const float t, const Context* ctx){
    const float eps = PRECISION_FACTOR * t * 0.57f;
    float3 e = (float3)(eps, 0.f, 0.f);
    return normalize((float3)(
		DE_RAY(p+e.xyy,ctx)-DE_RAY(p-e.xyy,ctx),
		DE_RAY(p+e.yxy,ctx)-DE_RAY(p-e.yxy,ctx),
		DE_RAY(p+e.yyx,ctx)-DE_RAY(p-e.yyx,ctx)));
}
*/




// Cameras 
#ifndef WITH_DEPTH_OF_FIELD

float3 RD(const float3 ro, const float3 ww,  const float3 uu, const float x, const float y, const int2 res, const float fov) {
    const float3 
//        ww = normalize(ta - ro),
//        uu = normalize(cross(ww, (float3)(0.f,0.f,1.f))), // up
        vv = normalize(cross(uu,ww));

    const float2 resF = convert_float2(res);
    const float px = (2.f * (x/resF.x) - 1.f) * resF.x/resF.y, 
                py = (2.f * (y/resF.y) - 1.f);  

    float3 er = normalize( (float3)( px, py, fov) );
    return normalize( er.x*uu + er.y*vv + er.z*ww );
    //return normalize( px*uu + py*vv + fov*ww );
}

#else


float3 RD_DOF(float3* ro, const float3 ww, const float3 uu, const float x, const float y, const int2 res, const float fov, const float focusDist, const float2 rv2) {
    const float3
   //     ww = normalize(ta - *ro),
   //     uu = normalize(cross(ww, (float3)(0.f,0.f,1.f))), // up
        vv = normalize(cross(uu,ww));
   
 //float2 rv2 = hash2(24.4316544311f*seed);  
    float2 resF = convert_float2(res);
    float2 q = ((float2)(x,y))/resF;

    float2 p = 2.f*q - 1.f;
    p.x *= resF.x/resF.y;

    float2 pt = p + rv2/resF;   
 

    float3 er = normalize( (float3)( pt.xy, fov*2.f ) );
    float3 rd = er.x*uu + er.y*vv + er.z*ww;

    float3 go = FOCUSBLUR*(float3)(2.f*rv2-1.f, 0.f);
    float3 gd = normalize( er*focusDist - go );

    *ro += go.x*uu + go.y*vv;
    rd += gd.x*uu + gd.y*vv;
    rd = normalize(rd);
    return rd;
}

#endif

float3 renderScene(float3 ro, float3 rd, const Context* ctx, float* dist) {
    float3 col = (float3)(0.f);
    float3 res = trace( ro, rd, ctx);
    float t = res.x;
    *dist = t;
    if (t>0.f) {
        float3 pos = ro + t*rd;
        float3 nor = calcNormal( pos, t, ctx);
        float3 ref = reflect( rd, nor);
 // Color
	col = .9f; //0.5f + 0.5f*cos( 6.2831f*res.y + (float3)(0.f,1.f,2.f) ); 

 // lighting        
        float3 lig = LIGHT_VEC; 
        float3 hal = normalize( lig-rd);

#ifdef WITH_AO
        float occ = calcAO( pos, nor, ctx);
#else
        float occ = 1.f;
#endif

#ifdef WITH_SHADOWS
       float sh = .5f+.5f*shadow( pos, lig, ctx );
#else
        float sh = 1.f;
#endif

#ifdef ONLY_AO
	col = (float3)occ*(.5f+.5f*sh);
#else

        float amb = .3;//clamp(0.5+0.5*nor.y, 0.f, 1.f );

        float dif = clamp( dot( nor, lig ), 0.f, 1.f );
        float bac = clamp( dot( nor, normalize((float3)(-lig.x,0.f,-lig.z))), 0.f, 1.f )*clamp( 1.f-pos.y,0.f,1.f);
        float dom = smoothstep( -0.1f, 0.1f, ref.y );
        float fre = clamp(1.0+dot(nor,rd),0.f,1.f);
        fre *= fre;
        float spe = pow(clamp( dot( ref, lig ), 0.f, 1.f ),16.f);

       // dom *= softshadow( pos, ref, 0.02f, 2.5f );

	float3 lin = (float3)(.5f) + 
            + 1.3f*sh*dif*(float3)(1.f,0.8f,0.55f)
            + 2.f*spe*(float3)(1.f,0.9f,0.7f)*dif
            + .5f*occ*( .4f*amb*(float3)(0.4f,0.6f,1.f) +
                    .5f*sh*(float3)(0.4f,0.6f,1.f) +
                   // .5f*bac*(float3)(0.25f,0.25f,0.25f) +
                    .25f*fre*(float3)(1.f,1.f,1.f));

	col = col*lin;

   // 	col = mix( col, (float3)(0.8f,0.9f,1.f), 1.f-exp( -0.0002f*t*t*t ) );
        // Light attenuation, based on the distances above.

#endif
	    // Shading.
       float lDist = t;
       float atten = 1.f/(1.f + lDist*.2f + lDist*0.1f); // + distlpsp*distlpsp*0.02
       col *= atten*col*occ;
       col = mix(col, BACK_COLOR, smoothstep(0.f, .95f, t/MAX_DIST_RAYMARCHING)); // exp(-.002*t*t), etc.

    } else {
	col = BACK_COLOR;// (float3)(.08f, .16f, .34f);  	
    }
	
    return sqrt(col);
}

__kernel void render(float3 ro, float3 ww, float3 uu, const float4 sliderMins, const float4 sliderMaxs, 
    write_only image2d_t outputImage, const float4 deltaPix, read_only image2d_t demPalette,
    global write_only float* zBuffer, const float4 camera) {
    
    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int2 outSize = get_image_dim(outputImage);
    
    if (x>=outSize.x || y>=outSize.y) return;
	Context ctx = {sliderMins, sliderMaxs};
    // create ray with depth of field
    const float fov = camera.x; // 3.f;
       	
  //  ro.y -= 3.f;	
  //  ta.y -= 3.f;

    float2 res = convert_float2(outSize);
    const float2 q = ((float2)(x,y)+deltaPix.xy)/res;

#ifdef WITH_DEPTH_OF_FIELD
    float FOCUSDISTANCE = camera.y; //.05f;//length(ro-ta)*.75f;
    const float3 rd = RD_DOF(&ro, ww, uu, (float)(x)+deltaPix.x, (float)(y)+deltaPix.y, outSize, fov, FOCUSDISTANCE, deltaPix.zw); 
#else
    const float3 rd = RD(ro, ww, uu, (float)(x)+deltaPix.x, (float)(y)+deltaPix.y, outSize, fov);
#endif	

    const float3 cback = (float3)(.1*(1.-length(q-.5)));
    float dist = 0.f;
    float3 col = renderScene(ro, rd, &ctx, &dist);

#ifdef WITH_VIGNETING
    col *= pow(16.*q.x*q.y*(1.-q.x)*(1.-q.y), .3); // vigneting
#endif

    float4 out = (float4)(clamp(col,(float3)(0.f), (float3)(1.f)),1.f);

#ifdef TEXTURE_OUT_FLOAT
    write_imagef(outputImage, (int2)(x, y), out); 
#else
    uint4 rgba = (uint4)((int)(out.z*256.f),
                         (int)(out.y*256.f),
                         (int)(out.x*256.f), 256);
    write_imageui(outputImage, (int2)(x, y), rgba); 
#endif
    zBuffer[x+y*outSize.x] = dist;
}


//__kernel void dummy(const float3 roo) {
//    
//}


