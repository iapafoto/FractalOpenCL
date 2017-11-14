#define TEXTURE_IN_FLOAT
//#define TEXTURE_OUT_FLOAT

#define DE_RAY pseudoKleinianSlider
#define DE_COLOR pseudoKleinianSliderColor

// #define DE kleinian
// #define DE deSpace
// #define DE apollonian

//#define WITH_SHADOWS
//#define WITH_SUN
//#define WITH_AO
//#define ONLY_AO
//#define WITH_VIGNETING
//#define WITH_DEPTH_OF_FIELD
//#define WITH_ISOLINE

#ifdef ONLY_AO
    #define BACK_COLOR (float3)(.8f, .8f, .8f) 
#else 
    #define BACK_COLOR (float3)(.08f, .16f, .34f) 
//    #define BACK_COLOR (float3)(.8f, .65f, .54f) 
#endif

#define FOCUSBLUR .005f

#define PRECISION_FACTOR 3e-4
#define MIN_DIST_RAYMARCHING .01f
#define MAX_DIST_RAYMARCHING 15.f

#define MIN_DIST_SHADOW 10.f*PRECISION_FACTOR
#define MAX_DIST_SHADOW 3.f
#define PRECISION_FACTOR_SHADOW PRECISION_FACTOR

#define MIN_DIST_AO 10.f*PRECISION_FACTOR
#define MAX_DIST_AO .2f
#define PRECISION_FACTOR_AO PRECISION_FACTOR

#define LIGHT_VEC normalize((float3)(.2,.7, 1.6) )

#define DIFF 0
#define REFR 1
#define SPEC 2
#define CHECK 3

#define PATH_TRACING_DEPTH 5
#define EYEPATHLENGTH 4


typedef struct Material {
	float3 emission; 
	float3 colour; 
	int type;
} Material;


const sampler_t sampler_linear = CLK_NORMALIZED_COORDS_FALSE |  CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

__constant float3 COLOR_WATER = (float3)(0.3f, 0.13f, 0.08f);
__constant float3 COLOR_BACK = (float3)(.42f,.46f,.48f);    

__constant int iter = 100;
__constant float eps = 0.001f, far = 3.f;


typedef struct matrix3 {
	float3 m[3];
}
matrix3;


typedef struct Context {
    const float4 mins;
    const float4 maxs;
    const matrix3 rot;
} Context;


inline float3 reflect(const float3 i, const float3 n){
  return i - 2.f * n * dot(n,i);
}

inline float3 refract(float3 v, float3 normal, float n){
	// assumes that n is already normalized
	float ct1 = dot(normal, -1 * v);
	float ct2 = sqrt(1 - n*n*(1 - ct1*ct1));
	return n * v + (n * ct1 - ct2) * normal;
}

// -------------------------------------------------------------------

inline float hash1(float* seed) {
    float fractptr;
    return fract(sin((*seed+=.1f))*43758.5453123f, &fractptr);
}
inline float2 hash2(float* seed) {
    float2 fractptr;
    return fract(sin((float2)((*seed+=.1)*43758.5453123f,(*seed+=.1)*22578.1459123f)), &fractptr);
}

inline float3 hash3(const float seed) {
    float3 fractptr;
    return fract(sin((float3)(seed,seed+.1f,seed+.2f))*(float3)(43758.5453123f,22578.1459123f,19642.3490423f), &fractptr);
}

float rand(float2 co){
	float fractptr;
	return fract(sin(dot(co.xy ,(float2)(12.9898f,78.233f))) * 43758.5453f, &fractptr);
}

matrix3  rotationMatrix3(const float3 v, const float angle) {
	const float c = cos(angle), s = sin(angle);
	matrix3 m;
	m.m[0] = (float3)(c + (1.f - c) * v.x * v.x, (1.f - c) * v.x * v.y - s * v.z, (1.f - c) * v.x * v.z + s * v.y);
	m.m[1] = (float3)((1.f - c) * v.x * v.y + s * v.z, c + (1.f - c) * v.y * v.y, (1.f - c) * v.y * v.z - s * v.x);
	m.m[2] = (float3)((1.f - c) * v.x * v.z - s * v.y, (1.f - c) * v.y * v.z + s * v.x, c + (1.f - c) * v.z * v.z);
	return m;
}




float3 cosWeightedRandomHemisphereDirection( const float3 n, float* seed ) {
  	float2 r = hash2(seed);
    
	float3  uu = normalize( cross( n, (float3)(0.f,1.f,1.f) ) );
	float3  vv = cross( uu, n );
	
	float ra = sqrt(r.y);
	float rx = ra*cos(6.2831f*r.x); 
	float ry = ra*sin(6.2831f*r.x);
	float rz = sqrt( 1.f-r.y );
	float3  rr = (float3)( rx*uu + ry*vv + rz*n );
    
    return normalize( rr );
}

float3 randomSphereDirection(float* seed) {
    float2 r = hash2(seed)*6.2831f;
	float3 dr = (float3)(sin(r.x)*(float2)(sin(r.y),cos(r.y)),cos(r.x));
	return dr;
}

float3 randomHemisphereDirection( const float3 n, float* seed ) {
	float3 dr = randomSphereDirection(seed);
	return dot(dr,n) * dr;
}

//-----------------------------------------------------
// light
//-----------------------------------------------------
/*
float4 lightSphere;

void initLightSphere( float time ) {
	lightSphere = (float4)( 3.f+2.f*sin(time),2.8f+2.f*sin(time*0.9f),3.f+4.f*cos(time*.7f), .5f );
}
*/

float3 sampleLight(float3 ro, float* seed, float4 lightSphere) {
    float3 n = randomSphereDirection( seed ) * lightSphere.w;
    return lightSphere.xyz + n;
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

float cylUnion(const float3 p){
    float xy = dot(p.xy,p.xy);
    float xz = dot(p.xz,p.xz);
    float yz = dot(p.yz,p.yz);
    return sqrt(min(xy,min(xz,yz))) - 1.f;
}

float cylIntersection(const float3 p){
    float xy = dot(p.xy,p.xy);
    float xz = dot(p.xz,p.xz);
    float yz = dot(p.yz,p.yz);
    return sqrt(max(xy,max(xz,yz))) - 1.f;
}



#define R(p, a) p=cos(a)*p+sin(a)*(float2)(p.y, -p.x)


float pseudoKleinianSlider(float3 p, const Context* ctx) {//knighty's pseudo kleinian
    float k,r2, scale=1.f, orb = 1.f;
    for(int i=0;i<7;i++) {
       // if(i==rotater)p.xy=p.xy*rmx;
        p = 2.f*clamp(p, ctx->mins.xyz, ctx->maxs.xyz)-p;
        r2 = dot(p,p);
        k = max(ctx->mins.w/dot(p,p),1.f);
        p *= k;
        scale *= k;
    }
    float rxy=length(p.xy);
    return .5f*max(rxy-ctx->maxs.w, /*fabs*/(rxy*p.z) / length(p))/scale;
}

float pseudoKleinianSliderSymetrique(float3 p, const Context* ctx) {//knighty's pseudo kleinian
    float k, scale = 1.f;
    for(int i=0;i<8;i++) {
        p = 2.f*clamp(p, ctx->mins.xyz, ctx->maxs.xyz)-p;
        k = max(ctx->mins.w / dot(p,p),1.f);
        p *= k;
        scale *= k;
    }
    return .5f*max(1.f-ctx->maxs.w/length(p), /*fabs*/(max(p.x, max(p.y, p.z))))/scale;
}

float3 pseudoKleinianSliderColor(float3 p, const Context* ctx) {//knighty's pseudo kleinian
    float k,r2, scale = 1.f, orb = 1.f;
    for(int i=0;i<8;i++) {
        p = 2.f*clamp(p, ctx->mins.xyz, ctx->maxs.xyz)-p;
        r2 = dot(p,p);
        orb = min(orb, r2);
        k = max(ctx->mins.w/r2,1.f);
        p *= k;
        scale *= k;
    }
    return (float3)(0.f, .25f+sqrt(orb), orb);
}



inline float rayIntersect(const float3 ro, const float3 rd, const Context* ctx, const float precision, const float mind, const float maxd) {
    float h, t = mind;
    for(int i=0; i<192; i++ ) {
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

float iSphere( float3 ro, float3 rd, float4 sph ) {
    float3 oc = ro - sph.xyz;
    float b = dot(oc, rd);
    float c = dot(oc, oc) - sph.w * sph.w;
    float h = b * b - c;
    if (h < 0.f) return -1.f;

	float s = sqrt(h);
	float t1 = -b - s;
	float t2 = -b + s;
	
	return t1 < 0.f ? t2 : t1;
}
float3 nSphere( float3 pos, float4 sph ) {
    return (pos-sph.xyz)/sph.w;
}

float2 intersect( float3 ro, float3 rd, const Context* ctx,  float3* normal, float4 lightSphere) {
	float2 res = (float2)( 1e20, -1.f );

	float3 intersectFractal = trace(ro,rd, ctx);
	res.x = intersectFractal.x;
	res.y = 1.f;
	float t;
	
   // t = iSphere( ro, rd, (float4)( 1.5,1.0, 2.7, 1.0) ); if( t>eps && t<res.x ) { res = (float2)( t, 1. ); normal = nSphere( ro+t*rd, (float4)( 1.5,1.0, 2.7,1.0) ); }
   // t = iSphere( ro, rd, (float4)( 4.0,1.0, 4.0, 1.0) ); if( t>eps && t<res.x ) { res = (float2)( t, 6. ); normal = nSphere( ro+t*rd, (float4)( 4.0f,1.0f, 4.0f,1.0f) ); }
    t = iSphere( ro, rd, lightSphere ); 
	if( t>eps && t<res.x ) { 
		res = (float2)( t, 0.0f );  
		*normal = nSphere( ro+t*rd, lightSphere ); 
	}
	else {
 		*normal = calcNormal(ro+rd*intersectFractal.x, intersectFractal.x, ctx);
	}
					  
    return res;					  
}

bool intersectShadow( float3 ro, float3 rd, float dist, const Context* ctx) {
    float t = rayIntersect(ro, rd, ctx, PRECISION_FACTOR, MIN_DIST_RAYMARCHING, dist/*MAX_DIST_RAYMARCHING*/);
    return (t>0.f && t<dist /*MAX_DIST_RAYMARCHING*/);
}



//-----------------------------------------------------
// materials
//-----------------------------------------------------

#define LIGHTCOLOR (float3)(16.86f, 10.76f, 8.2f)*1.3f
#define WHITECOLOR (float3)(.7295f, .7355f, .729f)*0.7f
#define GREENCOLOR (float3)(.117f, .4125f, .115f)*0.7f
#define REDCOLOR (float3)(.611f, .0555f, .062f)*0.7f

float3 matColor( float mat ) {
    float3 nor = (float3)(0.f, 0.95f, 0.f);
    if ( mat<3.5f ) nor = REDCOLOR;
    if ( mat<2.5f ) nor = GREENCOLOR;
    if ( mat<1.5f ) nor = WHITECOLOR;
    if ( mat<0.5f ) nor = LIGHTCOLOR;				  
    return nor;					  
}

bool matIsSpecular( float mat ) {
    return mat > 4.5f;
}

bool matIsLight( float mat ) {
    return mat < 0.5f;
}

//-----------------------------------------------------
// brdf
//-----------------------------------------------------

float3 getBRDFRay( float3 n, float3 rd,float m, bool* specularBounce, float* seed ) {
    *specularBounce = false;
    
    float3 r = cosWeightedRandomHemisphereDirection( n, seed );
    if(  !matIsSpecular( m ) ) {
        return r;
    } else {
        *specularBounce = true;
        
        float n1, n2, ndotr = dot(rd,n);
        
        if( ndotr > 0. ) {
            n1 = 1.f/1.5f; n2 = 1.f;
            n = -n;
        } else {
            n2 = 1.f/1.5f; n1 = 1.f;
        }
                
        float r0 = (n1-n2)/(n1+n2); r0 *= r0;
		float fresnel = r0 + (1.f-r0) * pow(1.f-fabs(ndotr),5.f);
        
        float3 ref;
        
        if( hash1(seed) < fresnel ) {
            ref = reflect( rd, n );
        } else {
            ref = refract( rd, n, n2/n1 );
        }
        
        return ref; // normalize( ref + 0.1f * r );
	}
}

//-----------------------------------------------------
// eyepath
//-----------------------------------------------------

float3 traceEyePath( float3 ro, float3 rd, Context* ctx, bool directLightSampling, float* seed ) {
    float3 tcol = (float3)(0.f);
    float3 fcol  = (float3)(1.f);
    
    bool specularBounce = true;
    
    float4 lightSphere;

    for( int j=0; j<EYEPATHLENGTH; ++j ) {
        float3 normal;
        
        float2 res = intersect( ro, rd, ctx, &normal, lightSphere );
        if( res.y < -0.5f ) {
            return tcol;
        }
        
        if( matIsLight( res.y ) ) {
            if( directLightSampling ) {
            	if( specularBounce ) tcol += fcol*LIGHTCOLOR;
            } else {
                tcol += fcol*LIGHTCOLOR;
            }
         //   basecol = (float3)(0.f);	// the light has no diffuse component, therefore we can return col
            return tcol;
        }
        
        ro = ro + res.x * rd;
        rd = getBRDFRay( normal, rd, res.y, specularBounce, seed );
        
        fcol *= matColor( res.y );

        float3 ld = sampleLight( ro, &seed, lightSphere) - ro;
        
        if( directLightSampling ) {
			float3 nld = normalize(ld);
            if( !specularBounce && j < EYEPATHLENGTH-1 && !intersectShadow( ro, nld, length(ld), ctx) ) {

                float cos_a_max = sqrt(1.f - clamp(lightSphere.w * lightSphere.w / dot(lightSphere.xyz-ro, lightSphere.xyz-ro), 0.f, 1.f));
                float weight = 2.f * (1.f - cos_a_max);

                tcol += (fcol * LIGHTCOLOR) * (weight * clamp(dot( nld, normal ), 0., 1.));
            }
        }
    }    
    return tcol;
}

float3 getRay(float3 ro, float3 look, float2 uv){
    float3 f = normalize(look - ro);
    float3 r = normalize((float3)(f.z,0.f,-f.x));
    float3 u = cross (f,r);
    return normalize(f + uv.x * r + uv.y * u);
}



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



__kernel void render(float3 ro, float3 ww, float3 uu, const float4 sliderMins, const float4 sliderMaxs, 
    write_only image2d_t outputImage, const float4 deltaPix, read_only image2d_t demPalette,
    global write_only float* zBuffer, const float4 camera) {
    
    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int2 outSize = get_image_dim(outputImage);
    


    if (x>=outSize.x || y>=outSize.y) return;
    
    const float3 RotVector = sliderMins.xyz;//(float3)(0.5f,-0.05f,-0.5f);
    matrix3 matRot = rotationMatrix3(normalize(RotVector), 99.f);

	Context ctx = {sliderMins, sliderMaxs, matRot};
    // create ray with depth of field
    const float fov = camera.x; // 3.f;
       	
   // ro.y -= 40.f;	
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
	float seed = deltaPix.x+deltaPix.y; //
	float3 col = traceEyePath( ro, rd,&ctx, true, &seed );
  //  float3 col = renderScene(ro, rd, &ctx, &dist);

#ifdef WITH_VIGNETING
    col *= pow(16.f*q.x*q.y*(1.f-q.x)*(1.f-q.y), .3f); // vigneting
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


