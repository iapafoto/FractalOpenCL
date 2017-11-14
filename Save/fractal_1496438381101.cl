#define TEXTURE_IN_FLOAT
//#define TEXTURE_OUT_FLOAT

#define DE_RAY testDE//pseudoKleinianSlider
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
#define MIN_DIST_RAYMARCHING .1f
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

float2 testDE(float3 p, const Context* ctx) {//knighty's pseudo kleinian
	float d1 = length(p)-1.5f;
	float d2 = length(p-(float3)(1.f))-1.f;
	if (d1 <d2) 
		return (float2)(d1,0.f);
	else  return (float2)(d2,1.f);
}

float2 pseudoKleinianSlider(float3 p, const Context* ctx) {//knighty's pseudo kleinian
    float k,r2, scale=1.f, orb = 1.f;
    for(int i=0;i<7;i++) {
       
        p = 2.f*clamp(p, ctx->mins.xyz, ctx->maxs.xyz)-p;
        r2 = dot(p,p);
        k = max(ctx->mins.w/dot(p,p),1.f);
        p *= k;
        scale *= k;
    }
    float rxy=length(p.xy);
    float d1 = .5f*max(rxy-ctx->maxs.w, /*fabs*/(rxy*p.z) / length(p))/scale;
	float d2 = length(p-(float3)(1.f))-1.f;
	if (d1 <d2) 
		return (float2)(d1,0.f);
	else  return (float2)(d2,1.f);

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



inline float2 rayIntersect(const float3 ro, const float3 rd, const Context* ctx, const float precision, const float mind, const float maxd) {
    float2 h;
    float t = mind;
    for(int i=0; i<192; i++ ) {
        h = DE_RAY(ro+rd*t, ctx);
        if (h.x<precision*t || t>maxd) 
            return (float2)(t,h.y);
        t += h.x;
    }
    return (float2)(-1.f,-1.f);
}

float3 trace(const float3 ro, const float3 rd, const Context* ctx ) {
    const float2 d = rayIntersect(ro, rd, ctx, PRECISION_FACTOR, MIN_DIST_RAYMARCHING, MAX_DIST_RAYMARCHING);
    if (d.x>0.f) {
        return (float3)(d.x, d.y, DE_COLOR(ro+rd*d.x, ctx).y);
    }
    return (float3)(-1.f, 1.f, 0.f);
}



float3 getRay(float3 ro, float3 look, float2 uv){
    float3 f = normalize(look - ro);
    float3 r = normalize((float3)(f.z,0.f,-f.x));
    float3 u = cross (f,r);
    return normalize(f + uv.x * r + uv.y * u);
}

float3 calcNormal( float3 pos, float t,const Context* ctx ){
    const float precis = PRECISION_FACTOR * t * 0.57f;
    const float3 e = (float3)(precis, -precis, 0.f);

    return normalize(e.xyy*DE_RAY(pos + e.xyy, ctx).x + 
		     e.yyx*DE_RAY(pos + e.yyx, ctx).x + 
		     e.yxy*DE_RAY(pos + e.yxy, ctx).x + 
                     e.xxx*DE_RAY(pos + e.xxx, ctx).x );
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


float3 renderScene( float3 ro, float3 rd, const Context* ctx, float* dist) {

		float2 r, rng = fract(ro.xy, &r);

		float3 finalCol = (float3)(0.0f,0.0f,0.0f);
        float3 fCum = (float3)(1.0f,1.0f,1.0f);

		Material obj, obj0;
		obj0.emission = (float3)(0.0f,0.f,0.0f);
		obj0.colour = (float3)(1.f,.5f,.5f);
		obj0.type = DIFF;

		Material objLight;
		objLight.emission = (float3)(1.0f,1.0f,1.0f)*2.f;
		objLight.colour = (float3)(1.f,1.f,1.f)*1.f;
		objLight.type = DIFF;

		float iFrame = hash1(rd.x*rd.y*ro.x);
    
		for (int depth = 0; depth < PATH_TRACING_DEPTH; depth++)
        {
              //  float t = 0.0;                            // distance to intersection
               // int id = 0;                               // id of intersected object
             //   Sphere obj;
				float3 res = trace( ro, rd, ctx);
				float t = res.x;

				if (t<0.f) {
//t = MAX_DIST_RAYMARCHING; 
//obj = objLight;

                    break;
				} else {
if (res.y>.5f) obj = objLight; else
					{
						obj = obj0;
//obj.colour =  0.5f + 0.5f*cos( 6.2831f*res.y + (float3)(0.f,1.f,2.f) );
					} 
				}


                float3 x = ro + rd * t;
                float3 n = calcNormal(x, t, ctx);
                float3 nl = dot(n,rd) < 0.0f ? n : n * -1.0f;
                float3 f = obj.colour;

             /*   float p = max(max(f.x,f.y),f.z);

				rng.x = rand( rng );
				
                if ( rng.x < p)
                    f = f / p;
                else
                    break; //R.R.
*/
                fCum = f * fCum;

                if (obj.type == DIFF || obj.type == CHECK) // Ideal DIFFUSE reflection
                {
					if( obj.type == CHECK )
						{
							if( (fmod(x.x,80.0f) < 40.0f && fmod(x.z,80.0) < 40.0f) || 
								(fmod(x.x,80.0f) > 40.0f && fmod(x.z,80.0f) > 40.0) )
								fCum *= 0.3f;
								
						}
				
                    float r1 = 2.f * 3.1415926536f * rand( rng ); rng.x = sin(r1 - iFrame);
                    float r2 = rand( rng ); 
		rng.y = sin(r2 + iFrame);
                    float r2s = sqrt(r2);
                    float3 w = nl;
					float3 u = normalize(cross( (fabs(w.x) > .1f ? (float3)(0.f, 1.0, 0.f) : (float3)(1.0,1.0,1.0)) , w));
					float3 v = cross(w,u);
                    float3 d = normalize(u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1.f - r2));
                    finalCol = finalCol + fCum * obj.emission;
                    ro = x;
		rd = d;
                    continue;
                }
                else
                    if (obj.type == SPEC)            // Ideal SPECULAR reflection
                    {
                        finalCol = finalCol + fCum * obj.emission;
                        ro = x;
		   rd = rd - n * 2.f * dot(n,rd);
                        continue;
                    }

                // Ideal dielectric REFRACTION
                float3 reflRayRo = x;
	      float3 reflRayRd = rd - n * 2.f * dot(n,rd);     
	      bool into = dot(n,nl) > 0.f;                // Ray from outside going in?

                float nc = 1.0;   // IOR of air
                float nt = 1.5f; // IOR of solid
                float nnt = into ? nc / nt : nt / nc;
                float ddn = dot(rd , nl);
                float cos2t = 1.f - nnt * nnt * (1.f - ddn * ddn);

                if (cos2t < 0.f)    // Total internal reflection
                {
                    finalCol = finalCol + fCum * obj.emission;
                    ro = reflRayRo;
		rd = reflRayRd;
                    continue;
                }

                float3 tdir = normalize(rd * nnt - n * ((into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t))));

                float a = nt - nc;
                float b = nt + nc;
                float R0 = a * a / (b * b);
                float c = 1.f - (into ? -ddn : dot(tdir,n));
                float Re = R0 + (1.f - R0) * c * c * c * c * c;
                float Tr = 1.f - Re;
                float P = .25f + .5 * Re;
                float RP = Re / P;
                float TP = Tr / (1.f - P);

		rng.y = rand(rng);
				
                if( rng.y < P )
                {
                    ro = reflRayRo;
                    rd = reflRayRd;
                    fCum = fCum * RP;
                    finalCol = finalCol + fCum * obj.emission;
                }
                else
                {
                    ro = x;
		rd = tdir;
                    fCum = fCum * TP;
                    finalCol = finalCol + fCum * obj.emission;
                }
				
				// we reached something bright, don't spawn any more rays
                if (length( obj.emission ) > 100.f)
                    break;

            }

            return finalCol;
	}

__kernel void render(float3 ro, float3 ww, float3 uu, const float4 sliderMins, const float4 sliderMaxs, 
    write_only image2d_t outputImage, const float4 deltaPix, read_only image2d_t demPalette,
    global write_only float* zBuffer, const float4 camera) {
    
    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int2 outSize = get_image_dim(outputImage);
    

//ro = (float3)(-6.f);
    if (x>=outSize.x || y>=outSize.y) return;
    
   // const float3 RotVector = sliderMins.xyz;//(float3)(0.5f,-0.05f,-0.5f);
    //matrix3 matRot = rotationMatrix3(normalize(RotVector), 99.f);

	Context ctx = {sliderMins, sliderMaxs};
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


