
// https://github.com/Syntopia/Fragmentarium/blob/master/Fragmentarium-Source/Examples/Include/Sky-Pathtracer.frag

#define Reflectivity .5f
#define Albedo .5f

#define DirectLight true
#define Stratify false



float3 ortho(float3 v) {
	//  See : http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts
	return fabs(v.x) > fabs(v.z) ? (float3)(-v.y, v.x, 0.f)  : (float3)(0.f, -v.z, v.y);
}

float3 getSampleBiased(float3  dir, float power, float2 seed) {
	dir = normalize(dir);
	// create orthogonal vector
	float3 o1 = normalize(ortho(dir));
	float3 o2 = normalize(cross(dir, o1));
	
	// Convert to spherical coords aligned to dir;
	float2 r = rand(seed+=(float2)(-1.f,1.f));
 
	if (Stratify) { r*=0.1; r+= cx;  cx = mod(cx + float2(0.1,0.9),1.0);}
	r.x=r.x*2.*PI;
	
	// This is  cosine^n weighted.
	// See, e.g. http://people.cs.kuleuven.be/~philip.dutre/GI/TotalCompendium.pdf
	// Item 36
	r.y = fpow(r.y,1.0/(power+1.0));
	
	float oneminus = sqrt(1.0-r.y*r.y);
	float3 sdir = cos(r.x)*oneminus*o1+
	sin(r.x)*oneminus*o2+
	r.y*dir;
	
	return sdir;
}

float3 getConeSample(float3 dir, float extent, float2 seed) {
	// Create orthogonal vector (fails for z,y = 0)
	dir = normalize(dir);
	float3 o1 = normalize(ortho(dir));
	float3 o2 = normalize(cross(dir, o1));
	
	// Convert to spherical coords aligned to dir
	float2 r =  rand2n(seed+=(float2)(-1.f,1.f));
	
	if (Stratify) {r*=0.1f; r+= cx;}
	r.x=r.x*2.f*PI;
	r.y=1.f-r.y*extent;
	
	float oneminus = sqrt(1.f-r.y*r.y);
	return cos(r.x)*oneminus*o1+sin(r.x)*oneminus*o2+r.y*dir;
}

float3 renderScene(float3 from, float3 dir, Context* ctx, float* dist)
{
	bool BiasedSampling = false;


	float3 hit = (float3)(0.0);
	float3 hitNormal = (float3)(0.0);
	
	float3 color = (float3)(1.0);
	float3 direct = (float3)(0.0);

	for( int i=0; i <RayDepth; i++ ) {

		float3 res = trace( from, dir, ctx);
		float hit = res.x;

		if (hit<=0.f|| hit >= MAX_DIST_RAYMARCHING) {

			hitNormal = calcNormal(from, hit, ctx);

			// We hit something
			if (rand(from+dir) > Reflectivity ) {
				#ifdef providesColor
				color *= baseColor(hit, hitNormal);
				#else
				color *= PALETTE(sqrt(res.y)); //getColor();
				#endif
				
				//color *= (1.0-Reflectivity);
				if  (!BiasedSampling) {
					// Unbiased sampling:
					// PDF = 1/(2*PI), BRDF = Albedo/PI
					dir = getConeSample( hitNormal,1.0);
					// modulate color with: BRDF*CosAngle/PDF
					color *= 2.0*Albedo*max(0.0,dot(dir,hitNormal));
				}
				else  {
					// Biased sampling (cosine weighted):
					// PDF = CosAngle / PI, BRDF = Albedo/PI
					dir =getSampleBiased( hitNormal, 1.0 );
					
					// modulate color	 with: BRDF*CosAngle/PDF
					color *= Albedo;
				}
				
				// Direct
				if (DirectLight) {
					float3 a;
					float3 b;
					float3 sunSampleDir = getConeSample(sunDirection, 1.f-sunAngularDiameterCos);
					float sunLight = dot(hitNormal, sunSampleDir);
					if (sunLight>0.0 && !trace(hit+ hitNormal*3.0*minDist,sunSampleDir,a,b)) {
						direct += color*sun(sunSampleDir)*sunLight *1E-5;
					}
				}
				
			}
			else {
				//color *=Reflectivity;
				dir=reflect(dir,hitNormal);
				color *= max(0.0, dot(dir, hitNormal));
			}
			
			// Choose new starting point for ray
			from =  hit + hitNormal*minDist*8.0;
		} else {
			if (DebugLast && i!=RayDepth-1) {
				return (float3)(0.0);
			}
			if (!DirectLight) return color * sunsky(dir);
			return direct + color * (i>0 ? sky(dir) : sunsky(dir));
		}
	}
	return direct;
}