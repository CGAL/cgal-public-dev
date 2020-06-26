#include <iostream>
#include <limits>
#include <assert.h>
#include <cmath>
#include <mm_malloc.h>

#include <embree3/rtcore.h> 

RTCDevice g_device;

struct vec3f { 
public:
    float x,y,z; 

    vec3f(){}
    vec3f(float _x, float _y, float _z): x(_x), y(_y), z(_z){}

    vec3f operator -(const vec3f& vec){
        vec3f result;
        result.x = x-vec.x;
        result.y = y-vec.y;
        result.z = z-vec.z;
        return result;
    }

    vec3f operator +(const vec3f& vec){
        vec3f result;
        result.x = x+vec.x;
        result.y = y+vec.y;
        result.z = z+vec.z;
        return result;
    }

    vec3f operator *(const float f){
        vec3f result;
        result.x = f*x;
        result.y = f*y;
        result.z = f*z;
        return result;
    }
}; 

inline float dot (const vec3f& a, const vec3f& b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

inline float sqr ( const float x ) { return x*x; }
double sqrt ( const double x ) { return ::sqrt (x); }
double rcp  ( const double x ) { return 1.0/x; }
template<typename T>  T clamp(const T& x, const T& lower = T(zero), const T& upper = T(one)) { return max(min(x,upper),lower); }




struct Ray
{
    Ray() {}
    Ray(const vec3f& org, 
        vec3f& dir, 
        float tnear = 0, 
        float tfar = std::numeric_limits<float>::infinity(), 
        float time = 0, 
        int mask = -1,
        unsigned int geomID = RTC_INVALID_GEOMETRY_ID, 
        unsigned int primID = RTC_INVALID_GEOMETRY_ID)
      : org(org), dir(dir), tnear(tnear), tfar(tfar), mask(mask), primID(primID), geomID(geomID)
    {
      instID[0] = RTC_INVALID_GEOMETRY_ID;
    }

public:
    vec3f org;
    vec3f dir;
    float tfar;
    float tnear;
    unsigned int mask;
    unsigned int id;
    unsigned int flags;

    vec3f Ng;
    float u;
    float v;
    unsigned int primID;           //!< primitive ID
    unsigned int geomID;           //!< geometry ID
    unsigned int instID[RTC_MAX_INSTANCE_LEVEL_COUNT]; 
};


struct Sphere{
    vec3f p;
    float r;
    RTCGeometry geometry;
    unsigned int geomID;
};

struct IntersectContext
{
  RTCIntersectContext context;
  void* userRayExt;               //!< can be used to pass extended ray data to callbacks
};

void sphereBoundsFunc(const struct RTCBoundsFunctionArguments* args)
{
  const Sphere* spheres = (const Sphere*) args->geometryUserPtr;
  RTCBounds* bounds_o = args->bounds_o;
  const Sphere& sphere = spheres[args->primID];
  bounds_o->lower_x = sphere.p.x-sphere.r;
  bounds_o->lower_y = sphere.p.y-sphere.r;
  bounds_o->lower_z = sphere.p.z-sphere.r;
  bounds_o->upper_x = sphere.p.x+sphere.r;
  bounds_o->upper_y = sphere.p.y+sphere.r;
  bounds_o->upper_z = sphere.p.z+sphere.r;
}

void sphereIntersectFunc(const RTCIntersectFunctionNArguments* args)
{
  int* valid = args->valid;
  void* ptr  = args->geometryUserPtr;
  Ray *ray = (Ray*)args->rayhit;
  RTCHit* hit = (RTCHit *)&ray->Ng.x;
  unsigned int primID = args->primID;
  
  assert(args->N == 1);
  const Sphere* spheres = (const Sphere*)ptr;
  const Sphere& sphere = spheres[primID];
  
  if (!valid[0]) return;

  const vec3f v = ray->org-sphere.p;
  const float A = dot(ray->dir,ray->dir);
  const float B = 2.0f*dot(v,ray->dir);
  const float C = dot(v,v) - sqr(sphere.r);
  const float D = B*B - 4.0f*A*C;
  if (D < 0.0f) return;
  const float Q = sqrt(D);
  const float rcpA = rcp(A);
  const float t0 = 0.5f*rcpA*(-B-Q);
  const float t1 = 0.5f*rcpA*(-B+Q);

  RTCHit potentialHit;
  potentialHit.u = 0.0f;
  potentialHit.v = 0.0f;
//   copyInstanceIdStack(args->context, potentialHit.instID);
  potentialHit.geomID = sphere.geomID;
  potentialHit.primID = primID;
  if ((ray->tnear < t0) & (t0 < ray->tfar))
  {
    int imask;
    bool mask = 1;
    {
      imask = mask ? -1 : 0;
    }
    
    const vec3f Ng = ray->org+t0*ray->dir-sphere.p;
    potentialHit.Ng_x = Ng.x;
    potentialHit.Ng_y = Ng.y;
    potentialHit.Ng_z = Ng.z;

    RTCFilterFunctionNArguments fargs;
    fargs.valid = (int*)&imask;
    fargs.geometryUserPtr = ptr;
    fargs.context = args->context;
    fargs.ray = (RTCRayN *)args->rayhit;
    fargs.hit = (RTCHitN*)&potentialHit;
    fargs.N = 1;

    const float old_t = ray->tfar;
    ray->tfar = t0;
    rtcFilterIntersection(args,&fargs);

    if (imask == -1)
      *hit = potentialHit;
    else
      ray->tfar = old_t;
  }

  if ((ray->tnear < t1) & (t1 < ray->tfar))
  {
    int imask;
    bool mask = 1;
    {
      imask = mask ? -1 : 0;
    }
    
    const vec3f Ng = ray->org+t1*ray->dir-sphere.p;
    potentialHit.Ng_x = Ng.x;
    potentialHit.Ng_y = Ng.y;
    potentialHit.Ng_z = Ng.z;

    RTCFilterFunctionNArguments fargs;
    fargs.valid = (int*)&imask;
    fargs.geometryUserPtr = ptr;
    fargs.context = args->context;
    fargs.ray = (RTCRayN *)args->rayhit;
    fargs.hit = (RTCHitN*)&potentialHit;
    fargs.N = 1;

    const float old_t = ray->tfar;
    ray->tfar = t1;
    rtcFilterIntersection(args,&fargs);

    if (imask == -1)
      *hit = potentialHit;
    else
      ray->tfar = old_t;
  }
}

void sphereFilterFunction(const RTCFilterFunctionNArguments* args)
{
  int* valid = args->valid;
  const IntersectContext* context = (const IntersectContext*) args->context;
  struct Ray* ray    = (struct Ray*)args->ray;
  //struct RTCHit* hit = (struct RTCHit*)args->hit;
  // const unsigned int N = args->N;
  assert(args->N == 1);


  /* avoid crashing when debug visualizations are used */
  if (context == nullptr)
    return;

  /* ignore inactive rays */
  if (valid[0] != -1) return;
  
  /* carve out parts of the sphere */
  const vec3f h = ray->org+ray->dir*ray->tfar;
  float v = abs(sin(10.0f*h.x)*cos(10.0f*h.y)*sin(10.0f*h.z));
  float T = clamp((v-0.1f)*3.0f,0.0f,1.0f);

  /* reject some hits */
  if (T < 0.5f) valid[0] = 0;
}

void* alignedMalloc(size_t size, size_t align) 
{
    if (size == 0)
        return nullptr;

    assert((align & (align-1)) == 0);
    void* ptr = _mm_malloc(size,align);

    if (size != 0 && ptr == nullptr)
        throw std::bad_alloc();

    return ptr;
}

Sphere* createAnalyticalSphere (RTCScene scene, const vec3f& p, float r)
{
    RTCGeometry geom = rtcNewGeometry(g_device, RTC_GEOMETRY_TYPE_USER);
    Sphere* sphere = (Sphere*) alignedMalloc(sizeof(Sphere),16);
    sphere->p = p;
    sphere->r = r;
    sphere->geometry = geom;
    sphere->geomID = rtcAttachGeometry(scene,geom);
    rtcSetGeometryUserPrimitiveCount(geom,1);
    rtcSetGeometryUserData(geom,sphere);
    rtcSetGeometryBoundsFunction(geom,sphereBoundsFunc,nullptr);
    rtcSetGeometryIntersectFunction(geom,sphereIntersectFunc);
    // rtcSetGeometryOccludedFunction (geom,sphereOccludedFunc);
    rtcCommitGeometry(geom);
    rtcReleaseGeometry(geom);
    return sphere;
}

Sphere* createAnalyticalSpheres (RTCScene scene, unsigned int N)
{
    RTCGeometry geom = rtcNewGeometry(g_device, RTC_GEOMETRY_TYPE_USER);
    Sphere* spheres = (Sphere*) alignedMalloc(N*sizeof(Sphere),16);
    unsigned int geomID = rtcAttachGeometry(scene,geom);
    for (unsigned int i=0; i<N; i++) {
        spheres[i].geometry = geom;
        spheres[i].geomID = geomID;
    }
    rtcSetGeometryUserPrimitiveCount(geom,N);
    rtcSetGeometryUserData(geom,spheres);
    rtcSetGeometryBoundsFunction(geom,sphereBoundsFunc,nullptr);
  
    rtcSetGeometryIntersectFunction(geom,sphereIntersectFunc);
    rtcSetGeometryIntersectFilterFunction(geom,sphereFilterFunction);
    
    rtcCommitGeometry(geom);
    rtcReleaseGeometry(geom);
    return spheres;
}

int main(int argc, char const *argv[])
{
    RTCDevice g_device = rtcNewDevice(NULL);

    RTCScene g_scene = rtcNewScene(g_device);
    Sphere* g_spheres = createAnalyticalSpheres(g_scene,5);
    g_spheres[0].p = vec3f( 0, 0,+1); g_spheres[0].r = 0.5f;
    g_spheres[1].p = vec3f(+1, 0, 0); g_spheres[1].r = 0.5f;
    g_spheres[2].p = vec3f( 0, 0,-1); g_spheres[2].r = 0.5f;
    g_spheres[3].p = vec3f(-1, 0, 0); g_spheres[3].r = 0.5f;
    g_spheres[4].p = vec3f( 0,+1, 0); g_spheres[4].r = 0.5f;
    rtcCommitScene(g_scene);

    struct RTCIntersectContext context;
    rtcInitIntersectContext(&context);

      /* initialize ray */
    // Ray ray(Vec3fa(camera.xfm.p), 
    //                     Vec3fa(normalize(x*camera.xfm.l.vx + y*camera.xfm.l.vy + camera.xfm.l.vz)), 
    //                     0.0f, inf, 0.0f, -1,
    //                     RTC_INVALID_GEOMETRY_ID, RTC_INVALID_GEOMETRY_ID);

    /* intersect ray with scene */
    // rtcIntersect1(g_scene,&context,RTCRayHit_(ray));

    rtcReleaseScene(g_scene);
    rtcReleaseDevice(g_device);
    return 0;
}
