#include <stdio.h>
#include <embree3/rtcore.h>
#include <limits>
#include <iostream>
#include <math.h>

struct Vertex   { float x,y,z; }; 
struct Triangle { unsigned int v0, v1, v2; };

void errorFunction(void* userPtr, enum RTCError error, const char* str)
{
  printf("error %d: %s\n", error, str);
}

RTCDevice initializeDevice()
{
  RTCDevice device = rtcNewDevice(NULL);

  if (!device)
    printf("error %d: cannot create device\n", rtcGetDeviceError(NULL));

  rtcSetDeviceErrorFunction(device, errorFunction, NULL);
  return device;
}

RTCScene initializeScene(RTCDevice device)
{
  RTCScene scene = rtcNewScene(device);

  RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
//   float* vertices = (float*) rtcSetNewGeometryBuffer(geom,
//                                                      RTC_BUFFER_TYPE_VERTEX,
//                                                      0,
//                                                      RTC_FORMAT_FLOAT3,
//                                                      3*sizeof(float),
//                                                      4);

  Vertex* vertices = (Vertex*) rtcSetNewGeometryBuffer(geom,RTC_BUFFER_TYPE_VERTEX,0,RTC_FORMAT_FLOAT3,sizeof(Vertex),3);


//   unsigned* indices = (unsigned*) rtcSetNewGeometryBuffer(geom,
//                                                           RTC_BUFFER_TYPE_INDEX,
//                                                           0,
//                                                           RTC_FORMAT_UINT3,
//                                                           3*sizeof(unsigned),
//                                                           2);
  Triangle* triangles = (Triangle*) rtcSetNewGeometryBuffer(geom,RTC_BUFFER_TYPE_INDEX,0,RTC_FORMAT_UINT3,sizeof(Triangle),1);

  if (vertices && triangles)
  {
    vertices[0].x = 0.f; vertices[0].y = 0.f; vertices[0].z = 0.f;
    vertices[1].x = 1.f; vertices[1].y = 0.f; vertices[1].z = 0.f;
    vertices[2].x = 0.f; vertices[2].y = 1.f; vertices[2].z = 0.f;
    // vertices[3].x = 1.f; vertices[3].y = 1.f; vertices[3].z = 0.f;
    // vertices[0] = 0.f; vertices[1] = 0.f; vertices[2] = 0.f;
    // vertices[3] = 1.f; vertices[4] = 0.f; vertices[5] = 0.f;
    // vertices[6] = 0.f; vertices[7] = 1.f; vertices[8] = 0.f;
    // vertices[9] = 1.f; vertices[10] = 1.f; vertices[11] = 0.f;

    triangles[0].v0 = 0; triangles[0].v1 = 1; triangles[0].v2 = 2;
    triangles[1].v0 = 3; triangles[1].v1 = 1; triangles[1].v2 = 2;

    // indices[0] = 0; indices[1] = 1; indices[2] = 2;
    // indices[3] = 3; indices[4] = 1; indices[5] = 2;
  }

  rtcCommitGeometry(geom);

  rtcAttachGeometry(scene, geom);
  rtcReleaseGeometry(geom);

  rtcCommitScene(scene);

  return scene;
}

void castRay(RTCScene scene, 
             float ox, float oy, float oz,
             float dx, float dy, float dz)
{
  struct RTCIntersectContext context;
  rtcInitIntersectContext(&context);

  struct RTCRayHit rayhit;
  rayhit.ray.org_x = ox;
  rayhit.ray.org_y = oy;
  rayhit.ray.org_z = oz;
  rayhit.ray.dir_x = dx;
  rayhit.ray.dir_y = dy;
  rayhit.ray.dir_z = dz;
  rayhit.ray.tnear = 0;
  rayhit.ray.tfar = std::numeric_limits<float>::infinity();
  rayhit.ray.mask = 0;
  rayhit.ray.flags = 0;
  rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
  rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

  rtcIntersect1(scene, &context, &rayhit);

  // printf("%f, %f, %f: ", ox, oy, oz);
  if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
  {
    printf("Found intersection on geometry %d, primitive %d at tfar=%f\n", 
           rayhit.hit.geomID,
           rayhit.hit.primID,
           rayhit.ray.tfar);

    float factor = rayhit.ray.tfar/ sqrt(pow(rayhit.ray.dir_x, 2)+ pow(rayhit.ray.dir_y, 2)+ pow(rayhit.ray.dir_z, 2));
    float x = rayhit.ray.org_x + factor * rayhit.ray.dir_x;
    float y = rayhit.ray.org_y + factor * rayhit.ray.dir_y;
    double z = rayhit.ray.org_z + (factor * rayhit.ray.dir_z);
    std::cout<<"Intersection Point : "<<x<<" "<<y<<" "<<z<<std::endl; 
    
  }
  else
    printf("Did not find any intersection.\n");
}

int main(int argc, char *argv[])
{
    RTCDevice device = initializeDevice();
    RTCScene scene = initializeScene(device);

    castRay(scene, 0.1f, 0.1f, -1, 0, 0, 1);
    /*the directions need to be normalised.*/

    rtcReleaseScene(scene);
    rtcReleaseDevice(device);

    return 0;
}
