#include <stdio.h>
#include <embree3/rtcore.h>
#include <limits>
#include <iostream>
#include <fstream>
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

  std::ifstream in("in.off");
  std::string off;
  int v, f, e;
  in >> off >> v >> f >> e;

  RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
//   float* vertices = (float*) rtcSetNewGeometryBuffer(geom,
//                                                      RTC_BUFFER_TYPE_VERTEX,
//                                                      0,
//                                                      RTC_FORMAT_FLOAT3,
//                                                      3*sizeof(float),
//                                                      v);

  Vertex* vertices = (Vertex*) rtcSetNewGeometryBuffer(geom,RTC_BUFFER_TYPE_VERTEX,0,RTC_FORMAT_FLOAT3,sizeof(Vertex),v);


//   unsigned* indices = (unsigned*) rtcSetNewGeometryBuffer(geom,
//                                                           RTC_BUFFER_TYPE_INDEX,
//                                                           0,
//                                                           RTC_FORMAT_UINT3,
//                                                           3*sizeof(unsigned),
//                                                           2);
  Triangle* triangles = (Triangle*) rtcSetNewGeometryBuffer(geom,RTC_BUFFER_TYPE_INDEX,0,RTC_FORMAT_UINT3,sizeof(Triangle),f);

  if (vertices && triangles)
  {
    for(int i = 0; i < v; i++){
      in >> vertices[i].x;
      in >> vertices[i].y;
      in >> vertices[i].z;
    }
    for(int i = 0; i < f; i++){
      in >> e; // 3
      in >> triangles[i].v0;
      in >> triangles[i].v1;
      in >> triangles[i].v2;
    }
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
    rtcSetSceneFlags(scene, RTC_SCENE_FLAG_ROBUST);


    castRay(scene, 0.1f, 0.1f, 0, 0, 0, 1);
    /*the directions need to be normalised.*/

    rtcReleaseScene(scene);
    rtcReleaseDevice(device);

    return 0;
}
