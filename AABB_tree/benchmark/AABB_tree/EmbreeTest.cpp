#include <iostream>
#include <fstream>
#include <limits>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <embree3/rtcore.h> 

#include <CGAL/Real_timer.h>

#include "RaysGenerate.h"

typedef CGAL::Simple_cartesian<float> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;

struct Vertex   { float x,y,z; }; 
struct Triangle { unsigned int v0, v1, v2; };

void errorFunction(void* userPtr, enum RTCError error, const char* str)
{
  std::cout<<"error "<<error<<": "<<str<<std::endl;
}

RTCDevice initializeDevice()
{
  RTCDevice device = rtcNewDevice(NULL);

  if (!device)
    printf("error %d: cannot create device\n", rtcGetDeviceError(NULL));

  rtcSetDeviceErrorFunction(device, errorFunction, NULL);
  return device;
}

void visualisation(const RTCRayHit& rayhit,std::ofstream& output){
        if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID){
            float factor = rayhit.ray.tfar/ sqrt(pow(rayhit.ray.dir_x, 2)+ pow(rayhit.ray.dir_y, 2)+ pow(rayhit.ray.dir_z, 2));
            float outX = rayhit.ray.org_x + factor * rayhit.ray.dir_x;
            float outY = rayhit.ray.org_y + factor * rayhit.ray.dir_y;
            float outZ = rayhit.ray.org_z + factor * rayhit.ray.dir_z;
            // std::cout<<rayhit.hit.primID<<std::endl;
            output<<outX<<" "<<outY<<" "<<outZ<<std::endl; 
        }
}

int main(int argc, char *argv[])
{   
    bool help = false;
    bool offFile = false;
    bool visual = false;

    for (int i = 1; i < argc; i++) {
        if ( (strcmp( "-h", argv[i]) == 0) || (strcmp( "-help", argv[i]) == 0)) 
            help = true;
        else if ( strcmp( "-o", argv[i]) == 0)
            offFile = true; 
        else if ( strcmp( "-v", argv[i]) == 0)
            visual = true; 
    }
    if(argc == 1 || help){
        std::cerr << "Usage: " << argv[0] << " <infile> <NumberOfRays> <XPoint> <YPoint> <ZPoint> <-o>[if the input file is .off]"<< std::endl;
        return 0;
    }
    else if(argc<5){
        std::cerr << "Too less arguments."<<std::endl;
        return 0;
    }

    const char* filename =  argv[1];
    std::ifstream input(filename);

    std::ofstream output;
    if(visual) output.open("EmbreeOut.xyz");

    std::stringstream ss(argv[2]);
    int _numberOfRays = 0;
    ss >> _numberOfRays ;

    float _xPoint, _yPoint, _zPoint;
    ss = std::stringstream(argv[3]);
    ss >> _xPoint;

    ss = std::stringstream(argv[4]);
    ss >> _yPoint;

    ss = std::stringstream(argv[5]);
    ss >> _zPoint;

    Mesh surfaceMesh;
        if(offFile)
        input >> surfaceMesh;
    else
        CGAL::read_ply(input, surfaceMesh);

    RTCDevice device = initializeDevice();

    RTCScene scene = rtcNewScene(device);
    // rtcSetSceneFlags(scene, RTC_SCENE_FLAG_ROBUST);

    RTCGeometry mesh = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

    Vertex* vertices = (Vertex*) rtcSetNewGeometryBuffer(mesh,
                                                         RTC_BUFFER_TYPE_VERTEX,
                                                         0,
                                                         RTC_FORMAT_FLOAT3,
                                                         sizeof(Vertex),
                                                         surfaceMesh.number_of_vertices());

    Triangle* triangles = (Triangle*) rtcSetNewGeometryBuffer(mesh,
                                                              RTC_BUFFER_TYPE_INDEX,
                                                              0,
                                                              RTC_FORMAT_UINT3,
                                                              sizeof(Triangle),
                                                              surfaceMesh.number_of_faces());
    
    if ( vertices && triangles) ;
    else
        std::cout<<"Error in creating Buffer objects."<<std::endl;

    for(vertex_descriptor vd : surfaceMesh.vertices()){
        Point data = surfaceMesh.point(vd);
        vertices[vd.idx()].x = data.x();       
        vertices[vd.idx()].y = data.y();       
        vertices[vd.idx()].z = data.z();       
    }

    
    for (face_descriptor fd : surfaceMesh.faces()){
        Mesh::Halfedge_index hf = surfaceMesh.halfedge(fd);
        int temp[3]; int i=0;
        for(Mesh::Halfedge_index hi : halfedges_around_face(hf, surfaceMesh)){
            Mesh::Vertex_index vi = target(hi, surfaceMesh);
            temp[i] = vi.idx();
            i++;
        }
        triangles[fd.idx()].v0 = temp[0];
        triangles[fd.idx()].v1 = temp[1];
        triangles[fd.idx()].v2 = temp[2];
        // std::cout<<"Index: "<<fd.idx()<<std::endl;
        // std::cout<<triangles[fd.idx()].v0<<" "<<triangles[fd.idx()].v1<<" "<<triangles[fd.idx()].v2<<std::endl;
    }

    CGAL::Real_timer time;
    time.start();

    rtcCommitGeometry(mesh);
    
    rtcAttachGeometry(scene, mesh);
    rtcReleaseGeometry(mesh);

    rtcCommitScene(scene);

    time.stop();
    std::cout << "  Construction time: " << time.time() << std::endl;
    time.reset();


    int numberOfRays = _numberOfRays; /*NUMBER OF RAY QUERIES*/
    RaysGenerate rg(numberOfRays); 

    time.start();
    
    struct RTCIntersectContext context;
    rtcInitIntersectContext(&context);
    
    for(size_t n=0; n!=numberOfRays; ++n){
    

        struct RTCRayHit rayhit;
        rayhit.ray.org_x =  _xPoint; /*POINT.X*/ 
        rayhit.ray.org_y =  _yPoint; /*POINT.Y*/
        rayhit.ray.org_z =  _zPoint; /*POINT.Z*/

        rayhit.ray.dir_x = rg.normalisedRayDirections[n]._x;
        rayhit.ray.dir_y = rg.normalisedRayDirections[n]._y;
        rayhit.ray.dir_z = rg.normalisedRayDirections[n]._z;

        rayhit.ray.tnear = 0;
        rayhit.ray.tfar = std::numeric_limits<float>::infinity();
        rayhit.ray.mask = 0;
        rayhit.ray.flags = 0;

        rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
        rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID;

        rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

        rtcIntersect1(scene, &context, &rayhit);
        if (visual) visualisation(rayhit, output);

    }
    time.stop();
    std::cout << "  Function() time: " << time.time() << std::endl;

    rtcReleaseScene(scene);    
    rtcReleaseDevice(device);
    return 0;
}
