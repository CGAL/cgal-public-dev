#include <iostream>
#include <fstream>
#include <limits>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <embree3/rtcore.h> 

#include <CGAL/Timer.h>

#include "RaysGenerate.h"

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;

struct Vertex   { float x,y,z,r; }; 
struct Triangle { int v0, v1, v2; };

int main(int argc, char *argv[])
{   
    bool help = false;
    for (int i = 1; i < argc; i++) { 
        if ( (strcmp( "-h", argv[i]) == 0) || (strcmp( "-help", argv[i]) == 0)) 
            help = true;
    }
    if(argc == 1 || help){
        std::cerr << "Usage: " << argv[0] << " <infile> <NumberOfRays> <XPoint> <YPoint> <ZPoint>"<< std::endl;
        return 0;
    }
    else if(argc<5){
        std::cerr << "Too less arguments."<<std::endl;
        return 0;
    }

    const char* filename =  argv[1];
    std::ifstream input(filename);

    std::stringstream ss(argv[2]);
    int _numberOfRays = 0;
    ss >> _numberOfRays ;

    double _xPoint, _yPoint, _zPoint;
    ss = std::stringstream(argv[3]);
    ss >> _xPoint;

    ss = std::stringstream(argv[4]);
    ss >> _yPoint;

    ss = std::stringstream(argv[5]);
    ss >> _zPoint;

    Mesh surfaceMesh;
    CGAL::read_ply(input, surfaceMesh);

    RTCDevice device = rtcNewDevice("verbose=0");
    RTCScene scene = rtcNewScene(device);

    RTCGeometry mesh = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

    Vertex* vertices = (Vertex*) rtcSetNewGeometryBuffer(mesh,RTC_BUFFER_TYPE_VERTEX,0,RTC_FORMAT_FLOAT3,sizeof(Vertex),surfaceMesh.number_of_vertices());

    for(vertex_descriptor vd : surfaceMesh.vertices()){
        Point data = surfaceMesh.point(vd);
        vertices[vd.idx()].x = data.x();       
        vertices[vd.idx()].y = data.y();       
        vertices[vd.idx()].z = data.z();       
    }

    Triangle* triangles = (Triangle*) rtcSetNewGeometryBuffer(mesh,RTC_BUFFER_TYPE_INDEX,0,RTC_FORMAT_UINT3,sizeof(Triangle),surfaceMesh.number_of_faces());
    
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
    }
    CGAL::Timer time;
    time.start();

    rtcCommitGeometry(mesh);
    unsigned int geomID = rtcAttachGeometry(scene, mesh);
    rtcReleaseGeometry(mesh);

    rtcCommitScene(scene);

    time.stop();
    std::cout << "  Construction time: " << time.time() << std::endl;
    time.reset();

    RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    RTCRayHit rayhit;
    rayhit.ray.org_x =  _xPoint; /*POINT.X*/ 
    rayhit.ray.org_y =  _yPoint; /*POINT.Y*/
    rayhit.ray.org_z =  _zPoint; /*POINT.Z*/

    rayhit.ray.tnear = 0.0;
    rayhit.ray.tfar = std::numeric_limits<double>::infinity();
    rayhit.ray.flags = 0;

    // rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID;
    // rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    
    int numberOfRays = _numberOfRays; /*NUMBER OF RAY QUERIES*/
    RaysGenerate rg(numberOfRays); 
    time.start();
    for(size_t n=0; n!=numberOfRays; ++n){
        rayhit.ray.dir_x = rg.rayDirections[n]._x;
        rayhit.ray.dir_x = rg.rayDirections[n]._y;
        rayhit.ray.dir_x = rg.rayDirections[n]._z;

        rtcIntersect1(scene, &context, &rayhit);

    }
    time.stop();
    std::cout << "  Function() time: " << time.time() << std::endl;
        
    rtcReleaseDevice(device);
    return 0;
}
