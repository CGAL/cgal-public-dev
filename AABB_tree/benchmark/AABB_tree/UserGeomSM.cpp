#include <iostream>
#include <limits>
#include <fstream>

#include <embree3/rtcore.h> 

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include "RaysGenerate.h"

typedef CGAL::Simple_cartesian<float> K;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef K::Ray_3 Ray;
typedef K::Vector_3 Vector;
typedef CGAL::Surface_mesh<Point> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;
 
struct SM
{   
    Mesh surfaceMesh;
    RTCGeometry geometry;
    unsigned int geomID;
};


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
        output<<outX<<" "<<outY<<" "<<outZ<<std::endl; 
    }
}

void SmBoundFunction(const struct RTCBoundsFunctionArguments* args){
    SM* sm = (SM*) args->geometryUserPtr;
    RTCBounds* bounds_o = args->bounds_o;
    unsigned int primID = args->primID;

    std::vector<Point> FacePoints;

    Mesh::Face_index fd(primID); 
    Mesh::Halfedge_index hf = sm->surfaceMesh.halfedge(fd);
    for(Mesh::Halfedge_index hi : halfedges_around_face(hf, sm->surfaceMesh)){
        Mesh::Vertex_index vi = target(hi, sm->surfaceMesh  );
        Point data = sm->surfaceMesh.point(vi);
        FacePoints.push_back(data);
    }
    bounds_o->lower_x = std::min({FacePoints[0].x(), FacePoints[1].x(), FacePoints[2].x()});
    bounds_o->lower_y = std::min({FacePoints[0].y(), FacePoints[1].y(), FacePoints[2].y()});
    bounds_o->lower_z = std::min({FacePoints[0].z(), FacePoints[1].z(), FacePoints[2].z()});
    bounds_o->upper_x = std::max({FacePoints[0].x(), FacePoints[1].x(), FacePoints[2].x()});
    bounds_o->upper_y = std::max({FacePoints[0].y(), FacePoints[1].y(), FacePoints[2].y()});
    bounds_o->upper_z = std::max({FacePoints[0].z(), FacePoints[1].z(), FacePoints[2].z()});
}

void SmIntersectionFunction(const RTCIntersectFunctionNArguments* args){
    
    int* valid = args->valid;
    SM* sm = (SM*) args->geometryUserPtr;
    struct RTCRayHit* rayhit = (RTCRayHit*)args->rayhit;
    unsigned int primID = args->primID;

    assert(args->N == 1);
    
    std::vector<Point> FacePoints;
    
    if (!valid[0]) return;

    Mesh::Face_index fd(primID); 
    Mesh::Halfedge_index hf = sm->surfaceMesh.halfedge(fd);
    for(Mesh::Halfedge_index hi : halfedges_around_face(hf, sm->surfaceMesh)){
        Mesh::Vertex_index vi = target(hi, sm->surfaceMesh  );
        Point data = sm->surfaceMesh.point(vi);
        FacePoints.push_back(data);
    }
    Triangle face(FacePoints[0], FacePoints[1], FacePoints[2]);

    Vector rayDirection(rayhit->ray.dir_x, rayhit->ray.dir_y, rayhit->ray.dir_z);
    Point rayOrgin(rayhit->ray.org_x, rayhit->ray.org_y, rayhit->ray.org_z);
    Ray ray(rayOrgin, rayDirection);

    auto v = CGAL::intersection(ray, face);
    if(v){
        rayhit->hit.geomID = sm->geomID;
        rayhit->hit.primID = primID;
        if (const Point *intersectionPoint = boost::get<Point>(&*v) ){
            float _distance = sqrt(CGAL::squared_distance(rayOrgin, *intersectionPoint));
            rayhit->ray.tfar = _distance;
        }
    }
}

int main(int argc, char  *argv[])
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
        std::cerr << "Usage: " << argv[0] << " <infile> <NumberOfRays> <XPoint> <YPoint> <ZPoint> <-o>[if the input file is .off]  <-v>[for visualisation]"<< std::endl;
        return 0;
    }
    else if(argc<5){
        std::cerr << "Too less arguments."<<std::endl;
        return 0;
    }

    const char* filename =  argv[1];
    std::ifstream input(filename);

    std::ofstream output;
    if(visual) output.open("UserGeomOut.xyz");

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

    /*Embree initialisation*/
    RTCDevice device = initializeDevice();
    RTCScene scene = rtcNewScene(device);

    RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
    SM sm;
    sm.surfaceMesh = surfaceMesh;
    sm.geometry = geom;
    sm.geomID = rtcAttachGeometry(scene, geom);

    rtcSetGeometryUserPrimitiveCount(geom, surfaceMesh.number_of_faces());
    rtcSetGeometryUserData(geom, &sm);
    rtcSetGeometryBoundsFunction(geom, SmBoundFunction, nullptr);
    rtcSetGeometryIntersectFunction(geom, SmIntersectionFunction);
    rtcCommitGeometry(geom);
    rtcReleaseGeometry(geom);

    rtcCommitScene(scene);

    int numberOfRays = _numberOfRays; /*NUMBER OF RAY QUERIES*/
    RaysGenerate rg(numberOfRays);

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

    rtcReleaseScene(scene);    
    rtcReleaseDevice(device);
    return 0;
}
