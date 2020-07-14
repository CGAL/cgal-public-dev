// Copyright (c) 2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri, Ayush Saraswat


#ifndef CGAL_EMBREE_AABB_TREE_H
#define CGAL_EMBREE_AABB_TREE_H

// #include <CGAL/license/Embree.h>
#include <CGAL/intersections.h>
#include <CGAL/boost/graph/helpers.h>

#include <embree3/rtcore.h>

#include <boost/optional.hpp>

#include <vector>
#include <limits>
#include <unordered_map>

namespace CGAL {
namespace Embree {


template <typename TriangleMesh, bool b>
struct Id2descriptor
{};
 
template <typename TriangleMesh>
struct Id2descriptor<TriangleMesh,true>
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  std::vector<face_descriptor> faceDescriptors;

  Id2descriptor(){}
  Id2descriptor(const TriangleMesh& tm)
  {
    for(face_descriptor fd : faces(tm))
    {
      faceDescriptors.push_back(fd);
    }
  }
 
  face_descriptor operator()(unsigned int i) const
  {
     return faceDescriptors[i];
  }
};
 
template <typename TriangleMesh>
struct Id2descriptor<TriangleMesh,false>
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  
  Id2descriptor(){}
  Id2descriptor(const TriangleMesh& tm)
  {}
 
 
  face_descriptor operator()(unsigned int i) const
  {
     return face_descriptor(i);
  }
};


  // AF: This is what you had called SM

template <typename TriangleMesh, typename GeomTraits, bool b>
struct Triangle_mesh_geometry {

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  typedef typename boost::property_map<TriangleMesh, vertex_point_t>::type Vertex_point_map;

  typedef std::pair<face_descriptor, TriangleMesh*> Primitive_id;
  typedef typename GeomTraits::Point_3 Point;
  typedef typename GeomTraits::Triangle_3 Triangle;
  typedef typename GeomTraits::Ray_3 Ray;
  typedef typename GeomTraits::Vector_3 Vector;

  const TriangleMesh* surfaceMesh;
  RTCGeometry rtc_geometry;
  unsigned int rtc_geomID;
  Vertex_point_map vpm;

  Id2descriptor<TriangleMesh, b> id2desc;

  Triangle_mesh_geometry()
  {}

  Triangle_mesh_geometry(const TriangleMesh& tm)
    : surfaceMesh(&tm), vpm(get(CGAL::vertex_point, const_cast<TriangleMesh&>(tm)))
  {
    Id2descriptor<TriangleMesh, b> _id2desc(tm);
    id2desc = _id2desc;
  }

  static void bound_function(const struct RTCBoundsFunctionArguments* args)
  {
    Triangle_mesh_geometry* self = (Triangle_mesh_geometry*) args->geometryUserPtr;
    RTCBounds* bounds_o = args->bounds_o;
    unsigned int primID = args->primID;

    std::vector<Point> FacePoints;
    face_descriptor fd = self->id2desc(primID);
    halfedge_descriptor hf = halfedge(fd, *self->surfaceMesh);
    for(halfedge_descriptor hi : halfedges_around_face(hf, *(self->surfaceMesh))){
        vertex_descriptor vi = target(hi, *(self->surfaceMesh));
        Point data = get(self->vpm,vi);
        FacePoints.push_back(data);
    }
    bounds_o->lower_x = std::min({FacePoints[0].x(), FacePoints[1].x(), FacePoints[2].x()});
    bounds_o->lower_y = std::min({FacePoints[0].y(), FacePoints[1].y(), FacePoints[2].y()});
    bounds_o->lower_z = std::min({FacePoints[0].z(), FacePoints[1].z(), FacePoints[2].z()});
    bounds_o->upper_x = std::max({FacePoints[0].x(), FacePoints[1].x(), FacePoints[2].x()});
    bounds_o->upper_y = std::max({FacePoints[0].y(), FacePoints[1].y(), FacePoints[2].y()});
    bounds_o->upper_z = std::max({FacePoints[0].z(), FacePoints[1].z(), FacePoints[2].z()});

  }

  static void intersection_function(const RTCIntersectFunctionNArguments* args)
  {
    Triangle_mesh_geometry* self = (Triangle_mesh_geometry*) args->geometryUserPtr;
    int* valid = args->valid;
    struct RTCRayHit* rayhit = (RTCRayHit*)args->rayhit;
    unsigned int primID = args->primID;

    assert(args->N == 1);
    std::vector<Point> FacePoints;
    if (!valid[0]) return;

    // face_descriptor fd(primID);
    face_descriptor fd = self->id2desc(primID);
    halfedge_descriptor hf = halfedge(fd, *self->surfaceMesh);
    for(halfedge_descriptor hi : halfedges_around_face(hf, *(self->surfaceMesh))){
        vertex_descriptor vi = target(hi, *(self->surfaceMesh));
        Point data = get(self->vpm,vi);
        FacePoints.push_back(data);
    }
    Triangle face(FacePoints[0], FacePoints[1], FacePoints[2]);

    Vector rayDirection(rayhit->ray.dir_x, rayhit->ray.dir_y, rayhit->ray.dir_z);
    Point rayOrgin(rayhit->ray.org_x, rayhit->ray.org_y, rayhit->ray.org_z);
    Ray ray(rayOrgin, rayDirection);

    auto v = CGAL::intersection(ray, face);
    if(v){
        rayhit->hit.geomID = self->rtc_geomID;
        rayhit->hit.primID = primID;
        if (const Point *intersectionPoint = boost::get<Point>(&*v) ){
            float _distance = sqrt(CGAL::squared_distance(rayOrgin, *intersectionPoint));
            rayhit->ray.tfar = _distance;
        }
    }
  }

  void insert_primitives()
  {
    rtcSetGeometryUserPrimitiveCount(rtc_geometry, num_faces(*surfaceMesh));
    rtcSetGeometryUserData(rtc_geometry, this);

    // AF: For the next two you have to find out how to write
    // the function pointer for a static member function

    // Ayush: for a static member function, we can directly pass a pointer, so the below should work fine.
    // https://isocpp.org/wiki/faq/pointers-to-members

    rtcSetGeometryBoundsFunction(rtc_geometry, bound_function, nullptr);
    rtcSetGeometryIntersectFunction(rtc_geometry, intersection_function);
    rtcCommitGeometry(rtc_geometry);

    rtcReleaseGeometry(rtc_geometry);
  }

  Primitive_id primitive_id(unsigned int primID) const
  {
    return std::make_pair(id2desc(primID), const_cast<TriangleMesh*>(surfaceMesh));
  }
};

/**
 * \ingroup PkgEmbreeRef
 * This class...
 */


//  AF:  Geometry is the above class
// AF: For GeomTraits you take a kernel, that is Simple_cartesian
template <typename Geometry, typename GeomTraits>
class AABB_tree {

  typedef typename Geometry::Primitive_id Primitive_id;
public:
  typedef std::pair<typename Geometry::Point, Primitive_id> Intersection_and_primitive_id;

  RTCDevice device;
  RTCScene scene;

  std::unordered_map<unsigned int, Geometry> id2geometry;
  std::list<Geometry> geometries;

public:
  AABB_tree()
  {
    device = rtcNewDevice(NULL);
    scene = rtcNewScene(device);
  }



  /// T is the surface mesh
  template<typename T>
  void insert (const T& t)
  {
    geometries.push_back(Geometry(t));
    Geometry& geometry = geometries.back();

    geometry.rtc_geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
    geometry.rtc_geomID = rtcAttachGeometry(scene, geometry.rtc_geometry);
    id2geometry.insert({geometry.rtc_geomID, geometry});
    geometry.insert_primitives();
    rtcCommitScene(scene);
  }

  template<typename Ray>
  boost::optional<Intersection_and_primitive_id> first_intersection(const Ray& query) const
  {
    struct RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    struct RTCRayHit rayhit;

    rayhit.ray.org_x =  query.source().x(); /*POINT.X*/
    rayhit.ray.org_y =  query.source().y(); /*POINT.Y*/
    rayhit.ray.org_z =  query.source().z(); /*POINT.Z*/

    rayhit.ray.dir_x = query.direction().dx()/ sqrt(pow(query.direction().dx(), 2) + pow(query.direction().dy(), 2) + pow(query.direction().dz(), 2));
    rayhit.ray.dir_y = query.direction().dy()/ sqrt(pow(query.direction().dx(), 2) + pow(query.direction().dy(), 2) + pow(query.direction().dz(), 2));
    rayhit.ray.dir_z = query.direction().dz()/ sqrt(pow(query.direction().dx(), 2) + pow(query.direction().dy(), 2) + pow(query.direction().dz(), 2));

    rayhit.ray.tnear = 0;
    rayhit.ray.tfar = std::numeric_limits<float>::infinity();
    rayhit.ray.mask = 0;
    rayhit.ray.flags = 0;

    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID;

    rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

    rtcIntersect1(scene, &context, &rayhit);

    unsigned int rtc_geomID = rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return boost::none;
    }

    float factor = rayhit.ray.tfar/ sqrt(pow(rayhit.ray.dir_x, 2)+ pow(rayhit.ray.dir_y, 2)+ pow(rayhit.ray.dir_z, 2));
    float outX = rayhit.ray.org_x + factor * rayhit.ray.dir_x;
    float outY = rayhit.ray.org_y + factor * rayhit.ray.dir_y;
    float outZ = rayhit.ray.org_z + factor * rayhit.ray.dir_z;
    typename Geometry::Point p(outX, outY, outZ);

    Geometry geometry = id2geometry.at(rtc_geomID);
    return boost::make_optional(std::make_pair(p, geometry.primitive_id(rayhit.hit.primID)));
  }


  template<typename Ray>
  boost::optional<Primitive_id> first_intersected_primitive(const Ray& query) const
  {
    struct RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    struct RTCRayHit rayhit;

    rayhit.ray.org_x =  query.source().x(); /*POINT.X*/
    rayhit.ray.org_y =  query.source().y(); /*POINT.Y*/
    rayhit.ray.org_z =  query.source().z(); /*POINT.Z*/

    rayhit.ray.dir_x = query.direction().dx()/ sqrt(pow(query.direction().dx(), 2) + pow(query.direction().dy(), 2) + pow(query.direction().dz(), 2));
    rayhit.ray.dir_y = query.direction().dy()/ sqrt(pow(query.direction().dx(), 2) + pow(query.direction().dy(), 2) + pow(query.direction().dz(), 2));
    rayhit.ray.dir_z = query.direction().dz()/ sqrt(pow(query.direction().dx(), 2) + pow(query.direction().dy(), 2) + pow(query.direction().dz(), 2));

    rayhit.ray.tnear = 0;
    rayhit.ray.tfar = std::numeric_limits<float>::infinity();
    rayhit.ray.mask = 0;
    rayhit.ray.flags = 0;

    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID;

    rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

    rtcIntersect1(scene, &context, &rayhit);

    unsigned int rtc_geomID = rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return boost::none;
    }

    Geometry geometry = id2geometry.at(rtc_geomID);

    return boost::make_optional(geometry.primitive_id(rayhit.hit.primID));
  }


};

} // namespace Embree
} // namespace CGAL
#endif // CGAL_EMBREE_AABB_TREE_H
