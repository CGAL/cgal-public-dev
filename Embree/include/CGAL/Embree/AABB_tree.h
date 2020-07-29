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
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <embree3/rtcore.h>

#include <boost/optional.hpp>

#include <type_traits>
#include <vector>
#include <limits>
#include <unordered_map>

namespace CGAL {
namespace Embree {

template<typename Ray, typename Segment>
struct Intersect_context : public RTCIntersectContext{
public:
  Ray ray;
  Segment segment;
  struct RTCRayHit rayhit;

  enum IntersectionType{
    FIRST = 0, ANY, ALL
  };

  enum QueryType{
    RAY_QUERY = 0, SEGMENT_QUERY
  };

  IntersectionType intersectionType;
  QueryType queryType;

  Intersect_context(IntersectionType i_type, const Ray& _ray)
  : intersectionType(i_type), ray(_ray)
  {
    queryType = RAY_QUERY;
    init_context();
    init(_ray);
  }

  Intersect_context(IntersectionType i_type, const Segment& _segment)
  : intersectionType(i_type), segment(_segment)
  {
    queryType = SEGMENT_QUERY;
    init_context();
    init(_segment);
  }

  void init_context(){
    rtcInitIntersectContext(this);
  }

  void init(const Ray& _ray){
    init_rayhit(_ray);
  }

  void init(const Segment& _segment){
    float segmentLength = sqrt(_segment.squared_length());
    init_rayhit(_segment);
    rayhit.ray.tfar = segmentLength;
  }

  template<typename T>
  void init_rayhit(const T& query){

    rayhit.ray.org_x =  query.source().x(); /*POINT.X*/
    rayhit.ray.org_y =  query.source().y(); /*POINT.Y*/
    rayhit.ray.org_z =  query.source().z(); /*POINT.Z*/

    rayhit.ray.dir_x = query.direction().dx()/ sqrt(square(query.direction().dx()) + square(query.direction().dy()) + square(query.direction().dz()));
    rayhit.ray.dir_y = query.direction().dy()/ sqrt(square(query.direction().dx()) + square(query.direction().dy()) + square(query.direction().dz()));
    rayhit.ray.dir_z = query.direction().dz()/ sqrt(square(query.direction().dx()) + square(query.direction().dy()) + square(query.direction().dz()));

    rayhit.ray.tnear = 0;
    rayhit.ray.tfar = std::numeric_limits<float>::infinity();

    rayhit.ray.mask = 0;
    rayhit.ray.flags = 0;

    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID;

    rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
  }

};

template <typename TriangleMesh, typename ConstructibleFromId>
struct Id2descriptor
{};

template <typename TriangleMesh>
struct Id2descriptor<TriangleMesh,Tag_false>
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
struct Id2descriptor<TriangleMesh,Tag_true>
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

template <typename Geometry, typename GeomTraits>
class AABB_tree;

template <typename TriangleMesh, typename GeomTraits, typename ConstructibleFromId = CGAL::Boolean_tag<std::is_constructible<typename boost::graph_traits<TriangleMesh>::vertex_descriptor, unsigned int>::value> >
class Triangle_mesh_geometry {

  typedef Triangle_mesh_geometry<TriangleMesh, GeomTraits, ConstructibleFromId> Self;
  friend AABB_tree<Self, GeomTraits>;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  typedef typename boost::property_map<TriangleMesh, vertex_point_t>::const_type Vertex_point_map;

public:
  typedef std::pair<face_descriptor, TriangleMesh*> Primitive_id;
  typedef typename GeomTraits::Point_3 Point;
  typedef typename GeomTraits::Triangle_3 Triangle;
  typedef typename GeomTraits::Ray_3 Ray;
  typedef typename GeomTraits::Vector_3 Vector;
  typedef typename GeomTraits::Segment_3 Segment;

private:
  const TriangleMesh* surface_mesh;
  RTCGeometry rtc_geometry;
  unsigned int rtc_geomID;
  Vertex_point_map vpm;
  Id2descriptor<TriangleMesh, ConstructibleFromId> id2desc;
  std::vector<std::pair<float, unsigned int>> allIntersections;

public:
  Triangle_mesh_geometry()
  {}

  Triangle_mesh_geometry(const TriangleMesh& tm)
    : surface_mesh(&tm), vpm(get(CGAL::vertex_point, tm)), id2desc(tm)
  {}

  static void bound_function(const struct RTCBoundsFunctionArguments* args)
  {
    Triangle_mesh_geometry* self = (Triangle_mesh_geometry*) args->geometryUserPtr;
    RTCBounds* bounds_o = args->bounds_o;
    unsigned int primID = args->primID;

    Bbox_3 bb;

    face_descriptor fd = self->id2desc(primID);
    halfedge_descriptor hf = halfedge(fd, *self->surface_mesh);
    for(halfedge_descriptor hi : halfedges_around_face(hf, *(self->surface_mesh))){
        vertex_descriptor vi = target(hi, *(self->surface_mesh));
        bb += get(self->vpm,vi).bbox();
    }
    bounds_o->lower_x = bb.xmin();
    bounds_o->lower_y = bb.ymin();
    bounds_o->lower_z = bb.zmin();
    bounds_o->upper_x = bb.xmax();
    bounds_o->upper_y = bb.ymax();
    bounds_o->upper_z = bb.zmax();
  }

  static void intersection_function(const RTCIntersectFunctionNArguments* args)
  {
    typedef Intersect_context<Ray, Segment> Intersect_context;
    Triangle_mesh_geometry* self = (Triangle_mesh_geometry*) args->geometryUserPtr;
    int* valid = args->valid;
    Intersect_context* context = (Intersect_context*)args->context;
    if (!valid[0]) {
      return;
    }

    struct RTCRayHit* rayhit = (RTCRayHit*)args->rayhit;
    unsigned int primID = args->primID;
    assert(args->N == 1);
    std::vector<Point> face_points;
    face_points.reserve(3);

    face_descriptor fd = self->id2desc(primID);
    halfedge_descriptor hf = halfedge(fd, *self->surface_mesh);
    for(halfedge_descriptor hi : halfedges_around_face(hf, *(self->surface_mesh))){
        vertex_descriptor vi = target(hi, *(self->surface_mesh));
        Point data = get(self->vpm,vi);
        face_points.push_back(data);
    }
    Triangle face(face_points[0], face_points[1], face_points[2]);
    
    auto v = context->queryType==Intersect_context::QueryType::RAY_QUERY 
    ? CGAL::intersection(context->ray, face) 
    : CGAL::intersection(context->segment, face);
    if(v){
        rayhit->hit.geomID = self->rtc_geomID;
        rayhit->hit.primID = primID;
        if (const Point *intersection_point = boost::get<Point>(&*v) ){
            float _distance = context->queryType==Intersect_context::QueryType::RAY_QUERY
            ? sqrt(CGAL::squared_distance(context->ray.source(), *intersection_point))
            : sqrt(CGAL::squared_distance(context->segment.source(), *intersection_point));
            if(context->intersectionType == Intersect_context::IntersectionType::FIRST)
              rayhit->ray.tfar = _distance;
            else if (context->intersectionType == Intersect_context::IntersectionType::ALL)
              self->allIntersections.push_back(std::make_pair(_distance, primID));
            else {
              rayhit->ray.tfar = _distance;
              // Makes the ray invalid, so there is no further traversal
              rayhit->ray.tnear = rayhit->ray.tfar + 1.0f;
            }
        }
    }
  }


  // static void closest_point(RTCPointQueryFunctionArguments* args)
  // {
  //   Triangle_mesh_geometry* self = (Triangle_mesh_geometry*) args->geometryUserPtr;
  //   unsigned int primID = args->primID;
  //   face_descriptor fd = self->id2desc(primID)

  //   // query position in world space
  //   Vec3fa q(args->query->x, args->query->y, args->query->z);
  // }


  void insert_primitives()
  {
    rtcSetGeometryUserPrimitiveCount(rtc_geometry, num_faces(*surface_mesh));
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
    return std::make_pair(id2desc(primID), const_cast<TriangleMesh*>(surface_mesh));
  }

  inline std::vector<std::pair<float, unsigned int>> getIntersections(){ return allIntersections; }
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

  typedef Bbox_3   Bounding_box;
  typedef std::size_t size_type;
  typedef typename Geometry::Ray Ray;
  typedef typename Geometry::Segment Segment;

private:
  RTCDevice device;
  RTCScene scene;

  std::unordered_map<unsigned int, Geometry*> id2geometry;
  std::vector<Geometry> geometries;

public:
  AABB_tree()
  {
    rtc_bind();
  }

    AABB_tree(bool robust)
    :AABB_tree()
  {
    if(robust)
      rtcSetSceneFlags(scene, RTC_SCENE_FLAG_ROBUST);
  }

  ~AABB_tree()
  {
    clear();
  }

  void rtc_bind()
  {
    device = rtcNewDevice(NULL);
    scene = rtcNewScene(device);
  }

  void rtc_unbind()
  {
    rtcReleaseScene(scene);
    rtcReleaseDevice(device);
  }

  bool empty() const
  {
    if (!geometries.size()) return false;
    return true;
  }


  void clear()
  {
    rtc_unbind();
    if (this->empty()) return;
    geometries.clear();
    id2geometry.clear();
  }


  Bounding_box bbox() const
  {
    Bounding_box bb;

    for (size_type i =0; i!=geometries.size(); ++i){
      Geometry g = geometries[i];
      bb += Polygon_mesh_processing::bbox(*(g->surface_mesh));
    }

    return bb;
  }


  size_type size() const
  {
    size_type number_of_primitives = 0;
    for (size_type i =0; i!=geometries.size(); ++i){
      Geometry g = geometries[i];
      number_of_primitives+= num_faces(*g.surface_mesh);
    }
    return number_of_primitives;
  }

  /// T is the surface mesh
  template<typename T>
  void insert (const T& t)
  {
    geometries.push_back(Geometry(t));
    Geometry* geometry = &(geometries.back());

    geometry->rtc_geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
    geometry->rtc_geomID = rtcAttachGeometry(scene, geometry->rtc_geometry);
    id2geometry.insert({geometry->rtc_geomID, geometry});
    geometry->insert_primitives();
    rtcCommitScene(scene);
  }

  /**
   *
   */

  template<typename Query>
  bool do_intersect(const Query& query) const
  {
    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::IntersectionType::ANY, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return false;
    }
    return true;
  }

  template<typename Query>
  size_type number_of_intersected_primitives(const Query& query) const
  {
    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::IntersectionType::ANY, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    Geometry* geometry = id2geometry.at(rtc_geomID);
    return (geometry->getIntersections()).size();
  }

  template<typename Query>
  boost::optional<Intersection_and_primitive_id> first_intersection(const Query& query) const
  {
    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::IntersectionType::ANY, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return boost::none;
    }

    float factor = context.rayhit.ray.tfar/ sqrt(square(context.rayhit.ray.dir_x)+ square(context.rayhit.ray.dir_y)+ square(context.rayhit.ray.dir_z));
    float outX = context.rayhit.ray.org_x + factor * context.rayhit.ray.dir_x;
    float outY = context.rayhit.ray.org_y + factor * context.rayhit.ray.dir_y;
    float outZ = context.rayhit.ray.org_z + factor * context.rayhit.ray.dir_z;
    typename Geometry::Point p(outX, outY, outZ);

    Geometry* geometry = id2geometry.at(rtc_geomID);
    return boost::make_optional(std::make_pair(p, geometry->primitive_id(context.rayhit.hit.primID)));
  }


  template<typename Query>
  boost::optional<Primitive_id> first_intersected_primitive(const Query& query) const
  {
    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::IntersectionType::ANY, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return boost::none;
    }

    Geometry* geometry = id2geometry.at(rtc_geomID);

    return boost::make_optional(geometry->primitive_id(context.rayhit.hit.primID));
  }

  template<typename Query, typename OutputIterator>
  OutputIterator all_intersections(const Query& query, OutputIterator out) const
  {
    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::IntersectionType::ANY, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    Geometry* geometry = id2geometry.at(rtc_geomID);
    std::vector<std::pair<float, unsigned int>> intersectionDistance = geometry->getIntersections();

    for(int i=0; i<intersectionDistance.size();i++){
      float factor = intersectionDistance[i].first/ sqrt(square(context.rayhit.ray.dir_x)+ square(context.rayhit.ray.dir_y)+ square(context.rayhit.ray.dir_z));
      float outX = context.rayhit.ray.org_x + factor * context.rayhit.ray.dir_x;
      float outY = context.rayhit.ray.org_y + factor * context.rayhit.ray.dir_y;
      float outZ = context.rayhit.ray.org_z + factor * context.rayhit.ray.dir_z;
      typename Geometry::Point p(outX, outY, outZ);

      *out++ = boost::make_optional(std::make_pair(p, geometry->primitive_id(intersectionDistance[i].second)));
    }
    // out stores the following type  ----->  boost::optional<Intersection_and_primitive_id>
    return out;

  }

  template<typename Query>
  boost::optional<Intersection_and_primitive_id> any_intersection(const Query& query) const
  {
    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::IntersectionType::ANY, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return boost::none;
    }

    float factor = context.rayhit.ray.tfar/ sqrt(square(context.rayhit.ray.dir_x)+ square(context.rayhit.ray.dir_y)+ square(context.rayhit.ray.dir_z));
    float outX = context.rayhit.ray.org_x + factor * context.rayhit.ray.dir_x;
    float outY = context.rayhit.ray.org_y + factor * context.rayhit.ray.dir_y;
    float outZ = context.rayhit.ray.org_z + factor * context.rayhit.ray.dir_z;
    typename Geometry::Point p(outX, outY, outZ);

    Geometry* geometry = id2geometry.at(rtc_geomID);
    return boost::make_optional(std::make_pair(p, geometry->primitive_id(context.rayhit.hit.primID)));
  }

};

} // namespace Embree
} // namespace CGAL
#endif // CGAL_EMBREE_AABB_TREE_H
