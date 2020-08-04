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
#include <deque>
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

  enum Intersection_type{
    FIRST = 0, ANY, ALL
  };

  enum Query_type{
    RAY_QUERY = 0, SEGMENT_QUERY
  };

  Intersection_type intersection_type;
  Query_type query_type;

  Intersect_context(Intersection_type i_type, const Ray& _ray)
  : intersection_type(i_type), ray(_ray)
  {
    query_type = RAY_QUERY;
    init_context();
    init(_ray);
  }

  Intersect_context(Intersection_type i_type, const Segment& _segment)
  : intersection_type(i_type), segment(_segment)
  {
    query_type = SEGMENT_QUERY;
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
  const TriangleMesh* triangle_mesh;
  RTCGeometry rtc_geometry;
  unsigned int rtc_geomID;
  Vertex_point_map vpm;
  Id2descriptor<TriangleMesh, ConstructibleFromId> id2desc;
  std::vector<std::pair<float, unsigned int>> all_intersections;

public:
  Triangle_mesh_geometry()
  {}

  Triangle_mesh_geometry(const TriangleMesh& tm)
    : triangle_mesh(&tm), vpm(get(CGAL::vertex_point, tm)), id2desc(tm)
  {}

  static void bound_function(const struct RTCBoundsFunctionArguments* args)
  {
    Triangle_mesh_geometry* self = (Triangle_mesh_geometry*) args->geometryUserPtr;
    RTCBounds* bounds_o = args->bounds_o;
    unsigned int primID = args->primID;

    Bbox_3 bb;

    face_descriptor fd = self->id2desc(primID);
    halfedge_descriptor hf = halfedge(fd, *self->triangle_mesh);
    for(halfedge_descriptor hi : halfedges_around_face(hf, *(self->triangle_mesh))){
        vertex_descriptor vi = target(hi, *(self->triangle_mesh));
        bb += get(self->vpm,vi).bbox();
    }
    bounds_o->lower_x = bb.xmin();
    bounds_o->lower_y = bb.ymin();
    bounds_o->lower_z = bb.zmin();
    bounds_o->upper_x = bb.xmax();
    bounds_o->upper_y = bb.ymax();
    bounds_o->upper_z = bb.zmax();
  }


  static Triangle triangle(const Self* tmg, face_descriptor fd)
  {
    const TriangleMesh& tm = *(tmg->triangle_mesh);
    halfedge_descriptor hd = halfedge(fd, tm);
    return Triangle(get(tmg->vpm, target(hd, tm)),
                    get(tmg->vpm, target(next(hd, tm), tm)),
                    get(tmg->vpm, source(hd, tm)));
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

    Triangle face = triangle(self, fd);

    auto v = context->query_type==Intersect_context::Query_type::RAY_QUERY
    ? CGAL::intersection(context->ray, face)
    : CGAL::intersection(context->segment, face);
    if(v){
        rayhit->hit.geomID = self->rtc_geomID;
        rayhit->hit.primID = primID;
        if (const Point *intersection_point = boost::get<Point>(&*v) ){
            auto _distance = context->query_type==Intersect_context::Query_type::RAY_QUERY
            ? sqrt(CGAL::squared_distance(context->ray.source(), *intersection_point))
            : sqrt(CGAL::squared_distance(context->segment.source(), *intersection_point));
            if(context->intersection_type == Intersect_context::Intersection_type::FIRST)
              rayhit->ray.tfar = _distance;
            else if (context->intersection_type == Intersect_context::Intersection_type::ALL)
              self->all_intersections.push_back(std::make_pair(_distance, primID));
            else {
              rayhit->ray.tfar = _distance;
              // Makes the ray invalid, so there is no further traversal
              rayhit->ray.tnear = rayhit->ray.tfar + 1.0f;
            }
        }
    }
  }


  template <typename ClosestPointResult>
  bool closest_point_function(RTCPointQueryFunctionArguments* args, ClosestPointResult* result) const
  {
    unsigned int primID = args->primID;
    unsigned int geomID = args->geomID;
    face_descriptor fd = id2desc(primID);

    Triangle t = triangle(this, fd);

    GeomTraits gt;
    Point cp = gt.construct_projected_point_3_object()(t, result->query);
    double d = sqrt(squared_distance(result->query, cp));
    if(d < args->query->radius){
      args->query->radius = d;
      result->result = cp;
      result->primID = primID;
      result->geomID = geomID;
      return true;
    }
    return false;
  }


  void insert_primitives()
  {
    rtcSetGeometryUserPrimitiveCount(rtc_geometry, num_faces(*triangle_mesh));
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
    return std::make_pair(id2desc(primID), const_cast<TriangleMesh*>(triangle_mesh));
  }

  inline const std::vector<std::pair<float, unsigned int>>& intersections() const
  {
    return all_intersections;
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

  typedef AABB_tree<Geometry,GeomTraits> Self;
  typedef typename Geometry::Primitive_id Primitive_id;
public:
  typedef std::pair<typename Geometry::Point, Primitive_id> Intersection_and_primitive_id;
  typedef std::pair<typename Geometry::Point, Primitive_id> Point_and_primitive_id;

  typedef Bbox_3   Bounding_box;
  typedef std::size_t size_type;
  typedef typename GeomTraits::Point_3 Point;
  typedef typename GeomTraits::Ray_3 Ray;
  typedef typename GeomTraits::Segment_3 Segment;

private:



struct Closest_point_result
{
  Closest_point_result(const Self* aabb_tree, const Point& query)
    : aabb_tree(aabb_tree)
    , query(query)
    , primID(RTC_INVALID_GEOMETRY_ID)
    , geomID(RTC_INVALID_GEOMETRY_ID)
  {}

  const Self* aabb_tree;
  Point query, result;
  unsigned int primID;
  unsigned int geomID;
};

  static bool closest_point_function(RTCPointQueryFunctionArguments* args)
  {
    unsigned int primID = args->primID;
    unsigned int geomID = args->geomID;
    Closest_point_result* result = (Closest_point_result*)args->userPtr;
    const Self* aabb_tree = result->aabb_tree;
    const Geometry& geometry = aabb_tree->geometries[geomID];
    return geometry.closest_point_function(args, result);
  }



  RTCDevice device;
  RTCScene scene;

  std::deque<Geometry> geometries;

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
    rtc_unbind();
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


  /// returns \c true, iff the tree contains no primitive.
  bool empty() const
  {
    if (geometries.empty()) return true;
    return false;
  }


  // clears the tree.
  void clear()
  {
    if (this->empty()) return;
    geometries.clear();
  }


  /// returns the axis-aligned bounding box of the whole tree.
  /// \pre `!empty()`
  Bounding_box bbox() const
  {
    Bounding_box bb;

    for (size_type i =0; i!=geometries.size(); ++i){
      const Geometry& g = geometries[i];
      bb += Polygon_mesh_processing::bbox(*(g.triangle_mesh));
    }

    return bb;
  }


  /// returns the number of primitives in the tree.
  size_type size() const
  {
    size_type number_of_primitives = 0;
    for (size_type i =0; i!=geometries.size(); ++i){
      const Geometry& g = geometries[i];
      number_of_primitives+= num_faces(*g.triangle_mesh);
    }
    return number_of_primitives;
  }

  /// T is the surface mesh
  template<typename T>
  void insert (const T& t)
  {
    unsigned int geomID = geometries.size();
    geometries.push_back(Geometry(t));
    Geometry& geometry = geometries.back();

    geometry.rtc_geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
    geometry.rtc_geomID =geomID;
    rtcAttachGeometryByID(scene, geometry.rtc_geometry, geomID);
    geometry.insert_primitives();

    rtcSetGeometryPointQueryFunction(geometry.rtc_geometry, closest_point_function);

    rtcCommitScene(scene);
  }

  /// \name Intersection Tests
  ///@{
  /**
   *
   */


  /// returns `true`, iff the query intersects at least one of
  /// the input primitives.
  /// \tparam Query may be `GeomTraits::Ray_3` or `GeomTraits::Segment_3`.
  template<typename Query>
  bool do_intersect(const Query& query) const
  {
    if (this->empty()) return false;

    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::ANY, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return false;
    }
    return true;
  }


  /// returns the number of primitives intersected by the query.
  /// \tparam Query may be `GeomTraits::Ray_3` or `GeomTraits::Segment_3`.
  template<typename Query>
  size_type number_of_intersected_primitives(const Query& query) const
  {
    if (this->empty()) return 0;

    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::ALL, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    const Geometry& geometry = geometries[rtc_geomID];
    return (geometry.intersections()).size();
  }


  /// puts in `out` the ids of all intersected primitives.
  /// This function does not compute the intersection points
  /// and is hence faster than the function `all_intersections()`
  /// function below.
  template<typename Query, typename OutputIterator>
  OutputIterator all_intersected_primitives (const Query& query, OutputIterator out) const
  {
    if (this->empty()) return out;

    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::ALL, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    const Geometry& geometry = geometries[rtc_geomID];
    const std::vector<std::pair<float, unsigned int>>& intersectionDistance = geometry.intersections();

    for(int i=0; i<intersectionDistance.size();i++){
      *out++ = boost::make_optional(geometry.primitive_id(intersectionDistance[i].second));
    }

    return out;
  }


  /// returns the id of the intersected primitive that is encountered first
  /// in the tree traversal, iff
  /// the query intersects at least one of the input primitives. No
  /// particular order is guaranteed over the tree traversal, such
  /// that, e.g, the primitive returned is not necessarily the
  /// closest from the source point of a ray query.
  /// \tparam Query may be `GeomTraits::Ray_3` or `GeomTraits::Segment_3`.
  template<typename Query>
  boost::optional<Primitive_id> any_intersected_primitive(const Query& query) const
  {
    if (this->empty()) return boost::none;

    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::ANY, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return boost::none;
    }

    const Geometry& geometry = geometries[rtc_geomID];

    return boost::make_optional(geometry.primitive_id(context.rayhit.hit.primID));
  }
  ///@}

  /// \name Intersections
  ///@{
  template<typename Query>
  boost::optional<Intersection_and_primitive_id> first_intersection(const Query& query) const
  {
    if (this->empty()) return boost::none;

    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::FIRST, query);

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

    const Geometry& geometry = geometries[rtc_geomID];
    return boost::make_optional(std::make_pair(p, geometry.primitive_id(context.rayhit.hit.primID)));
  }


  template<typename Query>
  boost::optional<Primitive_id> first_intersected_primitive(const Query& query) const
  {
    if (this->empty()) return boost::none;

    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::FIRST, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return boost::none;
    }

    const Geometry& geometry = geometries[rtc_geomID];

    return boost::make_optional(geometry.primitive_id(context.rayhit.hit.primID));
  }

  template<typename Query, typename OutputIterator>
  OutputIterator all_intersections(const Query& query, OutputIterator out) const
  {
    if (this->empty()) return out;

    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::ALL, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    const Geometry& geometry = geometries[rtc_geomID];
    const std::vector<std::pair<float, unsigned int>>& intersectionDistance = geometry.intersections();

    for(int i=0; i<intersectionDistance.size();i++){
      float factor = intersectionDistance[i].first/ sqrt(square(context.rayhit.ray.dir_x)+ square(context.rayhit.ray.dir_y)+ square(context.rayhit.ray.dir_z));
      float outX = context.rayhit.ray.org_x + factor * context.rayhit.ray.dir_x;
      float outY = context.rayhit.ray.org_y + factor * context.rayhit.ray.dir_y;
      float outZ = context.rayhit.ray.org_z + factor * context.rayhit.ray.dir_z;
      typename Geometry::Point p(outX, outY, outZ);

      *out++ = boost::make_optional(std::make_pair(p, geometry.primitive_id(intersectionDistance[i].second)));
    }
    // out stores the following type  ----->  boost::optional<Intersection_and_primitive_id>
    return out;

  }

  /// returns if any the intersection that is encountered first
  /// in the tree traversal. No particular
  /// order is guaranteed over the tree traversal, e.g, the
  /// primitive returned is not necessarily the closest from the query.
  ///
  /// \tparam Query may be `GeomTraits::Ray_3` or `GeomTraits::Segment_3`.
  template<typename Query>
  boost::optional<Intersection_and_primitive_id> any_intersection(const Query& query) const
  {
    if (this->empty()) return boost::none;

    typedef Intersect_context<Ray, Segment> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::ANY, query);

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

    const Geometry& geometry = geometries[rtc_geomID];
    return boost::make_optional(std::make_pair(p, geometry.primitive_id(context.rayhit.hit.primID)));
  }
 ///@}


  /// \name Distance Queries
  ///@{

  /// returns the point in the union of all input primitives which
  /// is closest to the query. In case there are several closest
  /// points, one arbitrarily chosen closest point is
  /// returned.
  /// \pre `!empty()`
  Point closest_point(const typename Geometry::Point &query) const
  {
    RTCPointQuery rtc_query;
    rtc_query.x = query.x();
    rtc_query.y = query.y();
    rtc_query.z = query.z();
    rtc_query.radius = std::numeric_limits<float>::infinity();
    rtc_query.time = 0.f;

    Closest_point_result result(this, query);
    RTCPointQueryContext context;
    rtcInitPointQueryContext(&context);
    rtcPointQuery(scene, &rtc_query, &context, nullptr, (void*)&result);

    return result.result;
  }


  /// returns a `Point_and_primitive_id` which realizes the
  /// smallest distance between the query point and all input
  /// primitives.
  /// \pre `!empty()`
  Point_and_primitive_id closest_point_and_primitive(const typename Geometry::Point &query) const
  {
    RTCPointQuery rtc_query;
    rtc_query.x = query.x();
    rtc_query.y = query.y();
    rtc_query.z = query.z();
    rtc_query.radius = std::numeric_limits<float>::infinity();
    rtc_query.time = 0.f;

    Closest_point_result result(this, query);
    RTCPointQueryContext context;
    rtcInitPointQueryContext(&context);
    rtcPointQuery(scene, &rtc_query, &context, nullptr, (void*)&result);

    const Geometry& geometry = geometries[result.geomID];
    return std::make_pair(result.result, geometry.primitive_id(result.geomID));
  }
 ///@}

};

} // namespace Embree
} // namespace CGAL
#endif // CGAL_EMBREE_AABB_TREE_H
