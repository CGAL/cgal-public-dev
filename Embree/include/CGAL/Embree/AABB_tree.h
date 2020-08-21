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

#define BOX_EXT FLT_EPSILON

namespace CGAL {
namespace Embree {

/**
 * \ingroup PkgEmbreeRef
 * Data structure created to store 
 * computations and return values of intersection and distance queries.
 * It takes care of the type of intersection and conveys information from one structure to another
 * in the pipeline.
 *
 * \tparam Geometry Surface_mesh or Polyhedron_3
 */
template<typename Geometry>
struct Intersect_context : public RTCIntersectContext{
/// \name Types
///@{
  typedef typename Geometry::Ray Ray;
  typedef typename Geometry::Segment Segment;
  typedef typename Geometry::Point Point;
  typedef std::tuple<Point, float, unsigned int> IntersectionData;
///@}

private:
  std::vector<IntersectionData> all_intersections;

public:
/// \name Attributes
///@{
  Ray ray;
  Segment segment;
  unsigned int counter =0;
  struct RTCRayHit rayhit;
///@}

/// \name Enum Types
///@{
  /// Enum type to supply the program with the type of information queries
  enum Intersection_type{
    FIRST = 0, ANY, ALL, COUNTER
  };
  /// as the name suggests this suplies the program with the object of querying.
  enum Query_type{
    RAY_QUERY = 0, SEGMENT_QUERY
  };
///@}

/// \name Attributes
///@{
  Intersection_type intersection_type;
  Query_type query_type;
///@}


/// \name Constructors
///@{
  /// Constructor to initialise when query type in a ray.
  Intersect_context(Intersection_type i_type, const Ray& _ray)
  : intersection_type(i_type), ray(_ray)
  {
    query_type = RAY_QUERY;
    init_context();
    init(_ray);
  }

  /// Constructor to initialise when query type in a segment.
  Intersect_context(Intersection_type i_type, const Segment& _segment)
  : intersection_type(i_type), segment(_segment)
  {
    query_type = SEGMENT_QUERY;
    init_context();
    init(_segment);
  }
///@}

/// \name Set-up
///@{
  /// initialise an EmbreeAPI context.
  void init_context(){
    rtcInitIntersectContext(this);
  }
  /// makes a call to init_rayhit for a ray.
  void init(const Ray& _ray){
    init_rayhit(_ray);
  }

  /// makes a call to init_rayhit for a segment.
  void init(const Segment& _segment){
    float segmentLength = sqrt(to_double(_segment.squared_length()));
    init_rayhit(_segment);
    rayhit.ray.tfar = segmentLength;
  }

  ///initialises the rayhit object used by the EmbreeAPI.
  /// \tparam T ray or segment
  template<typename T>
  void init_rayhit(const T& query){

    rayhit.ray.org_x =  to_double(query.source().x());
    rayhit.ray.org_y =  to_double(query.source().y());
    rayhit.ray.org_z =  to_double(query.source().z());

    double len = sqrt(square(to_double(query.direction().dx())) + square(to_double(query.direction().dy())) + square(to_double(query.direction().dz())));
    rayhit.ray.dir_x = to_double(query.direction().dx()) / len;
    rayhit.ray.dir_y = to_double(query.direction().dy()) / len;
    rayhit.ray.dir_z = to_double(query.direction().dz()) / len;

    rayhit.ray.tnear = 0;
    rayhit.ray.tfar = std::numeric_limits<float>::infinity();

    rayhit.ray.mask = 0;
    rayhit.ray.flags = 0;

    rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID;

    rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
  }
///@}

/// \name Helpers
///@{
  ///Getter to return the array which contain the intersection objects.
  inline std::vector<IntersectionData>& intersections()
  {
    return all_intersections;
  }
///@}

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


/**
 * \ingroup PkgEmbreeRef
 * Support structure for that stores a pointer to the userMesh, namely a surface_mesh or a polyhedron.
 * It also stores the RTCGeometry object used by the EmbreeAPI and a geomID assigned to it which 
 * can be used to query a geometry object. 
 *
 * \tparam TriangleMesh Surface_mesh or Polyhedron_3
 * \tparam GeomTraits a kernel
 */
template <typename TriangleMesh, typename GeomTraits, typename ConstructibleFromId = CGAL::Boolean_tag<std::is_constructible<typename boost::graph_traits<TriangleMesh>::vertex_descriptor, unsigned int>::value> >
class Triangle_mesh_geometry {

  typedef Triangle_mesh_geometry<TriangleMesh, GeomTraits, ConstructibleFromId> Self;
  friend AABB_tree<Self, GeomTraits>;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

  typedef typename boost::property_map<TriangleMesh, vertex_point_t>::const_type Vertex_point_map;

public:
/// \name Types
///@{
  typedef std::pair<face_descriptor, TriangleMesh*> Primitive_id;
  typedef typename GeomTraits::Point_3 Point;
  typedef typename GeomTraits::Triangle_3 Triangle;
  typedef typename GeomTraits::Ray_3 Ray;
  typedef typename GeomTraits::Vector_3 Vector;
  typedef typename GeomTraits::Segment_3 Segment;
  typedef typename GeomTraits::FT FT;
///@}

private:
  const TriangleMesh* triangle_mesh;
  RTCGeometry rtc_geometry;
  unsigned int rtc_geomID;
  Vertex_point_map vpm;
  Id2descriptor<TriangleMesh, ConstructibleFromId> id2desc;

public:
/// \name Constructors
///@{
  Triangle_mesh_geometry()
  {}

  Triangle_mesh_geometry(const TriangleMesh& tm)
    : triangle_mesh(&tm), vpm(get(CGAL::vertex_point, tm)), id2desc(tm)
  {}
///@}

/// \name Call-backs
///@{
  static void bound_function(const struct RTCBoundsFunctionArguments* args)
  {
    Triangle_mesh_geometry* self = (Triangle_mesh_geometry*) args->geometryUserPtr;
    RTCBounds* bounds_o = args->bounds_o;
    unsigned int primID = args->primID;

    const TriangleMesh& tm = *(self->triangle_mesh);
    face_descriptor fd = self->id2desc(primID);
    halfedge_descriptor hd = halfedge(fd, tm);

    Bbox_3 bb = get(self->vpm, target(hd, tm)).bbox();
    bb += get(self->vpm, target(next(hd, tm), tm)).bbox();
    bb += get(self->vpm, source(hd, tm)).bbox();

    bounds_o->lower_x = bb.xmin() == bb.xmax() ? (bb.xmin()-BOX_EXT) : bb.xmin();
    bounds_o->upper_x = bb.xmin() == bb.xmax() ? (bb.xmax()+BOX_EXT) : bb.xmax();

    bounds_o->lower_y = bb.ymin() == bb.ymax() ? (bb.ymin()-BOX_EXT) : bb.ymin();
    bounds_o->upper_y = bb.ymin() == bb.ymax() ? (bb.ymax()+BOX_EXT) : bb.ymax();

    bounds_o->lower_z = bb.zmin() == bb.zmax() ? (bb.zmin()-BOX_EXT) : bb.zmin();
    bounds_o->upper_z = bb.zmin() == bb.zmax() ? (bb.zmax()+BOX_EXT) : bb.zmax();
  }
///@}

/// \name Utility
///@{
  static Triangle triangle(const Self* tmg, face_descriptor fd)
  {
    const TriangleMesh& tm = *(tmg->triangle_mesh);
    halfedge_descriptor hd = halfedge(fd, tm);
    return Triangle(get(tmg->vpm, target(hd, tm)),
                    get(tmg->vpm, target(next(hd, tm), tm)),
                    get(tmg->vpm, source(hd, tm)));
  }
///@}


/// \name Call-backs
///@{
  static void intersection_function(const RTCIntersectFunctionNArguments* args)
  {
    typedef Intersect_context<Self> Intersect_context;
    Triangle_mesh_geometry* self = (Triangle_mesh_geometry*) args->geometryUserPtr;
    int* valid = args->valid;
    Intersect_context* context = (Intersect_context*)args->context;
    if (!valid[0]) {
      return;
    }

    struct RTCRayHit* rayhit = (RTCRayHit*)args->rayhit;
    unsigned int primID = args->primID;
    assert(args->N == 1);

    face_descriptor fd = self->id2desc(primID);

    Triangle face = triangle(self, fd);

    auto v = context->query_type==Intersect_context::Query_type::RAY_QUERY
      ? CGAL::intersection(context->ray, face)
      : CGAL::intersection(context->segment, face);
    if(v){
      rayhit->hit.geomID = self->rtc_geomID;
      rayhit->hit.primID = primID;

      const Point& source = (context->query_type==Intersect_context::Query_type::RAY_QUERY)?
        context->ray.source() : context->segment.source();

      float _distance;
      Point intersection_point;
      if (const Point *ip = boost::get<Point>(&*v) ){
        intersection_point = *ip;
        _distance = sqrt(to_double(CGAL::squared_distance(source, intersection_point)));
      }else{
        const Segment *intersection_segment = boost::get<Segment>(&*v);
        intersection_point = intersection_segment->source();
        _distance = sqrt(to_double(CGAL::squared_distance(source, intersection_point)));
        double distance_to_target = sqrt(to_double(CGAL::squared_distance(source, intersection_segment->target())));
        if(distance_to_target < _distance){
          _distance = distance_to_target;
          intersection_point = intersection_segment->target();
        }
      }

      if(context->intersection_type == Intersect_context::Intersection_type::COUNTER){
        context->counter++;
      } else if(context->intersection_type == Intersect_context::Intersection_type::FIRST){
        rayhit->ray.tfar = _distance;
        if(! context->intersections().empty()){
          context->intersections()[0] = std::make_tuple(intersection_point, _distance, primID);
        } else {
          context->intersections().push_back(std::make_tuple(intersection_point, _distance, primID));
        }
      } else if(context->intersection_type == Intersect_context::Intersection_type::ALL){
        context->intersections().push_back(std::make_tuple(intersection_point, _distance, primID));
      } else {
        rayhit->ray.tfar = _distance;
        context->intersections().push_back(std::make_tuple(intersection_point, _distance, primID));
        // Makes the ray invalid, so there is no further traversal
        rayhit->ray.tnear = rayhit->ray.tfar + 1.0f;
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
    double d = sqrt(to_double(squared_distance(result->query, cp)));
    if(d < args->query->radius){
      args->query->radius = d;
      result->result = cp;
      result->primID = primID;
      result->geomID = geomID;
      return true;
    }
    return false;
  }
///@}


/// \name Utility
///@{
  void insert_primitives()
  {
    rtcSetGeometryUserPrimitiveCount(rtc_geometry, num_faces(*triangle_mesh));
    rtcSetGeometryUserData(rtc_geometry, this);
    rtcSetGeometryBoundsFunction(rtc_geometry, bound_function, nullptr);
    rtcSetGeometryIntersectFunction(rtc_geometry, intersection_function);
    rtcCommitGeometry(rtc_geometry);
    rtcReleaseGeometry(rtc_geometry);
  }


  Primitive_id primitive_id(unsigned int primID) const
  {
    return std::make_pair(id2desc(primID), const_cast<TriangleMesh*>(triangle_mesh));
  }
///@}

};

/**
 * \ingroup PkgEmbreeRef
 * Data structure for efficient
 * intersection and distance computations in 3D. It uses Embree's user defined method to calculate faster
 * query on ray geometry intersections and distance computations
 * for 3D geometric objects. It uses another structure which contains the data of the mesh and commits the geometry to a scene of EmbreeAPI.
 *
 * \tparam Geometry Surface_mesh or Polyhedron_3
 * \tparam GeomTraits a kernel
 */
template <typename Geometry, typename GeomTraits>
class AABB_tree {

  typedef AABB_tree<Geometry,GeomTraits> Self;
public:
/// \name Types
///@{
  /// Identifier for a primitive in the tree.
  typedef typename Geometry::Primitive_id Primitive_id;
  /// Intersection Point and Primitive ID type.
  typedef std::pair<typename Geometry::Point, Primitive_id> Intersection_and_primitive_id;
  /// 3D Point and Primitive Id type
  typedef std::pair<typename Geometry::Point, Primitive_id> Point_and_primitive_id;
  /// Type of bounding box.
  typedef Bbox_3   Bounding_box;
  /// Unsigned integral size type.
  typedef std::size_t size_type;
  /// Type of 3D point.
  typedef typename GeomTraits::Point_3 Point;
  /// Type of ray query.
  typedef typename GeomTraits::Ray_3 Ray;
  /// Type of segment query.
  typedef typename GeomTraits::Segment_3 Segment;
///@}
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
  /// \name Creation
  ///@{

  /// constructs an empty tree, and initialises the underlying objects needed by the EmbreeAPI.
  /// Makes a call to Embree::AABB_tree::rtc_bind.
  AABB_tree()
  {
    rtc_bind();
  }

  /// constructs an empty tree, and initialises the underlying objects needed by the EmbreeAPI.
  /// user can supply a bool to make the underlying EmbreeAPI robust for calculations
  AABB_tree(bool robust)
    : AABB_tree()
  {
    if(robust)
      rtcSetSceneFlags(scene, RTC_SCENE_FLAG_ROBUST);
  }
  /// Destroys the Tree object
  /// makes a call to Embree::AABB_tree::rtc_unbind.
  ~AABB_tree()
  {
    rtc_unbind();
  }

  /// initializes the EmbreeAPI objects, namely the RTCDevice and the RTCScene.
  void rtc_bind()
  {
    device = rtcNewDevice(NULL);
    scene = rtcNewScene(device);
  }

  /// Can be called to free the EmbreeAPI objects, namely the Device and the Scene.
  /// the API can not be used further if this is called.
  void rtc_unbind()
  {
    rtcReleaseScene(scene);
    rtcReleaseDevice(device);
  }
  ///@}

  /// \name Operations
  ///@{

  /// returns \c true, iff the tree contains no primitive.
  bool empty() const
  {
    if (geometries.empty()) return true;
    return false;
  }


  /// clears the tree.
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

  /// T is the Surface_mesh/ Polyhedron_3
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

  ///@}

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

    typedef Intersect_context<Geometry> Intersect_context;
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

    typedef Intersect_context<Geometry> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::COUNTER, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return 0;
    }
    // const Geometry& geometry = geometries[rtc_geomID];
    return context.counter;
  }

  /// returns the primitive id of the face closest to the source point of the ray
  /// query.
  /// \tparam Query may be `GeomTraits::Ray_3` or `GeomTraits::Segment_3`.
  template<typename Query>
  boost::optional<Primitive_id> first_intersected_primitive(const Query& query) const
  {
    if (this->empty()) return boost::none;

    typedef Intersect_context<Geometry> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::FIRST, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return boost::none;
    }

    const Geometry& geometry = geometries[rtc_geomID];

    return boost::make_optional(geometry.primitive_id(context.rayhit.hit.primID));
  }

  /// puts in `out` the ids of all intersected primitives.
  /// This function does not compute the intersection points
  /// and is hence faster than the function `all_intersections()`
  /// function below.
  template<typename Query, typename OutputIterator>
  OutputIterator all_intersected_primitives (const Query& query, OutputIterator out) const
  {
    if (this->empty()) return out;

    typedef Intersect_context<Geometry> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::ALL, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return out;
    }

    const Geometry& geometry = geometries[rtc_geomID];

    for(int i=0; i<context.intersections().size();i++){
      *out++ = boost::make_optional(geometry.primitive_id(std::get<2>(context.intersections()[i])));
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

    typedef Intersect_context<Geometry> Intersect_context;
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

  /// returns the intersection and  primitive id closest to the source point of the ray
  /// query.
  /// \tparam Query may be `GeomTraits::Ray_3` or `GeomTraits::Segment_3`.
  template<typename Query>
  boost::optional<Intersection_and_primitive_id> first_intersection(const Query& query) const
  {
    if (this->empty()) return boost::none;

    typedef Intersect_context<Geometry> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::FIRST, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return boost::none;
    }

    typename Geometry::Point p = std::get<0>(context.intersections()[0]) ;

    const Geometry& geometry = geometries[rtc_geomID];
    return boost::make_optional(std::make_pair(p, geometry.primitive_id(context.rayhit.hit.primID)));
  }

  /// puts in `out` all intersections, as objects of
  /// `Intersection_and_primitive_id<Query>::%Type`,
  /// between the query and the input data to
  /// the iterator.
  /// \tparam Query may be `GeomTraits::Ray_3` or `GeomTraits::Segment_3`.
  template<typename Query, typename OutputIterator>
  OutputIterator all_intersections(const Query& query, OutputIterator out) const
  {
    if (this->empty()) return out;


    typedef Intersect_context<Geometry> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::ALL, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return out;
    }

    const Geometry& geometry = geometries[rtc_geomID];
    for(int i=0; i<context.intersections().size();i++){

      *out++ = boost::make_optional(std::make_pair(std::get<0>(context.intersections()[i]), geometry.primitive_id(std::get<2>(context.intersections()[i]))));
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

    typedef Intersect_context<Geometry> Intersect_context;
    Intersect_context context(Intersect_context::Intersection_type::ANY, query);

    rtcIntersect1(scene, &context, &(context.rayhit));

    unsigned int rtc_geomID = context.rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return boost::none;
    }

    typename Geometry::Point p = std::get<0>(context.intersections()[0]);

    const Geometry& geometry = geometries[rtc_geomID];
    return boost::make_optional(std::make_pair(p, geometry.primitive_id(context.rayhit.hit.primID)));
  }
 ///@}


  /// \name Distance Queries
  ///@{

  /// Returns the minimum squared distance between the query 
  /// point and all input primitives.
  /// \pre `!empty()`
  typename GeomTraits::FT squared_distance(const Point &query) const
  {
    return CGAL::squared_distance(query, this->closest_point(query));
  }

  /// returns the point in the union of all input primitives which
  /// is closest to the query. In case there are several closest
  /// points, one arbitrarily chosen closest point is
  /// returned.
  /// \pre `!empty()`
  Point closest_point(const Point &query) const
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
  Point_and_primitive_id closest_point_and_primitive(const Point &query) const
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
