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

#include <CGAL/license/Embree.h>

namespace CGAL {
namespae Embree {


  // AF: This is what you had called SM

template <typename TriangleMesh>
struct Triangle_mesh_geometry {
  typedef typename TriangleMesh:Face_index Face_index;
  typedef std::pair<Face_index,TriangleMesh*> Primitive_id;

  const TriangleMesh* surfaceMesh;
  RTCGeometry rtc_geometry;
  unsigned int rtc_geomID;

  Triangle_mesh_geometry(const TriangleMesh& tm)
    : surfaceMesh(&tm)
  {}

  static void bound_function(const struct RTCBoundsFunctionArguments* args)
  {
    // AF: move your code
  }

  static void intersection_function(const RTCIntersectFunctionNArguments* args)
  {
    // AF: move your code
  }

  void insert_primitives()
  {
    rtcSetGeometryUserPrimitiveCount(rtc_geometry, surfaceMesh.number_of_faces());
    rtcSetGeometryUserData(rtc_geometry, this);

    // AF: For the next two you have to find out how to write
    // the function pointer for a static member function
    rtcSetGeometryBoundsFunction(rtc_geometry, bound_function, nullptr);
    rtcSetGeometryIntersectFunction(rtc_geometry, intersection_function);
    rtcCommitGeometry(rtc_geometry);

    rtcReleaseGeometry(rtc_geometry);

    rtcCommitScene(rtc_scene);
  }

  Primitive_id primitive_id(unsigned int primID) const
  {
    return std::make_pair(Face_index(primID),surfaceMesh);
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

  RTCDevice device;
  RTCScene scene;

  // AF: As the ray intersection returns a geomID we need a mapping
  std::unordered_map<unsigned int, Geometry> id2geometry;
  std::list<Geometry> geometries;

  AABB_tree()
  {
    device = initializeDevice();
    scene = rtcNewScene(device);
  }

  /// T is the surface mesh
  template<typename T>
  void insert (cost T& t)
  {
    geometries.push_back(Geometry(t));
    const Geometry geometry = geometries.back();
    geometry.rtc_geometry = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
    geometry.rtc_geomID = rtcAttachGeometry(scene, geometry.rtc_geometry);
    geometry.insert_primitives();
  }


  template<typename Ray>
  boost::optional<Primitive_id>
  first_intersected_primitive(const Ray& query) const
  {
    // AF: no idea where the next two lines should go
    struct RTCIntersectContext context;
    rtcInitIntersectContext(&context);

    RTCRayHit rayhit;
    // AF initialize rayhit
    rtcIntersect1(scene, &context, &rayhit);

    unsigned int rtc_geomID = rayhit.hit.geomID;
    if(rtc_geomID == RTC_INVALID_GEOMETRY_ID){
      return boost::none;
    }
    Geometry geometry = id2geometry[rtc_geomID];

    return boost::make_optional(geometry.primitive(rayhit.hit.primID));
  }


};

} // namespace Embree
} // namespace CGAL
#endif // CGAL_EMBREE_AABB_TREE_H
