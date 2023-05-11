// Copyright (c) 2020 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_PMP_REMESHING_SIZING_FIELD_H
#define CGAL_PMP_REMESHING_SIZING_FIELD_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Kernel_traits.h>
#include <boost/optional.hpp>

namespace CGAL
{
/*!
* Sizing field virtual class
*/
template <class PolygonMesh>
class Sizing_field
{
private:
  typedef PolygonMesh                     PM;
  typedef typename boost::property_map<PM, boost::vertex_point_t>::const_type VPMap;
  typedef typename boost::property_traits<VPMap>::value_type Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel K;

public:
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PM>::vertex_descriptor   vertex_descriptor;
  typedef Point                                                 Point_3;
  typedef typename K::FT                                        FT;

public:
  virtual boost::optional<FT> is_too_long(const halfedge_descriptor& h) const = 0;
  virtual boost::optional<FT> is_too_long(const vertex_descriptor& va,
                                          const vertex_descriptor& vb) const = 0;
  virtual boost::optional<FT> is_too_short(const halfedge_descriptor& h) const = 0;
  virtual Point_3 split_placement(const halfedge_descriptor& h) const = 0;

};

}//end namespace CGAL

#endif //CGAL_PMP_REMESHING_SIZING_FIELD_H
