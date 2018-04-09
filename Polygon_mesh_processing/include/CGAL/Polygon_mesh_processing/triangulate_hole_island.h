// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Konstantinos Katrioplas

#ifndef CGAL_PMP_TRIANGULATE_HOLE_ISLAND_H
#define CGAL_PMP_TRIANGULATE_HOLE_ISLAND_H


#include <CGAL/Polygon_mesh_processing/internal/Hole_filling/island_triangulate_hole_polyline.h>
#include <CGAL/Polygon_mesh_processing/internal/Hole_filling/Triangulate_hole_polyline.h>


namespace CGAL {
namespace Polygon_mesh_processing{


/*!
\ingroup hole_filling_grp
triangulates a hole defined by points in a boundary and points in islands.
The area between islands is filled with triangles and the result
is saved as a polygon mesh in `mesh`.
The triangulated region does not contain any non-manifold edges or degenerate triangles.

@tparam PolygonMesh a model of `MutableFaceGraph`.
@tparam PointRange range of points, model of `Range`.

@param boundary the range of input points.
@param islands a std::vector of range of points, each containing points on an island.
@param mesh the resulting polygon mesh.
@param use_DT boolean to determine whether to use Delaunay triangulation as an optimization step.
@param correct_orientation boolean to determine whether the sequence of all islands are oriented
                           in the opposite direction to that of the boundary.

@return `mesh`

@todo use a polygon mesh as input to a wrapper around this function.
@todo the wrapper function should identify a boundary and islands.
*/
template <typename PointRange, typename PolygonMesh>
void triangulate_hole_islands(const PointRange& boundary,
                              const std::vector<PointRange>& islands,
                              PolygonMesh& mesh,
                              const bool& use_DT,
                              const bool& correct_orientation)
{
  // create domain from the boundary
  std::size_t boundary_size = boundary.size();
  if( *boost::begin(boundary) == *cpp11::prev(boost::end(boundary)) )
    --boundary_size;

  // indices
  std::vector<int> b_indices(boundary_size);
  for(std::size_t i = 0; i < boundary_size; ++i)
    b_indices[i] = i;

  CGAL::internal::Domain domain(b_indices);

  std::size_t number_of_vertices_on_islands = 0;
  // add islands if there are
  if(!islands.empty())
  {
    std::size_t ind = b_indices.size(); // start
    for(auto island : islands)
    {
      std::vector<int> h_ids;
      std::size_t island_size = island.size();
      // dont use the last one
      if( *boost::begin(island) == *cpp11::prev(boost::end(island)) )
        --island_size;

      number_of_vertices_on_islands += island_size;

      h_ids.reserve(island_size);
      // generate their ids
      for(std::size_t j = 0; j < island_size; ++j)
      {
        h_ids.push_back(ind);
        ++ind;
      }
      // push it to the domain
      domain.add_island(h_ids);
    }
  }

  std::cout << "Number of boundary vertices: " << b_indices.size() << std::endl;
  std::cout << "Number of vertices on islands: " << number_of_vertices_on_islands << std::endl;

  // assign access edge on the boundary
  int i =  0;
  int k = static_cast<int>(b_indices.size()) - 1;
  //std::size_t count = 0;

  // weight calculator
  typedef CGAL::internal::Weight_min_max_dihedral_and_area      Weight;
  typedef CGAL::internal::Weight_calculator<Weight,
                CGAL::internal::Is_not_degenerate_triangle>  WC;

  // lookup tables
  //typedef CGAL::internal::Lookup_table_map<Weight> WeightTable;
  //WeightTable W(n, Weight::DEFAULT());
  typedef CGAL::internal::Lookup_table_map<int> LambdaTable;
  int n = static_cast<int>(b_indices.size() + number_of_vertices_on_islands);
  LambdaTable lambda(n, -1);

  // put together points list
  PointRange points;
  points.reserve(boundary_size + number_of_vertices_on_islands);
  points.insert(points.end(), boundary.begin(), boundary.end() - 1); // boundary has not been reduced
  if(!islands.empty())
  {
    for(auto island : islands)
      points.insert(points.end(), island.begin(), island.end() - 1);
  }

  // output triangulation
  std::vector<std::vector<int> > triplets;

  if(use_DT)
  {
    CGAL::internal::Triangulate_hole_with_islands<PointRange, LambdaTable, WC>
        triangulation(domain, points, lambda, WC(), n, true, correct_orientation);
    triangulation.do_triangulation(i, k, triplets);
    triangulation.test_corectness(triplets);
    triangulation.visualize(points, triplets, mesh);
  }
  else
  {
    CGAL::internal::Triangulate_hole_with_islands<PointRange, LambdaTable, WC>
        triangulation(domain, points, lambda, WC(), n, false, correct_orientation);
    triangulation.do_triangulation(i, k, triplets);
    triangulation.visualize(points, triplets, mesh);
  }

}


} //namespace Polygon_mesh_processing
} //namespace CGAL

#endif // CGAL_PMP_TRIANGULATE_HOLE_ISLAND_H
