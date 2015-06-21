// Copyright (c) 2009-2013 GeometryFactory (France).
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
//
//
// Author(s)     : Laurent Rineau and Sebastien Loriot


#ifndef CGAL_ORIENT_POLYGON_SOUP
#define CGAL_ORIENT_POLYGON_SOUP

#include <boost/foreach.hpp>

#include <set>
#include <map>
#include <stack>
#include <algorithm>
#include <iostream>

namespace CGAL {

namespace internal {

template<class Point_3, class Polygon_3>
class Polygon_soup_orienter
{
/// Index types
  typedef typename std::iterator_traits<
            typename Polygon_3::iterator >::value_type                     V_ID;
  typedef typename std::vector<Polygon_3>::size_type                       P_ID;
//  typedef int                                                             CC_ID;
  typedef std::pair<V_ID, V_ID>                                       V_ID_pair;
/// Container types
  typedef std::vector<Point_3>                                           Points;
  typedef std::vector<Polygon_3>                                       Polygons;
  typedef std::map<V_ID_pair, std::set<P_ID> >                         Edge_map;
  typedef typename Edge_map::iterator                         Edge_map_iterator;
  typedef std::set<V_ID_pair>                                      Marked_edges;

/// Data members
  Points& points;             //< the set of input points
  Polygons& polygons;         //< the set of input polygons
  Edge_map edges;             //< the set of edges of the input polygons
  Marked_edges marked_edges;  //< the set of singular edges or edges incident
                              //<   to non-compatible orientation polygons

  /// for each polygon referenced by its position in `polygons`, indicates
  /// the connected component it belongs too after orientation.
//  std::vector< CC_ID > polygon_cc_id;
  /// for each vertex, indicates the list of polygon containing it
  std::vector< std::vector<P_ID> > incident_polygons_per_vertex;

/// Utility functions
  V_ID_pair canonical_edge(V_ID i, V_ID j)
  {
    return i<j ? V_ID_pair(i,j):V_ID_pair(j,i);
  }

  bool is_edge_marked(V_ID i, V_ID j)
  {
    return marked_edges.count(canonical_edge(i,j));
  }

  void set_edge_marked(V_ID i, V_ID j)
  {
    marked_edges.insert(canonical_edge(i,j));
  }

  cpp11::array<V_ID,3>
  get_neighbor_vertices(V_ID v_id, P_ID polygon_index)
  {
    std::size_t nbv = polygons[polygon_index].size(), pvid=0;
    for (; pvid!=nbv; ++pvid)
      if (v_id==polygons[polygon_index][pvid]) break;
    CGAL_assertion( pvid!=nbv );
    V_ID prev = polygons[polygon_index][ (pvid+nbv-1)%nbv ];
    V_ID next = polygons[polygon_index][ (pvid+1)%nbv ];
    return make_array(prev,v_id,next);
  }

  std::pair<V_ID,P_ID>
  next_cw_vertex_around_source(V_ID src, V_ID tgt)
  {
    typedef std::pair<V_ID,P_ID> VID_and_PID;
    if ( is_edge_marked(src,tgt) ) return VID_and_PID(src,300612);
    Edge_map_iterator em_it=edges.find(V_ID_pair(tgt, src));
    if ( em_it==edges.end() ) return VID_and_PID(src,300612);// the vertex is on the border
    CGAL_assertion(em_it->second.size()==1);
    P_ID p_id = *(em_it->second.begin());
    return VID_and_PID(get_neighbor_vertices(src, p_id)[2], p_id);
  }

  std::pair<V_ID,P_ID>
  next_ccw_vertex_around_target(V_ID src, V_ID tgt)
  {
    typedef std::pair<V_ID,P_ID> VID_and_PID;
    if ( is_edge_marked(src,tgt) ) return VID_and_PID(tgt,300612);
    Edge_map_iterator em_it=edges.find(V_ID_pair(tgt, src));
    if ( em_it==edges.end() ) return VID_and_PID(tgt,300612);// the vertex is on the border
    CGAL_assertion(em_it->second.size()==1);
    P_ID p_id = *(em_it->second.begin());
    return VID_and_PID(get_neighbor_vertices(tgt, p_id)[0], p_id);
  }

  void inverse_orientation(const std::size_t index) {
    std::reverse(polygons[index].begin(), polygons[index].end());
  }

  void replace_vertex_index_in_polygon(
    std::size_t polygon_id,
    V_ID old_index,
    V_ID new_index)
  {
    BOOST_FOREACH(V_ID& i, polygons[polygon_id])
      if( i==old_index )
        i=new_index;
  }

/// Functions filling containers
  void fill_edge_map() {
    // Fill edges
    edges.clear();
    for(P_ID i = 0; i < polygons.size(); ++i)
    {
      const P_ID size = polygons[i].size();
      for(P_ID j = 0; j < size; ++j) {
        V_ID i0 = polygons[i][j];
        V_ID i1 = polygons[i][ (j+1)%size];
        edges[V_ID_pair(i0, i1)].insert(i);
      }
    }

    // Fill non-manifold edges
    marked_edges.clear();
    for(P_ID i = 0; i < polygons.size(); ++i)
    {
      const P_ID size = polygons[i].size();
      for(P_ID j = 0; j < size; ++j) {
        V_ID i0 = polygons[i][j];
        V_ID i1 = polygons[i][ (j+1)%size ];

        std::size_t nb_edges = 0;
        Edge_map_iterator em_it = edges.find( V_ID_pair(i0, i1) );
        if ( em_it!=edges.end() ) nb_edges += em_it->second.size();
        em_it = edges.find( V_ID_pair(i1, i0) );
        if ( em_it!=edges.end() ) nb_edges += em_it->second.size();

        if( nb_edges > 2 ) set_edge_marked(i0,i1);
      }
    }
  }

  void fill_incident_polygons_per_vertex()
  {
    incident_polygons_per_vertex.resize(points.size());

    P_ID nb_polygons=polygons.size();
    for(P_ID ip=0; ip<nb_polygons; ++ip)
    {
      BOOST_FOREACH(V_ID iv, polygons[ip])
        incident_polygons_per_vertex[iv].push_back(ip);
    }
  }

public:

  Polygon_soup_orienter(Points& points, Polygons& polygons)
    : points(points), polygons(polygons)
  {
    fill_edge_map();
  }

  /// We try to orient polygon consistently by walking in the dual graph, from
  /// a not yet re-oriented polygon.
  /// We have an edge between two polygons if they share an edge, and this edge
  /// is shared by exactly two polygons. While walking along an edge, we reorient
  /// the polygon we walked in if its orientation is not compatible with the one
  /// we come from.
  /// If the polygon was already marked as oriented, then we cut the dual edge
  /// in the graph and the primal edge is marked.
  /// At the same time, we assign an id to each polygon in the same connected
  /// componenet of the dual graph.
  void orient()
  {
    std::vector<bool> oriented;
    std::stack<std::size_t> stack;
//    polygon_cc_id.resize(polygons.size(), -1);

    // We first consider all polygons as non-oriented
    oriented.resize(polygons.size());

    P_ID polygon_index = 0;

//    CC_ID current_cc_index=-1;
    while (polygon_index != polygons.size())
    {
      // We look for the first polygon not already oriented
      while ( polygon_index != polygons.size() && oriented[polygon_index] ) {
        ++polygon_index;
      }
      if(polygon_index == polygons.size()) break;

//      ++ current_cc_index; // visit a new connected component

      // we visit the connected component by crossing edges manifold edges
      oriented[polygon_index] = true;
      stack.push(polygon_index);
      while(! stack.empty() )
      {
        const P_ID to_be_oriented_index = stack.top();
        stack.pop();

//        CGAL_assertion(polygon_cc_id[to_be_oriented_index]==-1);
//        polygon_cc_id[to_be_oriented_index]=current_cc_index;

        const P_ID size = polygons[to_be_oriented_index].size();
        for(P_ID ih = 0 ; ih < size ; ++ih) {
          P_ID ihp1 = (ih+1)%size;
          const V_ID i1 = polygons[to_be_oriented_index][ih];
          const V_ID i2 = polygons[to_be_oriented_index][ihp1];

          if( is_edge_marked(i1,i2) ) continue;

          // edge (i1,i2)
          Edge_map_iterator it_same_orient = edges.find(V_ID_pair(i1, i2));
          // edges (i2,i1)
          Edge_map_iterator it_other_orient = edges.find(V_ID_pair(i2, i1));

          CGAL_assertion(it_same_orient != edges.end());
          CGAL_assertion(it_other_orient == edges.end() ||
                         it_other_orient->second.size()==1);

          if (it_same_orient->second.size() > 1)
          {
            CGAL_assertion(it_other_orient == edges.end());
            // one neighbor but with the same orientation
            P_ID index = *(it_same_orient->second.begin());
            if(index == to_be_oriented_index)
              index = *(++it_same_orient->second.begin());
            if(oriented[index])
            {
              // polygon already oriented but its orientation is not compatible ---> mark the edge and continue
              set_edge_marked(i1,i2);
              continue; // next edge
            }

            // reverse the orientation
            const P_ID size = polygons[index].size();
            for(P_ID j = 0; j < size; ++j) {
              V_ID i0 = polygons[index][j];
              V_ID i1 = polygons[index][(j+1)%size];
              Edge_map_iterator em_it = edges.find(V_ID_pair(i0, i1));
              CGAL_assertion_code(const bool r = )
                em_it->second.erase(index)
              CGAL_assertion_code(!= 0);
              CGAL_assertion(r);
              if ( em_it->second.empty() ) edges.erase(em_it);
            }
            inverse_orientation(index);
            for(P_ID j = 0; j < size; ++j) {
              V_ID i0 = polygons[index][j];
              V_ID i1 = polygons[index][(j+1)%size];
              edges[V_ID_pair(i0, i1)].insert(index);
            }
            // "inverse the orientation of polygon #index
            oriented[index] = true;
            stack.push(index);
          }
          else{
            if( it_other_orient != edges.end() ){
              CGAL_assertion(it_same_orient->second.size() == 1);
              CGAL_assertion(it_other_orient->second.size() == 1);
              // one polygon, same orientation
              const P_ID index = *(it_other_orient->second.begin());
              if(oriented[index]) continue; //nothing todo already processed and correctly oriented
              oriented[index] = true;
              // "keep the orientation of polygon #index
              stack.push(index);
            }
          }
        } // end for on all edges of one
      } // end while loop on the polygons of the connected component
    } // end while loop on all non-oriented polygons remaining
  }

  /// A vertex is said to be singular if its link is neither a cycle nor a chain,
  /// but several cycles and chains.
  /// For each such vertex v, we consider each set of polygons incident to v
  /// and sharing a non-marked edge incident to v. A copy of v is assigned to
  /// each but one set of incident polygons.
  void duplicate_singular_vertices()
  {
    fill_incident_polygons_per_vertex();
    std::vector< std::pair<V_ID, std::vector<P_ID> > > vertices_to_duplicate;

    V_ID nbv = static_cast<V_ID>( points.size() );
    for (V_ID v_id = 0; v_id < nbv; ++v_id)
    {
      const std::vector< P_ID >& incident_polygons = incident_polygons_per_vertex[v_id];

      if ( incident_polygons.empty() ) continue; //isolated vertex
      std::set<P_ID> visited_polygons;

      bool first_pass = true;
      BOOST_FOREACH(P_ID p_id, incident_polygons)
      {
        if ( !visited_polygons.insert(p_id).second ) continue; // already visited

        if (!first_pass)
        {
          vertices_to_duplicate.push_back(std::pair<V_ID, std::vector<P_ID> >());
          vertices_to_duplicate.back().first=v_id;
        }

        const cpp11::array<V_ID,3>& neighbors = get_neighbor_vertices(v_id,p_id);

        V_ID next = neighbors[2];

        if( !first_pass)
          vertices_to_duplicate.back().second.push_back(p_id);

        do{
          P_ID other_p_id;
          cpp11::tie(next, other_p_id) = next_cw_vertex_around_source(v_id, next);
          if (next==v_id) break;
          visited_polygons.insert(other_p_id);
          if( !first_pass)
            vertices_to_duplicate.back().second.push_back(other_p_id);
        }
        while(next!=neighbors[0]);

        if (next==v_id){
          /// turn the otherway round
          next = neighbors[0];
          do{
            P_ID other_p_id;
            cpp11::tie(next, other_p_id) = next_ccw_vertex_around_target(next, v_id);
            if (next==v_id) break;
            visited_polygons.insert(other_p_id);
            if( !first_pass)
              vertices_to_duplicate.back().second.push_back(other_p_id);
          }
          while(true);
        }
        first_pass=false;
      }
    }

    /// now duplicate the vertices
    typedef std::pair<V_ID, std::vector<P_ID> > V_ID_and_Polygon_ids;
    BOOST_FOREACH(const V_ID_and_Polygon_ids& vid_and_pids, vertices_to_duplicate)
    {
      V_ID new_index = static_cast<V_ID>(points.size());
      points.push_back( points[vid_and_pids.first] );
      BOOST_FOREACH(P_ID polygon_id, vid_and_pids.second)
        replace_vertex_index_in_polygon(polygon_id, vid_and_pids.first, new_index);
    }
  }
};
} // namespace internal

/**
 * Tries to consistently orient a soup of polygons in 3D space.
 * If it is not possible to produce a combinatorial manifold surface, some points will be
 * duplicated. These points are either an endpoint of edges incident to more than
 * two polygons, or an endpoint of an edge between two polygons with incompatible orientations
 * (during the re-orientation process), or a point shared by at least two polygons that do not
 * share an edge this point is incident to.
 * @tparam Point_3 the point type
 * @tparam Polygon_3 the Polygon type, being a container of indices
 *
 * @param[in,out] points points of the soup of polygons. Some points might be pushed back to resolve
 *                non-manifold or non-orientability issues.
 * @param[in, out] polygons each element in the vector describes a polygon using the index of the points in `points`.
 *
 * @return return false if some points where duplicated, thus producing a self-intersecting polyhedron
 *
 */
template <class Point_3, class Polygon_3>
bool orient_polygon_soup(std::vector<Point_3>& points,
                         std::vector< Polygon_3 >& polygons)
{
  std::size_t inital_nb_pts = points.size();
  internal::Polygon_soup_orienter<Point_3, Polygon_3> orienter(points, polygons);
  orienter.orient();
  orienter.duplicate_singular_vertices();

  return inital_nb_pts==points.size();
}

}// namespace CGAL

#endif // CGAL_ORIENT_POLYGON_SOUP
