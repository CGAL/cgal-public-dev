// Copyright (c) 2014  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s) : Andread Fabri

#ifndef CGAL_BOOST_GRAPH_INTERNAL_HELPERS_H
#define CGAL_BOOST_GRAPH_INTERNAL_HELPERS_H

#include <boost/graph/graph_traits.hpp>
#include <boost/optional.hpp>
#include <CGAL/boost/graph/iterator.h>

namespace CGAL {

// breaks a dependency loop between <CGAL/boost/graph/helpers.h>
// and <CGAL/boost/graph/iterator.h>
template <typename Graph> class Halfedge_around_target_iterator;

namespace internal {

template <typename Graph>
void
set_border(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
          , Graph& g)
{
  set_face(h, boost::graph_traits<Graph>::null_face(), g);
}

template <typename Graph>
typename boost::graph_traits<Graph>::halfedge_descriptor
copy(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
                    , Graph& g)
{
  typename boost::graph_traits<Graph>::edge_descriptor e = add_edge(g);
  typename boost::graph_traits<Graph>::halfedge_descriptor res = halfedge(e,g);
  typename boost::graph_traits<Graph>::halfedge_descriptor ropp = opposite(res, g);
  typename boost::graph_traits<Graph>::halfedge_descriptor hopp = opposite(h, g);
  set_target(res, target(h, g), g);
  set_target(hopp, target(hopp, g), g);
  set_face(res, face(h, g), g);
  set_face(ropp, face(hopp, g), g);
  // note that we cannot call set_next as it then would call set_prev on the  original
  return res;
 }


template <typename Graph>
void
set_vertex_halfedge(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
                    , Graph& g)
{ set_halfedge(target(h, g), h, g); }


template <typename Graph>
void
close_tip(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
          , Graph& g)
{
  // makes `opposite(h,g)' the successor of h.
  set_next( h, opposite(h, g), g);
}


template <typename Graph>
void
close_tip(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
          , typename boost::graph_traits<Graph>::vertex_descriptor const& v
          , Graph& g)
{
  // makes `h->opposite()' the successor of h and sets the incident
  // vertex of h to v.
  set_next(h, opposite(h, g), g);
  set_target(h, v, g);
  set_halfedge(v, h, g);
}

template <typename Graph>
void
insert_tip(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
           , typename boost::graph_traits<Graph>::halfedge_descriptor const& h2
           , Graph& g)
{
  set_next(h, next(h2,g), g);
  set_next(h2, opposite(h, g), g);
  set_target(h, target(h2, g), g);
}


template <typename Graph>
void
remove_tip(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
           , Graph& g)
{ 
  set_next(h, next(opposite(next(h, g), g), g), g);
}


template <typename Graph>
void 
set_face_in_face_loop(typename boost::graph_traits<Graph>::halfedge_descriptor h, 
                      typename boost::graph_traits<Graph>::face_descriptor f, 
                      Graph& g) 
{
  typename boost::graph_traits<Graph>::halfedge_descriptor end = h;
  do {
    set_face(h, f, g);
    h = next(h, g);
  } while ( h != end);
}
    

template <typename Graph>
void insert_halfedge(typename boost::graph_traits<Graph>::halfedge_descriptor const& h
                     , typename boost::graph_traits<Graph>::halfedge_descriptor const& f
                     , Graph& g)
{
  set_next(h, next(f, g), g);
  set_next(f, h, g);
  set_face(h, face(f, g), g);
}

template <typename Graph>
std::size_t
exact_num_vertices(const Graph& g)
{ 
  typename boost::graph_traits<Graph>::vertex_iterator beg, end;
  boost::tie(beg,end) = vertices(g);
  return std::distance(beg,end);
 }

template <typename Graph>
std::size_t
exact_num_halfedges(const Graph& g)
{ 
  typename boost::graph_traits<Graph>::halfedge_iterator beg, end;
  boost::tie(beg,end) = halfedges(g);
  return std::distance(beg,end);
 }

template <typename Graph>
std::size_t
exact_num_edges(const Graph& g)
{ 
  typename boost::graph_traits<Graph>::edge_iterator beg, end;
  boost::tie(beg,end) = edges(g);
  return std::distance(beg,end);
 }

template <typename Graph>
std::size_t
exact_num_faces(const Graph& g)
{ 
  typename boost::graph_traits<Graph>::face_iterator beg, end;
  boost::tie(beg,end) = faces(g);
  return std::distance(beg,end);
 }




} // internal
} // CGAL


#endif // CGAL_BOOST_GRAPH_INTERNAL_HELPERS_H
