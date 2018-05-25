// Copyright (c) 2017, 2018 GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_NODE_DICTIONARY_H
#define CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_NODE_DICTIONARY_H

#include <CGAL/Polyline_tracing/Motorcycle_graph_node.h>

#include <CGAL/assertions.h>
#include <CGAL/functional.h>
#include <CGAL/Polygon_mesh_processing/locate.h>

#include <boost/bimap.hpp>
#include <boost/container/slist.hpp>
#include <boost/unordered_set.hpp>

#include <algorithm>
#include <iostream>
#include <list>
#include <ostream>
#include <utility>

namespace CGAL {

namespace Polyline_tracing {

// -----------------------------------------------------------------------------
//                           Dictionary class
// The points involved in the motorcycle graph.
// -----------------------------------------------------------------------------
template<typename MotorcycleGraphTraits>
class Motorcycle_graph_node_dictionary
{
  typedef Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>           Self;

public:
  typedef MotorcycleGraphTraits                                             Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                               Triangle_mesh;

  typedef typename Geom_traits::FT                                          FT;
  typedef typename Geom_traits::Point_d                                     Point;

  typedef Motorcycle_graph_node<Geom_traits>                                Node;
  typedef Motorcycle_graph_node_base<Geom_traits, Node>                     Node_base;

  typedef boost::container::slist<Node_base>                                NB_container;
  typedef typename NB_container::iterator                                   NBC_it;

  // @todo We need a stable container for 'Node_container' because of siblings, but
  // unordered set is _not_ stable. Switching to unordered here is doable by
  // using another underlying data structure (a list, not a vector (not stable)) and
  // making an unordered set of pointers to that structure
  typedef std::set<Node>                                                    Node_container;
  typedef typename internal::Node_ptr_type<Node>::type                      Node_ptr;

  typedef typename Node_container::iterator                                 iterator;
  typedef typename Node_container::const_iterator                           const_iterator;

  typedef typename Node::Face_location                                      Face_location;
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;
  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                                   descriptor_variant;

  // Access
  Node_container& container() { return nodes_; }
  const Node_container& container() const { return nodes_; }

  inline iterator begin() { return nodes_.begin(); }
  inline const_iterator begin() const { return nodes_.begin(); }
  inline iterator end() { return nodes_.end(); }
  inline const_iterator end() const { return nodes_.end(); }
  inline std::size_t size() const { return nodes_.size(); }

  // Constructor
  Motorcycle_graph_node_dictionary() : node_bases_(), nodes_() { }

  // Functions
  void erase(Node_ptr pos, const bool erase_siblings = true);

  std::pair<Node_ptr, bool> find(const Face_location& loc) const;
  std::pair<Node_ptr, bool> find(const Node& e) const;

  Node_ptr get_sibling(Node_ptr e, const face_descriptor fd) const;

  std::pair<Node_ptr, bool> insert(const Face_location& loc, const Point& p, const Triangle_mesh& mesh);
  std::pair<Node_ptr, bool> insert(const Face_location& loc, const Point& p, const std::size_t i, const FT time, const Triangle_mesh& mesh);
  std::pair<Node_ptr, bool> insert(const Face_location& loc, const std::size_t i, const FT time, const Triangle_mesh& mesh);
  std::pair<Node_ptr, bool> insert(const Face_location& loc, const Triangle_mesh& mesh);

  // Ouput
  friend std::ostream& operator<<(std::ostream& out, const Self& d)
  {
    Node_ptr it = d.container().begin(), end = d.container().end();
    for(; it!=end; ++it)
      out << &*it << " " << *it << std::endl;

    return out;
  }

private:
  NB_container node_bases_;
  Node_container nodes_;
};

template<typename MotorcycleGraphTraits>
void
Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::
erase(Node_ptr pos, const bool erase_siblings /*= true*/)
{
  // WARNING: this function also erase (by default) the potential node's siblings !!

  face_descriptor fd = pos->face();
  NBC_it common_base = pos->base();
  CGAL_assertion(common_base != NBC_it());

  if(erase_siblings)
  {
    // remove the other siblings, if any
    typename Node_base::Siblings_container::iterator smit = common_base->siblings().begin(),
                                                     end = common_base->siblings().end();
    for(; smit!=end; ++smit)
    {
      // Location for that face
      face_descriptor adj_fd = (*smit)->face();
      CGAL_assertion(adj_fd != boost::graph_traits<Triangle_mesh>::null_face());

      if(adj_fd == fd)
        continue;

      nodes_.erase(*smit);
    }
  }

  // erase the base
  node_bases_.erase(common_base);

  // erase the main one
  nodes_.erase(pos);
}

template<typename MotorcycleGraphTraits>
std::pair<typename Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::Node_ptr, bool>
Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::
find(const Node& e) const
{
  Node_ptr res = container().find(e);
  return std::make_pair(res, (res != end()));
}

template<typename MotorcycleGraphTraits>
std::pair<typename Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::Node_ptr, bool>
Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::
find(const Face_location& loc) const
{
  Node dummy(loc);
  return find(dummy);
}

template<typename MotorcycleGraphTraits>
typename Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::Node_ptr
Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::
get_sibling(const Node_ptr e, const face_descriptor fd) const
{
  if(e->face() == fd)
    return e;

  return e->sibling(fd);
}

template<typename MotorcycleGraphTraits>
std::pair<typename Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::Node_ptr, bool>
Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::
insert(const Face_location& location, const Point& p, const Triangle_mesh& mesh)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  std::pair<Node_ptr, bool> is_insert_successful = container().emplace(location);

  if(!is_insert_successful.second)
  {
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cerr << "Warning: point already exists in the dictionary: "
              << &*(is_insert_successful.first) << std::endl
              << *(is_insert_successful.first) << std::endl;
#endif
  }
  else // Insertion is successful
  {
    // Initialize the base with all motorcycle information
    Node_base deb(p);
    node_bases_.push_front(deb);
    NBC_it common_base = node_bases_.begin();
    is_insert_successful.first->set_base(common_base);

    // If the point is on border: initialize and insert its siblings
    descriptor_variant dv = PMP::get_descriptor_from_location(location, mesh);
    face_descriptor initial_fd = location.first;

    CGAL_precondition(common_base->siblings().empty());

    if(dv.which() != 2) // point is not strictly on a face -> it's on the border
    {
      // fill for the initial face
      common_base->siblings().insert(is_insert_successful.first);

      if(dv.which() == 0) // on a vertex
      {
        vertex_descriptor vd = boost::get<vertex_descriptor>(dv);
        halfedge_descriptor hd = halfedge(vd, mesh);

        CGAL::Face_around_target_iterator<Triangle_mesh> fit, fend;
        boost::tie(fit, fend) = CGAL::faces_around_target(hd, mesh);
        std::size_t number_of_incident_faces = std::distance(fit, fend);
        for(; fit!=fend; ++fit)
        {
          face_descriptor fd = *fit;
          if(fd == initial_fd || fd == boost::graph_traits<Triangle_mesh>::null_face())
            continue;

          const Face_location location_in_fd = PMP::locate_in_adjacent_face(location, fd, mesh);

          std::pair<Node_ptr, bool> is_insert_successful_in_fd = container().emplace(location_in_fd);

          // Insertion must be successful: since we insert all the possible locations
          // when we insert a location on a border, we can't have had successful
          // insertion for the main node but not for its equivalent locations.
          CGAL_assertion(is_insert_successful_in_fd.second);

          // Set up the sibling information
          common_base->siblings().insert(is_insert_successful_in_fd.first);

          // Set the base
          is_insert_successful_in_fd.first->set_base(common_base);
        }
        CGAL_postcondition_code(if(PMP::is_on_mesh_border(location, mesh)))
        CGAL_postcondition(common_base->siblings().size() == number_of_incident_faces - 1 /*the null face*/);
        CGAL_postcondition_code(else)
        CGAL_postcondition(common_base->siblings().size() == number_of_incident_faces);
      }
      else if(dv.which() == 1) // on a halfedge
      {
        halfedge_descriptor hd = boost::get<halfedge_descriptor>(dv);
        if(!is_border(edge(hd, mesh), mesh))
        {
          face_descriptor fd = face(opposite(hd, mesh), mesh);
          const Face_location location_in_fd = PMP::locate_in_adjacent_face(location, fd, mesh);
          std::pair<Node_ptr, bool> is_insert_successful_in_fd = container().emplace(location_in_fd);
          CGAL_assertion(is_insert_successful_in_fd.second);

          common_base->siblings().insert(is_insert_successful_in_fd.first);
          CGAL_postcondition(common_base->siblings().size() == 2);

          is_insert_successful_in_fd.first->set_base(common_base);
        }
      }
    }

#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
    std::cout << "New point in the dictionary: "
              << &*(is_insert_successful.first) << std::endl
              << *(is_insert_successful.first) << std::endl;
#endif
  }

  return is_insert_successful;
}

template<typename MotorcycleGraphTraits>
std::pair<typename Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::Node_ptr, bool>
Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::
insert(const Face_location& loc, const Point& p, const std::size_t i, const FT time, const Triangle_mesh& mesh)
{
  std::pair<Node_ptr, bool> node = insert(loc, p, mesh);
  node.first->add_motorcycle(i, time);

  return node;
}

template<typename MotorcycleGraphTraits>
std::pair<typename Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::Node_ptr, bool>
Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::
insert(const Face_location& loc, const std::size_t i, const FT time, const Triangle_mesh& mesh)
{
  Point p = CGAL::Polygon_mesh_processing::location_to_point(loc, mesh);
  return insert(loc, p, i, time, mesh);
}

template<typename MotorcycleGraphTraits>
std::pair<typename Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::Node_ptr, bool>
Motorcycle_graph_node_dictionary<MotorcycleGraphTraits>::
insert(const Face_location& loc, const Triangle_mesh& mesh)
{
  Point p = CGAL::Polygon_mesh_processing::location_to_point(loc, mesh);
  return insert(loc, p, mesh);
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_NODE_DICTIONARY_H
