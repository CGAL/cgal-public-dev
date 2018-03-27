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

#ifndef CGAL_POLYLINE_TRACING_DICTIONARY_H
#define CGAL_POLYLINE_TRACING_DICTIONARY_H

#include <CGAL/Polyline_tracing/Dictionary_entry.h>

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
template<typename MotorcycleGraphTraits_>
class Dictionary
{
  typedef Dictionary<MotorcycleGraphTraits_>                                Self;

public:
  typedef MotorcycleGraphTraits_                                            Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                               Triangle_mesh;

  typedef typename Geom_traits::FT                                          FT;
  typedef typename Geom_traits::Point_d                                     Point;

  typedef Dictionary_entry<MotorcycleGraphTraits_>                          Dictionary_entry;
  typedef typename internal::DE_ptr_type<Dictionary_entry>::type            Dictionary_entry_ptr;
  typedef Dictionary_entry_base<MotorcycleGraphTraits_, Dictionary_entry>   Dictionary_entry_base;

  typedef boost::container::slist<Dictionary_entry_base>                    DEB_container;
  typedef typename DEB_container::iterator                                  DEBC_it;

  // @todo We want a stable container for 'Dictionary_entry_container' because of siblings, but
  // unordered set is _not_ stable. Switching to unordered here is doable by
  // using another underlying data structure (a list, not a vector (not stable)) and
  // making an unordered set of pointers to that structure
  typedef std::set<Dictionary_entry>                                        Dictionary_entry_container;
  typedef typename Dictionary_entry_container::iterator                     DEC_it;

  typedef typename Dictionary_entry::Face_location                          Face_location;

  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;
  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                                   descriptor_variant;

  // Access
  Dictionary_entry_container& entries() { return entries_; }
  const Dictionary_entry_container& entries() const { return entries_; }
  std::size_t size() const { return entries_.size(); }

  // Constructor
  Dictionary() : entries_bases_(), entries_() { }

  // Functions
  void erase(DEC_it pos, const bool erase_siblings = true);

  std::pair<DEC_it, bool> find(const Dictionary_entry& e) const;
  std::pair<DEC_it, bool> find(const Face_location& loc) const;

  Dictionary_entry_ptr get_sibling(DEC_it e, const face_descriptor fd) const;

  std::pair<DEC_it, bool> insert(const Face_location& loc, const Point& p, const Triangle_mesh& mesh);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const Point& p, const std::size_t i, const FT time, const Triangle_mesh& mesh);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const std::size_t i, const FT time, const Triangle_mesh& mesh);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const Triangle_mesh& mesh);

  // Ouput
  friend std::ostream& operator<<(std::ostream& out, const Self& d)
  {
    DEC_it it = d.entries().begin(), end = d.entries().end();
    for(; it!=end; ++it)
      out << &*it << " " << *it << std::endl;

    return out;
  }

private:
  DEB_container entries_bases_;
  Dictionary_entry_container entries_;
};

template<typename MotorcycleGraphTraits_>
void
Dictionary<MotorcycleGraphTraits_>::
erase(DEC_it pos, const bool erase_siblings /*= true*/)
{
  // WARNING: this function also erase (by default) the potential entry's siblings !!

  face_descriptor fd = pos->location().first;
  DEBC_it common_base = pos->base();
  CGAL_assertion(common_base != DEBC_it());

  if(erase_siblings)
  {
    // remove the other siblings, if any
    typename Dictionary_entry_base::Siblings_container::iterator smit = common_base->siblings().begin(),
                                                                 end = common_base->siblings().end();
    for(; smit!=end; ++smit)
    {
      // Location for that face
      face_descriptor adj_fd = (*smit)->location().first;
      CGAL_assertion(adj_fd != boost::graph_traits<Triangle_mesh>::null_face());

      if(adj_fd == fd)
        continue;

      entries_.erase(*smit);
    }
  }

  // erase the base
  entries_bases_.erase(common_base);

  // erase the main one
  entries_.erase(pos);
}

template<typename MotorcycleGraphTraits_>
std::pair<typename Dictionary<MotorcycleGraphTraits_>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits_>::
find(const Dictionary_entry& e) const
{
  DEC_it res = entries().find(e);
  return std::make_pair(res, (res != entries().end()));
}

template<typename MotorcycleGraphTraits_>
std::pair<typename Dictionary<MotorcycleGraphTraits_>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits_>::
find(const Face_location& loc) const
{
  Dictionary_entry dummy(loc);
  return find(dummy);
}

template<typename MotorcycleGraphTraits_>
typename Dictionary<MotorcycleGraphTraits_>::Dictionary_entry_ptr
Dictionary<MotorcycleGraphTraits_>::
get_sibling(const DEC_it e, const face_descriptor fd) const
{
  if(e->location().first == fd)
    return e;

  return e->sibling(fd);
}

template<typename MotorcycleGraphTraits_>
std::pair<typename Dictionary<MotorcycleGraphTraits_>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits_>::
insert(const Face_location& location, const Point& p, const Triangle_mesh& mesh)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  Dictionary_entry entry(location); // only the location is set up so far
  std::pair<DEC_it, bool> is_insert_successful = entries().insert(entry);

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
    Dictionary_entry_base deb(p);
    entries_bases_.push_front(deb);
    DEBC_it common_base = entries_bases_.begin();
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

          std::pair<DEC_it, bool> is_insert_successful_in_fd = entries().insert(location_in_fd);

          // Insertion must be successful: since we insert all the possible locations
          // when we insert a location on a border, we can't have had successful
          // insertion for the main entry but not for equivalent locations.
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
          std::pair<DEC_it, bool> is_insert_successful_in_fd = entries().insert(location_in_fd);
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

template<typename MotorcycleGraphTraits_>
std::pair<typename Dictionary<MotorcycleGraphTraits_>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits_>::
insert(const Face_location& loc, const Point& p, const std::size_t i, const FT time, const Triangle_mesh& mesh)
{
  std::pair<DEC_it, bool> entry = insert(loc, p, mesh);
  entry.first->add_motorcycle(i, time);

  return entry;
}

template<typename MotorcycleGraphTraits>
std::pair<typename Dictionary<MotorcycleGraphTraits>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits>::
insert(const Face_location& loc, const std::size_t i, const FT time, const Triangle_mesh& mesh)
{
  Point p = CGAL::Polygon_mesh_processing::location_to_point(loc, mesh);
  return insert(loc, p, i, time, mesh);
}

template<typename MotorcycleGraphTraits>
std::pair<typename Dictionary<MotorcycleGraphTraits>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits>::
insert(const Face_location& loc, const Triangle_mesh& mesh)
{
  Point p = CGAL::Polygon_mesh_processing::location_to_point(loc, mesh);
  return insert(loc, p, mesh);
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_DICTIONARY_H
