// Copyright (c) 2017 GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYLINE_TRACING_DICTIONARY_H
#define CGAL_POLYLINE_TRACING_DICTIONARY_H

#include <CGAL/array.h>
#include <CGAL/assertions.h>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/unordered_set.hpp>
#include <boost/variant.hpp>

#include <limits>
#include <map>
#include <utility>

namespace CGAL {

namespace Polyline_tracing {

// @todo put all that stuff in traits
namespace internal {

// Given a location, returns whether a point is on the border or not.
template<typename Face_location, typename PolygonMesh>
bool is_on_face_border(const Face_location& loc, const PolygonMesh& mesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor     face_descriptor;

  typedef typename Face_location::second_type                            Barycentric_coordinates;

  const face_descriptor fd = loc.first;
  const Barycentric_coordinates& baryc = loc.second;

  // the first barycentric coordinate corresponds to source(halfedge(fd, mesh), mesh)
  halfedge_descriptor hd = prev(halfedge(fd, mesh), mesh);

  // check if the point is a vertex
  for(int i=0; i<3; ++i)
  {
    if(baryc[i] == 1.) // coordinate at target(hd, mesh)
      return CGAL::is_border(target(hd, mesh), mesh);
    hd = next(hd, mesh);
  }
  CGAL_assertion(hd == prev(halfedge(fd, mesh), mesh));

  // check if the point is on an edge
  for(int i=0; i<3; ++i)
  {
    if(baryc[i] == 0) // coordinate at target(hd, mesh)
      return CGAL::is_border(prev(hd, mesh), mesh);
    hd = next(hd, mesh);
  }

  // point is strictly within the face, so it's not on the border
  return false;
}

// Given a Face_location (face + barycentric coordinates) and an adjacent face,
// compute the barycentric coordinates of the point in the other face.
template<typename Face_location, typename PolygonMesh>
Face_location compute_barycentric_coordinates(const Face_location& loc,
                                              const typename boost::graph_traits<PolygonMesh>::face_descriptor new_face,
                                              const PolygonMesh& mesh)
{
  CGAL_precondition(internal::is_on_face_border(loc, mesh));

  CGAL_assertion(false);
  return Face_location(); // @todo
}

// Compute the barycentric coordinates of 'p' that belongs in the face 'fd'
template<typename Point, typename Face_location, typename PolygonMesh>
Face_location compute_location(const Point& p,
                               const typename boost::graph_traits<PolygonMesh>::face_descriptor fd,
                               const PolygonMesh& mesh)
{
  CGAL_assertion(false);
  return Face_location(); // @todo
}

template<typename Point, typename Face_location, typename PolygonMesh>
Face_location compute_location(const Point& p, const PolygonMesh& mesh)
{
  CGAL_assertion(false);
  return Face_location(); // @todo
}

template<typename Face_location, typename PolygonMesh>
boost::variant<typename boost::graph_traits<PolygonMesh>::vertex_descriptor,
               typename boost::graph_traits<PolygonMesh>::halfedge_descriptor,
               typename boost::graph_traits<PolygonMesh>::face_descriptor>
get_descriptor_from_location(const Face_location& loc, const PolygonMesh& mesh)
{
  typedef boost::variant<
    typename boost::graph_traits<PolygonMesh>::vertex_descriptor,
    typename boost::graph_traits<PolygonMesh>::halfedge_descriptor,
    typename boost::graph_traits<PolygonMesh>::face_descriptor>   descriptor_variant;

  CGAL_assertion(false);
  return descriptor_variant();
}

template<typename Point, typename Face_location, typename PolygonMesh>
Point loc_to_point(const Face_location& loc, const PolygonMesh& mesh)
{
  CGAL_assertion(false);
  return Point(); // @todo
}

} // namespace internal

// Information associated to any point that is involved in the motorcycle graph algorithm
template<typename K, typename PolygonMesh>
class Dictionary_entry
{
  typedef Dictionary_entry<K, PolygonMesh>                         Self;

public:
  typedef typename K::FT                                           FT;
  typedef typename K::Point_2                                      Point;

  // A container of motorcycles that reach this point. We need to efficiently know:
  // - if a motorcycle is in the container
  // - the visiting times (when do they reach this point) to detect simultaneous collisions
  // Note that we need a multiset for the second container because motorcycles
  // are unique, but times are not
  typedef boost::bimap<
            boost::bimaps::set_of<int>, // set of motorcycles
            boost::bimaps::multiset_of<FT> > // multi-set of visiting times
                                                                   Visiting_motorcycles_container;
  typedef typename Visiting_motorcycles_container::value_type      value_type;
  typedef typename Visiting_motorcycles_container::size_type       size_type;

  typedef typename Visiting_motorcycles_container::left_iterator         VMC_left_it;
  typedef typename Visiting_motorcycles_container::left_const_iterator   VMC_left_cit;
  typedef typename Visiting_motorcycles_container::right_const_iterator  VMC_right_cit;

  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor     face_descriptor;

  // \brief An ordered pair specifying a location on the surface of the `Triangle_mesh`.
  // \details If `tm` is the input graph and given the pair (`f`, `bc`)
  // such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
  // and the vertices of the face `f` is the following:
  // - `w0 = source(halfedge(f, tm), tm)`
  // - `w1 = target(halfedge(f, tm), tm)`
  // - `w2 = target(next(halfedge(f, tm), tm), tm)`
  typedef typename CGAL::cpp11::array<FT, 3>                       Barycentric_coordinates;
  typedef std::pair<face_descriptor, Barycentric_coordinates>      Face_location;

  // Access
  const Point& point() const { return p; }
  void set_location(const Face_location& l) const { loc = l; }
  const Face_location& location() const { return loc; }
  bool is_blocked() const { return blocked; }
  void block() const { blocked = true; }
  const Visiting_motorcycles_container& visiting_motorcycles() const { return visiting_mcs; }

  // Constructor
  Dictionary_entry(const Point& p);
  Dictionary_entry(const Point& p, const Face_location& loc);

  // Most of these functions are not actually 'const' but the members they modify
  // are mutable. See next comment.
  void add_motorcycle(const int id, const FT time) const;
  bool has_motorcycle(const int id) const;
  bool has_simultaneous_collision() const;
  VMC_left_it find_motorcycle(const int id) const;
  size_type remove_motorcycle(const int id) const;

  friend bool operator<(const Self& lhs, const Self& rhs) {
    if(lhs.point() == rhs.point())
      return lhs.location().first < rhs.location().first;

    return lhs.point() < rhs.point();
  }

  // Output
  friend std::ostream& operator<<(std::ostream& out, const Self& dec) {
    out << "Point: (" << dec.point() << ") blocked: " << dec.is_blocked() << std::endl;
    out << "  Visiting motorcycles: " << std::endl;
    VMC_left_cit vmc_it = dec.visiting_motorcycles().left.begin();
    VMC_left_cit end = dec.visiting_motorcycles().left.end();
    for(; vmc_it!=end; ++vmc_it)
      out << "\t motorcycle: " << vmc_it->first << " time: " << vmc_it->second << std::endl;

    return out;
  }

private:
  // The position of the point
  const Point p;

  // Its location in the mesh
  mutable Face_location loc;

  // This class is meant to be an element of a set entry, which is const. However,
  // The members below must still be modified so they are made mutable.
  // It is safe because the set is only ordered with 'p', which is never modified.
  mutable Visiting_motorcycles_container visiting_mcs;
  mutable bool blocked;
};

template<typename K, typename PolygonMesh>
Dictionary_entry<K, PolygonMesh>::
Dictionary_entry(const Point& p)
  : p(p), loc(), visiting_mcs(), blocked(false)
{ }

template<typename K, typename PolygonMesh>
Dictionary_entry<K, PolygonMesh>::
Dictionary_entry(const Point& p, const Face_location& loc)
  : p(p), loc(loc), visiting_mcs(), blocked(false)
{ }

template<typename K, typename PolygonMesh>
void
Dictionary_entry<K, PolygonMesh>::
add_motorcycle(const int id, const FT time) const
{
  // the motorcycle `i` should not already exists in the list of motorcycles
  // that (might) reach this point
  // -- disabled to handle motorcycles with identical source and target
//  CGAL_precondition(visiting_mcs.left.find(id) == visiting_mcs.left.end());
  visiting_mcs.insert(value_type(id, time));
}

template<typename K, typename PolygonMesh>
bool
Dictionary_entry<K, PolygonMesh>::
has_motorcycle(const int id) const
{
  return (find_motorcycle(id) != visiting_mcs.left.end());
}

template<typename K, typename PolygonMesh>
bool
Dictionary_entry<K, PolygonMesh>::
has_simultaneous_collision() const
{
  CGAL_precondition(!visiting_mcs.empty());
  if(visiting_mcs.size() == 1)
    return false;

  // accessing 'right' of the bimap gives ordered time values
  // the first and second entries are thus the times for the two closest
  // points. If the times are equal, there is a simultaneous collision.
  VMC_right_cit first_mc_it = visiting_mcs.right.begin();
  VMC_right_cit second_mc_it = ++(visiting_mcs.right.begin());
  const FT first_time = first_mc_it->first;
  const FT second_time = second_mc_it->first;
  CGAL_assertion(first_time <= second_time);
  return (first_time == second_time);
}

template<typename K, typename PolygonMesh>
typename Dictionary_entry<K, PolygonMesh>::VMC_left_it
Dictionary_entry<K, PolygonMesh>::
find_motorcycle(const int id) const
{
  return visiting_mcs.left.find(id);
}

template<typename K, typename PolygonMesh>
typename Dictionary_entry<K, PolygonMesh>::size_type
Dictionary_entry<K, PolygonMesh>::
remove_motorcycle(const int id) const
{
  CGAL_precondition(find_motorcycle(id) != visiting_mcs.left.end());
  return visiting_mcs.left.erase(id);
}

// -----------------------------------------------------------------------------
//                           Dictionary class

template<typename K, typename PolygonMesh>
class Dictionary
{
public:
  typedef typename K::FT                                  FT;
  typedef typename K::Point_2                             Point;

  typedef Dictionary_entry<K, PolygonMesh>                Dictionary_entry;
  typedef typename Dictionary_entry::Face_location        Face_location;

  // @todo doesn't need to be an (ordered) set? (find out a good hash function...)
  typedef std::set<Dictionary_entry>                      Dictionary_entry_container;
  typedef typename Dictionary_entry_container::iterator   DEC_it;

  const Dictionary_entry_container& all_entries() const { return entries; }

  // Constructor
  Dictionary() : entries() { }

  // Functions
  DEC_it insert(const Point& p, const Face_location& loc = Face_location());
  DEC_it insert(const Point& p, const Face_location& loc, const int i, const FT time);
  DEC_it insert(const Point& p, const int i, const FT time);

private:
  Dictionary_entry_container entries;
};

template<typename K, typename PolygonMesh>
typename Dictionary<K, PolygonMesh>::DEC_it
Dictionary<K, PolygonMesh>::
insert(const Point& p, const Face_location& loc)
{
  Dictionary_entry e(p, loc);
  std::pair<DEC_it, bool> is_insert_successful = entries.insert(e);

  if(!is_insert_successful.second)
  {
    std::cerr << "Warning: point (" << p << ") already exists in the dictionary: "
              << *(is_insert_successful.first) << std::endl;
  }

  return is_insert_successful.first;
}

template<typename K, typename PolygonMesh>
typename Dictionary<K, PolygonMesh>::DEC_it
Dictionary<K, PolygonMesh>::
insert(const Point& p, const Face_location& loc, const int i, const FT time)
{
  DEC_it it = insert(p, loc);
  it->add_motorcycle(i, time);

  std::cout << "(New) point in the dictionary : " << *it << std::endl;

  return it;
}

template<typename K, typename PolygonMesh>
typename Dictionary<K, PolygonMesh>::DEC_it
Dictionary<K, PolygonMesh>::
insert(const Point& p, const int i, const FT time)
{
  return insert(p, Face_location(), i, time);
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_DICTIONARY_H
