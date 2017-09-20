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
#include <CGAL/Polygon_mesh_processing/locate.h>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>

#include <utility>

namespace CGAL {

namespace Polyline_tracing {

// This class represents a point that is involved in the motorcycle graph algorithm
// It is useful for robustness to regroup them in a single dictionary-like structure.
template<typename K, typename PolygonMesh>
class Dictionary_entry
{
  typedef Dictionary_entry<K, PolygonMesh>                         Self;

public:
  typedef typename K::FT                                           FT;
  typedef typename K::Point_2                                      Point;

  // A container of motorcycles that reach this point. We need to efficiently know:
  // - if a motorcycle reaches this point
  // - the ordered reaching times to detect simultaneous collisions
  // Note that we need a multiset for the second container because motorcycles
  // are unique, but arrival times might not be.
  typedef boost::bimap<
            boost::bimaps::set_of<std::size_t>, // set of motorcycles
            boost::bimaps::multiset_of<FT> > // multi-set of visiting times
                                                                   Visiting_motorcycles_container;
  typedef typename Visiting_motorcycles_container::value_type      value_type;
  typedef typename Visiting_motorcycles_container::size_type       size_type;

  typedef typename Visiting_motorcycles_container::left_iterator         VMC_left_it;
  typedef typename Visiting_motorcycles_container::left_const_iterator   VMC_left_cit;
  typedef typename Visiting_motorcycles_container::right_const_iterator  VMC_right_cit;

  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor     face_descriptor;

  // Points are not described through a Point_d, but through an ordered pair
  // specifying a location on the surface of the `Triangle_mesh`.
  //
  // If `tm` is the input graph and given the pair (`f`, `bc`)
  // such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
  // and the vertices of the face `f` is the following:
  // - `w0 = source(halfedge(f, tm), tm)`
  // - `w1 = target(halfedge(f, tm), tm)`
  // - `w2 = target(next(halfedge(f, tm), tm), tm)`
  typedef typename CGAL::cpp11::array<FT, 3>                       Barycentric_coordinates;
  typedef std::pair<face_descriptor, Barycentric_coordinates>      Face_location;

  // Access
  const Point& point() const { return p; }
  const Face_location& location() const { return loc; }
  bool is_blocked() const { return blocked; }
  void block() const { blocked = true; }
  const Visiting_motorcycles_container& visiting_motorcycles() const { return visiting_mcs; }

  // Constructor
  Dictionary_entry(const Face_location& loc, const Point& p);
  Dictionary_entry(const Face_location& loc);

  // Most of these functions are not actually 'const' but the members they modify
  // are mutable. See next comment.
  void add_motorcycle(const std::size_t id, const FT time) const;
  VMC_left_it find_motorcycle(const std::size_t id) const;
  bool has_motorcycle(const std::size_t id) const;
  bool has_simultaneous_collision() const;
  std::pair<bool, FT> motorcycle_visiting_time(const std::size_t id) const;
  size_type remove_motorcycle(const std::size_t id) const;

  // need to build a set of Dictionary_entry items
  friend bool operator<(const Self& lhs, const Self& rhs) {
    if(lhs.location().first == rhs.location().first)
    {
      return std::lexicographical_compare(lhs.location().second.begin(), lhs.location().second.end(),
                                          rhs.location().second.begin(), rhs.location().second.end());
    }

    return lhs.location().first < rhs.location().first;
  }

  // Output
  friend std::ostream& operator<<(std::ostream& out, const Self& dec) {
    out << "Point: (" << dec.point() << ") blocked: " << dec.is_blocked() << std::endl;
    out << "Location -- face: " << dec.location().first << " barycentric coordinates: { "
        << dec.location().second[0] << "; " << dec.location().second[1]
        << "; " << dec.location().second[2] << "}" << std::endl;
    out << "  Visiting motorcycles: " << std::endl;
    VMC_left_cit vmc_it = dec.visiting_motorcycles().left.begin();
    VMC_left_cit end = dec.visiting_motorcycles().left.end();
    for(; vmc_it!=end; ++vmc_it)
      out << "\t motorcycle: " << vmc_it->first << " time: " << vmc_it->second << std::endl;

    return out;
  }

private:
  // Its location in the mesh
  const Face_location loc;

  // The position of the point (only for information)
  const Point p;

  // This class is meant to be an element of a set entry, which is const. However,
  // the members below must still be modifiable so they are made mutable.
  // It is safe because the set is only ordered with the location, which is never modified.
  mutable bool blocked;
  mutable Visiting_motorcycles_container visiting_mcs;
};

template<typename K, typename PolygonMesh>
Dictionary_entry<K, PolygonMesh>::
Dictionary_entry(const Face_location& loc)
  : loc(loc), p(), blocked(false), visiting_mcs()
{
  CGAL_assertion(false);
  // @todo compute p with locate.h
}

template<typename K, typename PolygonMesh>
Dictionary_entry<K, PolygonMesh>::
Dictionary_entry(const Face_location& loc, const Point& p)
  : loc(loc), p(p), blocked(false), visiting_mcs()
{ }

template<typename K, typename PolygonMesh>
void
Dictionary_entry<K, PolygonMesh>::
add_motorcycle(const std::size_t id, const FT time) const
{
  // the motorcycle `i` should not already exists in the list of motorcycles
  // that (might) reach this point
  // -- disabled to handle motorcycles with identical source and target
//  CGAL_precondition(visiting_mcs.left.find(id) == visiting_mcs.left.end());
  visiting_mcs.insert(value_type(id, time));
}

template<typename K, typename PolygonMesh>
typename Dictionary_entry<K, PolygonMesh>::VMC_left_it
Dictionary_entry<K, PolygonMesh>::
find_motorcycle(const std::size_t id) const
{
  return visiting_mcs.left.find(id);
}

template<typename K, typename PolygonMesh>
bool
Dictionary_entry<K, PolygonMesh>::
has_motorcycle(const std::size_t id) const
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
std::pair<bool, typename Dictionary_entry<K, PolygonMesh>::FT>
Dictionary_entry<K, PolygonMesh>::
motorcycle_visiting_time(const std::size_t id) const
{
  VMC_left_it mit = find_motorcycle(id);

  if(mit == visiting_mcs.left.end())
    return std::make_pair(false, 0.);
  else
    return std::make_pair(true, mit->second);
}

template<typename K, typename PolygonMesh>
typename Dictionary_entry<K, PolygonMesh>::size_type
Dictionary_entry<K, PolygonMesh>::
remove_motorcycle(const std::size_t id) const
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
  void erase(DEC_it pos);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const Point& p);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const Point& p, const std::size_t i, const FT time);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const std::size_t i, const FT time, const PolygonMesh& mesh);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const PolygonMesh& mesh);

private:
  Dictionary_entry_container entries;
};

template<typename K, typename PolygonMesh>
void
Dictionary<K, PolygonMesh>::
erase(DEC_it pos)
{
  return entries.erase(pos);
}

template<typename K, typename PolygonMesh>
std::pair<typename Dictionary<K, PolygonMesh>::DEC_it, bool>
Dictionary<K, PolygonMesh>::
insert(const Face_location& loc, const Point& p)
{
  Dictionary_entry e(loc, p);
  std::pair<DEC_it, bool> is_insert_successful = entries.insert(e);

  if(!is_insert_successful.second)
  {
    std::cerr << "Warning: point (" << p << ") already exists in the dictionary: " << std::endl
              << *(is_insert_successful.first) << std::endl;
  }

  return is_insert_successful;
}

template<typename K, typename PolygonMesh>
std::pair<typename Dictionary<K, PolygonMesh>::DEC_it, bool>
Dictionary<K, PolygonMesh>::
insert(const Face_location& loc, const Point& p, const std::size_t i, const FT time)
{
  std::pair<DEC_it, bool> entry = insert(loc, p);
  entry.first->add_motorcycle(i, time);

  std::cout << "(New) point in the dictionary : " << *(entry.first) << std::endl;

  return entry;
}

template<typename K, typename PolygonMesh>
std::pair<typename Dictionary<K, PolygonMesh>::DEC_it, bool>
Dictionary<K, PolygonMesh>::
insert(const Face_location& loc, const std::size_t i, const FT time, const PolygonMesh& mesh)
{
  Point p = CGAL::Polygon_mesh_processing::internal::loc_to_point(loc, mesh);
  return insert(loc, p, i, time);
}

template<typename K, typename PolygonMesh>
std::pair<typename Dictionary<K, PolygonMesh>::DEC_it, bool>
Dictionary<K, PolygonMesh>::
insert(const Face_location& loc, const PolygonMesh& mesh)
{
  Point p = CGAL::Polygon_mesh_processing::internal::loc_to_point(loc, mesh);
  return insert(loc, p);
}

} // namespace Polyline_tracing
} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_DICTIONARY_H
