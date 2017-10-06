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
#include <boost/functional/hash.hpp>
#include <boost/unordered_set.hpp>

#include <utility>

namespace CGAL {

namespace Polyline_tracing {

// This class represents a point that is involved in the motorcycle graph algorithm
// It is useful for robustness to regroup them in a single dictionary-like structure.
template<typename MotorcycleGraphTraits>
class Dictionary_entry
{
  typedef Dictionary_entry<MotorcycleGraphTraits>                        Self;

public:
  typedef MotorcycleGraphTraits                                          Geom_traits;

  typedef typename Geom_traits::FT                                       FT;
  typedef typename Geom_traits::Point_d                                  Point;

  // A container of motorcycles that visit this point. We need to efficiently know:
  // - if a motorcycle visits this point,
  // - the ordered visiting times to detect simultaneous collisions.
  //
  // Note that we need multisets:
  // - for the first container, because a motorcycle might visit a point twice
  //   (e.g. if the trajectory makes a loop and self-intersects at the point);
  // - for the second container, because arrival times might not be unique.
  typedef boost::bimap<
            boost::bimaps::multiset_of<std::size_t>, // set of motorcycles
            boost::bimaps::multiset_of<FT> > // multi-set of visiting times
                                                                         Visiting_motorcycles_container;
  typedef typename Visiting_motorcycles_container::iterator              iterator;
  typedef typename Visiting_motorcycles_container::value_type            value_type;
  typedef typename Visiting_motorcycles_container::size_type             size_type;

  typedef typename Visiting_motorcycles_container::left_iterator         VMC_left_it;
  typedef typename Visiting_motorcycles_container::left_const_iterator   VMC_left_cit;
  typedef typename Visiting_motorcycles_container::right_const_iterator  VMC_right_cit;

  typedef typename Geom_traits::face_descriptor                          face_descriptor;

  typedef typename Geom_traits::Face_location                            Face_location;

  // Access
  const Point& point() const { return p; }
  const Face_location& location() const { return loc; }
  bool is_blocked() const { return blocked; }
  void block() const { blocked = true; }
  const Visiting_motorcycles_container& visiting_motorcycles() const { return visiting_mcs; }

  // Constructor
  Dictionary_entry(const Face_location& loc, const Point& p);

  // The following function is not actually 'const' but the members it modifies
  // are mutable.
  std::pair<iterator, bool> add_motorcycle(const std::size_t id, const FT time) const;

  // second bool indicates whether we found a motorcycle with id 'id' or not
  std::pair<VMC_left_it, bool> find_motorcycle(const std::size_t id) const;
  // check if there is a motorcycle visiting at time 'visiting_time'
  bool has_motorcycle(const std::size_t id, const FT visiting_time) const;
  // check if there is a motorcycle visiting at time between 'min_' and 'max_visiting_time'
  // the last parameter is optional and can be used to grab the visiting time
  // 'strictly_at_X' to include the boundary of the interval or not
  bool has_motorcycle(const std::size_t id, const FT min_visiting_time,
                      const FT max_visiting_time, FT& visiting_time,
                      const bool strictly_at_min = false,
                      const bool strictly_at_max = false) const;
  bool has_motorcycle(const std::size_t id, const FT min_visiting_time,
                      const FT max_visiting_time, const bool strictly_at_min = false,
                      const bool strictly_at_max = false) const;
  // check if the two earliest motorcycles meet at the same time
  bool has_simultaneous_collision() const;
  size_type remove_motorcycle(const std::size_t id) const;

  // to build a set<Dictionary_entry>
  friend bool operator<(const Self& lhs, const Self& rhs)
  {
    if(lhs.location().first == rhs.location().first)
    {
      // lexicographical compare over the barycentric coordinates
      return std::lexicographical_compare(lhs.location().second.begin(), lhs.location().second.end(),
                                          rhs.location().second.begin(), rhs.location().second.end());
    }

    return lhs.location().first < rhs.location().first;
  }

  friend bool operator==(const Self& lhs, const Self& rhs)
  {
    return (lhs.location().first == rhs.location().first &&
            lhs.location().second == rhs.location().second);
  }

  friend std::size_t hash_value(const Self& dec)
  {
    boost::hash<face_descriptor> face_hasher;
    std::size_t seed = 0;
    boost::hash_combine(seed, face_hasher(dec.location().first));
    boost::hash_combine(seed, boost::hash_range(dec.location().second.begin(),
                                                dec.location().second.end()));
    return seed;
  }

  // Output
  friend std::ostream& operator<<(std::ostream& out, const Self& dec)
  {
    out << "Point: (" << dec.point() << ") blocked: " << dec.is_blocked() << std::endl;
    out << "  Location -- face: " << dec.location().first << " barycentric coordinates: { "
        << dec.location().second[0] << "; " << dec.location().second[1]
        << "; " << dec.location().second[2] << "}" << std::endl;
    out << "  Visiting motorcycles: " << std::endl;
    VMC_left_cit vmc_it = dec.visiting_motorcycles().left.begin();
    VMC_left_cit end = dec.visiting_motorcycles().left.end();
    for(; vmc_it!=end; ++vmc_it)
      out << "\t motorcycle #" << vmc_it->first << " time: " << vmc_it->second << std::endl;

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

template<typename MotorcycleGraphTraits>
Dictionary_entry<MotorcycleGraphTraits>::
Dictionary_entry(const Face_location& loc, const Point& p)
  : loc(loc), p(p), blocked(false), visiting_mcs()
{ }

template<typename MotorcycleGraphTraits>
std::pair<typename Dictionary_entry<MotorcycleGraphTraits>::iterator, bool>
Dictionary_entry<MotorcycleGraphTraits>::
add_motorcycle(const std::size_t id, const FT time) const
{
  std::cout << " Point " << this << " is visited by motorcycle #" << id << " at time: " << time << std::endl;
  return visiting_mcs.insert(value_type(id, time));
}

template<typename MotorcycleGraphTraits>
std::pair<typename Dictionary_entry<MotorcycleGraphTraits>::VMC_left_it, bool>
Dictionary_entry<MotorcycleGraphTraits>::
find_motorcycle(const std::size_t id) const
{
  // Since 'lower_bound' is used, it returns the first motorcycle fitting this
  VMC_left_it it = visiting_mcs.left.lower_bound(id);
  bool found_motorcycle = (it != visiting_mcs.left.end() && it->first == id);

  return std::make_pair(it, found_motorcycle);
}

template<typename MotorcycleGraphTraits>
bool
Dictionary_entry<MotorcycleGraphTraits>::
has_motorcycle(const std::size_t id, const FT visiting_time) const
{
  std::pair<VMC_left_it, bool> mres = find_motorcycle(id);
  VMC_left_it mit = mres.first;
  bool is_valid_iterator = mres.second;

  while(is_valid_iterator)
  {
    if(mit->second == visiting_time)
      return true;

    ++mit;
    is_valid_iterator = (mit != visiting_mcs.left.end() && mit->first == id);
  }

  return false;
}

template<typename MotorcycleGraphTraits>
bool
Dictionary_entry<MotorcycleGraphTraits>::
has_motorcycle(const std::size_t id, const FT min_visiting_time,
               const FT max_visiting_time, FT& visiting_time,
               const bool strictly_at_min, const bool strictly_at_max) const
{
  std::pair<VMC_left_cit, bool> mres = find_motorcycle(id);
  VMC_left_cit mit = mres.first;
  bool is_valid_iterator = mres.second; // = false if we couldn't find the motorcycle

  while(is_valid_iterator) // while still considering the motorcycle 'id'
  {
    CGAL_assertion(mit->first == id);
    visiting_time = mit->second;

    if((visiting_time > min_visiting_time ||
        (!strictly_at_min && visiting_time == min_visiting_time)) &&
       (visiting_time < max_visiting_time ||
        (!strictly_at_max && visiting_time == max_visiting_time)))
      return true;

    ++mit;
    is_valid_iterator = (mit != visiting_mcs.left.end() && mit->first == id);
  }

  return false;
}

template<typename MotorcycleGraphTraits>
bool
Dictionary_entry<MotorcycleGraphTraits>::
has_motorcycle(const std::size_t id, const FT min_visiting_time,
               const FT max_visiting_time, const bool strictly_at_min,
               const bool strictly_at_max) const
{
  FT useless;
  return has_motorcycle(id, min_visiting_time, max_visiting_time, useless,
                        strictly_at_min, strictly_at_max);
}

template<typename MotorcycleGraphTraits>
bool
Dictionary_entry<MotorcycleGraphTraits>::
has_simultaneous_collision() const
{
  CGAL_precondition(!visiting_mcs.empty());
  if(visiting_mcs.size() == 1)
    return false;

  // accessing 'right' of the bimap gives ordered time values
  // the first and second entries are thus the times for the two closest
  // points. If the times are equal (or almost), there is a simultaneous collision.
  VMC_right_cit first_mc_it = visiting_mcs.right.begin();
  VMC_right_cit second_mc_it = ++(visiting_mcs.right.begin());
  const FT first_time = first_mc_it->first;
  const FT second_time = second_mc_it->first;
  CGAL_assertion(first_time <= second_time);

#ifdef CGAL_MOTORCYCLE_GRAPH_ROBUSTNESS_CODE
  // Add a little bit of tolerance
  return (CGAL::abs(second_time - first_time) < (std::numeric_limits<FT>::epsilon() *
                                                 CGAL::abs(first_time + second_time)));
#else
  return (second_time == first_time);
#endif

}

template<typename MotorcycleGraphTraits>
typename Dictionary_entry<MotorcycleGraphTraits>::size_type
Dictionary_entry<MotorcycleGraphTraits>::
remove_motorcycle(const std::size_t id) const
{
  // Note that this will remove all instances of motorcycle 'id' in the multiset
  return visiting_mcs.left.erase(id);
}

// -----------------------------------------------------------------------------
//                           Dictionary class

template<typename MotorcycleGraphTraits>
class Dictionary
{
  typedef Dictionary<MotorcycleGraphTraits>               Self;

public:
  typedef MotorcycleGraphTraits                           Geom_traits;
  typedef typename Geom_traits::Triangle_mesh             Triangle_mesh;

  typedef typename Geom_traits::FT                        FT;
  typedef typename Geom_traits::Point_d                   Point;

  typedef Dictionary_entry<MotorcycleGraphTraits>         Dictionary_entry;
  typedef typename Dictionary_entry::Face_location        Face_location;

  typedef boost::unordered_set<Dictionary_entry>          Dictionary_entry_container;
  typedef typename Dictionary_entry_container::iterator   DEC_it;

  // Access
  const Dictionary_entry_container& all_entries() const { return entries; }

  // Constructor
  Dictionary() : entries() { }

  // Functions
  std::pair<DEC_it, bool> find(const Dictionary_entry& e) const;
  std::pair<DEC_it, bool> find(const Face_location& loc) const;
  DEC_it erase(DEC_it pos);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const Point& p);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const Point& p, const std::size_t i, const FT time);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const std::size_t i, const FT time, const Triangle_mesh& mesh);
  std::pair<DEC_it, bool> insert(const Face_location& loc, const Triangle_mesh& mesh);

  // Ouput
  friend std::ostream& operator<<(std::ostream& out, const Self& d)
  {
    DEC_it it = d.all_entries().begin(), end = d.all_entries().end();
    for(; it!=end; ++it)
      out << &*it << " " << *it << std::endl;

    return out;
  }

private:
  Dictionary_entry_container entries;
};

template<typename MotorcycleGraphTraits>
std::pair<typename Dictionary<MotorcycleGraphTraits>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits>::
find(const Dictionary_entry& e) const
{
  DEC_it res = entries.find(e);
  return std::make_pair(res, (res != entries.end()));
}

template<typename MotorcycleGraphTraits>
std::pair<typename Dictionary<MotorcycleGraphTraits>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits>::
find(const Face_location& loc) const
{
  Dictionary_entry dummy(loc, Point());
  return find(dummy);
}

template<typename MotorcycleGraphTraits>
typename Dictionary<MotorcycleGraphTraits>::DEC_it
Dictionary<MotorcycleGraphTraits>::
erase(DEC_it pos)
{
  return entries.erase(pos);
}

template<typename MotorcycleGraphTraits>
std::pair<typename Dictionary<MotorcycleGraphTraits>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits>::
insert(const Face_location& loc, const Point& p)
{
  Dictionary_entry e(loc, p);
  std::pair<DEC_it, bool> is_insert_successful = entries.insert(e);

  if(!is_insert_successful.second)
  {
    std::cerr << "Warning: point already exists in the dictionary: "
              << &*(is_insert_successful.first) << std::endl
              << *(is_insert_successful.first) << std::endl;
  }
  else
  {
    std::cout << "New point in the dictionary: "
              << &*(is_insert_successful.first) << std::endl
              << *(is_insert_successful.first) << std::endl;
  }

  return is_insert_successful;
}

template<typename MotorcycleGraphTraits>
std::pair<typename Dictionary<MotorcycleGraphTraits>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits>::
insert(const Face_location& loc, const Point& p, const std::size_t i, const FT time)
{
  std::pair<DEC_it, bool> entry = insert(loc, p);
  entry.first->add_motorcycle(i, time);

  return entry;
}

template<typename MotorcycleGraphTraits>
std::pair<typename Dictionary<MotorcycleGraphTraits>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits>::
insert(const Face_location& loc, const std::size_t i, const FT time, const Triangle_mesh& mesh)
{
  Point p = CGAL::Polygon_mesh_processing::location_to_point(loc, mesh);
  return insert(loc, p, i, time);
}

template<typename MotorcycleGraphTraits>
std::pair<typename Dictionary<MotorcycleGraphTraits>::DEC_it, bool>
Dictionary<MotorcycleGraphTraits>::
insert(const Face_location& loc, const Triangle_mesh& mesh)
{
  Point p = CGAL::Polygon_mesh_processing::location_to_point(loc, mesh);
  return insert(loc, p);
}

} // namespace Polyline_tracing
} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_DICTIONARY_H
