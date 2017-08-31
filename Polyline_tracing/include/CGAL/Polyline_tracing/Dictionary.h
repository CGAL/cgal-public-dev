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

#include <CGAL/assertions.h>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/unordered_set.hpp>

#include <limits>
#include <map>
#include <utility>

namespace CGAL {

namespace Polyline_tracing {

// Information associated to any point that is involved in the motorcycle graph algorithm
template<typename K>
class Dictionary_entry
{
  typedef Dictionary_entry<K>                                      Self;

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

  const Point& point() const { return p; }
  bool is_blocked() const { return blocked; }
  void block() const { blocked = true; }
  const Visiting_motorcycles_container& visiting_motorcycles() const { return visiting_mcs; }

  Dictionary_entry(const Point& p);

  // Most of these functions are not actually 'const' but the members they modify
  // are mutable. See next comment.
  void add_motorcycle(const int id, const FT time) const;
  bool has_motorcycle(const int id) const;
  bool has_simultaneous_collision() const;
  VMC_left_it find_motorcycle(const int id) const;
  size_type remove_motorcycle(const int id) const;

  friend bool operator<(const Self& lhs, const Self& rhs) {
    return lhs.point() < rhs.point();
  }

  // Output
  friend std::ostream& operator<<(std::ostream& out, const Self& dec) {
    out << "Point: " << dec.point() << " blocked: " << dec.is_blocked() << std::endl;
    out << " Visiting motorcycles: " << std::endl;
    VMC_left_cit vmc_it = dec.visiting_motorcycles().left.begin();
    VMC_left_cit end = dec.visiting_motorcycles().left.end();
    for(; vmc_it!=end; ++vmc_it)
      out << " motorcycle: " << vmc_it->first << " time: " << vmc_it->second << std::endl;

    return out;
  }

private:
  const Point p;

  // This class is meant to be an element of a set entry, which is const. However,
  // The members below must still be modified so they are made mutable.
  // It is safe because the set is only ordered with 'p', which is never modified.
  mutable Visiting_motorcycles_container visiting_mcs;
  mutable bool blocked;
};

template<typename K>
Dictionary_entry<K>::
Dictionary_entry(const Point& p)
  : p(p), visiting_mcs(), blocked(false)
{ }

template<typename K>
void
Dictionary_entry<K>::
add_motorcycle(const int id, const FT time) const
{
  // the motorcycle `i` should not already exists in the list of motorcycles
  // that (might) reach this point
  CGAL_precondition(visiting_mcs.left.find(i) == visiting_mcs.left.end());
  visiting_mcs.insert(value_type(id, time));
}

template<typename K>
bool
Dictionary_entry<K>::
has_motorcycle(const int id) const
{
  return (find_motorcycle(id) != visiting_mcs.left.end());
}

template<typename K>
bool
Dictionary_entry<K>::
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

template<typename K>
typename Dictionary_entry<K>::VMC_left_it
Dictionary_entry<K>::
find_motorcycle(const int id) const
{
  return visiting_mcs.left.find(id);
}

template<typename K>
typename Dictionary_entry<K>::size_type
Dictionary_entry<K>::
remove_motorcycle(const int id) const
{
  CGAL_precondition(find_motorcycle(id) != visiting_mcs.left.end());
  return visiting_mcs.left.erase(id);
}

// -----------------------------------------------------------------------------
//                           Dictionary class

template<typename K>
class Dictionary
{
public:
  typedef typename K::FT                                  FT;
  typedef typename K::Point_2                             Point;

  typedef Dictionary_entry<K>                             Dictionary_entry;
  typedef std::set<Dictionary_entry>                      Dictionary_entry_container;
  typedef typename Dictionary_entry_container::iterator   DEC_it;

  // Constructor
  Dictionary() : entries() { }

  // Functions
  DEC_it insert(const Point& p);
  DEC_it insert(const Point& p, const int i, const FT dist);

private:
  Dictionary_entry_container entries;
};

template<typename K>
typename Dictionary<K>::DEC_it
Dictionary<K>::
insert(const Point& p)
{
  Dictionary_entry e(p);
  std::pair<DEC_it, bool> is_insert_successful = entries.insert(e);

  if(!is_insert_successful.second)
    std::cerr << "Warning: point (" << p << ") already exists in the dictionary" << std::endl;

  return is_insert_successful.first;
}

template<typename K>
typename Dictionary<K>::DEC_it
Dictionary<K>::
insert(const Point& p, const int i, const FT dist)
{
  DEC_it it = insert(p);
  it->add_motorcycle(i, dist);

  return it;
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_DICTIONARY_H
