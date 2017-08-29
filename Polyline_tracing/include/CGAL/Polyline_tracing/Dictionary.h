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

#include <boost/unordered_set.hpp>

#include <limits>
#include <map>
#include <utility>

namespace CGAL {

namespace Polyline_tracing {

template<typename K>
struct Motorcycle_set_comparer
{
  typedef typename K::FT FT;

  bool operator()(const std::pair<int, FT>& lhs,
                  const std::pair<int, FT>& rhs) const {
    return lhs.first < rhs.first;
  }
};

template<typename K>
class Dictionary_entry
{
public:
  typedef typename K::FT                                  FT;
  typedef typename K::Point_2                             Point;

  bool is_blocked() const { return blocked; }
  void block() { blocked = true; }

  Dictionary_entry();
  void add_motorcycle(const int i, const FT distance);

private:
  boost::unordered_set<int> motorcycles;
  bool blocked;
  int i, j; // ids of the two motorcycles with smallest distance at p
  FT dist_at_i, dist_at_j;
};

template<typename K>
Dictionary_entry<K>::
Dictionary_entry()
  : motorcycles(), blocked(false), i(-1), j(-1),
    dist_at_i(std::numeric_limits<FT>::max()),
    dist_at_j(std::numeric_limits<FT>::max())
{ }

template<typename K>
void
Dictionary_entry<K>::
add_motorcycle(const int id, const FT distance)
{
  // the motorcycle `i` should not already exists in the list of motorcycle
  // that (might) reach this point
  CGAL_precondition(motorcycles.find(i) == motorcycles.end());

  motorcycles.insert(id);

  if(distance < dist_at_i)
  {
    i = id;
    dist_at_i = distance;
  }
  else if(distance < dist_at_j)
  {
    j = id;
    dist_at_j = distance;
  }
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
  typedef std::map<Point, Dictionary_entry>               Dictionary_entry_container;
  typedef typename Dictionary_entry_container::iterator   DEC_it;

  Dictionary() : entries() { }

  DEC_it insert(const Point& p, const int i, const FT dist);

private:
  Dictionary_entry_container entries;
};

template<typename K>
typename Dictionary<K>::DEC_it
Dictionary<K>::
insert(const Point& p, const int i, const FT dist)
{
  Dictionary_entry e;
  std::pair<DEC_it, bool> is_insert_successful = entries.insert(std::make_pair(p, e));

  if(!is_insert_successful.second)
    std::cerr << "Warning: point " << p << " already exists in the dictionary" << std::endl;

  DEC_it it = is_insert_successful.first;
  it->second.add_motorcycle(i, dist);

  return it;
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_DICTIONARY_H
