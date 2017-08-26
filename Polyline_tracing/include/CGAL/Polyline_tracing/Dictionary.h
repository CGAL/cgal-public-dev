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

#ifndef CGAL_DICTIONARY_H
#define CGAL_DICTIONARY_H

#include <CGAL/assertions.h>

#include <limits>
#include <set>

namespace CGAL {

namespace Polyline_tracing {

template<typename K>
class Dictionary_entry
{
public:
  typedef typename K::FT                                  FT;
  typedef typename K::Point_2                             Point;

  bool is_blocked() const { return blocked; }
  void block() { blocked = true; }

  Dictionary_entry(const Point& p = Point());
  void add_motorcycle(const int i, const FT distance);
  bool operator<(const Dictionary_entry& e);

private:
  Point p;
  std::set<int> motorcycles;
  bool blocked;
  int i, j; // ids of the two motorcycles with smallest distance at p
  FT dist_at_i, dist_at_j;
};

template<typename K>
Dictionary_entry<K>::
Dictionary_entry(const Point& p)
  : p(p), motorcycles(), blocked(false), i(-1), j(-1),
    dist_at_i(std::numeric_limits<FT>::max()),
    dist_at_j(std::numeric_limits<FT>::max())
{ }

template<typename K>
void
Dictionary_entry<K>::
add_motorcycle(const int id, const FT distance)
{
  CGAL_precondition(motorcycles.find(i) == motorcycles.end());

  motorcycles.insert(id);

  // what about three intersections @fixme
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

template<typename K>
bool
Dictionary_entry<K>::
operator<(const Dictionary_entry<K>& e)
{
  return p < e.p;
}

// -----------------------------------------------------------------------------

template<typename K>
class Dictionary
{
public:
  typedef typename K::FT                                  FT;
  typedef typename K::Point_2                             Point;

  typedef Dictionary_entry<K>                             Dictionary_entry;
  typedef std::set<Dictionary_entry>                      Dictionary_entry_container;
  typedef Dictionary_entry_container::iterator            DEC_it;

  Dictionary() : entries() { }

  void insert(const Point& p, const int i, const FT dist);

private:
  Dictionary_entry_container entries;
};

template<typename K>
void
Dictionary<K>::
insert(const Point& p, const int i, const FT dist)
{
  Dictionary_entry e(p);

  std::pair<DEC_it, bool> is_insert_successful = entries.insert(e);

  DEC_it it = is_insert_successful.first;
  it->add_motorcycle(i, dist);
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_DICTIONARY_H
