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

#ifndef CGAL_MOTORCYCLE_H
#define CGAL_MOTORCYCLE_H

#include <list>
#include <set>

namespace CGAL {

namespace Polyline_tracing {

template<typename K>
class Motorcycle
{
public:
  typedef typename K::FT                                        FT;
  typedef typename K::Point_2                                   Point;
  typedef typename K::Vector_2                                  Vector;

  int id() const { return i; }
  void set_id(int id) { i = id; }
  bool is_crashed() const { return crashed; }
  void crash() { crashed = true; }

  const Point& position() const { return pos; }
  void position(const Point& p) { pos = p; }
  const Point& destination() const { return dest; }
  void destination(const Point& p) { dest = p; }

  Motorcycle(const int id, const Point& p, const Point& v);

  bool has_reached_end();
  void trace();

private:
  int i;
  bool crashed;

  Point pos;
  Point source;
  Point dest;
  FT distance;
  Point confirmed;
  Point tentative;
  std::list<Point> path;
  std::set<Point> target_points;
};

template<typename K>
Motorcycle<K>::
Motorcycle(const int id, const Point& p, const Point& d)
  : i(id), crashed(false),
    pos(p), source(p), dest(d), distance(0.),
    confirmed(p), tentative(p),
    path(), target_points()
{
  target_points.insert(pos);
  target_points.insert(dest);
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_MOTORCYCLE_H
