// Copyright (c) 2010  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
//
// Author(s)     : Asaf Porat          <asafpor1@post.tau.ac.il>

#ifndef LINE_THROUGH_SEGMENTS_ISOLATED_POINTS_H
#define LINE_THROUGH_SEGMENTS_ISOLATED_POINTS_H

/*************************************************************
 * This class represent a data base of isolated 2 dimensional points.
 * Each entry at the list holds a point and the line that created the point.
 * The data base is sorted according to the x coordinate of the points if two
 * x coordinates equal than according to the y coordinates.
 *
 *************************************************************/

#include <list>

namespace CGAL {

template <typename Point_and_segment_pair, typename Compare>
class Lines_through_segments_isolated_points
{
public:
  typedef typename std::list<Point_and_segment_pair >::iterator iterator;
  typedef typename std::list<Point_and_segment_pair >::const_iterator
  const_iterator;
  typedef Point_and_segment_pair                         IP_point_and_line_pair;
         
private:
  std::list<Point_and_segment_pair > m_points_list;
  bool m_sorted;
  /*
    Comparison functor class - taking two values of the same type than those
    contained in the list object, returns true if the first argument goes
    before the second argument in the specific order.
    (i.e., if the first is less than the second), and false otherwise. 
  */
  Compare m_compare;
         
public:
         
  /* iterator support */
  iterator begin() { return m_points_list.begin(); }
  const_iterator begin() const { return m_points_list.begin(); }
  iterator end() { return m_points_list.end(); }
  const_iterator end() const { return m_points_list.end(); }
  iterator erase(iterator position) { return m_points_list.erase(position); }
         
  Lines_through_segments_isolated_points() {m_sorted = false;}
  Lines_through_segments_isolated_points(Compare& _compare) 
  {
    m_compare = _compare;m_sorted = false;
  }
         
  void add_element(Point_and_segment_pair pl)
  {
    m_points_list.push_back(pl);
    m_sorted = false;
  }
         
  void sort_points()
  {
    m_points_list.sort(m_compare);
    m_sorted = true;
  }
      
  ~Lines_through_segments_isolated_points()
  {
  }
         
  /*
    The following function finds all duplicated points and removes them from
    the list.
    output:
    list of duplicated points, each element at the list is a point 
    and two line segment that created it.
  */
  template <typename Point_and_two_objs>
  void remove_duplicated_points(std::list<Point_and_two_objs>& ret_points_list)
  {
    iterator it, next_it;
            
    /* Delete duplicated points.*/
    it = m_points_list.begin();
    while (it != m_points_list.end())
    {
      next_it = it;
      next_it++;
      if (next_it == m_points_list.end())
        break;
               
#if ISOLATED_POINTS_DEBUG
      std::cout << "Remove Point  = (" << ((*it).get_point()) 
                << ",line = "
                << *((*it).get_segment()) << ")" << std::endl;
#endif /* ISOLATED_POINTS_DEBUG */
              
      if (m_compare(*it,*next_it) && m_compare(*next_it,*it))
      {
        Point_and_segment_pair curr = *it;
        ret_points_list.push_back(Point_and_two_objs(it->get_point(),
                                                     it->get_segment(),
                                                     next_it->get_segment()));
        it = m_points_list.erase(it);
                  
        /* remove all duplicated points. */
        while (it != m_points_list.end() &&
               /* Compare equals true if first <= second */
               (m_compare(curr,*it) && m_compare(*it,curr)))
        {
          it = m_points_list.erase(it);
        }
      }
      else
      {
        it++;
      }
    }
  }
         
  unsigned int size()
  {
    return this->m_points_list.size();
  }

  void clear()
  {
    return this->m_points_list.clear();
  }
         
  std::string to_string()
  {
    iterator it;
    std::ostringstream o;
    o << "this->points_list.size = " << this->m_points_list.size() 
      << std::endl;
            
    for (it = this->m_points_list.begin();
         it != this->m_points_list.end(); ++it)
    {
      o << "Point  = (" << (*it).get_point() << "," <<
         "Line = " << *((*it).get_segment()) << ")" << std::endl;
    }
            
    return o.str();
  }
};

template <typename Point, typename Segment_3>
class Point_and_segment_pair
{
public:
  typedef Point PL_Point;
    
private:
  Point m_point;
  const Segment_3* m_segment;
    
public:
  Point_and_segment_pair(Point _point,const Segment_3* _segment)
  {
    m_point = _point;
    m_segment = _segment;
  }
    
  const Point get_point() const 
  {
    return m_point;
  }
    
  const Segment_3* get_segment() const
  {
    return m_segment;
  }
};

template <typename Rational_point_3, typename Rational_segment_3>
class Compare_points_on_line
{
  typedef Point_and_segment_pair<Rational_point_3, Rational_segment_3>
  Point_on_line_and_segment_pair;

public:
  bool operator() (const Point_on_line_and_segment_pair& first,
                   const Point_on_line_and_segment_pair& second)
  {
    if (first.get_point().x() < second.get_point().x())
      return true;
    else if (first.get_point().x() == second.get_point().x() && 
             first.get_point().y() < second.get_point().y())
      return true;
    else if (first.get_point().x() == second.get_point().x() && 
             first.get_point().y() == second.get_point().y() && 
             first.get_point().z() <= second.get_point().z())
      return true;
    return false;
  }
};
  
template <typename Point_2, typename Rational_segment_3>
class Compare_points_on_plane
{
  typedef Point_and_segment_pair<Point_2,Rational_segment_3>
  Point_on_plane_and_line_pair;

public:
  bool operator()(const Point_on_plane_and_line_pair& first,
                  const Point_on_plane_and_line_pair& second)
  {
    if (first.get_point().x() < second.get_point().x())
      return true;
    else if (first.get_point().x() == second.get_point().x() &&
             first.get_point().y() <= second.get_point().y())
      return true;
    return false;
  }
};
  
template <typename Point_on_sphere_2, typename Rational_segment_3>
class Compare_points_on_sphere
{
  typedef Point_and_segment_pair<Point_on_sphere_2, Rational_segment_3>
  Point_on_sphere_and_line_pair;

public:
  bool operator()(const Point_on_sphere_and_line_pair& first,
                  const Point_on_sphere_and_line_pair& second)
  {
    typedef typename Point_on_sphere_2::Algebraic Algebraic;
    Algebraic first_len =
      CGAL::sqrt(Algebraic(first.get_point().dx() * first.get_point().dx() +
                           first.get_point().dy() * first.get_point().dy() +
                           first.get_point().dz() * first.get_point().dz()));
            
    Algebraic second_len =
      CGAL::sqrt(Algebraic(second.get_point().dx() * second.get_point().dx() +
                           second.get_point().dy() * second.get_point().dy() +
                           second.get_point().dz() * second.get_point().dz()));
            
    Algebraic first_dx = 
      first.get_point().dx() / first_len;

    Algebraic second_dx = 
      second.get_point().dx() / second_len;          
                        
    if (first_dx < second_dx)
      return true;
            
    typename Point_on_sphere_2::Algebraic first_dy = 
      first.get_point().dy() / first_len;

    typename Point_on_sphere_2::Algebraic second_dy = 
      second.get_point().dy() / second_len;          
            
    if (first_dx == second_dx && first_dy < second_dy)
      return true;

    typename Point_on_sphere_2::Algebraic first_dz = 
      first.get_point().dz() / first_len;

    typename Point_on_sphere_2::Algebraic second_dz = 
      second.get_point().dz() / second_len;          

    if (first_dx == second_dx && first_dy == second_dy && 
        first_dz <= second_dz)
      return true;
            
    return false;
  }
};
  
} //namespace CGAL

#endif /*ISOLATED_POINTS_H */
