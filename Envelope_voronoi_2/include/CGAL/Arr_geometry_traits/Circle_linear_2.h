// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL: 
// $Id: 
// 
//
// Author(s)     : Ron Wein        <wein@post.tau.ac.il>

#ifndef CGAL_CIRCULE_LINEAR_2_H
#define CGAL_CIRCULE_LINEAR_2_H

/*! \file
 * Header file for the _Circle_linear_2<Kernel, Filter> class.
 */
#include <CGAL/Arr_geometry_traits/Circle_segment_2.h>
#include <CGAL/Arr_geometry_traits/One_root_number.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Handle_for.h>

#include <list>
#include <map>
#include <ostream>

namespace CGAL {


/*! \class
 * Representation of a circle, a circular arc, or a linear segment, a linear
 * ray or a line segment.
 */
template <class Kernel_, bool Filter_>
class _Circle_linear_2 : protected _Circle_segment_2< Kernel_, Filter_ >
{
public:

  typedef Kernel_                                          Kernel;
  typedef typename Kernel::FT                              NT;
  typedef _One_root_point_2<NT, Filter_>                   Point_2;
  typedef typename Kernel::Circle_2                        Circle_2;
  typedef typename Kernel::Segment_2                       Segment_2;
  typedef typename Kernel::Ray_2                           Ray_2;
  typedef typename Kernel::Line_2                          Line_2;
  typedef _Circle_segment_2< Kernel_, Filter_>             Base;
  typedef _Circle_linear_2< Kernel_, Filter_>              Self;

private:

  typedef typename Point_2::CoordNT                        CoordNT;

  // Data members:
  bool          _has_source;  // True if we have source pt. (ray, segment, arc).
  bool          _has_target;  // True if we have source pt. (ray, segment, arc).
public:

  /*! Default constructor. */
  _Circle_linear_2 () : Base (),
    _has_source (false),
    _has_target (false)
  {}

  /*!
   * Constructor from a line segment.
   * \param seg The segment.
   */
  _Circle_linear_2 (const Segment_2& seg) : Base(seg),
    _has_source (true),
    _has_target (true)
  {}

  /* Constructor from a ray.
   * \param ray The ray.
   */
  _Circle_linear_2 (const Ray_2& ray) :
    Base (),
    _has_source (true),
    _has_target (false)
  {
    Kernel k;
    this->_line = k.construct_line_2_object() (ray);
    this->_is_full  = false;
    this->_has_radius = false;
    this->_source =  Point_2(ray.source().x(), ray.source().y());
    this->_orient = COLLINEAR;
  }

  /* Constructor from a line.
   * \param line The line.
   */
  _Circle_linear_2 (const Line_2& line) :
    Base (),
    _has_source (false),
    _has_target (false)
  {
    this->_line = line;
    this->_is_full = true;
    this->_has_radius = false;
    this->_orient = COLLINEAR;
  }


  /*!
  * Constructor from of a line segment.
   * \param ps The source point.
   * \param pt The target point.
   */
  _Circle_linear_2 (const typename Kernel::Point_2& ps,
                  const typename Kernel::Point_2& pt) : 
    Base (ps, pt),
    _has_source (true),
    _has_target (true)  
  {}

  /*!
   * Constructor of a segment, given a supporting line and two endpoints,
   * which need not necessarily have rational coordinates.
   * \param line The supporting line.
   * \param source The source point.
   * \param target The target point.
   * \pre Both endpoints lie on the supporting line.
   */
  _Circle_linear_2 (const Line_2& line,
                  const Point_2& source, const Point_2& target) :
    Base (line, source, target),
    _has_source (true),
    _has_target (true)
  {}

  /*!
   * Constructor from a circle.
   * \param circ The circle.
   */
  _Circle_linear_2 (const Circle_2& circ) :
    Base (circ),
    _has_source (false),
    _has_target (false)
  {}

  /*!
   * Constructor from a circle.
   * \param c The circle center.
   * \param r The radius.
   * \param orient The orientation of the circle.
   */
  _Circle_linear_2 (const typename Kernel::Point_2& c,
                  const NT& r,
                  Orientation orient = COUNTERCLOCKWISE) :
    Base (c, r, orient),
    _has_source (false),
    _has_target (false)
  {}

  /*!
   * Constructor of a circular arc, given a supporting circle and two
   * endpoints, which need not necessarily have rational coordinates.
   * The orientation of the circle determines the orientation of the arc.
   * \param circ The supporting circle.
   * \param source The source point.
   * \param target The target point.
   * \pre Both endpoints lie on the supporting circle.
   */
  _Circle_linear_2 (const Circle_2& circ,
                    const Point_2& source, const Point_2& target) :
    Base (circ, source, target),
    _has_source (true),
    _has_target (true)
  {}

  /*!
   * Constructor of a circular arc, given a supporting circle and two
   * endpoints, which need not necessarily have rational coordinates.
   * \param c The circle center.
   * \param r The radius.
   * \param orient The orientation of the circle.
   * \param source The source point.
   * \param target The target point.
   * \pre Both endpoints lie on the supporting circle.
   */
  _Circle_linear_2 (const typename Kernel::Point_2& c,
                    const NT& r, Orientation orient,
                    const Point_2& source, const Point_2& target) :
    Base (c, r, orient, source, target),
    _has_source (true),
    _has_target (true)
  {}

  /*!
   * Constructor of a circular arc, from the given three points, in case of
   * three collinear points, a segment will be constructed.            
   * \param p1 The arc source.
   * \param p2 A point in the interior of the arc.
   * \param p3 The arc target.
   * \pre p1 and p3 are not equal.
   */
   _Circle_linear_2 (const typename Kernel::Point_2& p1,
                     const typename Kernel::Point_2& p2,
                     const typename Kernel::Point_2& p3) :
    Base (p1, p2, p3),
    _has_source (true),
    _has_target (true)
  {}

  /*!
   * Get the orientation of the curve. 
   * \return COLLINEAR in case of a line segment,
   *         CLOCKWISE or COUNTERCLOCKWISE for circular curves.
   */
  inline Orientation orientation () const
  {
    return Base::orientation();
  }

  /*! Check if the curve is has a source point. */
  bool has_source () const
  {
    return _has_source;
  }

  /*! Check if the curve is has a target point. */
  bool has_target () const
  {
    return _has_target;
  }

  /*!
   * Get the source point.
   * \pre The curve has a source point.
   */
  const Point_2& source () const
  {
    CGAL_precondition (_has_source);
    return this->_source;
  }

  /*!
   * Get the target point.
   * \pre The curve is not a full circle.
   */
  const Point_2& target () const
  {
    CGAL_precondition (_has_target);
    return this->_target;
  }

  /*!
   * Get the supporting line.
   * \pre The curve orientation is COLLINEAR.
   */
  const Line_2& supporting_line () const
  {
    return Base::supporting_line();
  }

  /*!
   * Get the supporting circle.
   * \pre The curve orientation is not COLLINEAR.
   */
  const Circle_2& supporting_circle () const
  {
    return Base::supporting_circle();
  }

  /*! Check if the curve is a full circle. */
  bool is_full () const
  {
    return Base::is_full();
  }
  
  /*!
   * Get the vertical tangency points the arc contains.
   * \param vpts Output: The vertical tagnecy points.
   * \pre The curve is circular.
   * \return The number of points (0, 1, or 2).
   */
  unsigned int vertical_tangency_points (Point_2 *vpts) const
  {
    return Base::vertical_tangency_points(vpts);
  }
};


/*!
 * Exporter for line segments and circular arcs.
 */
template <class Kernel, bool Filter>
std::ostream& 
operator<< (std::ostream& os, 
            const _Circle_linear_2<Kernel, Filter>& c)
{
  if (c.orientation() == COLLINEAR)
  {
    if (c.is_full())
      os << "line: " << c.supporting_line() << std::endl;
    else if (c.has_source() && c.has_target())
      os<< "segment: " << c.source() << " -> " << c.target() << std::endl;
    else
    {
      if (c.has_source())
        os<< "ray: " << c.source() << " supporting line " 
          << c.supporting_line() << std::endl;
      else
        os<< "ray: " << c.target() << " supporting line " 
          << c.supporting_line() << std::endl;
    }
  }
  else
  {
    if(!c.is_full())
    {
      os << "circular arc: " << c.supporting_circle() << ' '
         << c.source() << " -> " << c.target() << std::endl;
    }
    else
    {
      os << "circular arc: " << c.supporting_circle()<<std::endl;
    }
  }

  return (os);
}

/*! \class
 * Representation of an x-monotone circular arc.
 */
template <class Kernel_, bool Filter_>
class _X_monotone_circle_linear_2 
  : protected _X_monotone_circle_segment_2< Kernel_, Filter_ >
{
public:

  typedef Kernel_                                          Kernel;
  typedef _X_monotone_circle_linear_2<Kernel, Filter_>    Self;
  typedef typename Kernel::FT                              NT;
  typedef _One_root_point_2<NT, Filter_>                   Point_2;
  typedef typename Kernel::Circle_2                        Circle_2;
  typedef typename Kernel::Ray_2                           Ray_2;
  typedef typename Kernel::Line_2                          Line_2;
  typedef typename Kernel::Direction_2                     Direction_2;
  typedef typename Point_2::CoordNT                        CoordNT;

  typedef _X_monotone_circle_segment_2< Kernel_, Filter_ > Base;

  // Type definition for the intersection points mapping.
  typedef typename Base::Curve_id_pair            Curve_id_pair;
  typedef typename Base::Intersection_point_2     Intersection_point_2;
  typedef typename Base::Intersection_list        Intersection_list;


  typedef typename Base::Intersection_map             Intersection_map;
  typedef typename Base::Intersection_map_entry       Intersection_map_entry;
  typedef typename Base::Intersection_map_iterator    Intersection_map_iterator;

protected:
  bool          _has_source;
  bool          _has_target;

public:

  /*!
   * Default constructor.
   */
  _X_monotone_circle_linear_2 () : Base()
  {}

  /*!
   * Construct an arc from a line segment.
   * \param line The supporting line.
   * \param source The source point.
   * \param target The target point.
   */
  _X_monotone_circle_linear_2 (const Line_2& line,
                                const Point_2& source, const Point_2& target,
                                unsigned int index = 0) :
  Base(line, source, target, index), 
  _has_source(true), 
  _has_target(true)
  {}

  /*!
   * Construct an arc from a ray.
   * \param line The supporting line. The ray is in the same direction of 
   *    the line.
   * \param source The source of the ray.
   */
  _X_monotone_circle_linear_2 (const Line_2& line,
                               const Point_2& source,
                               unsigned int index = 0) :
  Base(), _has_source(true), _has_target(false)
  {
    const Direction_2 pos_dir = Direction_2(0, 1);
    const Direction_2 neg_dir = Direction_2(0, -1);

    Kernel k;
    this->_first = line.a();
    this->_second = line.b();
    this->_third = line.c();
    this->_source = source;
    this->_info = (index << Base::INDEX_SHIFT_BITS);
    
    Direction_2 dir = k.construct_direction_2_object() (line);
    
    if (CGAL::sign(this->_second) == ZERO)
    {
      this->_info = (this->_info | Base::IS_VERTICAL_SEGMENT_MASK);
      if (k.equal_2_object()(dir, pos_dir))
      {
        this->_info = (this->_info | Base::IS_DIRECTED_RIGHT_MASK);
      }
      else
      {
        CGAL_assertion(k.equal_2_object()(dir, neg_dir));
      }
    }
    else
    {
      if(k.counterclockwise_in_between_2_object()(dir, neg_dir, pos_dir))
        this->_info = (this->_info | Base::IS_DIRECTED_RIGHT_MASK);
    }
  }

  /*!
   * Construct an arc from a line.
   * \param line The line.
   */
  _X_monotone_circle_linear_2 (const Line_2& line,
                               unsigned int index = 0) :
  Base(), _has_source(false), _has_target(false)
  {
    const Direction_2 pos_dir = Direction_2(0, 1);
    const Direction_2 neg_dir = Direction_2(0, -1);

    Kernel k;
    this->_first = line.a();
    this->_second = line.b();
    this->_third = line.c();
    this->_info = (index << Base::INDEX_SHIFT_BITS);
    
    Direction_2 dir = k.construct_direction_2_object() (line);
    
    if (CGAL::sign(this->_second) == ZERO)
    {
      this->_info = (this->_info | Base::IS_VERTICAL_SEGMENT_MASK);
      if (k.equal_2_object()(dir, pos_dir))
      {
        this->_info = (this->_info | Base::IS_DIRECTED_RIGHT_MASK);
      }
      else
      {
        CGAL_assertion(k.equal_2_object()(dir, neg_dir));
      }
    }
    else
    {
      if(k.counterclockwise_in_between_2_object()(dir, neg_dir, pos_dir))
        this->_info = (this->_info | Base::IS_DIRECTED_RIGHT_MASK);
    }

  }

  /*!
   * Construct a segment arc from two kernel points
   * \param source the source point.
   * \ param target the target point.
   * \pre source and target are not equal.
   */
  _X_monotone_circle_linear_2 (const typename Kernel::Point_2& source,
                                const typename Kernel::Point_2& target) :
  Base(source, target),
  _has_source(true),
  _has_target(true)
  {  }
     

  /*! 
   * Construct a circular arc.
   * \param line The supporting line.
   * \param source The source point.
   * \param target The target point.
   * \param orient The orientation of the arc.
   */
  _X_monotone_circle_linear_2 (const Circle_2& circ,
                                const Point_2& source, const Point_2& target,
                                Orientation orient,
                                unsigned int index = 0) :
    Base(circ, source, target, orient, index),
    _has_source(true),
    _has_target(true)
    {}

  /*! Check if the curve has source point. */
  inline bool has_source () const
  {
    return _has_source;
  }

  /*! Check if the curve has target point. */
  inline bool has_target () const
  {
    return _has_target;
  }

  /*! Get the source point. */
  inline const Point_2& source () const
  {
    CGAL_precondition(has_source());
    return (this->_source);
  }

  /*! Get the target point. */
  inline const Point_2& target () const
  {
    CGAL_precondition(this->has_target());
    return (this->_target);
  }

  /*! Check if the curve has left endpoint. */
  inline bool has_left () const
  {
    return (((this->_info & this->IS_DIRECTED_RIGHT_MASK) != 0) ? 
            this->_has_source : this->_has_target);
  }

  /*! Check if the curve has right endpoint. */
  inline bool has_right () const
  {
    return (((this->_info & this->IS_DIRECTED_RIGHT_MASK) != 0) ? 
            this->_has_target : this->_has_source);
  }

  /*! Get the left endpoint of the arc. */
  inline const Point_2& left () const
  {
    CGAL_precondition(has_left());
    return (((this->_info & this->IS_DIRECTED_RIGHT_MASK) != 0) ? 
            this->_source : this->_target);
  }

  /*! Get the right endpoint of the arc. */
  inline const Point_2& right () const
  {
    CGAL_precondition(has_right());
    return (((this->_info & this->IS_DIRECTED_RIGHT_MASK) != 0) ? 
            this->_target : this->_source);
  }

  /*! Set the left endpoint of the arc. */
  inline void set_left (const Point_2& p, bool has_left = true)
  {
    if ((this->_info & this->IS_DIRECTED_RIGHT_MASK) != 0)
    {
      this->_source = p;
      _has_source = has_left;
    }
    else
    { 
      this->_target = p;
      _has_target = has_left;
    }
  }

  /*! Set the left endpoint of the arc. */
  inline void set_right (const Point_2& p, bool has_right)
  {
    if ((this->_info & this->IS_DIRECTED_RIGHT_MASK) != 0)
    {
      this->_target = p;
      _has_target = has_right;
    }
    else
    { 
      this->_source = p;
      _has_source = has_right;
    }
  }

  /*! Check if the arc is linear. */
  inline bool is_linear () const
  {
    return Base::is_linear();
  }

  /*! Check if the arc is circular. */
  inline bool is_circular () const
  {
    return Base::is_circular();
  }

  /*! return true iff the arc is directed right lexicoraphically. */
  bool is_directed_right() const
  {
    return Base::is_directed_right();
  }

  
  /*!
   * Get the supporting line.
   * \pre The curve orientation is COLLINEAR.
   */
  const Line_2 supporting_line () const
  {
    return Base::supporting_line();
  }

  /*!
   * Get the supporting circle. 
   * \pre The arc is circular.
   */
  Circle_2 supporting_circle () const
  {
    return Base::supporting_circle();
  }

  /*! Get the coefficient of y in the equation of the supporting line. */
  inline const NT& a () const
  {
    return Base::a();
  }

  /*! Get the coefficient of y in the equation of the supporting line. */
  inline const NT& b () const
  {
    return Base::b();
  }

  /*! Get the coefficient of y in the equation of the supporting line. */
  inline const NT& c () const
  {
    return Base::c();
  }

  /*!
   * Check whether the given point is in the x-range of the arc.
   */
  bool is_in_x_range (const Point_2& p) const
  {
    // full line
    if (has_source() == false && this->has_target() == false)
    {
      if (is_vertical())
        return CGAL::compare (p.x(), right().x()) == EQUAL;

      // this is a full line
      return true;
    }
    
    if (has_left())
    {
      Comparison_result res = CGAL::compare (p.x(), left().x());
      if (res == SMALLER)
        return (false);
      else if (res == EQUAL)
        return (true);
    }

    if (has_right())
      return (CGAL::compare (p.x(), right().x()) != LARGER);
    
    return true;
  }

  /*! Check if the arc is a vertical segment. */
  inline bool is_vertical () const
  {
    return ((this->_info & this->IS_VERTICAL_SEGMENT_MASK) != 0);
  }

  /*! Get the orientation of the arc. */ 
  inline Orientation orientation() const
  {
    return Base::orientation();
  }

  /*!
   * Check the position of a given point with respect to the arc.
   */
  Comparison_result point_position (const Point_2& p) const
  {
    if (this->is_linear())
      return (_line_point_position (p));
    else
      return (_circ_point_position (p));
  }

  /*!
   * Compare the two arcs to the right of their intersection point.
   */
  Comparison_result compare_to_right (const Self& cv, const Point_2& p) const
  {
    return Base::compare_to_right(cv, p);
  }

  /*!
   * Compare the two arcs to the left of their intersection point.
   */
  Comparison_result compare_to_left (const Self& cv, const Point_2& p) const
  {
    return Base::compare_to_left(cv, p);
  }

  /*!
   * Check if the two curves are equal.
   */
  bool equals (const Self& cv) const
  {
    if (! this->has_same_supporting_curve (cv))
      return (false);
    
    if (has_left() != cv.has_left())
      return false;
    if (has_right() != cv.has_right())
      return false;

    if (this->is_linear())
    {
      bool res;
      if (has_left())
      {
        res = left().equals(cv.left());
        if (res == false)
          return false;
      }
      if (has_right())
        return right().equals(cv.right());
      return true;
    }

    // Once again, opposite circular arcs are considered to be equal:
    return ((this->orientation() == cv.orientation() &&
             this->_source.equals (cv._source) && 
             this->_target.equals (cv._target)) ||
            (this->orientation() != cv.orientation() &&
             this->_source.equals (cv._target) && 
             this->_target.equals (cv._source)));
  }

  /*!
   * Split the curve at a given point into two sub-arcs.
   */
  void split (const Point_2& p, Self& c1, Self& c2) const
  {
    // split doesn't set our variables, so we use our = first.
    c1 = *this;
    c2 = *this;
    Base::split(p, c1, c2);

    // Change the endpoint, such that c1 lies to the right of c2:
    if (this->is_directed_right())
    {
      c1._has_target = true;
      c2._has_source = true;
    }
    else
    {
      c1._has_source = true;
      c2._has_target = true;
    }
  }

  /*!
   * Compute the intersections between the two arcs or segments.
   */
  template <class OutputIterator>
  OutputIterator intersect (const Self& cv, OutputIterator oi,
                            Intersection_map *inter_map = NULL) const
  {
    // First check whether the two arcs have the same supporting curve.
    if (has_same_supporting_curve (cv))
    {
      // Check for overlaps between the two arcs.
      Self    overlap;

      if (_compute_overlap (cv, overlap))
      {
        // There can be just a single overlap between two x-monotone arcs:
        *oi = CGAL::make_object (overlap);
        ++oi;
        return (oi);
      }

      // In case there is not overlap and the supporting curves are the same,
      // there cannot be any intersection points, unless the two arcs share
      // a common end point.
      // Note that in this case we do not define the multiplicity of the
      // intersection points we report.
      unsigned int  mult = 0;
      if (source().equals (cv.source()))
      {
        *oi = CGAL::make_object (std::make_pair (left(), mult));
        ++oi;
      }

      if (source().equals (cv.target()))
      {
        *oi = CGAL::make_object (std::make_pair (right(), mult));
        ++oi;
      }
      if (target().equals (cv.source()))
      {
        *oi = CGAL::make_object (std::make_pair (left(), mult));
        ++oi;
      }

      if (target().equals (cv.target()))
      {
        *oi = CGAL::make_object (std::make_pair (right(), mult));
        ++oi;
      }

      return (oi);
    }

    // Before computing the intersection points between the two supporting
    // curves, check if their intersection has already been computed and
    // cached.
    Curve_id_pair                id_pair;
    Intersection_map_iterator    map_iter;
    Intersection_list            inter_list;
    bool                         invalid_ids = false;

    if (inter_map != NULL && this->_index() != 0 && cv._index() != 0)
    {
      if (this->_index() < cv._index())
        id_pair = Curve_id_pair (this->_index(), cv._index());
      else
        id_pair = Curve_id_pair (cv._index(), this->_index());
      
      map_iter = inter_map->find (id_pair);
    }
    else
    {
      // In case one of the IDs is invalid, we do not look in the map neither
      // we cache the results.
      if (inter_map != NULL)
        map_iter = inter_map->end();
      invalid_ids = true;
    }

    if (inter_map == NULL || map_iter == inter_map->end())
    {
      // Compute the intersections points between the two supporting curves.
      if (this->is_linear())
      {
        if (cv.is_linear())
          _lines_intersect (cv, inter_list);
        else
          cv._circ_line_intersect (*this, inter_list);
      }
      else
      {
        if (cv.is_linear())
          _circ_line_intersect (cv, inter_list);
        else
          _circs_intersect (cv, inter_list);
      }

      // Cache the result.
      if (! invalid_ids)
        (*inter_map)[id_pair] = inter_list;
    }
    else
    {
      // Obtain the precomputed intersection points from the map.
      inter_list = (*map_iter).second;
    }

    // Report only the intersection points that lie on both arcs.
    typename Intersection_list::const_iterator   iter;

    for (iter = inter_list.begin(); iter != inter_list.end(); ++iter)
    {
      if (this->_is_between_endpoints (iter->first) &&
          cv._is_between_endpoints (iter->first))
      {
        *oi = CGAL::make_object (*iter);
        ++oi;
      }
    }

    return (oi);
  }

  /*!
   * Check whether it is possible to merge our arc with the given arc.
   */
  bool can_merge_with (const Self& cv) const
  {
    // In order to merge the two arcs, they should have the same supporting
    // curve.
    if (! this->has_same_supporting_curve (cv))
      return (false);

    // Check if the left endpoint of one curve is the right endpoint of the
    // other.
    return (has_right() && cv.has_left() && right().equals (cv.left()) ||
            has_left() && cv.has_right() && left().equals (cv.right()));
  }

  /*!
   * Merge our arc with the given arc.
   * \pre The two arcs are mergeable.
   */
  void merge (const Self& cv)
  {
    CGAL_precondition (this->can_merge_with (cv));

    // Check if we should extend the arc to the left or to the right.
    if (has_right() && right().equals (cv.left()))
    {
      // Extend the arc to the right.
      this->set_right(cv.right(), cv.has_right());
    }
    else
    {
      CGAL_precondition (has_left() && left().equals (cv.right()));

      // Extend the arc to the left.
      this->set_left(cv.left(), cv.has_left());
    }

    return;
  }

  /*! construct an opposite arc. */
  Self construct_opposite() const
  {
    Self opp_cv = Base::construct_opposite();
    opp_cv._has_source = this->_has_target;
    opp_cv._has_target = this->_has_source;

    return (opp_cv);
  }

protected:

  /// \name Auxiliary functions for the point_position predicate.
  //@{

  /*!
   * Check the position of a given point with respect to a line segment.
   */
  Comparison_result _line_point_position (const Point_2& p) const
  {
    // Check if we have a vertical segment.

    CGAL_precondition (this->is_in_x_range(p));

    Comparison_result    res;

    if (this->is_vertical())
    {
      // left() is the lower endpoint:
      res = CGAL::compare (p.y(), left().y());

      if (res != LARGER)
        return (res);

      // left() is the upper endpoint:
      res = CGAL::compare (p.y(), right().y());

      if (res != SMALLER)
        return (res);

      // p lies in the interior of the vertical segment:
      return (EQUAL);
    }

    // Compute the y-coordinate of the vertical projection of p onto the
    // supporting line.
    const CoordNT        y_proj = (this->a() * p.x() + this->c()) / 
      (-this->b());

    return (CGAL::compare (p.y(), y_proj));
  }

  /*!
   * Check the position of a given point with respect to a circular arc.
   */
  Comparison_result _circ_point_position (const Point_2& p) const
  {

    Comparison_result   c_res = CGAL::compare (p.y(), this->y0());

    if (this->_is_upper())
    {
      // Check if p lies below the "equator" (while the arc lies above it):
      if (c_res == SMALLER)
        return (SMALLER);
    }
    else
    {
      // Check if p lies above the "equator" (while the arc lies below it):
      if (c_res == LARGER)
        return (LARGER);
    }

    // Check if p lies inside the supporting circle, namely we have to check
    // whether (p.x() - x0)^2 + (p.y() - y0)^2 < r^2:
    Comparison_result   res =
                         CGAL::compare (CGAL::square (p.x() - this->x0()),
                                        this->sqr_r() - 
                                        CGAL::square (p.y() - this->y0()));

    if (res == EQUAL)
      // p lies on the circle:
      return (EQUAL);

    if (this->_is_upper())
    {
      // If p is inside the circle, it lies below the upper arc:
      return (res);
    }
    else
    {
      // If p is inside the circle, it lies above the lower arc:
      return (res == SMALLER ? LARGER : SMALLER);
    }
  }
  //@}

  /// \name Auxiliary functions for computing intersections.
  //@{

  /*!
   * Check if the given point lies on the arc.
   * \pre p lies on the supporting curve.
   */
  bool _is_between_endpoints (const Point_2& p) const
  {
    if (this->is_linear())
    {
      if (this->is_vertical())
      {
        // Check if the point is in the y-range of the arc.
        // Note that left() is the lower endpoint and right() is the upper
        // endpoint of the segment in this case. 
        Comparison_result    res = CGAL::compare (p.y(), left().y());

        if (has_left() && res == SMALLER)
          return (false);
        else if (res == EQUAL)
          return (true);

        if (has_right())
          return (CGAL::compare (p.y(), right().y()) != LARGER);
        
        return true;
      }

      // For non-vertical segments, it is sufficient to check if the point
      // is in the x-range of the arc.
      return (this->is_in_x_range (p));
    }

    // The supporting curve is a circle:
    // Check whether p lies on the upper or on the lower part of the circle.
    Comparison_result   c_res = CGAL::compare (p.y(), this->y0());

    if ((this->_is_upper() && c_res == SMALLER) ||
        (! this->_is_upper() && c_res == LARGER))
    {
      // The point lies on the other half of the circle:
      return (false);
    }

    // Check if the point is in the x-range of the arc.
    return (this->is_in_x_range (p));
  }

  /*!
   * Check if the given point lies in the interior of the arc.
   * \pre p lies on the supporting curve.
   */
  bool _is_strictly_between_endpoints (const Point_2& p) const
  {
    if (this->_has_source && p.equals (this->_source))
      return (false);
    if (this->_has_target && p.equals (this->_target))
      return false;
    
    return (_is_between_endpoints (p));
  }

  /*!
   * Compute the overlap with a given arc having the same supporting curve.
   * \param cv The given arc.
   * \param overlap Output: The overlapping arc (if any).
   * \pre The given arc has the same supporting curve as this.
   * \return Whether we found an overlap.
   */
  bool _compute_overlap (const Self& cv, Self& overlap) const
  {
    // Check if the two arcs are identical.
    if (this->is_linear())
    {
      // In case of line segments we can swap the source and target:
      if (equals(cv))
      {
        overlap = cv;
        return (true);
      }
    }
    else
    {
      if ((this->orientation() == cv.orientation() &&
           this->_source.equals (cv._source) && 
           this->_target.equals (cv._target)) ||
          (this->orientation() != cv.orientation() &&
           this->_source.equals (cv._target) && 
           this->_target.equals (cv._source)))
      {
        overlap = cv;
        return (true);
      }
    }

    Kernel kernel;
    typename Kernel::Compare_xy_2  compare_xy = kernel.compare_xy_2_object();
    bool cv_left = cv.has_left() ? 
      _is_strictly_between_endpoints (cv.left()) : (!has_left());
    bool cv_right = cv.has_right() ? 
      _is_strictly_between_endpoints (cv.right()) : (!has_right());

    if (cv_left)
    {
      if (cv_right)
      {
        // Case 1 - *this:     +----------->
        //             cv:       +=====>
        overlap = cv;
        return (true);
      }
      else
      {
        // Case 2 - *this:     +----------->
        //             cv:               +=====>
        overlap = *this;
        overlap.set_left(cv.left(), cv.has_left());

        return (true);
      }
    }
    else if (cv_right)
    {
      // Case 3 - *this:     +----------->
      //             cv:   +=====>
      overlap = *this;
      overlap.set_right(cv.right(), cv.has_right());

      return (true);
    }
    else if (cv._is_between_endpoints (this->_source) &&
             cv._is_between_endpoints (this->_target) &&
             (cv._is_strictly_between_endpoints (this->_source) ||
              cv._is_strictly_between_endpoints (this->_target)))
    {
      // Case 4 - *this:     +----------->
      //             cv:   +================>
      overlap = *this;
      return (true);
    }

    // If we reached here, there are no overlaps:
    return (false);
  }

  public:
  template <class OutputIterator>
  void approximate(OutputIterator oi, unsigned int n) const
  {
    const double x_left = CGAL::to_double(this->source().x());
    const double y_left = CGAL::to_double(this->source().y());

    const double x_right = CGAL::to_double(this->target().x());
    const double y_right = CGAL::to_double(this->target().y());
    if(this->is_linear())
    {
      *oi = std::make_pair(x_left, y_left);
      ++oi;

      *oi = std::make_pair(x_right, y_right);
      ++oi;
      return;
    }

    // Otherwise, sample (n - 1) equally-spaced points in between.
    const double  app_xcenter = CGAL::to_double (this->_first);
    const double  app_ycenter = CGAL::to_double (this->_second);
    const double  app_sqr_rad = CGAL::to_double (this->_third);
   
    const double  x_jump = (x_right - x_left) / n;
    double        x, y;
    double        disc;
    unsigned int        i;

    const bool is_up = this->_is_upper();
    *oi = std::make_pair (x_left, y_left);   // The left point.
    ++oi;
    for (i = 1; i < n; i++)
    {
      x = x_left + x_jump*i;
      disc = app_sqr_rad - CGAL::square(x - app_xcenter);
      CGAL_precondition(disc >= 0);
      if(is_up)
        y = app_ycenter + std::sqrt(disc);
      else
        y = app_ycenter - std::sqrt(disc);

      *oi = std::make_pair(x, y);
      ++oi;
    }
    *oi = std::make_pair(x_right, y_right);   // The right point.
    ++oi;
  }

  //@}
};

/*!
 * Exporter for circular arcs (or line segments).
 */
template <class Kernel, bool Filter>
std::ostream& 
operator<< (std::ostream& os, 
            const _X_monotone_circle_linear_2<Kernel, Filter> & arc)
{
  typedef typename _X_monotone_circle_linear_2<Kernel, Filter>::Base Base;

  if (! arc.is_linear())
    os << "(" << arc.supporting_circle() << ") ";

  if (arc.has_source())
    os << "[";
  else
    os << "<";
  os << ((Base*)&arc)->source() << " --> " << ((Base*)&arc)->target();
  
  if (arc.has_target())
    os << "]";
  else
    os << ">";

  os << std::endl;
  return (os);
}

} //namespace CGAL

#endif
