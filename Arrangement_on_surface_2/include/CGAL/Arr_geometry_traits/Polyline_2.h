// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_POLYLINE_2_H
#define CGAL_ARR_POLYLINE_2_H

/*! \file
 * Header file for the polyline classes used by the
 * Arr_polycurve_basic_traits_2, Arr_polycurve_traits_2, and
 * Arr_polyline_traits_2 classes.
 */

#if (defined __GNUC__)
  #warning Polyline_2.h is DEPRECATED, please include Polycurve_2.h instead
#elif (defined _MSC_VER)
  #pragma message("Polyline_2.h is DEPRECATED, please include Polycurve_2.h instead")
#endif

#include <CGAL/Arr_geometry_traits/Polycurve_2.h>

namespace CGAL {

namespace polyline {

template <typename SubcurveType_2, typename PointType_2>
class Polyline_2 : public internal::Polycurve_2<SubcurveType_2, PointType_2> {
public:
  typedef SubcurveType_2                                Subcurve_type_2;
  typedef PointType_2                                   Point_type_2;

  typedef internal::Polycurve_2<Subcurve_type_2, Point_type_2> Base;

  typedef typename Base::Subcurve_type_2                Segment_type_2;
  typedef typename Base::size_type                      Segments_size_type;
  typedef typename Base::Subcurve_iterator              Segment_iterator;
  typedef typename Base::Subcurve_const_iterator        Segment_const_iterator;
  typedef typename Base::Subcurve_const_reverse_iterator
    Segment_const_reverse_iterator;

  /*! Construct default. */
  Polyline_2() : Base() {}

  /*! Construct from a subcurve. */
  Polyline_2(const Subcurve_type_2& subcurve) : Base(subcurve) {}

  /*! Construct from a range. */
  template <typename InputIterator>
  Polyline_2(InputIterator begin, InputIterator end) : Base(begin, end) {}

  /*! Obtain an iterator for the polycurve subcurves. */
  Segment_const_iterator begin_segments() const
  { return this->begin_subcurves(); }

  /*! Obtain a past-the-end iterator for the polycurve subcurves. */
  Segment_const_iterator end_segments() const
  { return this->end_subcurves(); }

  /*! Obtain a reverse iterator for the polycurve subcurves. */
  Segment_const_reverse_iterator rbegin_segments() const
  { return this->rbegin_subcurves(); }

  /*! Obtain a reverse past-the-end iterator for the polycurve points. */
  Segment_const_reverse_iterator rend_segments() const
  { this->returnrend_subcurves() ; }

  /*! Obtain the number of subcurves that comprise the poyline.
   * \return The number of subcurves.
   */
  Segments_size_type number_of_segments() const
  { return this->number_of_subcurves(); }
};

template <typename SubcurveType_2, typename PointType_2>
class X_monotone_polyline_2 :
    public internal::X_monotone_polycurve_2<SubcurveType_2, PointType_2> {
public:
  typedef SubcurveType_2                                Subcurve_type_2;
  typedef PointType_2                                   Point_type_2;

  typedef internal::X_monotone_polycurve_2<Subcurve_type_2, Point_type_2> Base;

  typedef typename Base::Subcurve_type_2                Segment_type_2;
  typedef typename Base::size_type                      Segments_size_type;
  typedef typename Base::Subcurve_iterator              Segment_iterator;
  typedef typename Base::Subcurve_const_iterator        Segment_const_iterator;
  typedef typename Base::Subcurve_const_reverse_iterator
    Segment_const_reverse_iterator;

  /*! Construct default. */
  X_monotone_polyline_2() : Base() {}

  /*! Construct from a subcurve. */
  X_monotone_polyline_2(Subcurve_type_2 seg) : Base(seg) {}

  /*! Construct from a range.
   * Similar to the constructor of a general polycurve.
   * Like in the case of general polycurve, for the sake of backwards
   * compatibility we have to keep an implementation of construction
   * from a range of points. DO NOT USE THIS CONSTRUCTION.
   */
  template <typename InputIterator>
  X_monotone_polyline_2(InputIterator begin, InputIterator end) :
    Base(begin, end)
  {}

  /*! Obtain an iterator for the polycurve subcurves. */
  Segment_const_iterator begin_segments() const
  { return this->begin_subcurves(); }

  /*! Obtain a past-the-end iterator for the polycurve subcurves. */
  Segment_const_iterator end_segments() const
  { return this->end_subcurves(); }

  /*! Obtain a reverse iterator for the polycurve subcurves. */
  Segment_const_reverse_iterator rbegin_segments() const
  { return this->rbegin_subcurves(); }

  /*! Obtain a reverse past-the-end iterator for the polycurve points. */
  Segment_const_reverse_iterator rend_segments() const
  { this->returnrend_subcurves() ; }

  /*! Obtain the number of subcurves that comprise the poyline.
   * \return The number of subcurves.
   */
  Segments_size_type number_of_segments() const
  { return this->number_of_subcurves(); }
};

} // namespace polyline
} //namespace CGAL

#endif
