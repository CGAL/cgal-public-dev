// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_SEGMENTS_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_SEGMENTS_2_H

#include <CGAL/license/Principal_component_analysis.h>


#include <CGAL/basic.h>
#include <CGAL/centroid.h>
#include <CGAL/pca_fitting_2.h>
#include <CGAL/compute_moment_2.h>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/Subiterator.h>

#include <iterator>
#include <list>
#include <cmath>

namespace CGAL {

namespace internal {

// Fits a line to a set of 2D segments
template < typename InputIterator, typename K, typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_2& line,   // best fit line
                               typename K::Point_2& c,     // centroid
                               const typename K::Segment_2*,// used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<1>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  // types
  typedef typename K::Segment_2 Segment;

  typename DiagonalizeTraits::Covariance_matrix covariance = { { 0., 0., 0. } };
  compute_centroid_and_covariance_2(first, beyond, c, covariance, (Segment*)nullptr, k, tag);

  return fitting_line_2(covariance, c, line, k, diagonalize_traits);
} // end linear_least_squares_fitting_2 for segment set with 1D tag


// Fits a line to the vertices of a set of 2D segments
template < typename InputIterator, typename K, typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_2& line,   // best fit line
                               typename K::Point_2& c,     // centroid
                               const typename K::Segment_2*,// used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<0>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  // types
  typedef typename K::Point_2  Point;
  typedef typename K::Segment_2 Segment;
  auto converter = [](const Segment& s, std::size_t idx) -> Point { return s[idx]; };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_2
    (make_subiterator<Point, 2> (first, converter),
     make_subiterator<Point, 2> (beyond),
     line,c,(Point*)nullptr,k,tag,
     diagonalize_traits);
} // end linear_least_squares_fitting_2 for segment set with 0D tag

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_SEGMENTS_2_H
