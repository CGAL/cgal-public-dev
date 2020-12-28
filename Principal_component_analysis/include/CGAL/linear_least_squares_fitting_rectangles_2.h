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

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_RECTANGLES_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_RECTANGLES_2_H

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

// Fits a line to a 2D rectangle set.
// Returns a fitting quality (1 - lambda_min/lambda_max):
//  1 is best  (zero variance orthogonally to the fitting line);
//  0 is worst (isotropic case, returns a line with horizontal
//              direction by default)

template < typename InputIterator, typename K, typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_2& line,   // best fit line
                               typename K::Point_2& c,     // centroid
                               const typename K::Iso_rectangle_2*,// used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<2>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  // types
  typedef typename K::Iso_rectangle_2 Iso_rectangle;

  typename DiagonalizeTraits::Covariance_matrix covariance = { { 0., 0., 0. } };
  compute_centroid_and_covariance_2(first, beyond, c, covariance, (Iso_rectangle*)nullptr, k, tag);

  return fitting_line_2(covariance, c, line, k, diagonalize_traits);
} // end linear_least_squares_fitting_2 for rectangle set with 2D tag

template < typename InputIterator, typename K, typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_2& line,   // best fit line
                               typename K::Point_2& c,     // centroid
                               const typename K::Iso_rectangle_2*,// used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<1>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  // types
  typedef typename K::Iso_rectangle_2 Iso_rectangle;
  typedef typename K::Segment_2         Segment_2;
  auto converter = [](const Iso_rectangle& r, std::size_t idx) -> Segment_2 { return Segment_2(r[idx], r[(idx+1)%4]); };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_2
    (make_subiterator<Segment_2, 4> (first, converter),
     make_subiterator<Segment_2, 4> (beyond),
     line,c,(Segment_2*)nullptr,K(),tag,
     diagonalize_traits);

} // end linear_least_squares_fitting_2 for rectangle set with 1D tag


template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_2& line,   // best fit line
                               typename K::Point_2& c,     // centroid
                               const typename K::Iso_rectangle_2*,// used for indirection
                               const K& k,                   // kernel
                               const CGAL::Dimension_tag<0>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  // types
  typedef typename K::Iso_rectangle_2 Iso_rectangle;
  typedef typename K::Point_2         Point_2;
  auto converter = [](const Iso_rectangle& r, std::size_t idx) -> Point_2 { return r[idx]; };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_2
    (make_subiterator<Point_2, 4> (first, converter),
     make_subiterator<Point_2, 4> (beyond),
     line,c,(Point_2*)nullptr,K(),tag,
     diagonalize_traits);

} // end linear_least_squares_fitting_2 for rectangle set with 0D tag

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_RECTANGLES_2_H
