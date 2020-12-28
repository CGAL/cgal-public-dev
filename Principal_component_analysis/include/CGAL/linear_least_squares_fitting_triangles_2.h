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

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_TRIANGLES_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_TRIANGLES_2_H

#include <CGAL/license/Principal_component_analysis.h>


#include <CGAL/basic.h>
#include <CGAL/centroid.h>
#include <CGAL/pca_fitting_2.h>
#include <CGAL/compute_moment_2.h>
#include <CGAL/Linear_algebraCd.h>

#include <iterator>
#include <vector>
#include <cmath>

namespace CGAL {

<<<<<<< HEAD
namespace internal {
// Fits a line to a 2D triangle set.
// Returns a fitting quality (1 - lambda_min/lambda_max):
//  1 is best  (zero variance orthogonally to the fitting line);
//  0 is worst (isotropic case, returns a line with horizontal
//              direction by default)

template < typename InputIterator,
           typename Kernel, typename DiagonalizeTraits >
typename Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename Kernel::Line_2& line,   // best fit line
                               typename Kernel::Point_2& c,     // centroid
                               const typename Kernel::Triangle_2*,// used for indirection
                               const Kernel&,                   // kernel
                               const CGAL::Dimension_tag<2>& tag,
                               const DiagonalizeTraits&)
{
  // types
  typedef typename K::Triangle_2 Triangle;

  typename DiagonalizeTraits::Covariance_matrix covariance = { { 0., 0., 0. } };
  compute_centroid_and_covariance_2(first, beyond, c, covariance, (Triangle*)nullptr, k, tag);

  return fitting_line_2(covariance, c, line, k, diagonalize_traits);
} // end linear_least_squares_fitting_2 for triangle set with 2D tag

template < typename InputIterator,
           typename Kernel,
           typename DiagonalizeTraits >
typename Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename Kernel::Line_2& line,   // best fit line
                               typename Kernel::Point_2& c,     // centroid
                               const typename Kernel::Triangle_2*,// used for indirection
                               const Kernel&,                   // kernel
                               const CGAL::Dimension_tag<1>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  // types
  typedef typename Kernel::Triangle_2 Triangle;
  typedef typename Kernel::Segment_2  Segment;
  auto converter = [](const Triangle& t, std::size_t idx) -> Segment { return Segment(t[idx], t[(idx+1)%3]); };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_2
    (make_subiterator<Segment, 3> (first, converter),
     make_subiterator<Segment, 3> (beyond),
     line,c,(Segment*)nullptr,Kernel(),tag,
     diagonalize_traits);

} // end linear_least_squares_fitting_2 for triangle set with 1D tag

template < typename InputIterator,
           typename Kernel,
           typename DiagonalizeTraits >
typename Kernel::FT
linear_least_squares_fitting_2(InputIterator first,
                               InputIterator beyond,
                               typename Kernel::Line_2& line,   // best fit line
                               typename Kernel::Point_2& c,     // centroid
                               const typename Kernel::Triangle_2*,// used for indirection
                               const Kernel&,                   // kernel
                               const CGAL::Dimension_tag<0>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  // types

  typedef typename Kernel::Triangle_2 Triangle;
  typedef typename Kernel::Point_2 Point;
  auto converter = [](const Triangle& t, std::size_t idx) -> Point { return t[idx]; };

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  return linear_least_squares_fitting_2
    (make_subiterator<Point, 3> (first, converter),
     make_subiterator<Point, 3> (beyond),
     line,c,(Point*)nullptr,Kernel(),tag,
     diagonalize_traits);
} // end linear_least_squares_fitting_2 for triangle set with 0D tag

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_TRIANGLES_2_H
