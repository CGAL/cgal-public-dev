// Copyright (c) 2020  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez and Simon Giraudot

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_COMPUTE_MOMENT_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_COMPUTE_MOMENT_2_H

#include <CGAL/license/Principal_component_analysis.h>

#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/Default_diagonalize_traits.h>
#include <CGAL/Dimension.h>

namespace CGAL {

namespace internal {

// Initialize a matrix in n dimension by an array or numbers
template <typename FT>
typename CGAL::Linear_algebraCd<FT>::Matrix
init_matrix(const int n,
            FT entries[])
{
  CGAL_assertion(n > 1); // dimension > 1
  typedef typename CGAL::Linear_algebraCd<FT>::Matrix Matrix;

  Matrix m(n);
  int i,j;
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++)
      m[i][j] = entries[i*n+j];

  return m;
} // end initialization of matrix

// Computes closed form order-2 moment matrix of a 2D point set
// with respect to a given reference point.
template < typename InputIterator, typename Moment, typename K>
void
compute_moment_2(InputIterator first,
                 InputIterator beyond,
                 Moment& moment, // moment matrix
                 const typename K::Point_2& reference, // given reference point
                 const K&, // kernel
                 const typename K::Point_2*,   // used for indirection
                 const CGAL::Dimension_tag<0>& tag)
{
  // types
  typedef typename K::FT       FT;
  typedef typename K::Point_2  Point;
  typedef typename K::Vector_2 Vector;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute moment matrix as a 2D symmetric matrix.
  for (InputIterator it = first;
       it != beyond;
       it++)
  {
    const Point& p = *it;
    const Vector d = p - reference; // centered data point wrt reference point
    moment[0] += d.x() * d.x();
    moment[1] += d.x() * d.y();
    moment[2] += d.y() * d.y();
  }
} // end compute_moment_2 for a 2D point set


// Computes closed form order-2 moment matrix of 2D segments
// with respect to a given reference point.
template < typename InputIterator, typename Moment, typename K>
void
compute_moment_2(InputIterator first,
                 InputIterator beyond,
                 Moment& moment, // moment matrix
                 const typename K::Point_2& reference, // given reference point
                 const K&,                             // kernel
                 const typename K::Segment_2*,   // used for indirection
                 const CGAL::Dimension_tag<1>& tag)
{
  // types
  typedef typename K::FT       FT;
  typedef typename K::Line_2   Line;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Segment_2 Segment;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  // assemble 2nd order moment about the origin.
  FT temp[4] = { 1.0, 0.5, 0.5, 1.0 };
  Matrix canonical_moment = FT(1.0 / 3.0) * init_matrix<FT>(2, temp);

  FT mass = 0.0; // accumulate segment lengths
  for (InputIterator it = first;
       it != beyond;
       it++)
  {
    // Now for each segment, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Segment& t = *it;

    // defined for convenience.
    // FT example = CGAL::to_double(t[0].x());
    FT delta[4] = { t[0].x(), t[1].x(),
                    t[0].y(), t[1].y() };
    Matrix transformation = init_matrix<FT>(2, delta);
    const FT length = CGAL::approximate_sqrt(t.squared_length());

    // Find the 2nd order moment for the segment wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = length * transformation * canonical_moment * LA::transpose(transformation);

    // add to moment matrix
    moment[0] += transformation[0][0];
    moment[1] += transformation[0][1];
    moment[2] += transformation[1][1];

    mass += length;
  }

  CGAL_assertion_msg (mass != FT(0), "Can't compute PCA of null measure.");

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the moment.
  moment[0] += mass * (-1.0 * reference.x() * reference.x());
  moment[1] += mass * (-1.0 * reference.x() * reference.y());
  moment[2] += mass * (-1.0 * reference.y() * reference.y());
}

// Computes closed form order-2 moment matrix of 2D iso rectangles
// with respect to a given reference point.
template < typename InputIterator, typename Moment, typename K>
void
compute_moment_2(InputIterator first,
                 InputIterator beyond,
                 Moment& moment, // moment matrix
                 const typename K::Point_2& reference, // given reference point
                 const K&,                             // kernel
                 const typename K::Iso_rectangle_2*,   // used for indirection
                 const CGAL::Dimension_tag<2>& tag)
{
  // types
  typedef typename K::FT       FT;
  typedef typename K::Iso_rectangle_2 Iso_rectangle;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // assemble order-2 matrix as a semi-definite matrix.
  FT mass = 0.0;
  // moment = { { 0., 0., 0. } };

  // assemble reference 2nd order moment about the origin.
  FT temp[4] = { 1.0 / 3.0, 0.25,
                 0.25,  1.0 / 3.0 };
  Matrix unit_moment = init_matrix<FT>(2, temp);

  for (InputIterator it = first;
       it != beyond;
       it++)
  {
    // Now for each rectangle, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Iso_rectangle& t = *it;

    // defined for convenience.
    const FT x0 = t.xmin();
    const FT y0 = t.ymin();
    const FT x1 = t.xmax();
    const FT y2 = t.ymax();

    FT delta[4] = { x1 - x0, 0.0,
                    0.0, y2 - y0 };

    // TODO: use Eigen
    Matrix transformation = init_matrix<FT>(2, delta);
    const FT area = (x1 - x0) * (y2 - y0);

    CGAL_assertion(area != 0.0);

    // Find the 2nd order moment for the rectangle wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = area * transformation * unit_moment * LA::transpose(transformation);

    // Translate the 2nd order moment to the center of the rectangle.
    const FT xav0 = (x1 - x0) / 2.0;
    const FT yav0 = (y2 - y0) / 2.0;

    // and add to moment matrix
    moment[0] += transformation[0][0] + area * (x0 * xav0 * 2 + x0 * x0);
    moment[1] += transformation[0][1] + area * (x0 * yav0 + xav0 * y0 + x0 * y0);
    moment[2] += transformation[1][1] + area * (y0 * yav0 * 2 + y0 * y0);

    mass += area;
  }

  CGAL_assertion_msg (mass != FT(0), "Can't compute PCA of null measure.");

  // Translate the 2nd order moment calculated about the given reference point
  moment[0] += mass * (-1.0 * reference.x() * reference.x());
  moment[1] += mass * (-1.0 * reference.x() * reference.y());
  moment[2] += mass * (-1.0 * reference.y() * reference.y());

} // end compute_moment_2 for a set of 2D iso rectangles


// Computes closed form order-2 moment matrix of 2D discs
// with respect to a given reference point.
template < typename InputIterator, typename Moment, typename K>
void
compute_moment_2(InputIterator first,
                 InputIterator beyond,
                 Moment& moment, // moment matrix
                 const typename K::Point_2& reference, // given reference point
                 const K&,                             // kernel
                 const typename K::Circle_2*,   // used for indirection
                 const CGAL::Dimension_tag<2>&)
{
  typedef typename K::FT FT;
  typedef typename K::Circle_2 Circle;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  FT mass = 0.0;

  // assemble 2nd order moment about the origin.
  FT temp[4] = { 0.25, 0.0,
                 0.0,  0.25 };
  Matrix canonical_moment = init_matrix<FT>(2, temp);
  //  Matrix moment = Matrix(2,true,PI);

  for (InputIterator it = first;
       it != beyond;
       it++)
  {
    // Now for each circle, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Circle& t = *it;

    // defined for convenience.
    // FT example = CGAL::to_double(t[0].x());
    FT radius = CGAL::approximate_sqrt(t.squared_radius());
    FT delta[4] = { radius, 0.0,
                    0.0, radius };
    Matrix transformation = init_matrix<FT>(2, delta);
    const FT area = t.squared_radius();
    // CGAL_assertion(area != 0.0); // not needed

    // Find the 2nd order moment for the circle wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = area * transformation * canonical_moment * LA::transpose(transformation);

    // Translate the 2nd order moment to the center of the circle.
    FT x0 = t.center().x();
    FT y0 = t.center().y();

    // and add to covariance matrix
    moment[0] += transformation[0][0] + area * x0 * x0;
    moment[1] += transformation[0][1] + area * x0 * y0;
    moment[2] += transformation[1][1] + area * y0 * y0;

    mass += area;
  }

  CGAL_assertion_msg (mass != FT(0), "Can't compute PCA of null measure.");

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  moment[0] += mass * (-1.0 * reference.x() * reference.x());
  moment[1] += mass * (-1.0 * reference.x() * reference.y());
  moment[2] += mass * (-1.0 * reference.y() * reference.y());
} // end compute_moment_2 for 2D discs


// Computes closed form order-2 moment matrix of 2D circles
// with respect to a given reference point.
template < typename InputIterator, typename Moment, typename K>
void
compute_moment_2(InputIterator first,
                 InputIterator beyond,
                 Moment& moment, // moment matrix
                 const typename K::Point_2& reference, // given reference point
                 const K&,                             // kernel
                 const typename K::Circle_2*,   // used for indirection
                 const CGAL::Dimension_tag<1>&)
{
  // types
  typedef typename K::FT       FT;
  typedef typename K::Line_2   Line;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Circle_2 Circle;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  FT mass = 0.0;

  // assemble 2nd order moment about the origin.
  FT temp[4] = { 1.0, 0.0,
                 0.0, 1.0 };
  Matrix canonical_moment = init_matrix<FT>(2, temp);

  for (InputIterator it = first;
       it != beyond;
       it++)
  {
    // Now for each circle, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Circle& t = *it;

    // defined for convenience.
    // FT example = CGAL::to_double(t[0].x());
    FT radius = CGAL::approximate_sqrt(t.squared_radius());
    FT delta[4] = { radius, 0.0,
                    0.0, radius };
    Matrix transformation = init_matrix<FT>(2, delta);
    const FT length = FT(2.0) * radius;

    // Find the 2nd order moment for the circle wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = FT(0.5) * length * transformation * canonical_moment * LA::transpose(transformation);

    // Translate the 2nd order moment to the center of the circle.
    FT x0 = t.center().x();
    FT y0 = t.center().y();

    // and add to covariance matrix
    moment[0] += transformation[0][0] + length * x0 * x0;
    moment[1] += transformation[0][1] + length * x0 * y0;
    moment[2] += transformation[1][1] + length * y0 * y0;

    mass += length;
  }

  CGAL_assertion_msg (mass != FT(0), "Can't compute PCA of null measure.");

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  moment[0] += mass * (-1.0 * reference.x() * reference.x());
  moment[1] += mass * (-1.0 * reference.x() * reference.y());
  moment[2] += mass * (-1.0 * reference.y() * reference.y());

}


// Computes closed form order-2 moment matrix of 2D circles
// with respect to a given reference point.
template < typename InputIterator, typename Moment, typename K>
void
compute_moment_2(InputIterator first,
                 InputIterator beyond,
                 Moment& moment, // moment matrix
                 const typename K::Point_2& reference, // given reference point
                 const K&,                             // kernel
                 const typename K::Triangle_2*,   // used for indirection
                 const CGAL::Dimension_tag<2>& tag)
{
  // types
  typedef typename K::FT       FT;
  typedef typename K::Line_2   Line;
  typedef typename K::Vector_2 Vector;
  typedef typename K::Triangle_2 Triangle;
  typedef typename CGAL::Linear_algebraCd<FT> LA;
  typedef typename LA::Matrix Matrix;

  // assemble the 2nd order moment about the origin.
  FT temp[4] = { 1.0 / 12.0, 1.0 / 24.0,
                 1.0 / 24.0, 1.0 / 12.0 };

  Matrix canonical_moment = init_matrix<FT>(2, temp);

  FT mass = 0.0; // accumulates triangle areas
  for (InputIterator it = first; it != beyond; it++)
  {
    // Now for each triangle, construct the 2nd order moment about the origin.
    // assemble the transformation matrix.
    const Triangle& t = *it;

    // defined for convenience.
    FT x0 = t[0].x();
    FT y0 = t[0].y();
    FT delta[4] = { t[1].x() - x0, t[2].x() - x0,
                    t[1].y() - y0, t[2].y() - y0 };

    Matrix transformation = init_matrix<FT>(2, delta);
    const FT area = FT(0.5) * CGAL::abs(LA::determinant(transformation));
    if (area == 0.0)
      continue;

    // Find the 2nd order moment for the triangle wrt to the origin by an affine transformation.

    // Transform the standard 2nd order moment using the transformation matrix
    transformation = 2 * area * transformation * canonical_moment * LA::transpose(transformation);

    // Translate the 2nd order moment to (x0, y0).
    const FT xav0 = (delta[0] + delta[1]) / 3.0;
    const FT yav0 = (delta[2] + delta[3]) / 3.0;

    // and add to the covariance matrix
    moment[0] += transformation[0][0] + area * (x0 * xav0 * 2 + x0 * x0);
    moment[1] += transformation[0][1] + area * (x0 * yav0 + xav0 * y0 + x0 * y0);
    moment[2] += transformation[1][1] + area * (y0 * yav0 * 2 + y0 * y0);

    mass += area;
  }

  CGAL_assertion_msg (mass != FT(0), "Can't compute PCA of null measure.");

  // Translate the 2nd order moment calculated about the origin to
  // the center of mass to get the covariance.
  moment[0] += mass * (-1.0 * reference.x() * reference.x());
  moment[1] += mass * (-1.0 * reference.x() * reference.y());
  moment[2] += mass * (-1.0 * reference.y() * reference.y());
}

// compute centroid and covariance matrix
template < typename InputIterator,
           typename K,
           typename Moment,
           typename Dimension,
           typename Primitive>
void
compute_centroid_and_covariance_2(InputIterator first,
                                  InputIterator beyond,
                                  typename K::Point_2& c,       // centroid
                                  Moment& covariance,
                                  const Primitive* primitive,  // used for indirection
                                  const K& k,                   // kernel
                                  const Dimension& tag)
{
  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first, beyond, k, tag);

  // assemble covariance matrix
  compute_moment_2(first, beyond, covariance, c, k, primitive, tag);

} // end compute_centroid_and_covariance_3

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_COMPUTE_MOMENT_2_H
