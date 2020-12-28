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

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_PCA_FITTING_3_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_PCA_FITTING_3_H

#include <CGAL/license/Principal_component_analysis.h>

#include <CGAL/Linear_algebraCd.h>
#include <CGAL/Dimension.h>

namespace CGAL {

namespace internal {

// compute the eigen values and vectors of the covariance
// matrix and deduces the best linear fitting plane.
// returns fitting quality
template < typename K, typename DiagonalizeTraits >
typename K::FT
fitting_plane_3(typename DiagonalizeTraits::Covariance_matrix& covariance, // covariance matrix
                const typename K::Point_3& c,       // centroid
                typename K::Plane_3& plane,         // best fit plane
                const K&,                           // kernel
                const DiagonalizeTraits& )          // Diagonalize traits
{
  typedef typename K::FT       FT;
  typedef typename K::Plane_3  Plane;
  typedef typename K::Vector_3 Vector;

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in ascending order,
  // eigen vectors are sorted in accordance.
  typename DiagonalizeTraits::Vector eigen_values = {{ 0. , 0., 0. }};
  typename DiagonalizeTraits::Matrix eigen_vectors = {{ 0., 0., 0.,
                                                        0., 0., 0.,
                                                        0., 0., 0. }};
  DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix
    (covariance, eigen_values, eigen_vectors);

  // degenerate case
  if(eigen_values[0] == eigen_values[1] &&
     eigen_values[1] == eigen_values[2])
  {
    // assemble a default horizontal plane that goes
    // through the centroid.
    plane = Plane(c, Vector(FT(0.0), FT(0.0), FT(1.0)));
    return FT(0.0);
  }
  else // regular and line case
  {
    Vector normal(eigen_vectors[0], eigen_vectors[1], eigen_vectors[2]);
    plane = Plane(c, normal);
    return FT(1.0) - eigen_values[0] / eigen_values[1];
  } // end regular case
}

// compute the eigen values and vectors of the covariance
// matrix and deduces the best linear fitting line
// (this is an internal function)
// returns fitting quality
template < typename K, typename DiagonalizeTraits >
typename K::FT
fitting_line_3(typename DiagonalizeTraits::Covariance_matrix& covariance, // covariance matrix
               const typename K::Point_3& c,       // centroid
               typename K::Line_3& line,           // best fit line
               const K&,                           // kernel
               const DiagonalizeTraits& )                 // Diagonalize traits
{
  typedef typename K::FT       FT;
  typedef typename K::Line_3   Line;
  typedef typename K::Vector_3 Vector;

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in ascending order,
  // eigen vectors are sorted in accordance.
  typename DiagonalizeTraits::Vector eigen_values = {{ 0. , 0., 0. }};
  typename DiagonalizeTraits::Matrix eigen_vectors = {{ 0., 0., 0.,
                                                        0., 0., 0.,
                                                        0., 0., 0. }};
  DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix
    (covariance, eigen_values, eigen_vectors);

  // isotropic case (infinite number of directions)
  if(eigen_values[0] == eigen_values[1] &&
     eigen_values[0] == eigen_values[2])
  {
    // assemble a default line along x axis which goes
    // through the centroid.
    line = Line(c, Vector(FT(1.0), FT(0.0), FT(0.0)));
    return (FT)0.0;
  }
  else
  {
    // regular case
    Vector direction(eigen_vectors[6], eigen_vectors[7], eigen_vectors[8]);
    line = Line(c, direction);
    return (FT)1.0 - eigen_values[1] / eigen_values[2];
  }
}

} // end namespace internal

} //namespace CGAL


#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_PCA_FITTING_3_H
