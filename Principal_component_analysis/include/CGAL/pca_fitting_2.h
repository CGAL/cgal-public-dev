// Copyright (c) 2020  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_PCA_FITTING_2_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_PCA_FITTING_2_H

#include <CGAL/license/Principal_component_analysis.h>

#include <CGAL/Linear_algebraCd.h>
#include <CGAL/Dimension.h>

namespace CGAL {

namespace internal {

// compute the eigen values and vectors of the covariance
// matrix and deduces the best linear fitting line
// (this is an internal function)
// Returns a fitting quality (1 - lambda_min/lambda_max):
//  1 is best  (zero variance orthogonally to the fitting line);
//  0 is worst (isotropic case, returns a line with horizontal
//              direction by default)

template < typename K, typename DiagonalizeTraits >
typename K::FT
fitting_line_2(typename DiagonalizeTraits::Covariance_matrix& covariance, // covariance matrix
               const typename K::Point_2& c,       // centroid
               typename K::Line_2& line,           // best fit line
               const K&,                           // kernel
               const DiagonalizeTraits&)          // Diagonalize traits
{
  typedef typename K::FT       FT;
  typedef typename K::Line_2   Line;
  typedef typename K::Vector_2 Vector;

  // solve for eigenvalues and eigenvectors.
  // eigen values are sorted in ascending order,
  // eigen vectors are sorted in accordance.
  typename DiagonalizeTraits::Vector eigen_values = { { 0. , 0. } };
  typename DiagonalizeTraits::Matrix eigen_vectors = { { 0., 0., 0. } }; // TOFIX ?????????
  DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix
    (covariance, eigen_values, eigen_vectors);

  // check unicity and build fitting line accordingly
  if (eigen_values[0] != eigen_values[1])
  {
    // regular case
    line = Line(c, Vector(eigen_vectors[2], eigen_vectors[3]));
    return (FT)1.0 - eigen_values[0] / eigen_values[1];
  }
  else
  {
    // isotropic case (infinite number of directions)
    // by default: assemble a line that goes through
    // the centroid and with a default horizontal vector.
    line = Line(c, Vector(1.0, 0.0));
    return (FT)0.0;
  }

}

} // end namespace internal

} //namespace CGAL


#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_PCA_FITTING_2_H
