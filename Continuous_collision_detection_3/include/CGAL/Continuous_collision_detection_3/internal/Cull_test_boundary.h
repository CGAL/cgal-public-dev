// Copyright (c) 2023 GeometryFactory (France).
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef CULL_TEST_BOUNDARY_H
#define CULL_TEST_BOUNDARY_H

#include <vector>

namespace CGAL {
namespace Collisions {
namespace internal {



template <class K>
bool all_on_one_side_of_plane(
  const std::vector<typename K::Point_3*>& points_,
  const typename K::Vector_3& plane_normal
) {

  auto desired_orientation = ::CGAL::sign(
    ::CGAL::scalar_product( plane_normal, *(points_.front()) - ::CGAL::ORIGIN )
  );

  // If some of the points are on the plane,
  // it's possible that a line or patch contains the origin
  if( desired_orientation == ::CGAL::ZERO ) { return false; }

  return std::all_of(
    std::next(points_.begin()), // Don't need to check the first one
    points_.end(),
    [&plane_normal, &desired_orientation](const auto& m_p) {
      // TODO: there is a construction here, we need to use a predicate instead ==> add a new predicate in the Kernel
      return (
            ::CGAL::sign(::CGAL::scalar_product( plane_normal, *m_p - ::CGAL::ORIGIN ))
        ==  desired_orientation // check same side of plane
      );
    }
  );
}

template <class K>
bool cull_test_boundary(const std::vector<typename K::Point_3*>& points_)
{
  // TODO: instantiate these normal vectors once at
  //       compile time
  for (int i = -1; i < 1; ++i) {
  for (int j = -1; j < 1; ++j) {
  for (int k = -1; k < 1; ++k) {
    if(
      all_on_one_side_of_plane<K>(
        points_,
        typename K::Vector_3(typename K::FT(i), typename K::FT(j), typename K::FT(k))
      )
    ) {
      return true;
    }
  }
  }
  }
  return false;
}



}
}
}

#endif