// Copyright (c) 2023
// INRIA Sophia-Antipolis (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef BILINEAR_PATCH_3_PREDICATES_H
#define BILINEAR_PATCH_3_PREDICATES_H

#include <CGAL/Named_function_parameters.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_kernel_selector.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/array.h>
#include <CGAL/intersections.h>

#include <iostream>


namespace CGAL {
namespace Bilinear_patch {
namespace internal {

  template <typename K>
  auto signed_scaled_planar_distance(
    const typename K::Point_3& x,
    const typename K::Point_3& p,
    const typename K::Point_3& q,
    const typename K::Point_3& r
  ) -> typename K::FT
  {
    return ::CGAL::scalar_product(x-p, ::CGAL::cross_product(q-p, r-p));
  }

  template <typename K>
  auto signed_scaled_patch_distance(
    const typename K::Point_3& x,
    const typename K::Point_3& v0,
    const typename K::Point_3& v1,
    const typename K::Point_3& v2,
    const typename K::Point_3& v3
  ) -> typename K::FT
  {
    return (
        (
            signed_scaled_planar_distance<K>(x, v0, v1, v3)
          * signed_scaled_planar_distance<K>(x, v1, v2, v3)
        )
      - (
            signed_scaled_planar_distance<K>(x, v0, v1, v2)
          * signed_scaled_planar_distance<K>(x, v0, v2, v3)
        )
    );
  }

  template <typename K, template <class Kernel> class Pred>
  struct Get_filtered_predicate_FT
  {
    typedef typename ::CGAL::Exact_kernel_selector<K>::Exact_kernel     Exact_kernel_ft;
    typedef typename ::CGAL::Exact_kernel_selector<K>::C2E              C2E_ft;
    typedef ::CGAL::Simple_cartesian<Interval_nt_advanced>              Approximate_kernel;
    typedef ::CGAL::Cartesian_converter<K, Approximate_kernel>          C2F;

    typedef ::CGAL::Filtered_predicate<
      Pred<Exact_kernel_ft>,
      Pred<Approximate_kernel>,
      C2E_ft, C2F
    > type;
  };

  template <typename K>
  struct has_on_pred_impl
  {
    typedef bool result_type;

    bool operator()(
      const typename K::Point_3& x,
      const typename K::Point_3& v0,
      const typename K::Point_3& v1,
      const typename K::Point_3& v2,
      const typename K::Point_3& v3
    ) const
    {
      using FT = typename K::FT; 
      return (
            ::CGAL::compare<FT, FT>(
              ::CGAL::Bilinear_patch::internal::signed_scaled_patch_distance<K>(
                x, v0, v1, v2, v3
              ), 
              FT{0}
            )
        ==  ::CGAL::EQUAL
      );
    }
  };

  template <typename K, bool has_filtered_predicate = K::Has_filtered_predicates>
  struct has_on_pred : public has_on_pred_impl<K> 
  { 
    using has_on_pred_impl<K>::operator(); 
  };

  template <typename K>
  struct has_on_pred<K, true> : public Get_filtered_predicate_FT<K, has_on_pred_impl>::type
  {
    using Get_filtered_predicate_FT<K, has_on_pred_impl>::type::operator();
  };

  // ORIENTATION PREDICATE
  template <typename K>
  struct orientation_pred_impl
  {
    typedef ::CGAL::Orientation result_type;

    result_type operator()(
      const typename K::Point_3& x,
      const typename K::Point_3& v0,
      const typename K::Point_3& v1,
      const typename K::Point_3& v2,
      const typename K::Point_3& v3
    ) const
    {
      using FT = typename K::FT; 
      FT dist = ::CGAL::Bilinear_patch::internal::signed_scaled_patch_distance<K>(
        x, v0, v1, v2, v3
      );

      return ::CGAL::enum_cast<result_type> (
        ::CGAL::compare<FT, FT>(dist, FT{0})
      );
    }
  };

  template <typename K, bool has_filtered_predicate = K::Has_filtered_predicates>
  struct orientation_pred : public orientation_pred_impl<K> 
  { 
    using orientation_pred_impl<K>::operator(); 
  };

  template <typename K>
  struct orientation_pred<K, true> : public Get_filtered_predicate_FT<K, orientation_pred_impl>::type
  {
    using Get_filtered_predicate_FT<K, orientation_pred_impl>::type::operator();
  };

}
}
}

#endif