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

#ifndef BILINEAR_PATCH_3_H
#define BILINEAR_PATCH_3_H

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
#include <variant>

namespace CGAL {


// Predicates
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



/*!
\ingroup PkgCollisions3Classes

The class `CGAL::Bilinear_patch_3` provides a bilinear-patch implementation for cartesian
and related kernels. The patch is defined by the bilinear interpolation between four points
in space, \f$ { x_0, x_1, x_2, x_3 \in \mathbb{K}^3 }\f$. Any point on the surface can be 
expressed as \f$ { p(u,v) = (1-u)(1-v)x_0 + u(1-v)x_1 + (1-u)vx_2 + uvx_3  }\f$. 

*/
template <class Kernel>
class Bilinear_patch_3
{
  typedef Bilinear_patch::internal::has_on_pred<Kernel, Kernel::Has_filtered_predicates> Has_on_predicate;
  typedef Bilinear_patch::internal::orientation_pred<Kernel, Kernel::Has_filtered_predicates> Orientation_predicate;

  typedef typename Kernel::FT           FT;
  typedef typename Kernel::Point_3      Point_3;
  typedef typename Kernel::Vector_3     Vector_3;
  typedef typename Kernel::Triangle_3   Triangle_3;
  typedef typename Kernel::Segment_3    Segment_3;

  typedef          std::array<Point_3, 4>             Rep;
  typedef typename Kernel::template Handle<Rep>::type Base;

  Base base;

  Origin origin = ::CGAL::ORIGIN;

  Has_on_predicate      has_on_pred      = Has_on_predicate();
  Orientation_predicate orientation_pred = Orientation_predicate();

  std::vector<Triangle_3> triangles_;
  std::vector<Segment_3>  boundary_;

  typename Kernel::Tetrahedron_3 bounding_tetrahedron;

  bool IS_PLANAR_{false};
  bool HAS_TRIANGULAR_DECOMPOSITION_{false};
  bool IS_DEGENERATE_{false};
  bool COLLINEAR_012_{false};
  bool COLLINEAR_230_{false};
  bool COLLINEAR_013_{false};
  bool COLLINEAR_123_{false};

public:
  /// \name Types
  /// @{

  /// @brief The kernel type
  using K = Kernel;

  /// @brief The type of tetrahedron returned by `Bilinear_patch_3::tetrahedron()`
  typedef typename K::Tetrahedron_3 Tetrahedron_3;

  /// @}

  /// \name Creation
  /// @{

  /// @brief Creates an empty bilinear patch
  Bilinear_patch_3() {}

  /// @brief Creates a bilinear patch with corners p, q, r, and s.
  /// @details The constructor assumes that the vertices are provided in an order
  /// such that pqrsp is a valid cycle around the edges of the bilinear patch.
  Bilinear_patch_3(const Point_3 &p, const Point_3 &q, const Point_3 &r, const Point_3 &s)
    : base(CGAL::make_array(p, q, r, s))
    {
      bounding_tetrahedron = Tetrahedron_3(p, q, r, s);

      IS_PLANAR_ = ::CGAL::coplanar(p,q,r,s);

      if(IS_PLANAR_) {
        COLLINEAR_012_ = ::CGAL::collinear(vertex(0), vertex(1), vertex(2));
        COLLINEAR_230_ = ::CGAL::collinear(vertex(2), vertex(3), vertex(0));
        COLLINEAR_013_ = ::CGAL::collinear(vertex(0), vertex(1), vertex(3));
        COLLINEAR_123_ = ::CGAL::collinear(vertex(1), vertex(2), vertex(3));
        HAS_TRIANGULAR_DECOMPOSITION_ = (!COLLINEAR_012_) || (!COLLINEAR_230_) || (!COLLINEAR_013_);
        IS_DEGENERATE_ = IS_PLANAR_ && (!HAS_TRIANGULAR_DECOMPOSITION_);
      }

      if(HAS_TRIANGULAR_DECOMPOSITION_)
      {
        populate_triangles_();
      }

      populate_boundary_();
    }

  /// @}

  /// \name Operators
  /// @{

  /// @brief Returns true if the bilinear patches are equivalent.
  /// @details The bilinear patches are deemed equivalent if each patch
  /// is defined by the same set of vertices and the order of the vertices 
  /// in each patch produces equivalent cycles.
  bool  operator==(const Bilinear_patch_3& bp) const;

  /// @brief Returns true if `operator==()` returns false.
  bool  operator!=(const Bilinear_patch_3& bp) const;

  friend std::ostream& operator<<(std::ostream& os, Bilinear_patch_3 const& bp)
  {
      return (os << bp.vertex(0) << "\n" << bp.vertex(1) << "\n"<< bp.vertex(2) << "\n"<< bp.vertex(3) << "\n" );
  }

  /// @}

  /// \name Predicates
  /// @{

  /// @brief Returns the orientation of the point with respect to the bilinear patch.
  /// @details The orientation is determined by evaluation of the function
  /// `signed_scaled_patch_distance()`.
  ::CGAL::Orientation orientation(const Point_3& p) const;

  /// @brief Returns true if the point lies on the bilinear patch
  /// @details If the bilinear patch is non-planar, then the point is determined
  /// to lie on the patch by evaluation of the function `signed_scaled_patch_distance()`.
  /// If the bilinear patch is planar but nondegenerate, the bilinear patch is decomposed into  
  /// triangles, and the `Triangle_3::has_on()` method is used. If the bilinear patch is
  /// degenerate, then it is decomposed into segments, and the `Segment_3::has_on()` method is
  /// used.
  bool  has_on(const Point_3& p) const;

  /// @brief Returns true if the bilinear patch is equivalent to a collection of segments. 
  bool  is_degenerate() const;

  /// @brief Returns true if the corners of the bilinear patch are coplanar. 
  bool  is_planar() const;

  /// @}

  /// \name Methods
  /// @{

  /// @brief Returns a point on the bilinear patch given its corresponding parametric values.
  Point_3 get_point_from_parametric_coordinates(const FT& u, const FT& v) const;

  /// @brief Returns the point corresponding to the ith vertex. 
  const Point_3 & vertex(int i) const;

  /// @brief Returns the point corresponding to the ith vertex.
  const Point_3 & operator[](int i) const;

  /// @brief Returns the tetrahedron characterized by the four vertices of the bilinear patch.
  const Tetrahedron_3 & tetrahedron() const;

  /// @brief Returns a signed, scaled patch distance between the bilinear patch and the provided point.
  /// @details The signed scaled patch distance is computed as \f$ { d_{013}(x)d_{123}(x) - d_{012}(x)d_{023}(x) } \f$,
  /// where \f$ { d_{pqr}(x) } \f$ is six-times the signed distance between the x and the plane defined by p, q, and r.
  FT signed_scaled_patch_distance( const Point_3& x) const;

  /// @brief Returns a reference to the triangles that compromise a planar bilinear patch.
  /// \pre The bilinear patch must be planar and non-degenerate for this method to be well-defined.
  const std::vector<Triangle_3>& triangles() const;

  /// @}

private:
  void populate_triangles_();
  void populate_boundary_();

};


template <class K>
auto Bilinear_patch_3<K>::get_point_from_parametric_coordinates(const FT& u, const FT& v) const -> Point_3
{
  FT ONE{1.};
  // Interpolates between the points of the bilnear patch
  Vector_3 interpolant = (
      (ONE-u) * (ONE-v) * (vertex(0) - origin)
    + (ONE-u) * (v)     * (vertex(1) - origin)
    + (u)     * (v)     * (vertex(2) - origin)
    + (u)     * (ONE-v) * (vertex(3) - origin)
  );
  return origin + interpolant;
}

template <class K>
auto Bilinear_patch_3<K>::signed_scaled_patch_distance(
  const Point_3 & x
) const -> FT
{
  return Bilinear_patch::internal::signed_scaled_patch_distance<K>(
    x,
    vertex(0),
    vertex(1),
    vertex(2),
    vertex(3));
}

template < class K >
bool
Bilinear_patch_3<K>::operator==(const Bilinear_patch_3<K> &bp) const
{
  if (::CGAL::identical(base, bp.base))
      return true;

  int i;
  for(i=0; i<4; i++)
    if ( vertex(0) == bp.vertex(i) )
      break;

  return (
        (i<4)
    &&  vertex(1) == bp.vertex(i+1)
    &&  vertex(2) == bp.vertex(i+2)
    &&  vertex(3) == bp.vertex(i+3)
  );
}

template < class K >
inline
bool
Bilinear_patch_3<K>::operator!=(const Bilinear_patch_3<K> &bp) const
{
  return !(*this == bp);
}

template < class K >
const typename Bilinear_patch_3<K>::Point_3 &
Bilinear_patch_3<K>::vertex(int i) const
{
  if (i<0) i=(i%4)+4;
  else if (i>3) i=i%4;
  return (
    (i==0) ?  get_pointee_or_identity(base)[0] :
    (i==1) ?  get_pointee_or_identity(base)[1] :
    (i==2) ?  get_pointee_or_identity(base)[2] :
              get_pointee_or_identity(base)[3]
  );
}

template < class K >
inline
const typename Bilinear_patch_3<K>::Point_3 &
Bilinear_patch_3<K>::operator[](int i) const
{
  return vertex(i);
}

template <class K>
::CGAL::Orientation
Bilinear_patch_3<K>::orientation(const Point_3 &p) const
{
  CGAL_precondition(!is_planar());
  return orientation_pred(
    p,
    vertex(0),
    vertex(1),
    vertex(2),
    vertex(3)
  );
}

template <class K>
bool
Bilinear_patch_3<K>::has_on(const Point_3 &p) const
{
  bool has_on_{false};
  if( HAS_TRIANGULAR_DECOMPOSITION_ )
  {
    // If the patch can be decomposed into triangles,
    // check those instead
    for( const auto& t : triangles_)
    {
      has_on_ = (has_on_ || t.has_on(p));
    }
  }
  else if( IS_DEGENERATE_ )
  {
    // If the patch is truly degenerate, then it
    // consists only of its boundary edges
    for( const auto& s : boundary_ )
    {
      has_on_ = (has_on_ || s.has_on(p));
    }
  }
  else if(
        bounding_tetrahedron.has_on_bounded_side(p)
    ||  bounding_tetrahedron.has_on_boundary(p)
  )
  {
    // If the patch is not degenerate, then we can evaluate
    // the signed distance function, which should return zero
    // for any points on the patch
    has_on_ = has_on_pred(
      p,
      vertex(0),
      vertex(1),
      vertex(2),
      vertex(3)
    );
  }
  return has_on_;
}

template < class K >
bool
Bilinear_patch_3<K>::is_degenerate() const
{
    return IS_DEGENERATE_;
}

template < class K >
bool
Bilinear_patch_3<K>::is_planar() const
{
  return IS_PLANAR_;
}

template < class K >
const typename Bilinear_patch_3<K>::Tetrahedron_3 &
Bilinear_patch_3<K>::tetrahedron() const
{
  return bounding_tetrahedron;
}

template <class K>
void Bilinear_patch_3<K>::populate_triangles_()
{

  if( !COLLINEAR_012_ && !COLLINEAR_230_ && !COLLINEAR_013_ && !COLLINEAR_123_ )
  {
    // Case 1
    // No collinear triplets and orientations agree
    Vector_3 normal_012 = normal(vertex(0), vertex(1), vertex(2));
    Vector_3 normal_230 = normal(vertex(2), vertex(3), vertex(0));
    if( scalar_product(normal_012, normal_230) > typename K::FT(0) )
    {
      triangles_.reserve(2);
      triangles_.push_back(Triangle_3(vertex(0), vertex(1), vertex(2)));
      triangles_.push_back(Triangle_3(vertex(2), vertex(3), vertex(0)));
      return;
    }

    // Case 2a
    // No collinear triplets, but e01 and e23 intersect
    const auto intersection_01_23 = intersection(Segment_3(vertex(0), vertex(1)), Segment_3(vertex(2), vertex(3)));
    if( const Point_3* intersection_point = std::get_if<Point_3>(&*intersection_01_23) ) {
      // This cannot be false or return a segment, since nothing is collinear
      // but an intersection has occurred.
      triangles_.reserve(2);
      triangles_.push_back(Triangle_3(vertex(0), *intersection_point, vertex(3)));
      triangles_.push_back(Triangle_3(vertex(1), vertex(2), *intersection_point));
      return;
    }

    // Case 2b
    // No collinear triplets, but e12 and e30 intersect
    const auto intersection_12_30 = intersection(Segment_3(vertex(1), vertex(2)), Segment_3(vertex(3), vertex(0)));
    if( const Point_3* intersection_point = std::get_if<Point_3>(&*intersection_12_30) ) {
      // This cannot be false or return a segment, since nothing is collinear
      // but an intersection has occurred.
      triangles_.reserve(2);
      triangles_.push_back(Triangle_3(vertex(1), *intersection_point, vertex(0)));
      triangles_.push_back(Triangle_3(vertex(2), vertex(3), *intersection_point));
      return;
    }
  }

  // Case 3
  // Exactly three of the four vertices are collinear. A single triangle covers the patch.
  // Use the largest triangle.
  triangles_.reserve(1);
  if( COLLINEAR_012_ )
  {
    size_t i{0}, j{1};
    if ( are_ordered_along_line(vertex(0), vertex(1), vertex(2)) ) { j = 2; }
    else if( are_ordered_along_line(vertex(2), vertex(0), vertex(1)) ) { i = 2; }
    triangles_.push_back(Triangle_3(vertex(i), vertex(j), vertex(3)));
    return;
  }
  else if( COLLINEAR_230_ )
  {
    size_t i{0}, j{3};
    if ( are_ordered_along_line(vertex(0), vertex(3), vertex(2)) ) { j = 2; }
    else if( are_ordered_along_line(vertex(2), vertex(0), vertex(3)) ) { i = 2; }
    triangles_.push_back(Triangle_3(vertex(i), vertex(j), vertex(1)));
    return;
  }
  else if( COLLINEAR_013_ )
  {
    size_t i{0}, j{3};
    if ( are_ordered_along_line(vertex(1), vertex(0), vertex(3)) ) { i = 1; }
    else if( are_ordered_along_line(vertex(0), vertex(3), vertex(1)) ) { j = 1; }
    triangles_.push_back(Triangle_3(vertex(i), vertex(j), vertex(2)));
    return;
  }
  else
  {
    size_t i{1}, j{3};
    if ( are_ordered_along_line(vertex(1), vertex(3), vertex(2)) ) { j = 2; }
    else if( are_ordered_along_line(vertex(2), vertex(1), vertex(3)) ) { i = 2; }
    triangles_.push_back(Triangle_3(vertex(i), vertex(j), vertex(0)));
    return;
  }
  // Unreachable
}

template <class K>
void Bilinear_patch_3<K>::populate_boundary_()
{
  boundary_.reserve(4);
  boundary_.push_back(Segment_3(vertex(0), vertex(1)));
  boundary_.push_back(Segment_3(vertex(1), vertex(2)));
  boundary_.push_back(Segment_3(vertex(2), vertex(3)));
  boundary_.push_back(Segment_3(vertex(3), vertex(0)));
}

template <class K>
auto Bilinear_patch_3<K>::triangles() const -> const std::vector<Triangle_3>& 
{
  CGAL_precondition(IS_PLANAR_ && !IS_DEGENERATE_); 
  return triangles_;
}



} //namespace CGAL

#endif // CGAL_CARTESIAN_BILINEAR_PATCH_3_H
