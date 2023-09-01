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

#include <CGAL/Continuous_collision_detection_3/internal/Bilinear_patch_3_predicates.h>

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



template <class R_>
class BilinearPatchC3
{
  typedef Bilinear_patch::internal::has_on_pred<R_, R_::Has_filtered_predicates> Has_on_predicate;
  typedef Bilinear_patch::internal::orientation_pred<R_, R_::Has_filtered_predicates> Orientation_predicate;

  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Triangle_3           Triangle_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Tetrahedron_3        Tetrahedron_3;
  typedef          Interval_nt_advanced     IA_NT;

  typedef          std::array<Point_3, 4>         Rep;
  typedef typename R_::template Handle<Rep>::type Base;

  Base base;
  Origin origin = ::CGAL::ORIGIN;
  Has_on_predicate has_on_pred = Has_on_predicate();
  Orientation_predicate orientation_pred = Orientation_predicate();

public:
  std::vector<Triangle_3> triangles_;
  std::vector<Segment_3>  boundary_;
  typedef R_                                R;
  typedef typename R::Tetrahedron_3         Tetrahedron;

  Tetrahedron bounding_tetrahedron;

  bool IS_PLANAR_{false};
  bool HAS_TRIANGULAR_DECOMPOSITION_{false};
  bool IS_DEGENERATE_{false};
  bool COLLINEAR_012_{false};
  bool COLLINEAR_230_{false};
  bool COLLINEAR_013_{false};
  bool COLLINEAR_123_{false};

  BilinearPatchC3() {}
  // NOTE:  specify the points as a cycle around the perimeter!
  //                           q00,              q01,              q11,              q10
  BilinearPatchC3(const Point_3 &p, const Point_3 &q, const Point_3 &r, const Point_3 &s)
    : base(CGAL::make_array(p, q, r, s))
    {
      bounding_tetrahedron = Tetrahedron(p, q, r, s);

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

  Point_3 operator()(const FT& u, const FT& v) const; // Returns a point given its corresponding parametric values
  bool  operator==(const BilinearPatchC3& bp) const;
  bool  operator!=(const BilinearPatchC3& bp) const;
  friend std::ostream& operator<<(std::ostream& os, BilinearPatchC3 const& bp)
  {
      return (os << bp.vertex(0) << "\n" << bp.vertex(1) << "\n"<< bp.vertex(2) << "\n"<< bp.vertex(3) << "\n" );
  }


  ::CGAL::Orientation orientation(const Point_3& p) const;
  bool  has_on(const Point_3& p) const;
  bool  is_degenerate() const;
  bool  is_planar() const;

  const Point_3 & vertex(int i) const;
  const Point_3 & operator[](int i) const;
  const Tetrahedron_3 & tetrahedron() const;
  FT signed_scaled_patch_distance( const Point_3& x) const;


private:
  void populate_triangles_();
  void populate_boundary_();

};

template <class R>
auto BilinearPatchC3<R>::operator()(const FT& u, const FT& v) const -> Point_3
{
  FT ONE{1.};
  // Interpolates between the points of the bilnear patch
  Vector_3 interpolant = (
      (ONE-u) * (ONE-v) * (vertex(0) - origin)
    + (ONE-u) * (v)     * (vertex(1) - origin)
    + (u)     * (ONE-v) * (vertex(2) - origin)
    + (u)     * (v)     * (vertex(3) - origin)
  );
  return origin + interpolant;
}

template <class R>
auto BilinearPatchC3<R>::signed_scaled_patch_distance(
  const Point_3 & x
) const -> FT
{
  return Bilinear_patch::internal::signed_scaled_patch_distance<R>(
    x,
    vertex(0),
    vertex(1),
    vertex(2),
    vertex(3));
}

template < class R >
bool
BilinearPatchC3<R>::operator==(const BilinearPatchC3<R> &bp) const
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

template < class R >
inline
bool
BilinearPatchC3<R>::operator!=(const BilinearPatchC3<R> &bp) const
{
  return !(*this == bp);
}

template < class R >
const typename BilinearPatchC3<R>::Point_3 &
BilinearPatchC3<R>::vertex(int i) const
{
  if (i<0) i=(i%4)+4;
  else if (i>3) i=i%4;
  return (i==0) ? get_pointee_or_identity(base)[0] :
         (i==1) ? get_pointee_or_identity(base)[1] :
         (i==2) ? get_pointee_or_identity(base)[2] :
                  get_pointee_or_identity(base)[3];
}

template < class R >
inline
const typename BilinearPatchC3<R>::Point_3 &
BilinearPatchC3<R>::operator[](int i) const
{
  return vertex(i);
}

template <class R>
::CGAL::Orientation
BilinearPatchC3<R>::orientation(const Point_3 &p) const
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

template <class R>
bool
BilinearPatchC3<R>::has_on(const Point_3 &p) const
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

template < class R >
bool
BilinearPatchC3<R>::is_degenerate() const
{
    return IS_DEGENERATE_;
}

template < class R >
bool
BilinearPatchC3<R>::is_planar() const
{
  return IS_PLANAR_;
}

template < class R >
const typename BilinearPatchC3<R>::Tetrahedron_3 &
BilinearPatchC3<R>::tetrahedron() const
{
  return bounding_tetrahedron;
}

template <class R>
void BilinearPatchC3<R>::populate_triangles_()
{

  if( !COLLINEAR_012_ && !COLLINEAR_230_ && !COLLINEAR_013_ && !COLLINEAR_123_ )
  {
    // Case 1
    // No collinear triplets and orientations agree
    Vector_3 normal_012 = normal(vertex(0), vertex(1), vertex(2));
    Vector_3 normal_230 = normal(vertex(2), vertex(3), vertex(0));
    if( scalar_product(normal_012, normal_230) > typename R::FT(0) )
    {
      triangles_.reserve(2);
      triangles_.push_back(Triangle_3(vertex(0), vertex(1), vertex(2)));
      triangles_.push_back(Triangle_3(vertex(2), vertex(3), vertex(0)));
      return;
    }

    // Case 2a
    // No collinear triplets, but e01 and e23 intersect
    const auto intersection_01_23 = intersection(Segment_3(vertex(0), vertex(1)), Segment_3(vertex(2), vertex(3)));
    if( intersection_01_23 ) {
      // This cannot be false or return a segment, since nothing is collinear
      // but an intersection has occurred.
      const Point_3* intersection_point = boost::get<Point_3>(&*intersection_01_23);
      triangles_.reserve(2);
      triangles_.push_back(Triangle_3(vertex(0), *intersection_point, vertex(3)));
      triangles_.push_back(Triangle_3(vertex(1), vertex(2), *intersection_point));
      return;
    }

    // Case 2b
    // No collinear triplets, but e12 and e30 intersect
    const auto intersection_12_30 = intersection(Segment_3(vertex(1), vertex(2)), Segment_3(vertex(3), vertex(0)));
    if( intersection_12_30 ) {
      // This cannot be false or return a segment, since nothing is collinear
      // but an intersection has occurred.
      const Point_3* intersection_point = boost::get<Point_3>(&*intersection_12_30);
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

template <class R>
void BilinearPatchC3<R>::populate_boundary_()
{
  boundary_.reserve(4);
  boundary_.push_back(Segment_3(vertex(0), vertex(1)));
  boundary_.push_back(Segment_3(vertex(1), vertex(2)));
  boundary_.push_back(Segment_3(vertex(2), vertex(3)));
  boundary_.push_back(Segment_3(vertex(3), vertex(0)));
}



} //namespace CGAL

#endif // CGAL_CARTESIAN_BILINEAR_PATCH_3_H
