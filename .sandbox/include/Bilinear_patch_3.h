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

#ifndef BILINEAR_PATCH_3_H
#define BILINEAR_PATCH_3_H

#include <CGAL/Handle_for.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/array.h>
#include <iostream>

namespace CGAL {

template <class R_>
class BilinearPatchC3
{

  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Triangle_3           Triangle_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Tetrahedron_3        Tetrahedron_3;
  typedef          Interval_nt_advanced     IA_NT;

  typedef std::array<Point_3, 4>                  Rep;
  typedef typename R_::template Handle<Rep>::type Base;

  Base base;
  
  Origin origin = ::CGAL::ORIGIN;
  bool PQR_COLLINEAR_{false};
  bool QRS_COLLINEAR_{false};
  bool IS_PLANAR_{false};
  bool IS_SEGMENTS_{false};
  bool IS_DEGENERATE_{false};

public:
  std::vector<Segment_3> segments_;
  std::vector<Triangle_3> triangles_;
  typedef R_                                R;
  typedef typename R::Tetrahedron_3         Tetrahedron;

  Tetrahedron bounding_tetrahedron;

  BilinearPatchC3() {}
  // NOTE:  specify the points as a cycle around the perimeter!
  //                           q00,              q01,              q11,              q10
  BilinearPatchC3(const Point_3 &p, const Point_3 &q, const Point_3 &r, const Point_3 &s)
    : base(CGAL::make_array(p, q, r, s)) 
    {
      bounding_tetrahedron = Tetrahedron(p, q, r, s);

      IS_SEGMENTS_ = collinear(p,q,r) && collinear(q,r,s);
      IS_PLANAR_ = bounding_tetrahedron.is_degenerate() && !IS_SEGMENTS_;
      IS_DEGENERATE_ = IS_PLANAR_ || IS_SEGMENTS_;

      if( IS_PLANAR_ ) {
        triangles_.reserve(2);
        triangles_.push_back(R::Triangle_3(vertex(0), vertex(1), vertex(2)));
        triangles_.push_back(R::Triangle_3(vertex(2), vertex(3), vertex(0)));
      }

      if( IS_SEGMENTS_ ) {
        segments_.reserve(4);
        segments_.push_back(R::Segment_3(vertex(0), vertex(1)));
        segments_.push_back(R::Segment_3(vertex(1), vertex(2)));
        segments_.push_back(R::Segment_3(vertex(2), vertex(3)));
        segments_.push_back(R::Segment_3(vertex(3), vertex(0)));
      }
    }

  bool  operator==(const BilinearPatchC3 &bp) const;
  bool  operator!=(const BilinearPatchC3 &bp) const;

  friend std::ostream& operator<<(std::ostream& os, BilinearPatchC3 const& bp)
  {
      return (os << bp.vertex(0) << "\n" << bp.vertex(1) << "\n"<< bp.vertex(2) << "\n"<< bp.vertex(3) << "\n" );
  } 

  bool  has_on(const Point_3 &p) const;
  bool  is_degenerate() const;
  bool  is_planar() const;

  const Point_3 & vertex(int i) const;
  const Point_3 & operator[](int i) const;
  const Tetrahedron_3 & tetrahedron() const;

  Point_3 operator()(const FT& u, const FT& v) const;

  IA_NT aux_phi(const Point_3 & x) const;

private:
  IA_NT aux_h12(const Point_3 & x) const;
  IA_NT aux_h03(const Point_3 & x) const;
  IA_NT aux_g( const Point_3 & x, const Point_3 & p, const Point_3 & q, const Point_3 & r ) const;
};

template <class R>
auto BilinearPatchC3<R>::operator()(const FT& u, const FT& v) const -> Point_3
{
  FT ONE{1.};

  Vector_3 q00 = vertex(0) - origin;
  Vector_3 q01 = vertex(1) - origin;
  Vector_3 q11 = vertex(2) - origin;
  Vector_3 q10 = vertex(3) - origin;

  Vector_3 interpolant = (
      (ONE-u)*(ONE-v)*q00
    +     (ONE-u)*(v)*q01
    +     (u)*(ONE-v)*q10
    +         (u)*(v)*q11
  );

  return origin + interpolant;
}


template < class R >
Interval_nt_advanced
BilinearPatchC3<R>::aux_g(const Point_3 & x, const Point_3 & p, const Point_3 & q, const Point_3 & r) const
{

  typename R::Vector_3 rp = R::Vector_3(p, r);
  typename R::Vector_3 qp = R::Vector_3(p, q);
  typename R::Vector_3 xp = R::Vector_3(p, x);

  Interval_nt_advanced xp_x(xp.x());
  Interval_nt_advanced xp_y(xp.y()); 
  Interval_nt_advanced xp_z(xp.z());

  Interval_nt_advanced qp_x(qp.x());
  Interval_nt_advanced qp_y(qp.y());
  Interval_nt_advanced qp_z(qp.z());

  Interval_nt_advanced rp_x(rp.x());
  Interval_nt_advanced rp_y(rp.y());
  Interval_nt_advanced rp_z(rp.z());

  Interval_nt_advanced cross_x = qp_y*rp_z - qp_z*rp_y;
  Interval_nt_advanced cross_y = qp_z*rp_x - qp_x*rp_z;
  Interval_nt_advanced cross_z = qp_x*rp_y - qp_y*rp_x;

  Interval_nt_advanced result = xp_x*cross_x + xp_y*cross_y + xp_z*cross_z;

  return result;
}

template < class R >
Interval_nt_advanced
BilinearPatchC3<R>::aux_h12(const Point_3 & x)  const
{
  return aux_g(x, vertex(0), vertex(1), vertex(3))*aux_g(x, vertex(1), vertex(2), vertex(3));
}

template < class R >
Interval_nt_advanced
BilinearPatchC3<R>::aux_h03(const Point_3 & x) const
{
  return aux_g(x, vertex(0), vertex(1), vertex(2))*aux_g(x, vertex(0), vertex(2), vertex(3));
}

template < class R >
Interval_nt_advanced
BilinearPatchC3<R>::aux_phi(const Point_3 & x) const
{
  return aux_h12(x) - aux_h03(x);
}

template < class R >
bool
BilinearPatchC3<R>::operator==(const BilinearPatchC3<R> &bp) const
{
  if (CGAL::identical(base, bp.base))
      return true;

  int i;
  for(i=0; i<4; i++)
    if ( vertex(0) == bp.vertex(i) )
       break;

  // I assumed this test was to ensure that the set of vertices is not only the same,
  // but that the orientation is the same as well
  return (i<4) && vertex(1) == bp.vertex(i+1) && vertex(2) == bp.vertex(i+2) && vertex(3) == bp.vertex(i+3);
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

template < class R >
inline
bool
BilinearPatchC3<R>::has_on(const Point_3 &p) const
{
  // The patch is a collection of segments...
  if( IS_SEGMENTS_ ) { 
    return (
          segments_[0].has_on(p)
      ||  segments_[1].has_on(p)
      ||  segments_[2].has_on(p)
      ||  segments_[3].has_on(p)
    );
  }
  
  // If three of the points are collinear, then the
  // function phi == 0 everywhere, and we have to partition
  // the surface into two triangles (it's necessarily planar)
  if( IS_PLANAR_ ) {
  // std::cout << "Planar degeneracy\n";
    return (
          triangles_[0].has_on(p)
      ||  triangles_[1].has_on(p)
    );
  }

  // The only remaining possibility is that it's a proper bilinear patch, then the function
  // phi == 0, if and only if p is on the patch.
  return abs(to_double(aux_phi(p))) < 1e-14; 
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

} //namespace CGAL

#endif // CGAL_CARTESIAN_BILINEAR_PATCH_3_H
