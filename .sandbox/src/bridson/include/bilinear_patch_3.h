// Copyright (c) 2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jeffrey Cochran

#ifndef CGAL_CARTESIAN_BILINEAR_PATCH_3_H
#define CGAL_CARTESIAN_BILINEAR_PATCH_3_H

#include <CGAL/Handle_for.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/array.h>

namespace CGAL {

template <class R_>
class BilinearPatchC3
{
  
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Triangle_3           Triangle_3;
  typedef typename R_::Tetrahedron_3        Tetrahedron_3;
  typedef typename Interval_nt_advanced     IA_NT;

  typedef std::array<Point_3, 4>                  Rep;
  typedef typename R_::template Handle<Rep>::type Base;

  Base base;

public:
  typedef R_                                R;
  typedef typename R::Tetrahedron_3         Tetrahedron;

  Tetrahedron bounding_tetrahedron;

  BilinearPatchC3() {}

  BilinearPatchC3(const Point_3 &p, const Point_3 &q, const Point_3 &r, const Point_3 &s)
    : base(CGAL::make_array(p, q, r, s)) {
      bounding_tetrahedron = Tetrahedron(p, q, r, s);
    }

  bool  operator==(const BilinearPatchC3 &bp) const;
  bool  operator!=(const BilinearPatchC3 &bp) const;

  bool  has_on(const Point_3 &p) const;
  bool  is_degenerate() const;
  bool  is_planar() const;

  const Point_3 & vertex(int i) const;
  const Point_3 & operator[](int i) const;
  const Tetrahedron_3 & tetrahedron() const;

  IA_NT aux_phi(const Point_3 & x) const;

private:
  IA_NT aux_h12(const Point_3 & x) const;
  IA_NT aux_h03(const Point_3 & x) const;
  IA_NT aux_g( const Point_3 & x, const Point_3 & p, const Point_3 & q, const Point_3 & r ) const;
};

template < class R >
Interval_nt_advanced
BilinearPatchC3<R>::aux_g(const Point_3 & x, const Point_3 & p, const Point_3 & q, const Point_3 & r) const
{

  R::Vector_3 rp = R::Vector_3(p, r);
  R::Vector_3 qp = R::Vector_3(p, q);
  R::Vector_3 xp = R::Vector_3(p, x);

  Interval_nt_advanced xp_x = Interval_nt_advanced(xp.x());
  Interval_nt_advanced xp_y = Interval_nt_advanced(xp.y());
  Interval_nt_advanced xp_z = Interval_nt_advanced(xp.z());

  Interval_nt_advanced qp_x = Interval_nt_advanced(qp.x());
  Interval_nt_advanced qp_y = Interval_nt_advanced(qp.y());
  Interval_nt_advanced qp_z = Interval_nt_advanced(qp.z());

  Interval_nt_advanced rp_x = Interval_nt_advanced(rp.x());
  Interval_nt_advanced rp_y = Interval_nt_advanced(rp.y());
  Interval_nt_advanced rp_z = Interval_nt_advanced(rp.z());

  Interval_nt_advanced cross_x = qp_y*rp_z - qp_z*rp_y;
  Interval_nt_advanced cross_y = qp_z*rp_x - qp_x*rp_z;
  Interval_nt_advanced cross_z = qp_x*rp_y - qp_y*rp_x;

  return xp_x*cross_x + xp_y*cross_y + xp_z*cross_z;
}

template < class R >
Interval_nt_advanced
BilinearPatchC3<R>::aux_h12(const Point_3 & x)  const
{
  return aux_g(x, vertex(0), vertex(1), vertex(2))*aux_g(x, vertex(1), vertex(3), vertex(2));
}

template < class R >
Interval_nt_advanced
BilinearPatchC3<R>::aux_h03(const Point_3 & x) const
{
  return aux_g(x, vertex(0), vertex(1), vertex(3))*aux_g(x, vertex(0), vertex(3), vertex(2));
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
  return aux_phi(p).do_overlap(0);
}

template < class R >
bool
BilinearPatchC3<R>::is_degenerate() const
{
  return (collinear(vertex(0),vertex(1),vertex(2)) && collinear(vertex(1),vertex(2),vertex(3)));
}

template < class R >
bool
BilinearPatchC3<R>::is_planar() const
{
  return bounding_tetrahedron.is_degenerate();
}

template < class R >
const typename BilinearPatchC3<R>::Tetrahedron_3 &
BilinearPatchC3<R>::tetrahedron() const
{
  return bounding_tetrahedron;
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_BILINEAR_PATCH_3_H
