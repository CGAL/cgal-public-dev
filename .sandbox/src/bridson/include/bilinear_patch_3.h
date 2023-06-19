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

  typedef std::array<Point_3, 4>          Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                        R;
  typedef typename R::Tetrahedron_3 Tetrahedron;

  Tetrahedron bounding_tetrahedron;

  BilinearPatchC3() {}

  BilinearPatchC3(const Point_3 &p, const Point_3 &q, const Point_3 &r, const Point_3 &s)
    : base(CGAL::make_array(p, q, r, s)) {
      bounding_tetrahedron = Tetrahedron(p, q, r, s);
    }

  bool  operator==(const BilinearPatchC3 &bp) const;
  bool  operator!=(const BilinearPatchC3 &bp) const;

  // bool  has_on(const Point_3 &p) const;
  bool  is_degenerate() const;
  bool  is_planar() const;

  const Point_3 & vertex(int i) const;
  const Point_3 & operator[](int i) const;
  const Tetrahedron_3 & tetrahedron() const;

  FT phi(const Point_3 & x);

private:
  FT h12(const Point_3 & x);
  FT h03(const Point_3 & x);
  FT g( const Point_3 & x, const Point_3 & p, const Point_3 & q, const Point_3 & r );
};

template < class R >
typename R::FT
BilinearPatchC3<R>::g(const Point_3 & x, const Point_3 & p, const Point_3 & q, const Point_3 & r) 
{
  R::Vector_3 rp = R::Vector_3(p, r);
  R::Vector_3 qp = R::Vector_3(p, q);
  R::Vector_3 xp = R::Vector_3(p, x);

  R::Point_3 cross = R::Point_3(
    qp.y()*rp.z() - qp.z()*rp.y(),
    qp.z()*rp.x() - qp.x()*rp.z(),
    qp.x()*rp.y() - qp.y()*rp.x()
  );

  return xp.x()*cross.x() + xp.y()*cross.y() + xp.z()*cross.z();
}

template < class R >
typename R::FT
BilinearPatchC3<R>::h12(const Point_3 & x) 
{
  return g(x, vertex(0), vertex(1), vertex(2))*g(x, vertex(1), vertex(3), vertex(2));
}

template < class R >
typename R::FT
BilinearPatchC3<R>::h03(const Point_3 & x) 
{
  return g(x, vertex(0), vertex(1), vertex(3))*g(x, vertex(0), vertex(3), vertex(2));
}

template < class R >
typename R::FT
BilinearPatchC3<R>::phi(const Point_3 & x) 
{
  return h12(x) - h03(x);
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

// template < class R >
// inline
// bool
// BilinearPatchC3<R>::
// has_on(const typename BilinearPatchC3<R>::Point_3 &p) const
// {
//   return R().has_on_3_object()
//                (static_cast<const typename R::BilinearPatchC3&>(*this), p);
// }

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
