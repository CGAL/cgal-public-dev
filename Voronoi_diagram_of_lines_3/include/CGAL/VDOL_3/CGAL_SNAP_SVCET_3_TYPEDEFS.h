// Copyright (c) 2005-2007 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id:
// 
//
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>
//            Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//

// macro for internal usage, lists every public type but the functors 

#ifndef CGAL_SNAP_SVCET_3_TYPEDEFS
#define CGAL_SNAP_SVCET_3_TYPEDEFS                                      \
  typedef typename SVCET_3::Surface_3             Surface_3;            \
  typedef typename SVCET_3::Xy_monotone_surface_3 Xy_monotone_surface_3; \
                                                                        \
  typedef typename SVCET_3::FT                    FT;                   \
  typedef typename SVCET_3::Integer               Integer;              \
  typedef typename SVCET_3::Rational              Rational;             \
                                                                        \
  typedef typename SVCET_3::Algebraic_kernel_d_2     Algebraic_kernel_d_2;  \
  typedef typename SVCET_3::Algebraic_kernel_d_1     Algebraic_kernel_d_1;  \
                                                                        \
  typedef typename SVCET_3::Linear_kernel          Linear_kernel;       \
  typedef typename SVCET_3::Line_3                 Line_3;              \
  typedef typename SVCET_3::Vector_3               Vector_3;            \
  typedef typename SVCET_3::Point_3                Point_3;             \
  typedef typename SVCET_3::Plane_3                Plane_3;             \
                                                                        \
  typedef typename SVCET_3::Curved_kernel_2        Curved_kernel_2;     \
  typedef typename SVCET_3::Point_2                Point_2;             \
  typedef typename SVCET_3::Curve_2                Curve_2;             \
  typedef typename SVCET_3::X_monotone_curve_2     X_monotone_curve_2;  \
  typedef typename SVCET_3::Multiplicity           Multiplicity;        \
  typedef typename SVCET_3::Coordinate_1           Coordinate_1;        \
                                                                        \
  typedef typename SVCET_3::Bisector               Bisector;            \
  typedef typename SVCET_3::Monomial               Monomial;            \
  typedef typename SVCET_3::Poly_int_3             Poly_int_3;          \
  typedef typename SVCET_3::PT_int_3               PT_int_3;            \
  typedef typename SVCET_3::Poly_int_2             Poly_int_2;          \
  typedef typename SVCET_3::PT_int_2               PT_int_2;            \
  typedef typename SVCET_3::Poly_int_1             Poly_int_1;          \
  typedef typename SVCET_3::PT_int_1               PT_int_1;            \
  typedef typename SVCET_3::Poly_rat_3             Poly_rat_3;          \
  typedef typename SVCET_3::PT_rat_3               PT_rat_3;            \
  typedef typename SVCET_3::Poly_rat_2             Poly_rat_2;          \
  typedef typename SVCET_3::PT_rat_2               PT_rat_2;            \
  typedef typename SVCET_3::Poly_rat_1             Poly_rat_1;          \
  typedef typename SVCET_3::PT_rat_1               PT_rat_1;            \
  typedef typename SVCET_3::Poly_lazy_rat_3        Poly_lazy_rat_3;     \
  typedef typename SVCET_3::PT_lazy_rat_3          PT_lazy_rat_3;       \
  typedef typename SVCET_3::Poly_lazy_rat_2        Poly_lazy_rat_2;     \
  typedef typename SVCET_3::PT_lazy_rat_2          PT_lazy_rat_2;       \
  typedef typename SVCET_3::Poly_lazy_rat_1        Poly_lazy_rat_1;     \
  typedef typename SVCET_3::PT_lazy_rat_1          PT_lazy_rat_1;       \
                                                                        \
  typedef typename SVCET_3::Obj_list               Obj_list;            \

#endif
