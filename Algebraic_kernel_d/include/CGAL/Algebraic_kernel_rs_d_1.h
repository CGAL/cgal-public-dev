// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.

// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_ALGEBRAIC_KERNEL_RS_D_1
#define CGAL_ALGEBRAIC_KERNEL_RS_D_1

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Polynomial.h>
#include "RS/rs2_isolator_1.h"
#ifdef CGAL_USE_RS3
#include "RS/rs23_k_isolator_1.h"
#include "RS/rs3_refiner_1.h"
#include "RS/rs3_k_refiner_1.h"
#else
#include "RS/bisection_refiner_1.h"
#endif
#include "RS/ak_z_1.h"

// This file defines the generic RS algebraic kernel. The Bound type is
// fixed to Gmpfr and all the template parameters, except the input
// polynomial type, are fixed. Please note that, for defining a kernel for
// the type Polynomial<Gmpz>, it is useless to use this header file; please
// refer to the file CGAL/Algebraic_kernel_rs_gmpz_d_1.h.

namespace CGAL{

// Choice of the z-isolator: RS default or RS-k.
#ifdef CGAL_RS_USE_K
typedef CGAL::RS23_k_isolator_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                ZIsolator;
#else
typedef CGAL::RS2::RS2_isolator_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                ZIsolator;
#endif

// Choice of the z-refiner: bisection, RS3 or RS3-k.
#ifdef CGAL_USE_RS3
#ifdef CGAL_RS_USE_K
typedef CGAL::RS3::RS3_k_refiner_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                ZRefiner;
#else
typedef CGAL::RS3::RS3_refiner_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                ZRefiner;
#endif
#else
typedef CGAL::Bisection_refiner_1<CGAL::Polynomial<CGAL::Gmpz>,CGAL::Gmpfr>
                                                ZRefiner;
#endif

template <class P_>
class Algebraic_kernel_rs_d_1:
public CGAL::RS_AK1::Algebraic_kernel_z_1<
                P_,
                CGAL::Polynomial<CGAL::Gmpz>,
                CGAL::RS_AK1::Polynomial_converter_1<
                        P_,
                        CGAL::Polynomial<CGAL::Gmpz> >,
                CGAL::Gmpfr,
                ZIsolator,
                ZRefiner>{
};

}

#endif // CGAL_ALGEBRAIC_KERNEL_RS_D_1
