// Copyright (c) 2011 National and Kapodistrian University of Athens (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
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

#ifndef CGAL_SIMPLE_AK_RS_1
#define CGAL_SIMPLE_AK_RS_1

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Polynomial.h>
#include <CGAL/RS/simple_functors_1.h>
#include <CGAL/RS/simple_rs2_isolator_1.h>
#include <CGAL/RS/simple_rs3_refiner_1.h>
#include <CGAL/RS/simple_ak_1.h>

namespace CGAL{

typedef CGAL::SimpleAK1::Simple_algebraic_kernel_1<
        CGAL::Polynomial<CGAL::Gmpz>,
        CGAL::Gmpfr,
        CGAL::RS2::Simple_rs2_isolator_1<CGAL::Polynomial<CGAL::Gmpz>,
                                         CGAL::Gmpfr>,
        CGAL::RS3::Simple_rs3_refiner_1<CGAL::Polynomial<CGAL::Gmpz>,
                                        CGAL::Gmpfr> >

                Simple_algebraic_kernel_rs_gmpz_d_1;

}

#endif // CGAL_SIMPLE_AK_RS_1
