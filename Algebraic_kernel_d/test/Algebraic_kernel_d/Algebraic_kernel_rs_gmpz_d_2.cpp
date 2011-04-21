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

#include <CGAL/basic.h>

#if defined(CGAL_USE_GMP) && defined(CGAL_USE_MPFI) && defined(CGAL_USE_RS3)

#include <CGAL/Algebraic_kernel_rs_gmpz_d_2.h>

int main(){
        typedef CGAL::Algebraic_kernel_rs_gmpz_d_2      AK;
        typedef AK::Polynomial_2                        Polynomial_2;
        typedef AK::Coefficient                         Coefficient;
        typedef AK::Algebraic_real_2                    Algebraic_real_2;
        typedef AK::Bound                               Bound;
        typedef AK::Multiplicity_type                   Multiplicity_type;
        typedef CGAL::Polynomial_traits_d<Polynomial_2> Polynomial_traits_2;

        AK ak; // an object of Algebraic_kernel_rs_gmpz_d_2

        // TODO: this is the entire test; when the kernel is complete, this
        // line should be uncommented and the rest of the file eliminated
        //CGAL::test_algebraic_kernel_2<AK>(ak);

        CGAL::set_pretty_mode(std::cout);

        // create a Solve_2 object
        AK::Solve_2 solve_2 = ak.solve_2_object();

        // create a polynomial system to solve
        Polynomial_traits_2::Construct_polynomial construct_polynomial_2;
        std::pair<CGAL::Exponent_vector,Coefficient> coeffs_x[1]
                = {std::make_pair(CGAL::Exponent_vector(1,0),Coefficient(1))};
        Polynomial_2 x=construct_polynomial_2(coeffs_x,coeffs_x+1);
        std::pair<CGAL::Exponent_vector,Coefficient> coeffs_y[1]
                = {std::make_pair(CGAL::Exponent_vector(0,1),Coefficient(1))};
        Polynomial_2 y=construct_polynomial_2(coeffs_y,coeffs_y+1);
        Polynomial_2 f=x-CGAL::ipower(y,2);     // f(x,y)=x-y^2
        Polynomial_2 g=y-1;                     // g(x,y)=y-1
        std::cout<<"f(x,y)="<<f<<std::endl;
        std::cout<<"g(x,y)="<<g<<std::endl;

        // solve it in R^2 and in the box [0,1]*[2,3]
        typedef std::vector<
                std::pair<Algebraic_real_2,Multiplicity_type> > rvec;
        rvec roots1,roots2;
        // this is the first flavor of the Solve_2 functor
        solve_2(f,g,std::back_inserter(roots1));
        // and this is the second one
        solve_2(f,g,
                Bound(0),Bound(1),Bound(2),Bound(3),
                std::back_inserter(roots2));

        // show the solutions
        std::cout<<"solutions are:"<<std::endl;
        for(rvec::iterator i=roots1.begin();i!=roots1.end();++i)
                std::cout<<i->first<<" mult "<<i->second<<std::endl;
        std::cout<<"solutions in the box are:"<<std::endl;
        for(rvec::iterator i=roots2.begin();i!=roots2.end();++i)
                std::cout<<i->first<<" mult "<<i->second<<std::endl;

        return 0;
}
#else
int main(){
        return 0;
}
#endif
