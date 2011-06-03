#include <CGAL/config.h>

#if !defined CGAL_USE_GMPXX && !defined CGAL_USE_TOPCOM
#error This needs GMPXX or TOPCOM
#else

#include <iostream>
#include <vector>
#include <CGAL/LinBox/mpq_class_field.h>
#include <CGAL/LinBox/LA_LinBox.h>
#ifdef CGAL_USE_GMPXX
#include <gmpxx.h>
#else
#include <Rational.h>
#endif

int main(){
#ifdef CGAL_USE_GMPXX
        typedef mpq_class                                       FT;
#else
        typedef Rational                                        FT;
#endif
        typedef CGAL::Linbox_rational_field<FT>                 Field;
        typedef Field::Element                                  Element;
        typedef CGAL::LA_LinBox<Field>                          LA;
        typedef LA::Matrix                                      Matrix;
        typedef LA::Vector                                      Vector;

        Field F;
        std::vector<Vector> A;
        size_t cols=0,rows=1,ker_dim;

        LA traits;

        Matrix M(3,3);
        traits.set_matrix_entry(M,0,0,13);
        traits.set_matrix_entry(M,0,1,10);
        traits.set_matrix_entry(M,0,2,1);
        traits.set_matrix_entry(M,1,0,1);
        traits.set_matrix_entry(M,1,1,8);
        traits.set_matrix_entry(M,1,2,12);
        traits.set_matrix_entry(M,2,0,12);
        traits.set_matrix_entry(M,2,1,2);
        traits.set_matrix_entry(M,2,2,-11);

        // compute rank
        //int r=traits.rank(M);

        // compute determinant
        Element D;
        D=traits.determinant(M);

        // compute nullspace
        Matrix spanning_vectors;
        ker_dim=traits.homogeneous_linear_solver(M,spanning_vectors);

        // show results
        std::cout<<"M = "<<M<<std::endl;
        //std::cout<<"rank is "<<r<<std::endl;
        std::cout<<"determinant is "<<D<<std::endl;
        std::cout<<"kernel is "<<spanning_vectors<<std::endl;

        return 0;
}
#endif // !defined CGAL_USE_GMPXX && !defined CGAL_USE_TOPCOM
