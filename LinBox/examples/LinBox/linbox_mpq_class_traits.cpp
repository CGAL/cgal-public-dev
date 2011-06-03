#if !defined CGAL_USE_GMPXX && !defined CGAL_USE_TOPCOM
#error This needs GMPXX or TOPCOM
#else

#include <iostream>
#include <vector>
#include <CGAL/LinBox/linbox_mpq_class_field.h>
#include <CGAL/LinBox/Linear_algebra_traits_linbox.h>
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
        typedef CGAL::Linear_algebra_traits_linbox<Field>       LA;
        typedef typename LA::Matrix                             Matrix;
        typedef typename LA::Vector                             Vector;

        Field F;
        std::vector<Vector> A;
        size_t cols=0,rows=1,ker_dim;

        LA traits;

        Vector v[3];
        std::vector<Vector> vectors;
        v[0].push_back(13);v[0].push_back(1);v[0].push_back(12);
        v[1].push_back(10);v[1].push_back(8);v[1].push_back(2);
        v[2].push_back(1);v[2].push_back(12);v[2].push_back(-11);
        vectors.push_back(v[0]);
        vectors.push_back(v[1]);
        vectors.push_back(v[2]);
        Matrix M(vectors);

        // compute rank
        int r=traits.rank(M);

        // compute determinant
        Element D;
        D=traits.determinant(M);

        // compute nullspace
        Matrix spanning_vectors;
        ker_dim=traits.homogeneous_linear_solver(M,spanning_vectors);

        // show results
        std::cout<<"M = "<<M<<std::endl;
        std::cout<<"rank is "<<r<<std::endl;
        std::cout<<"determinant is "<<D<<std::endl;
        std::cout<<"kernel is "<<spanning_vectors<<std::endl;

        return 0;
}
#endif // !defined CGAL_USE_GMPXX && !defined CGAL_USE_TOPCOM
