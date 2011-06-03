// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_LINBOX_LA_LINBOX_H
#define CGAL_LINBOX_LA_LINBOX_H

#include <CGAL/LinBox/methods.h>
#include <CGAL/LinBox/dense_matrix.h>

namespace CGAL{

        // _LS is a LinBox algebraic structure and _AL is an allocator,
        // which defaults to the CGAL allocator.
        template<class _LS,class _AL=CGAL_ALLOCATOR(_LS)>
        class LA_LinBox{
                public:
                typedef _AL                                     AL;
                typedef _LS                                     LS;
                typedef CGAL::Linbox_dense_matrix<LS>           Matrix;
                typedef typename Matrix::Vector                 Vector;
                typedef typename LS::Element                    FT;
                typedef FT                                      RT;
                typedef LA_LinBox<LS,AL>                        Self;
                typedef const FT*                               const_iterator;
                typedef FT*                                     iterator;
                private:
                typedef typename Matrix::Identity               Identity;
                public:
                LA_LinBox(){};

                static Matrix transpose(const Matrix&);
                static bool inverse(const Matrix&,Matrix&,FT&,Vector&);
                static Matrix inverse(const Matrix&,RT&);
                static FT determinant(const Matrix&,Matrix&,Matrix&,
                                      std::vector<int>&,Vector&);
                static bool verify_determinant(const Matrix&,const RT&,
                                               const Matrix&,const Matrix&,
                                               const std::vector<int>&,
                                               const Vector&);
                static FT determinant(const Matrix&,
                                      CGAL::Method=CGAL_DEFAULT);
                static Sign sign_of_determinant(const Matrix&,
                                                CGAL::Method=CGAL_DEFAULT);
                static bool linear_solver(const Matrix&,const Vector&,
                                          Vector&,FT&,Matrix&,Vector&);
                static bool linear_solver(const Matrix&,const Vector&,Vector&,
                                          FT&,Vector&);
                static bool linear_solver(const Matrix&,
                                          const Vector&,
                                          Vector&,
                                          FT&,
                                          CGAL::Method=CGAL_DEFAULT);
                static bool is_solvable(const Matrix&,const Vector&);
                static bool homogeneous_linear_solver(const Matrix&,Vector&);
                static int homogeneous_linear_solver(const Matrix&,Matrix&);
                static int independent_columns(const Matrix&,std::vector<int>&);
                static int rank(const Matrix&,CGAL::Method=CGAL_DEFAULT);
        };

} // namespace CGAL

#include <CGAL/LinBox/LA_LinBox_impl.h>

#endif // CGAL_LINBOX_LA_LINBOX_H
