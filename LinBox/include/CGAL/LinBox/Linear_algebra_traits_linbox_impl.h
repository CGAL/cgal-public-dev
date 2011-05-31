// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_LINBOX_LINEAR_ALGEBRA_TRAITS_LINBOX_IMPL_H
#define CGAL_LINBOX_LINEAR_ALGEBRA_TRAITS_LINBOX_IMPL_H

#include <linbox/blackbox/blas-blackbox.h>
#include <linbox/algorithms/matrix-inverse.h>
#include <linbox/solutions/det.h>
#include <linbox/solutions/methods.h>
#include <linbox/solutions/solve.h>
//#include <linbox/solutions/nullspace.h> // for LinBox < 1.2
#include <linbox/algorithms/dense-nullspace.h> // for LinBox 1.2
#include <linbox/algorithms/linbox-tags.h> // for LinBox 1.2
#include <linbox/solutions/rank.h>

namespace CGAL {

        template <class FT,class AL>
        typename Linear_algebra_traits_linbox<FT,AL>::Matrix
        Linear_algebra_traits_linbox<FT,AL>::
        transpose(const Matrix &M){
                std::vector<Vector> Mrows;
                for(int i=0;i<M.row_dimension();++i)
                        Mrows.push_back(M.row(i));
                return Matrix(Mrows);
        }

        template <class FT,class AL>
        bool
        Linear_algebra_traits_linbox<FT,AL>::
        inverse(const Matrix &M,Matrix &I,FT &D,Vector &c){
                // TODO
                CGAL_error_msg("not implemented");
        }

        // XXX for the moment, this works only with rational matrices
        template <class FT,class AL>
        typename Linear_algebra_traits_linbox<FT,AL>::Matrix
        Linear_algebra_traits_linbox<FT,AL>::
        inverse(const Matrix &M,FT &D){
                CGAL_assertion_msg(M.column_dimension()==M.row_dimension(),
                                   "matrix is not square");
                CGAL_assertion_msg(Self::determinant(M)!=FT(0),
                                   "determinant should not be zero");
                Matrix A(M);
                CGAL_assertion_code(int result=)
                LinBox::MatrixInverse::matrixInverseIn(A.get_field(),A);
                CGAL_assertion(result==0);
                CGAL_assertion(M*A==Matrix(M.row_dimension(),Identity()));
                D=FT(1);
                return A;
        }

        template <class FT,class AL>
        typename Linear_algebra_traits_linbox<FT,AL>::FT
        Linear_algebra_traits_linbox<FT,AL>::
        determinant(const Matrix &M,Matrix &L,Matrix &U,
                    std::vector<int> &q,Vector &c){
                // TODO
                CGAL_error_msg("not implemented");
        }

        template <class FT,class AL>
        bool
        Linear_algebra_traits_linbox<FT,AL>::
        verify_determinant(const Matrix &M,
                           const FT &D,
                           const Matrix &L,
                           const Matrix &U,
                           const std::vector<int> &q,
                           const Vector &c){
                // TODO
                CGAL_error_msg("not implemented");
        }

        template <class FT,class AL>
        typename Linear_algebra_traits_linbox<FT,AL>::FT
        Linear_algebra_traits_linbox<FT,AL>::
        determinant(const Matrix &M){
                CGAL_assertion_msg(M.column_dimension()==M.row_dimension(),
                                   "matrix is not square");
                FT d;
                // TODO: Choose a method, between LinBox::Method::Hybrid
                // (default), Blackbox, Elimination, Wiedemann,
                // BlasElimination, SparseElimination.
                LinBox::det(d,M,
                            typename LinBox::ClassifyRing<LS>::categoryTag(),
                            LinBox::Method::Wiedemann());
                return d;
        }

        template <class FT,class AL>
        Sign
        Linear_algebra_traits_linbox<FT,AL>::
        sign_of_determinant(const Matrix &M){
                CGAL_assertion_msg(M.column_dimension()==M.row_dimension(),
                                   "matrix is not square");
                return determinant(M).sign();
        }

        template <class FT,class AL>
        bool
        Linear_algebra_traits_linbox<FT,AL>::
        linear_solver(const Matrix &M,const Vector &b,
                      Vector &x,FT &D,Matrix &spanning_vectors,
                      Vector &c){
                // TODO
                CGAL_error_msg("not implemented");
        }

        template <class FT,class AL>
        bool
        Linear_algebra_traits_linbox<FT,AL>::
        linear_solver(const Matrix &M,const Vector &b,
                      Vector &x,FT &D,Vector &c){
                // TODO
                CGAL_error_msg("not implemented");
        }

        template <class FT,class AL>
        bool
        Linear_algebra_traits_linbox<FT,AL>::
        linear_solver(const Matrix &M,const Vector &b,Vector &x,FT &D){
                CGAL_assertion(M.row_dimension()==b.dimension());
                // TODO: choose method
                //D=FT(1);
                LinBox::solve(x,D,M,b,LinBox::BlasEliminationTraits());
                return !x.is_zero();
        }

        template <class FT,class AL>
        bool
        Linear_algebra_traits_linbox<FT,AL>::
        is_solvable(const Matrix &M,const Vector &b){
                // TODO
                CGAL_error_msg("not implemented");
        }

        template <class FT,class AL>
        bool
        Linear_algebra_traits_linbox<FT,AL>::
        homogeneous_linear_solver(const Matrix &M,Vector &x){
                // TODO
                CGAL_error_msg("not implemented");
        }

        // TODO: this doesn't work for integer fields (I think it is
        // because of the division)
        template <class FT,class AL>
        int
        Linear_algebra_traits_linbox<FT,AL>::
        homogeneous_linear_solver(const Matrix &M,Matrix &spanning_vectors){
                size_t ldk,ker_dim;
                FT *kernel;
                Matrix temp(M);
                LinBox::NullSpaceBasis(M.get_field(),
                                       // LinBox::FFLAS::FflasRight, // <1.2
                                       LinBox::LinBoxTag::Right, // LinBox 1.2
                                       M.row_dimension(), // rows
                                       M.column_dimension(), // columns
                                       temp.FullIterator(),
                                       M.column_dimension(),// leading dim of M
                                       kernel,
                                       ldk,
                                       ker_dim);
                spanning_vectors=Matrix(M.column_dimension(),ker_dim);
                for(size_t col=0;col<ker_dim;++col)
                        for(size_t row=0;row<M.column_dimension();++row)
                                spanning_vectors.setEntry(row,
                                                          col,
                                                          kernel[row*ldk+col]);
                return ker_dim;
        }

        template <class FT,class AL>
        int
        Linear_algebra_traits_linbox<FT,AL>::
        independent_columns(const Matrix &M,std::vector<int> &q){
                // TODO
                CGAL_error_msg("not implemented");
        }

        template <class FT,class AL>
        int
        Linear_algebra_traits_linbox<FT,AL>::
        rank(const Matrix &M){
                unsigned long result;
                // TODO: Choose a method, between Method::Wiedemann
                // (default), BlasElimination or SparseElimination.
                LinBox::rank(result,M);
                CGAL_assertion_msg((unsigned long)((int)result)==result,
                                   "the result does not fit in an int");
                return (int)result;
        }

} // namespace CGAL

#endif // CGAL_LINBOX_LINEAR_ALGEBRA_TRAITS_LINBOX_IMPL_H
