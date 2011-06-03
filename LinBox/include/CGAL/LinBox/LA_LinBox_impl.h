// $URL$
// $Id$
//
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_LINBOX_LA_LINBOX_IMPL_H
#define CGAL_LINBOX_LA_LINBOX_IMPL_H

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
        typename LA_LinBox<FT,AL>::Matrix
        LA_LinBox<FT,AL>::
        set_matrix_entry(Matrix &M,unsigned i,unsigned j,const FT &value){
                M.setEntry(i,j,value);
                return M;
        }

        template <class FT,class AL>
        typename LA_LinBox<FT,AL>::Matrix
        LA_LinBox<FT,AL>::
        transpose(const Matrix &M){
                std::vector<Vector> Mrows;
                for(int i=0;i<M.row_dimension();++i)
                        Mrows.push_back(M.row(i));
                return Matrix(Mrows);
        }

        template <class FT,class AL>
        bool
        LA_LinBox<FT,AL>::
        inverse(const Matrix &M,Matrix &I,FT &D,Vector &c){
                // TODO
                CGAL_error_msg("not implemented");
        }

        // XXX for the moment, this works only with rational matrices
        template <class FT,class AL>
        typename LA_LinBox<FT,AL>::Matrix
        LA_LinBox<FT,AL>::
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
        typename LA_LinBox<FT,AL>::FT
        LA_LinBox<FT,AL>::
        determinant(const Matrix &M,Matrix &L,Matrix &U,
                    std::vector<int> &q,Vector &c){
                // TODO
                CGAL_error_msg("not implemented");
        }

        template <class FT,class AL>
        bool
        LA_LinBox<FT,AL>::
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
        typename LA_LinBox<FT,AL>::FT
        LA_LinBox<FT,AL>::
        determinant(const Matrix &M,CGAL::Method method){
                CGAL_assertion_msg(M.column_dimension()==M.row_dimension(),
                                   "matrix is not square");
                FT d;
                switch(method){
                        case CGAL_HYBRID:
                                LinBox::det(d,M,
                                            typename LinBox::ClassifyRing<LS>::
                                                categoryTag(),
                                            LinBox::Method::Hybrid());
                                break;
                        case CGAL_BLACKBOX:
                                LinBox::det(d,M,
                                            typename LinBox::ClassifyRing<LS>::
                                                categoryTag(),
                                            LinBox::Method::Blackbox());
                                break;
                        case CGAL_WIEDEMANN:
                                LinBox::det(d,M,
                                            typename LinBox::ClassifyRing<LS>::
                                                categoryTag(),
                                            LinBox::Method::Wiedemann());
                                break;
                        case CGAL_ELIMINATION:
                                LinBox::det(d,M,
                                            typename LinBox::ClassifyRing<LS>::
                                                categoryTag(),
                                            LinBox::Method::Elimination());
                                break;
                        case CGAL_BLAS_ELIMINATION:
                                LinBox::det(d,M,
                                            typename LinBox::ClassifyRing<LS>::
                                                categoryTag(),
                                            LinBox::Method::BlasElimination());
                                break;
                        case CGAL_SPARSE_ELIMINATION:
                                LinBox::det(d,M,
                                            typename LinBox::ClassifyRing<LS>::
                                                categoryTag(),
                                            LinBox::Method::SparseElimination()
                                           );
                                break;
                        default:
                                CGAL_assertion(method==CGAL_DEFAULT);
                                // this comes from experimental results
                                if(M.row_dimension()<18)
                                        LinBox::det(d,M,
                                        typename LinBox::ClassifyRing<LS>::
                                                categoryTag(),
                                        LinBox::Method::SparseElimination());
                                else
                                        LinBox::det(d,M,
                                            typename LinBox::ClassifyRing<LS>::
                                                categoryTag(),
                                            LinBox::Method::Wiedemann());
                }
                return d;
        }

        template <class FT,class AL>
        Sign
        LA_LinBox<FT,AL>::
        sign_of_determinant(const Matrix &M,CGAL::Method method){
                CGAL_assertion_msg(M.column_dimension()==M.row_dimension(),
                                   "matrix is not square");
                return determinant(M,method).sign();
        }

        template <class FT,class AL>
        bool
        LA_LinBox<FT,AL>::
        linear_solver(const Matrix &M,const Vector &b,
                      Vector &x,FT &D,Matrix &spanning_vectors,
                      Vector &c){
                // TODO
                CGAL_error_msg("not implemented");
        }

        template <class FT,class AL>
        bool
        LA_LinBox<FT,AL>::
        linear_solver(const Matrix &M,const Vector &b,
                      Vector &x,FT &D,Vector &c){
                // TODO
                CGAL_error_msg("not implemented");
        }

        template <class FT,class AL>
        bool
        LA_LinBox<FT,AL>::
        linear_solver(const Matrix &M,
                      const Vector &b,
                      Vector &x,
                      FT &D,
                      CGAL::Method method){
                CGAL_assertion(M.row_dimension()==b.dimension());
                // TODO: choose method
                //D=FT(1);
                //LinBox::solve(x,D,M,b,LinBox::BlasEliminationTraits());
                std::cout<<"solve\nb="<<b<<"\nM="<<M<<std::endl;
                LinBox::solve(x,
                              M,
                              b,
                              LinBox::Method::Hybrid());
                std::cout<<"result is "<<x<<std::endl;
                return !x.is_zero();
        }

        template <class FT,class AL>
        bool
        LA_LinBox<FT,AL>::
        is_solvable(const Matrix &M,const Vector &b){
                // TODO
                CGAL_error_msg("not implemented");
        }

        template <class FT,class AL>
        bool
        LA_LinBox<FT,AL>::
        homogeneous_linear_solver(const Matrix &M,Vector &x){
                // TODO
                CGAL_error_msg("not implemented");
        }

        // TODO: this doesn't work for integer fields (I think it is
        // because of the division)
        template <class FT,class AL>
        int
        LA_LinBox<FT,AL>::
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
        LA_LinBox<FT,AL>::
        independent_columns(const Matrix &M,std::vector<int> &q){
                // TODO
                CGAL_error_msg("not implemented");
        }

        template <class FT,class AL>
        int
        LA_LinBox<FT,AL>::
        rank(const Matrix &M,CGAL::Method method){
                unsigned long result;
                switch(method){
                        case CGAL_DEFAULT:
                                LinBox::rank(result,M);
                                break;
                        case CGAL_WIEDEMANN:
                                LinBox::rank(result,
                                             M,
                                             typename LinBox::ClassifyRing<LS>::
                                                categoryTag(),
                                             LinBox::Method::Wiedemann());
                                break;
                        case CGAL_BLAS_ELIMINATION:
                                LinBox::rank(result,
                                             M,
                                             typename LinBox::ClassifyRing<LS>::
                                                categoryTag(),
                                             LinBox::Method::BlasElimination());
                                break;
                        case CGAL_SPARSE_ELIMINATION:
                                LinBox::rank(result,
                                             M,
                                             typename LinBox::ClassifyRing<LS>::
                                                categoryTag(),
                                             LinBox::Method::
                                                SparseElimination());
                                break;
                        default:
                                CGAL_error_msg("not implemented");
                }
                CGAL_assertion_msg((unsigned long)((int)result)==result,
                                   "the result does not fit in an int");
                return (int)result;
        }

} // namespace CGAL

#endif // CGAL_LINBOX_LA_LINBOX_IMPL_H
