// Copyright (c) 1997-2012  ETH Zurich (Switzerland).
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
// $URL:$
// $Id:$
// 
//
// Author(s)     :  Yves Brise
//                  Bernd Gaertner <gaertner@inf.ethz.ch>

#ifndef CGAL_QP_LU_FACTORIZATION_H
#define CGAL_QP_LU_FACTORIZATION_H

#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/tags.h>
#include <CGAL/QP_solver/QP_sparse_matrix.h>

#include <boost/shared_ptr.hpp>

#include <typeinfo>


namespace CGAL {

// =============================================================================
// helper objects
// =============================================================================
class Init_Natural_Numbers {
  public:
    Init_Natural_Numbers (int start = 1): current(start) {}
    int operator() () { return current++; }
  private:
    int current;
};

  template <typename ET>
  class Dummy_Matrix_Provider {
    
  public:
    Dummy_Matrix_Provider(CGAL::QP_sparse_matrix<ET>& matrix,
                          int csize,
                          int bosize):
    matrix_ptr_(new QP_sparse_matrix<ET>(matrix)),
    csize_(csize),
    bosize_(bosize){}    
    
    boost::shared_ptr<CGAL::QP_sparse_matrix<ET> >
    get_basis_matrix(int& csize, int& bosize, bool is_linear) {
      csize = csize_;
      bosize = bosize_;
      return boost::shared_ptr<CGAL::QP_sparse_matrix<ET> > (matrix_ptr_);
    }
    
  private:
    boost::shared_ptr<CGAL::QP_sparse_matrix<ET> > matrix_ptr_;
    int csize_;
    int bosize_;
  };
  
// Used to manipulate the permutations
struct Conditional_inc : public std::binary_function<int,int,int> {
  int operator() (int k, int a) const {return (a < k ? a : a+1);}
};

struct Conditional_dec : public std::binary_function<int,int,int> {
  int operator() (int k, int a) const {return (a < k ? a : a-1);}
};




// =============================================================================
// declarations
// =============================================================================
template <typename ET, typename Is_LP, typename Matrix_Provider>
class QP_LU_factorization;
template <typename ET, typename Is_LP, typename Matrix_Provider>
std::ostream& operator<<(std::ostream& out,
                   const QP_LU_factorization<ET, Is_LP, Matrix_Provider>& f);


// =============================================================================
// class interface
// =============================================================================
template <typename ET, typename Is_LP, typename Matrix_Provider>
class QP_LU_factorization {

  // ===========================================================================
  // friend functions
  // ===========================================================================
  friend
  std::ostream& operator<< <>(std::ostream& out,
                      const QP_LU_factorization<ET, Is_LP, Matrix_Provider>& f);


  private:
    QP_LU_factorization (); // in order to prevent default construction

  public:
  
    // =========================================================================
    // typedefs
    // =========================================================================
    typedef QP_sparse_matrix<ET>                          matrix_t;
    typedef QP_sparse_vector<ET>                          vector_t;
    typedef typename QP_sparse_matrix<ET>::index_pair_t   index_pair_t;
    typedef std::vector<ET>                               dense_vector_t;
    typedef std::vector<int>                              permutation_t;
    
    
  
  
    // =========================================================================
    // creation and initialization
    // =========================================================================
    //QP_LU_factorization (CGAL::Verbose_ostream&  verbose_ostream,
    //                     CGAL::QP_sparse_matrix<ET>& to_factor);
    QP_LU_factorization (CGAL::Verbose_ostream&  verbose_ostream,
                         Matrix_Provider* matrix_provider);
            
            
    // =========================================================================
    // Getter & setter methods
    // =========================================================================
    //std::auto_ptr<QP_sparse_matrix<ET> > get_matrix_L() const;
    
    const ET& get_denominator() const {return d_;}
    const dense_vector_t& get_denominators() const {return d_row_;}
    

    // =========================================================================
    // Debug methods, to be used with care, preconditions may apply
    // =========================================================================
    
    // PRE: This assumes gmpzf as exact instantiation type
    int report_largest_GMPZF() const {
      typename QP_sparse_matrix<ET>::value_iterator it_row, it_row_end;
      int ret = 0;
      
      // check U
      for (int i = 0; i < n_; ++i) {
        it_row = matrix_u_->begin_row_value(i);
        it_row_end = matrix_u_->end_row_value(i);
        while (it_row != it_row_end) {
          int tmp;
          if ((tmp = mpz_sizeinbase((*it_row->second).man(), 2)) > ret) {
            ret = tmp;
          }
          ++it_row;
        }
      }
      // check L
      for (int i = n_; i < n_; ++i) {
        it_row = matrix_l_->begin_row_value(i);
        it_row_end = matrix_l_->end_row_value(i);
        while (it_row != it_row_end) {
          int tmp;
          if ((tmp = mpz_sizeinbase((*it_row->second).man(), 2)) > ret) {
            ret = tmp;
          }
          ++it_row;
        }
      }
      
      // check denominators
      for (int i = 0; i < d_row_.size(); ++i) {
        int tmp;
        if ((tmp = mpz_sizeinbase((d_row_[i]).man(), 2)) > ret) {
          ret = tmp;
        }
      }
      return CGAL::QP_solver_debug::timer.bit_size = ret;
    }


    int report_matrix_size() const {
      return CGAL::QP_solver_debug::timer.matrix_size = n_;
    }
    
    double report_density() {
      return CGAL::QP_solver_debug::timer.density = ( (matrix_l_ != 0 && matrix_u_ != 0) ? (matrix_l_->get_density()+matrix_u_->get_density()) / 2.0 : 0.0);
    }
    
    
    void test();
    
    
    // =========================================================================
    // Operations (public)
    // =========================================================================
    
    
    
    void set_invalid();
    
    std::ostream& print(std::ostream& out) const;
    
    // TAG: TODO Templatize the solve routines.
    // FTRAN ( LU y = v )
    template < typename ForwardIterator, typename OutputIterator >
    void  solve( ForwardIterator v_l_it, ForwardIterator v_x_it,
                     OutputIterator y_l_it,  OutputIterator y_x_it,
                     bool is_linear, bool is_phase_I);
    
    
    
    // special matrix-vector multiplication functions for LPs
    template < class ForwardIterator, class OutputIterator >  inline
    void  solve_l( ForwardIterator v_x_it, OutputIterator y_l_it) {
      CGAL_qpe_assertion( is_linear_ || is_phase_I_);
      solve__l( v_x_it, y_l_it);
    }
    
    template < class ForwardIterator, class OutputIterator >  inline
    void  solve_x( ForwardIterator v_l_it, OutputIterator y_x_it) {
      CGAL_qpe_assertion( is_linear_ || is_phase_I_);
	    solve__x( v_l_it, y_x_it); 
	  }
    
    // Rank 1 update
    // A' = A + a y z^t
    // where a is a scalar, y and z are vectors of appropriate dimensions
    bool rank_1_update(ET alpha, QP_sparse_vector<ET> y, QP_sparse_vector<ET> z);
    
    // Enlarge matrix by unit row/column
    bool enlarge(unsigned int k, unsigned int a,  bool is_basic);

    
    // Shrink matrix one row/column
    bool shrink(unsigned int k, bool is_basic);
	  
    // To retreive the internal permutation with regards to LU pivoting
	  void get_perm_row(permutation_t& perm) const;
    void get_perm_col(permutation_t& perm) const;
    void get_perm_row_inv(permutation_t& perm) const;
    void get_perm_col_inv(permutation_t& perm) const;
    
    // access to internal permutation
    inline int perm_row(int i) const {
      return perm_row_[i];
    }
    inline int perm_row_inv(int i) const {
      return perm_row_inv_[i];
    }
    inline int perm_col(int j) const {
      return perm_column_[j];
    }
    inline int perm_col_inv(int j) const {
      return perm_column_inv_[j];
    }
    
    // permutation bookkeeping
    void swap_rows(int i, int j);
    void swap_columns(int i, int j);
    
    void swap_columns_physically(int i, int j);
    void swap_rows_physically(int i, int j);
    
    boost::shared_ptr<QP_sparse_matrix<ET> > recover_original_matrix();
	  
	  
    
    // =========================================================================
    // update functions
    // =========================================================================


  private:
  
    // =========================================================================
    // typedefs
    // =========================================================================
    /*
    struct l_component_t {
      l_component_t(int a = 0, int b = 0, ET c = static_cast<ET>(1)): i(a), j(b), value(c) {}
      int i, j;
      ET  value;
    };
    typedef std::vector<l_component_t> matrix_L_t;
    */
    
    typedef typename vector_t::sparse_iterator_t Sparse_vector_iterator_t;
    typedef typename vector_t::sparse_const_iterator_t Sparse_const_vector_iterator_t;
    typedef typename matrix_t::value_iterator Sparse_matrix_value_iterator_t;
    typedef typename matrix_t::value_const_iterator Sparse_matrix_value_const_iterator_t;
    
    // Number type category (either Field_tag or Unique_factorization_domain_tag)
    typedef  Algebraic_structure_traits<ET> Traits;
    typedef typename Traits::Algebraic_category Category;

    
    
    // =========================================================================
    // Operations (private)
    // =========================================================================
    
    // factorization
    index_pair_t getPivotElement(const int k) const;
    index_pair_t getPivotElement_no_pivoting(const int k) const;
    index_pair_t getPivotElement_partial_pivoting(const int k) const;
    index_pair_t getPivotElement_markowitz(const int k) const;
    index_pair_t getPivotElement_markowitz_extended(const int k) const;
  
    void compute_factorization(bool is_linear_);
    void compute_factorization_gauss();
  
    //std::auto_ptr<QP_sparse_matrix<ET> > invert_L() const;
    
    
    
    // FTRAN ( LU y = v ) private
    // matrix-vector multiplication
    template < class ForIt, class OutIt>      // QP case, integral
    void  solve_QP( ForIt v_l_it, ForIt v_x_it,
                 OutIt y_l_it, OutIt y_x_it);
  
    template < class ForIt, class OutIt>      // LP case
    void  solve_LP( ForIt v_l_it, ForIt v_x_it,
		                OutIt y_l_it, OutIt y_x_it);


    // special matrix-vector multiplication functions for LPs
    // BTRAN (y^T LU = v^T)
    template < class ForIt, class OutIt >
    void  solve__l( ForIt v_x_it, OutIt y_l_it);
    // FTRAN (LU y = v)
    template < class ForIt, class OutIt >
    void  solve__x( ForIt v_l_it, OutIt y_x_it);
  
  
  // TAG: TODO revert to private
  public:
    // computing least common multiple
    /*
    template < class ET >
    ET gcd(ET a, ET b); */
    typename Algebraic_structure_traits< typename Coercion_traits<ET,ET>::Type >
    ::Gcd::result_type
    lcm(ET a, ET b);
    typename Algebraic_structure_traits< typename Coercion_traits<ET,ET>::Type >
    ::Gcd::result_type lcm(ET a, ET b, Field_tag );
    typename Algebraic_structure_traits< typename Coercion_traits<ET,ET>::Type >
    ::Gcd::result_type lcm(ET a, ET b, Unique_factorization_domain_tag );
    

    
    
    // =========================================================================
    // data members
    // =========================================================================
  private:
    // TAG: TODO remove valid_ member and replace with update functions
    bool                        valid_;     // factorization is valid  
  
    const ET                    et0_;
    const ET                    et1_;
    const ET                    et2_;
    
    int                         n_;         // size of the matrix
    
    boost::shared_ptr<matrix_t>   matrix_u_;
    boost::shared_ptr<matrix_t>   matrix_l_;
    
    // TODO: use typedef for d_row_
    dense_vector_t              d_row_;     // row denominators
    ET                          d_;         // denominator (always positive)
    bool                        d_neg_;     // indicates that d_ should be negative

    unsigned int                l_;         // minimum of `n' and `m'
    int                         csize_;     // size of `C = E \cup S_N' // TAG: REMOVE "formerly s"
    int                         bosize_;    // size of `B_O' // TAG: REMOVE "formerly b"

    
    bool                        is_phase_I_; // flag indicating phase I
    bool                        is_phase_II_;// flag indicating phase II
    bool                        is_linear_;  // flag indicating a linear program

    CGAL::Verbose_ostream&      vout_;      // used for diagnostic output
    
    permutation_t               perm_row_, perm_column_;
    permutation_t               perm_row_inv_, perm_column_inv_;
    
    Matrix_Provider*            matrix_provider_; // pointer to a class
                                                  // that implements
                                                  // and get_basis_matrix()
         
         
    // =========================================================================
    // old stuff (to be abandoned) // TAG: TODO
    // =========================================================================
  public:
    // vector-vector multiplication ( y = u^T v )
    template < class InputIterator1, class InputIterator2 >  inline
    typename std::iterator_traits<InputIterator1>::value_type
    inner_product( InputIterator1 u_l_it, InputIterator1 u_x_it,
		   InputIterator2 v_l_it, InputIterator2 v_x_it)
        { 
          /*if (!valid_) compute_factorization(is_linear_ || is_phase_I_);*/
        return inner_product_l( u_l_it, v_l_it)
	       + inner_product_x( u_x_it, v_x_it); }
    
    template < class InputIterator1, class InputIterator2 >  inline
    typename std::iterator_traits<InputIterator1>::value_type
    inner_product_l( InputIterator1 u_l_it, InputIterator2 v_l_it)
        { /*if (!valid_) compute_factorization(is_linear_ || is_phase_I_);*/
        return inner_product( u_l_it, v_l_it, csize_); }
    
    template < class InputIterator1, class InputIterator2 >  inline
    typename std::iterator_traits<InputIterator1>::value_type
    inner_product_x( InputIterator1 u_x_it, InputIterator2 v_x_it)
        { /*if (!valid_) compute_factorization(is_linear_ || is_phase_I_);*/
        return inner_product( u_x_it, v_x_it, bosize_); }
    
    // vector-vector multiplication  
    template < class InIt1, class InIt2 > inline
    typename std::iterator_traits<InIt1>::value_type  
    inner_product( InIt1 u_it, InIt2 v_it, unsigned int n) const
    {
	  typedef  typename std::iterator_traits<InIt1>::value_type  NT;
    
        // compute u^T v
        NT sum = NT(0);
        for ( unsigned int count = 0; count < n; ++count, ++u_it, ++v_it) {
          sum += NT(*u_it) * NT(*v_it);
        }
        return sum;
    }                                                  
};


// =============================================================================
// Non-member functions
// =============================================================================

template <typename ET, typename Is_LP, typename Matrix_Provider>
std::ostream& operator<<( std::ostream& out,
                  const QP_LU_factorization<ET, Is_LP, Matrix_Provider>& f) {
  return f.print(out);
}



} //namespace CGAL

#include <CGAL/QP_solver/QP_LU_factorization_impl.h>

#endif //CGAL_QP_LU_FACTORIZATION_H

// ===== EOF ==================================================================
