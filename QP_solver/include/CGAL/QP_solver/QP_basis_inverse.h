// Copyright (c) 1997-2007  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
//
// Author(s)     : Sven Schoenherr
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp
//                 Kaspar Fischer 
                                                                               
#ifndef CGAL_QP_SOLVER_QP_BASIS_INVERSE_H
#define CGAL_QP_SOLVER_QP_BASIS_INVERSE_H

#include <CGAL/QP_solver/basic.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <vector>

namespace CGAL {
                    
// =================
// class declaration
// =================
template < class ET_, class Is_LP_ >
class QP_basis_inverse;

// ===============
// class interface
// ===============
template < class ET_, class Is_LP_ >
class QP_basis_inverse {
public:
  // self
  typedef  ET_                        ET;
  typedef  Is_LP_                     Is_LP;
  typedef  QP_basis_inverse<ET,Is_LP>
  Self;

private:
    
  // private types
  typedef std::vector<ET>             Row;
  typedef std::vector<Row>            Matrix;

public:

  /*
   * Note: Some member functions below are suffixed with '_'.
   * They are member templates and their declaration is "hidden",
   * because they are also implemented in the class interface.
   * This is a workaround for M$-VC++, which otherwise fails to
   * instantiate them correctly.
   * Look for "TAG: M$-VC++ Workaround"
   */

  // creation and initialization
  // ---------------------------
  QP_basis_inverse( CGAL::Verbose_ostream&  verbose_ostream);

  // set-up
  void  set( int n, int m, int nr_equalities);
    
  // init
    
  // TAG: M$-VC++ Workaround
  template < class InputIterator >                            // phase I
  void  init_( unsigned int art_size, InputIterator art_first);
    

  /*
    template < class InputIterator >                            // phase II
    void  init_( ...);
  */

  // transition to phase II
  template < class InputIterator >                            // QP case
  void  transition_( InputIterator twice_D_it);

  void  transition( );                                        // LP case

  // access
  // ------
  const ET&  denominator( ) const { return denominator_; }

  const ET&  entry( unsigned int row, unsigned int column) const
  { return entry( row, column, Is_LP()); }
  
  
  // TAG: DEBUG
  
  int report_largest_GMPZF() const {
    int tmp;
    int ret = 0;
    for (int i = 0; i < M.size(); ++i) {
      for (int j = 0; j < M[i].size(); ++j) {      
        if ((tmp = mpz_sizeinbase((M[i][j]).man(), 2)) > ret) {
          ret = tmp;
        }
      }
    }
    if  ((tmp = mpz_sizeinbase(denominator_.man(), 2)) > ret) {
      ret = tmp;
    }
    return CGAL::QP_solver_debug::timer.bit_size = ret;
  }
  
  // TAG: DEBUG
  int report_matrix_size() const {
      return CGAL::QP_solver_debug::timer.matrix_size = size_E_cup_SN_+size_BO_;
  }
  
  
  

  // multiplication functions
  // ------------------------
  // matrix-vector multiplication ( y = M v )
  template < class ForwardIterator, class OutputIterator >  inline
  void  solve( ForwardIterator v_l_it, ForwardIterator v_x_it,
               OutputIterator y_l_it,  OutputIterator y_x_it) const
  { solve( v_l_it, v_x_it, y_l_it, y_x_it, Is_LP(), Tag_true()); }
    
  // special matrix-vector multiplication functions for LPs
  template < class ForwardIterator, class OutputIterator >  inline
  void  solve_l( ForwardIterator v_x_it, OutputIterator y_l_it) const
  { CGAL_qpe_assertion( is_LP || is_phaseI);
    solve__l( v_x_it, y_l_it); }
    
  template < class ForwardIterator, class OutputIterator >  inline
  void  solve_x( ForwardIterator v_l_it, OutputIterator y_x_it) const
  { CGAL_qpe_assertion( is_LP || is_phaseI);
    solve__x( v_l_it, y_x_it); }
    
    
  /*
  // TAG: YVES delete
  // vector-matrix multiplication ( x^T = u^T M )
  template < class ForwardIterator, class OutputIterator >  inline
  void  multiply_transposed( ForwardIterator u_l_it, ForwardIterator u_x_it,
  OutputIterator x_l_it,  OutputIterator x_x_it)
  const
  { multiply( u_l_it, u_x_it, x_l_it, x_x_it); } // M_B^{-1} is symmetric
    
  // special vector-matrix multiplication functions for LPs
  template < class ForwardIterator, class OutputIterator >  inline
  void  multiply_transposed_l( ForwardIterator u_x_it,
  OutputIterator x_l_it) const
  { multiply_l( u_x_it, x_l_it); }        // Note: M_B^{-1} is symmetric
    
  template < class ForwardIterator, class OutputIterator >  inline
  void  multiply_transposed_x( ForwardIterator u_l_it,
  OutputIterator x_x_it) const
  { multiply_x( u_l_it, x_x_it); }        // Note: M_B^{-1} is symmetric
  */
    
  // vector-vector multiplication ( y = u^T v )
  template < class InputIterator1, class InputIterator2 >  inline
  typename std::iterator_traits<InputIterator1>::value_type
  inner_product( InputIterator1 u_l_it, InputIterator1 u_x_it,
                 InputIterator2 v_l_it, InputIterator2 v_x_it) const
  { return inner_product_l( u_l_it, v_l_it)
      + inner_product_x( u_x_it, v_x_it); }
    
  template < class InputIterator1, class InputIterator2 >  inline
  typename std::iterator_traits<InputIterator1>::value_type
  inner_product_l( InputIterator1 u_l_it, InputIterator2 v_l_it) const
  { return inner_product( u_l_it, v_l_it, size_E_cup_SN_); }
    
  template < class InputIterator1, class InputIterator2 >  inline
  typename std::iterator_traits<InputIterator1>::value_type
  inner_product_x( InputIterator1 u_x_it, InputIterator2 v_x_it) const
  { return inner_product( u_x_it, v_x_it, size_BO_); }
        
    
  // update functions
  // ----------------
  // entering of original variable (update type U1)
  template < class ForwardIterator >
  void  enter_original_( ForwardIterator y_l_it,
                         ForwardIterator y_x_it, const ET& z);
    
  // leaving of original variable (update type U2)
  void  leave_original( );
    
  // entering of slack variable (update type U3)
  void  enter_slack( );

  // leaving of slack variable (update type U4)
  template < class ForwardIterator >
  void  leave_slack_( ForwardIterator u_x_it);
    
  // replacing of original by original variable (update type U5)
  template < class ForwardIterator >
  void  enter_original_leave_original_( ForwardIterator y, unsigned int k);
    
  // replacing of slack by slack variable (update type U6)
  template < class ForwardIterator >
  void  enter_slack_leave_slack_( ForwardIterator u, unsigned int k);
    
  // replacing of slack by original variable (update type U7)
  template < class ForwardIterator1, class ForwardIterator2 >
  void  enter_original_leave_slack_( ForwardIterator1 y, ForwardIterator2 u);
    
  // replacing of original by slack variable (update type U8)
  void  enter_slack_leave_original( );
    
    
  // replacing of original by original variable with precondition in QP-case
  // for phaseII                               (update type UZ_1)
  template < class ForwardIterator >
  void  z_replace_original_by_original(ForwardIterator y_l_it,
                                       ForwardIterator y_x_it,
                                       const ET& s_delta, const ET& s_nu,
                                       unsigned int k_i);


  // replacing of original by slack variable with precondition in QP-case
  // for phaseII                               (update type UZ_2)
  void  z_replace_original_by_slack( );


  // replacing of slack by original variable with precondition in QP-case
  // for phaseII                               (update type UZ_3)
  template < class ForwardIterator >
  void  z_replace_slack_by_original(ForwardIterator y_l_it,
                                    ForwardIterator y_x_it,
                                    ForwardIterator u_x_it, const ET& hat_kappa,
                                    const ET& hat_nu);


  // replacing of slack by slack variable with precondition in QP-case
  // for phaseII                               (update type UZ_4)
  template < class ForwardIterator >
  void  z_replace_slack_by_slack(ForwardIterator u_x_it, unsigned int k_j);
    

  // copying of reduced basis inverse row in (upper) C-half
  template < class OutIt >
  void  copy_row_in_C(OutIt y_l_it, OutIt y_x_it, unsigned int k);
    
  // copying of reduced basis inverse row in (lower) B_O-half
  template < class OutIt >
  void  copy_row_in_B_O(OutIt y_l_it, OutIt y_x_it, unsigned int k);


  // swap functions
  void  swap_variable( unsigned int j)                // ``to the end'' of R
  { CGAL_qpe_assertion( j < size_BO_);
    swap_variable( j, Is_LP()); }
  void  swap_constraint( unsigned int i)              // ``to the end'' of P
  { CGAL_qpe_assertion( i < size_E_cup_SN_);
    swap_constraint( i, Is_LP()); }

private:
    
  // constants
  const ET                 et0, et1, et2;
                                        
  // data members
  Matrix                   M;         // basis inverse, stored row-wise
  ET                       denominator_;         // denominator

  unsigned int             min_N_M_;         // minimum of `n' and `m'
  unsigned int             size_E_cup_SN_;         // size of `E \cup S_N'
  unsigned int             size_BO_;         // size of `B_O'

  bool                     is_phaseI; // flag indicating phase I
  bool                     is_phaseII;// flag indicating phase II
  const bool               is_LP;     // flag indicating a linear    program
  const bool               is_QP;     // flag indicating a quadratic program

  CGAL::Verbose_ostream&   vout;      // used for diagnostic output

  Row                      x_l, tmp_l;  // used in the 
  Row                      x_x, tmp_x;  // update functions

  // private member functions
  // ------------------------
  // set-up
  void  set( Tag_false);        // QP case
  void  set( Tag_true );        // LP case

    
    
  // init
  // TAG: M$-VC++ Workaround
  template < class InIt >                                     // QP case
  void  init_( unsigned int art_size, InIt art_first, Tag_false);
  template < class InIt >                                     // LP case
  void  init_( unsigned int art_size, InIt art_first, Tag_true );
    

  // access
  const ET&  entry( unsigned int row,
                    unsigned int column, Tag_false) const;    // QP case
  const ET&  entry( unsigned int row,
                    unsigned int column, Tag_true ) const;    // LP case

  // matrix-vector multiplication
  template < class ForIt, class OutIt, class Use1stArg >      // QP case
  void  solve_( ForIt v_l_it, ForIt v_x_it,
                OutIt y_l_it, OutIt y_x_it, Tag_false, Use1stArg) const;
  template < class ForIt, class OutIt, class  DummyArg >      // LP case
  void  solve_( ForIt v_l_it, ForIt v_x_it,
                OutIt y_l_it, OutIt y_x_it, Tag_true,   DummyArg) const;

    
  // TAG: M$-VC++ Workaround
  // special matrix-vector multiplication functions for LPs
  template < class ForIt, class OutIt >
  void  solve__l_( ForIt v_x_it, OutIt y_l_it) const;
  template < class ForIt, class OutIt >
  void  solve__x_( ForIt v_l_it, OutIt y_x_it) const;
    

  // in-place update
  template < class ForIt >                                    // QP case
  void  update_inplace_QP_( ForIt  y_l_it, ForIt  y_x_it,
                            const ET&  d_new, const ET&  d_old);
  template < class ForIt1, class ForIt2 >                     // LP case
  void  update_inplace_LP_( ForIt1 x_x_it, ForIt2 y_x_it,
                            const ET&  d_new, const ET&  d_old);
                              
  template < class ForIt >                                  // QP case only
  void  z_update_inplace( ForIt psi1_l_it, ForIt psi1_x_it,
                          ForIt psi2_l_it, ForIt psi2_x_it,
                          const ET& omega0, const ET& omega1,
                          const ET& omega2, const ET& omega3); 

  void  update_entry( ET& entry,   const ET& d_new,
                      const ET& y, const ET& d_old) const;

  // swap functions
  void  swap_variable  ( unsigned int, Tag_true );            // LP case
  void  swap_variable  ( unsigned int, Tag_false);            // QP case
  void  swap_constraint( unsigned int, Tag_true );            // LP case
  void  swap_constraint( unsigned int, Tag_false);            // QP case      
    
  // diagnostic output
  void  print( );

  // ----------------------------------------------------------------------------

  // ===============================
  // class implementation (template)
  // ===============================

public:

  // creation and initialization
  // ---------------------------
  // init
  template < class InputIterator >
  void
  init( unsigned int art_size, InputIterator art_first)
  {
    CGAL_qpe_assertion_msg( art_size <= min_N_M_, \
                            "There are more equality constraints than original variables!");

    init( art_size, art_first, Is_LP());
    denominator_ = et1;
    CGAL_qpe_assertion( size_E_cup_SN_ == art_size);
    CGAL_qpe_assertion( size_BO_ == art_size);

    is_phaseI  = true;
    is_phaseII = false;

    if ( x_x.size() < art_size) {
      x_x.insert( x_x.end(), art_size-x_x.size(), et0);
    }
        
    if ( tmp_x.size() < art_size) {
      tmp_x.insert( tmp_x.end(), art_size-tmp_x.size(), et0);
    }
        
    CGAL_qpe_debug {
      if ( vout.verbose()) print();
    }
  }

  // transition to phase II
  template < class InputIterator >                            // QP case
  void  transition( InputIterator twice_D_it)
  {
    typename Matrix::iterator  m_it1, m_it2, p_begin, r_begin;
    typename Row   ::iterator  x_it;
    unsigned int      row, col;

    // fill missing rows
    for (row= 0; row< size_E_cup_SN_; ++ row) {
      CGAL_qpe_assertion(M[row].size()==0);
      M[row].insert(M[row].end(), row+1, et0);
    }

    // compute new basis inverse [ upper-left part: -(Q^T * 2 D_B * Q) ]
    // -----------------------------------------------------------------
    // compute 'Q^T * 2 D_B' ( Q = A_B^-1 ) 
    p_begin = M.begin();
    r_begin = M.begin();
    if (size_BO_ > 0) r_begin += min_N_M_+size_E_cup_SN_-1;        // initialize only if for loops 
    // are entered
    for ( col = 0; col < size_BO_; ++col, ++twice_D_it) {
      ++p_begin;

      // get column of D (D symmetric)
      std::copy( *twice_D_it, *twice_D_it+size_BO_, x_l.begin());

      // compute 'Q^T * 2 D_Bj'
      solve__l( x_l.begin(), x_x.begin());

      // store resulting column in 'P' and 'R'
      x_it  = x_x.begin();
      m_it2 = r_begin;
      for ( row = min_N_M_+col; row >= min_N_M_; --row, --m_it2, ++x_it) {
        (*m_it2)[ row] = *x_it;
      }
      m_it1 = p_begin;
      for ( row = col+1; row <  size_E_cup_SN_; ++row, ++m_it1, ++x_it) {
        (*m_it1)[ col] = *x_it;
      }
    }

    // compute '(Q^T * 2 D_B) * Q' ( Q = A_B^-1 )
    m_it1 = M.begin();
    m_it2 = r_begin;
    for ( row = 0; row < size_E_cup_SN_; ++row, ++m_it1, --m_it2) {

      // get row of '(Q^T * 2 D_B)' (stored in 'P' and 'R')
      std::copy(m_it1->begin()  ,m_it1->begin()+row,    x_l.begin());
      std::copy(m_it2->begin()+min_N_M_,m_it2->begin()+(min_N_M_+size_BO_-row),
                x_l.begin()+row);

      // compute '(Q^T * 2 D_B)_i * Q'
      solve__l( x_l.begin(), x_x.begin());

      // negate and store result in 'P'
      std::transform( x_x.begin(), x_x.begin()+row+1,
                      m_it1->begin(), std::negate<ET>());

      // clean up in 'R'
      std::fill_n( m_it2->begin()+min_N_M_, size_BO_-row, et0);
    }

    // scale A_B^-1
    m_it1 = M.begin()+min_N_M_;
    for ( row = 0; row < size_E_cup_SN_; ++row, ++m_it1) {
      std::transform( m_it1->begin(), m_it1->begin()+size_E_cup_SN_, m_it1->begin(),
                      std::bind2nd( std::multiplies<ET>(), denominator_));
    }

    // new denominator: |det(A_B)|^2
    denominator_ *= denominator_;

    // update status
    transition();
  }

  // update functions
  // ----------------
  // entering of original variable (update type U1)
  template < class ForwardIterator >
  void
  enter_original( ForwardIterator y_l_it,
                  ForwardIterator y_x_it, const ET& z)
  {
    
    // TAG: DEBUG
    ++CGAL::QP_solver_debug::timer.counter_Oin;
  
    // assert QP case
    Assert_compile_time_tag( Tag_false(), Is_LP());
    
    // update matrix in-place
    // ----------------------
    // handle sign of new denominator
    CGAL_qpe_assertion( z != et0);
    bool  z_neg = ( z < et0);

    // update matrix
    update_inplace_QP( y_l_it, y_x_it, z, ( z_neg ? -denominator_ : denominator_));
                                                                      
    // append new row and column ("after R")
    // -------------------------------------
    typename Row::iterator  row_it;
    ForwardIterator           y_it;
    unsigned int            col, k = min_N_M_+(++size_BO_);
    
    //      // allocate new row, if necessary
    //      // BG: replaced this by the ensure_physical_row construct below
    //      if ( k <= M.size()) {
    //           row_it = M[ k-1].begin();
    //      } else {
    //           row_it = M.insert( M.end(), Row( k, et0))->begin();
    //           x_x.push_back( et0);
    //           tmp_x.push_back( et0);
    //      }
    ensure_physical_row(k-1);
    row_it = M[ k-1].begin();
    
    // store entries in new row
    for (   col = 0,                              y_it = y_l_it;
            col < size_E_cup_SN_;
            ++col,     ++row_it,                  ++y_it         ) {
      *row_it = ( z_neg ? *y_it : -( *y_it));
    }
    for (   col = min_N_M_,     row_it += min_N_M_-size_E_cup_SN_,           y_it = y_x_it;
            col < k-1;
            ++col,       ++row_it,                ++y_it         ) {
      *row_it = ( z_neg ? *y_it : -( *y_it));
    }
    *row_it = ( z_neg ? -denominator_ : denominator_);
    
    // store new denominator
    // ---------------------
    denominator_ = ( z_neg ? -z : z);
    CGAL_qpe_assertion( denominator_ > et0);
    
    CGAL_qpe_debug {
      if ( vout.verbose()) print();
    }
  }
    
  // leaving of slack variable (update type U4)
  template < class ForwardIterator >
  void
  leave_slack( ForwardIterator u_x_it)
  {
    
    // TAG: DEBUG
    ++CGAL::QP_solver_debug::timer.counter_Sout;
  
    // assert QP case
    Assert_compile_time_tag( Tag_false(), Is_LP());
    
    // update matrix in-place
    // ----------------------
    // compute new row/column of basis inverse
    solve( u_x_it,                               // dummy (not used)
           u_x_it, x_l.begin(), x_x.begin(),
           Tag_false(),                          // QP
           Tag_false());                         // ignore 1st argument
    ET    z     = -inner_product_x( x_x.begin(), u_x_it);
    bool  z_neg = ( z < et0);
    CGAL_qpe_assertion( z != et0);
    
    // update matrix
    update_inplace_QP( x_l.begin(), x_x.begin(), z, ( z_neg ? -denominator_ : denominator_));
                                                                      
    // insert new row and column ("after P")
    // -------------------------------------
    typename Row   ::iterator  row_it, x_it;
    typename Matrix::iterator  col_it;
    unsigned int               col, k = min_N_M_+size_BO_;
    
    // store entries in new row
    if (M[size_E_cup_SN_].size()==0) {
      // row has to be filled first
      M[size_E_cup_SN_].insert(M[size_E_cup_SN_].end(), size_E_cup_SN_+1, et0);
    }
    for (   col = 0,   row_it = M[ size_E_cup_SN_].begin(),        x_it = x_l.begin();
            col < size_E_cup_SN_;
            ++col,     ++row_it,                      ++x_it              ) {
      *row_it = ( z_neg ? *x_it : -( *x_it));
    }
    *row_it = ( z_neg ? -denominator_ : denominator_);
    
    for (   col = min_N_M_,   col_it = M.begin()+min_N_M_,          x_it = x_x.begin();
            col < k;
            ++col,     ++col_it,                      ++x_it              ) {
      (*col_it)[ size_E_cup_SN_] = ( z_neg ? *x_it : -( *x_it));
    }
    ++size_E_cup_SN_;
    
    // store new denominator
    // ---------------------
    denominator_ = ( z_neg ? -z : z);
    CGAL_qpe_assertion( denominator_ > et0);
    
    CGAL_qpe_debug {
      if ( vout.verbose()) print();
    }
  }

  // replacing of original variable (update type U5) [ replace column ]
  template < class RandomAccessIterator >
  void
  enter_original_leave_original( RandomAccessIterator y_x_it, unsigned int k)
  {
  
    
    // TAG: DEBUG
    ++CGAL::QP_solver_debug::timer.counter_O_O;
  
    // assert LP case or phase I
    CGAL_qpe_assertion( is_LP || is_phaseI);
    CGAL_qpe_assertion( k < size_BO_);

    // update matrix in place
    // ----------------------
    typename Matrix::iterator  matrix_it;
    typename Row   ::iterator     row_it, row_k_it, row_k;

    // handle sign of new denominator
    ET    z     = y_x_it[ k];
    bool  z_neg = ( z < et0);
    if ( z_neg) denominator_ = -denominator_;

    // QP (in phase I)?
    matrix_it = M.begin();
    if ( is_QP) matrix_it += min_N_M_;
    row_k = ( matrix_it+k)->begin();

    // rows: 0..s-1 without k
    unsigned int  row, col;
    ET            minus_y;
    for (   row = 0;
            row < size_E_cup_SN_;
            ++row,     ++matrix_it, ++y_x_it) {
      if ( row != k) {

        // columns: 0..b-1
        minus_y = -( *y_x_it);
        for (   col = 0, row_it = matrix_it->begin(), row_k_it = row_k;
                col < size_BO_;
                ++col,   ++row_it,                    ++row_k_it       ){
        
          // update in place
          update_entry( *row_it, z, minus_y * *row_k_it, denominator_);
        }
      }
    }

    // rows: k (flip signs, if necessary)
    if ( z_neg) {
      for (   col = 0,   row_it = row_k;
              col < size_BO_;
              ++col,     ++row_it        ) {
        
        *row_it = -( *row_it);
      }
    }

    // store new denominator
    // ---------------------
    denominator_ = ( z_neg ? -z : z);
    CGAL_qpe_assertion( denominator_ > et0);

    // diagnostic output
    CGAL_qpe_debug {
      if ( vout.verbose()) print();
    }
  }
    
  // replacing of slack variable (update type U6) [ replace row ]
  template < class ForwardIterator >
  void
  enter_slack_leave_slack( ForwardIterator u_x_it, unsigned int k)
  {
  
    
    // TAG: DEBUG
    ++CGAL::QP_solver_debug::timer.counter_S_S;
    
    // assert LP case or phase I
    CGAL_qpe_assertion( is_LP || is_phaseI);
    CGAL_qpe_assertion( k < size_E_cup_SN_);

    // compute new row of basis inverse
    solve__l( u_x_it, x_x.begin());

    // update matrix in place
    // ----------------------
    typename Matrix::iterator  matrix_it;
    typename Row   ::iterator     row_it, x_it;

    // handle sign of new denominator
    ET    z     = x_x[ k];
    bool  z_neg = ( z < et0);
    if ( z_neg) denominator_ = -denominator_;

    // QP (in phase I)?
    matrix_it = M.begin();
    if ( is_QP) matrix_it += min_N_M_;

    // rows: 0..s-1
    unsigned int          row, col;
    ET            minus_m_row;
    for (   row = 0;
            row < size_E_cup_SN_;
            ++row,     ++matrix_it) {

      // columns: 0..b-1
      minus_m_row = -( *matrix_it)[ k];
      for (   col = 0,   row_it = matrix_it->begin(), x_it = x_x.begin();
              col < size_BO_;
              ++col,     ++row_it,                    ++x_it             ){

        if ( col != k) {                // all columns but k

          // update in place
          update_entry( *row_it, z, minus_m_row * *x_it, denominator_);

        } else {                        // column k

          // flip sign, if necessary
          if ( z_neg) *row_it = -( *row_it);

        }
      }
    }

    // store new denominator
    // ---------------------
    denominator_ = ( z_neg ? -z : z);
    CGAL_qpe_assertion( denominator_ > et0);

    // diagnostic output
    CGAL_qpe_debug {
      if ( vout.verbose()) print();
    }
  }
    
  // replacing of slack by original variable (update type U7) [ grow ]
  template < class ForwardIterator1, class ForwardIterator2 >
  void  enter_original_leave_slack( ForwardIterator1 y_x_it,
                                    ForwardIterator2 u_x_it)
  {
    
    // TAG: DEBUG
    ++CGAL::QP_solver_debug::timer.counter_O_S;
  
    // assert LP case or phase I
    CGAL_qpe_assertion( is_LP || is_phaseI);

    // update matrix in-place
    // ----------------------
    // compute new row of basis inverse
    solve__l( u_x_it, x_x.begin());
    ET    z     = denominator_*u_x_it[ size_BO_] - inner_product_x( y_x_it, u_x_it);
    bool  z_neg = ( z < et0);
    CGAL_qpe_assertion( z != et0);
        
    // update matrix
    update_inplace_LP( x_x.begin(), y_x_it, z, ( z_neg ? -denominator_ : denominator_));
                                                  
    // append new row and column
    // -------------------------
    typename Matrix::iterator  matrix_it;
    typename Row   ::iterator     row_it, x_it;
    unsigned int                  row, col;

    // QP (in phase I)?
    if ( is_QP) {
      ensure_physical_row(min_N_M_+size_BO_);
      row_it = M[min_N_M_+size_BO_].begin();
      matrix_it = M.begin() + min_N_M_;
    } else {
      row_it = M[size_E_cup_SN_].begin();
      matrix_it = M.begin();
    }       
        
    // store 'x_x' in new row 
    x_it = x_x.begin();
    for ( col = 0; col < size_BO_; ++col, ++row_it, ++x_it) {
      *row_it = ( z_neg ? *x_it : -( *x_it));
    }
    *row_it = ( z_neg ? -denominator_ : denominator_);
        
    // store 'y_x' in new col
    for ( row = 0; row < size_E_cup_SN_; ++row, ++matrix_it, ++y_x_it) {
      (*matrix_it)[size_BO_] = ( z_neg ? *y_x_it : -( *y_x_it));
    }
    ++size_E_cup_SN_; ++size_BO_;
    
    // store new denominator
    // ---------------------
    denominator_ = ( z_neg ? -z : z);
    CGAL_qpe_assertion( denominator_ > et0);
    
    CGAL_qpe_debug {
      if ( vout.verbose()) print();
    }
  }
  // due to basis_matrix_stays_regular fix, needs reconsideration
  //private:

  // private member functions
  // ------------------------
  // init (QP case)
  template < class InIt >  inline
  void
  init( unsigned int art_size, InIt art_first, Tag_false)
  {
    // only Q is used to store A_B^-1 in phase I
    for ( size_E_cup_SN_ = min_N_M_, size_BO_ = 0; size_BO_ < art_size; ++size_E_cup_SN_, ++size_BO_, ++art_first) {
      ensure_physical_row(size_E_cup_SN_);
      M[ size_E_cup_SN_][ size_BO_] = ( art_first->second ? -et1 : et1);
    }

    size_E_cup_SN_ = art_size;
  }

  // init (LP case)
  template < class InIt >  inline
  void
  init( unsigned int art_size, InIt art_first, Tag_true)
  {
    for ( size_E_cup_SN_ = 0; size_E_cup_SN_ < art_size; ++size_E_cup_SN_, ++art_first) {
      std::fill_n( M[ size_E_cup_SN_].begin(), art_size, et0);
      M[ size_E_cup_SN_][ size_E_cup_SN_] = ( art_first->second ? -et1 : et1);
    }
    size_BO_ = art_size;
  }
    
  // append row in Q if no allocated row available
  void ensure_physical_row (unsigned int row) {
    unsigned int rows = static_cast<unsigned int>(M.size());
    CGAL_qpe_assertion(rows >= row);
    if (rows == row) {
      M.push_back(Row(row+1, et0));
            
      // do we have to grow x_x?
      CGAL_qpe_assertion(x_x.size() >= row-min_N_M_);
      if (x_x.size() == row-min_N_M_)
        x_x.push_back(et0);
            
      // do we have to grow tmp_x?
      CGAL_qpe_assertion(tmp_x.size() >= row-min_N_M_);
      if (tmp_x.size() == row-min_N_M_)
        tmp_x.push_back(et0);
            
      CGAL_qpe_assertion(M[row].size()==row+1);
      CGAL_qpe_debug {
        if ( vout.verbose()) {
          vout << "physical row " << (row) << " appended in Q\n";
        }
      }
    }
  }
    
  // matrix-vector multiplication (QP case)
  template < class ForIt, class OutIt, class Use1stArg >
  void
  solve( ForIt v_l_it, ForIt v_x_it,
         OutIt y_l_it, OutIt y_x_it, Tag_false,
         Use1stArg) const
  {
    // use 'LP' functions in phase I
    if ( is_phaseI) {
      solve__l( v_x_it, y_l_it);
      solve__x( v_l_it, y_x_it);
      return;
    }

    // phase II
    typename Matrix::const_iterator  matrix_it;
    typename Row   ::const_iterator     row_it;     // left  of diagonal
    typename Matrix::const_iterator  column_it;     // right of diagonal
    ForIt                                 v_it;
    
    unsigned int  row, count, k = min_N_M_+size_BO_;
    ET            sum;
    
    // compute  P v_l + Q^T v_x     
    for (   row = 0,   matrix_it = M.begin();
            row < size_E_cup_SN_;
            ++row,                                                ++y_l_it) {
      sum = et0;

      // P v_l
      if ( Use1stArg::value ) {

        // P: left of diagonal (including)
        for (   row_it =  matrix_it->begin(),            v_it = v_l_it;
                row_it != matrix_it->end();
                ++row_it,                                ++v_it) {
          sum += *row_it * *v_it;
        }

        // P: right of diagonal (excluding)
        for (   count = row+1,   column_it  = ++matrix_it;
                count < size_E_cup_SN_;
                ++count,         ++column_it,                ++v_it) {
          sum += (*column_it)[ row] * *v_it;
        }
      }
    
      // Q^T:
      for (   count = 0,       column_it  = M.begin()+min_N_M_,   v_it = v_x_it;
              count < size_BO_;
              ++count,         ++column_it,                ++v_it) {
        sum += (*column_it)[ row] * *v_it;
      }
    
      // store result
      *y_l_it = sum;
    }
    
    // compute  Q v_l + R v_x
    for (   row = min_N_M_,   matrix_it = M.begin()+min_N_M_;
            row < k;
            ++row,                                                ++y_x_it) {
              
      sum = et0;

      // Q v_l
      if ( Use1stArg::value ) {

        // Q:
        for (   count = 0,  row_it = matrix_it->begin(), v_it = v_l_it;
                count < size_E_cup_SN_;
                ++count,    ++row_it,                    ++v_it) {
          sum += *row_it * *v_it;
        }
      }
    
      // R: left of diagonal (including)
      for ( row_it =  matrix_it->begin()+min_N_M_, v_it = v_x_it;
            row_it != matrix_it->end();
            ++row_it,                       ++v_it) {
        sum += *row_it * *v_it;
      }
    
      // R: right of diagonal (excluding)
      for (   count = row+1,   column_it = ++matrix_it;
              count < k;
              ++count,         ++column_it,                ++v_it) {
        sum += (*column_it)[ row] * *v_it;
      }
    
      // store result
      *y_x_it = sum;
    }
  }
    
  // matrix-vector multiplication (LP case)
  template < class ForIt, class OutIt, class Dummy > inline
  void
  solve( ForIt v_l_it, ForIt v_x_it,
         OutIt y_l_it, OutIt y_x_it, Tag_true, Dummy) const
  {
    solve__l( v_x_it, y_l_it);
    solve__x( v_l_it, y_x_it);
  }
    
  // special matrix-vector multiplication functions for LPs
  template < class ForIt, class OutIt > inline
  void
  solve__l( ForIt v_x_it, OutIt y_l_it) const
  {
    typename Matrix::const_iterator  matrix_it = M.begin();
    typename Matrix::const_iterator  column_it;
    ForIt                                 v_it;
    
    unsigned int  row, count;
    ET            sum;
    
    // QP?
    if ( is_QP) matrix_it += min_N_M_;

    // compute  Q^T v_x
    for ( row = 0; row < size_E_cup_SN_; ++row,                              ++y_l_it) {
      sum = et0;
    
      for (   count = 0,   column_it = matrix_it,   v_it = v_x_it;
              count < size_BO_;
              ++count,     ++column_it,             ++v_it) {
        sum += (*column_it)[ row] * *v_it;
      }
    
      *y_l_it = sum;
    }
  }
    
  template < class ForIt, class OutIt > inline
  void
  solve__x( ForIt v_l_it, OutIt y_x_it) const
  {
    typename Matrix::const_iterator  matrix_it = M.begin();
    unsigned int  row;

    // QP?
    if ( is_QP) matrix_it += min_N_M_;

    // compute  Q v_l
    for (   row = 0;
            row < size_BO_;
            ++row,     ++matrix_it, ++y_x_it) {

      *y_x_it = inner_product( matrix_it->begin(), v_l_it, size_E_cup_SN_);
    }
  }
    
  // vector-vector multiplication  
  template < class InIt1, class InIt2 > inline
  typename std::iterator_traits<InIt1>::value_type  
  inner_product( InIt1 u_it, InIt2 v_it, unsigned int n) const
  {
    typedef  typename std::iterator_traits<InIt1>::value_type  NT;
    
    // compute u^T v
    NT sum = NT( 0);
        
    // TAG: DEBUG
    /*
    std::cout << "INV inner_product_l" << std::endl;
    std::cout << "n: " << n << std::endl;
    */
            
    for ( unsigned int count = 0; count < n; ++count, ++u_it, ++v_it) {
      sum += NT(*u_it) * NT(*v_it);
            
      // TAG: DEBUG
      //std::cout << "(" << NT(*u_it) << "," << NT(*v_it) << "), ";
    }
        
    // TAG: DEBUG
    //std::cout << std::endl;
    
    return sum;
  }
    
  // in-place update
  template < class ForIt >                                    // QP case
  void  update_inplace_QP( ForIt y_l_it, ForIt y_x_it,
                           const ET& d_new, const ET& d_old)
  {
    typename Matrix::      iterator  matrix_it;
    typename Row   ::      iterator     row_it;
    typename Row   ::const_iterator      y_it1, y_it2;
    
    unsigned int  row, col, k = min_N_M_+size_BO_;
    
    // rows: 0..s-1  ( P )
    for (   row = 0,   y_it1 = y_l_it,   matrix_it = M.begin();
            row < size_E_cup_SN_;
            ++row,     ++y_it1,          ++matrix_it            ) {
    
      // columns: 0..row  ( P )
      for (   row_it =  matrix_it->begin(),   y_it2 = y_l_it;
              row_it != matrix_it->end();
              ++row_it,                       ++y_it2         ) {
    
        update_entry( *row_it, d_new, *y_it1 * *y_it2, d_old);
      }
    }
    
    // rows: l..k-1  ( Q R )
    for (   row = min_N_M_,   y_it1 = y_x_it,   matrix_it += min_N_M_-size_E_cup_SN_;
            row < k;
            ++row,     ++y_it1,          ++matrix_it       ) {
    
      // columns: 0..s-1  ( Q )
      for (   col = 0,   row_it =  matrix_it->begin(),   y_it2 = y_l_it;
              col < size_E_cup_SN_;
              ++col,     ++row_it,                       ++y_it2         ){
    
        update_entry( *row_it, d_new, *y_it1 * *y_it2, d_old);
      }
    
      // columns: l..k-1  ( R )
      for (              row_it += min_N_M_-size_E_cup_SN_,                  y_it2 = y_x_it;
                         row_it != matrix_it->end();
                         ++row_it,                       ++y_it2         ){
    
        update_entry( *row_it, d_new, *y_it1 * *y_it2, d_old);
      }
    }
  }
    
  template < class ForIt1, class ForIt2 >                     // LP case
  void  update_inplace_LP( ForIt1 x_x_it, ForIt2 y_x_it,
                           const ET& d_new, const ET& d_old)
  {
    typename Matrix::      iterator  matrix_it;
    typename Row   ::      iterator     row_it;
    ForIt1                                x_it;
    
    unsigned int  row, col;
    ET            y;

    // QP (in phase I)?
    matrix_it = M.begin();
    if ( is_QP) matrix_it += min_N_M_;

    // rows: 0..s-1  ( Q )
    for (   row = 0;
            row < size_E_cup_SN_;
            ++row,     ++y_x_it, ++matrix_it) {
    
      // columns: 0..b-1  ( Q )
      y = *y_x_it;
      for (   col = 0,   row_it =  matrix_it->begin(),   x_it = x_x_it;
              col < size_BO_;
              ++col,     ++row_it,                       ++x_it         ){
    
        update_entry( *row_it, d_new, y * *x_it, d_old);
      }
    }
  }
    
    
  template < class RandomAccessIterator >
  typename std::iterator_traits<RandomAccessIterator>::value_type 
  inv_M_B_row_dot_col( int row, RandomAccessIterator u_l_it) const
  {
    typename Row::const_iterator row_it;
    if ( is_QP) {
      row_it = M[min_N_M_ + row].begin();
    } else {
      row_it = M[row].begin();
    }
    return inner_product(row_it, u_l_it, size_BO_);        
  }

}; // QP_basis_inverse<ET, Is_LP>

// ----------------------------------------------------------------------------

// =============================
// class implementation (inline)
// =============================

// creation
template < class ET_, class Is_LP_ >  inline
QP_basis_inverse<ET_,Is_LP_>::
QP_basis_inverse( CGAL::Verbose_ostream&  verbose_ostream)
  : et0( 0), et1( 1), et2( 2),
    is_LP( Is_LP::value ), is_QP( ! is_LP),
    vout( verbose_ostream)
{ }

// transition (LP case)
template < class ET_, class Is_LP_ >  inline
void  QP_basis_inverse<ET_,Is_LP_>::
transition( )
{
  is_phaseI  = false;
  is_phaseII = true;

  CGAL_qpe_debug {
    if ( vout.verbose()) print();
  }
}

// set-up (QP case)
template < class ET_, class Is_LP_ >  inline
void  QP_basis_inverse<ET_,Is_LP_>::
set( Tag_false)
{
  M.reserve( min_N_M_);
  // only allocate empty rows
  for ( unsigned int i = 0; i < min_N_M_; ++i )
    M.push_back(Row(0, et0)); 
}
    
// set-up (LP case)
template < class ET_, class Is_LP_ >  inline
void  QP_basis_inverse<ET_,Is_LP_>::
set( Tag_true)
{
  M.reserve( min_N_M_);
  for ( unsigned int i = 0; i < min_N_M_; ++i) M.push_back( Row( min_N_M_, et0));
}

// access (QP case)
template < class ET_, class Is_LP_ >  inline
const ET_&  QP_basis_inverse<ET_,Is_LP_>::
entry( unsigned int r, unsigned int c, Tag_false) const
{
  CGAL_qpe_assertion( ( r < size_E_cup_SN_) || ( ( r >= min_N_M_) && ( r < min_N_M_+size_BO_)));
  CGAL_qpe_assertion( ( c < size_E_cup_SN_) || ( ( c >= min_N_M_) && ( c < min_N_M_+size_BO_)));
  return ( c < r ? M[ r][ c] : M[ c][ r]);
}

// access (LP case)
template < class ET_, class Is_LP_ >  inline
const ET_&  QP_basis_inverse<ET_,Is_LP_>::
entry( unsigned int r, unsigned int c, Tag_true) const
{
  CGAL_qpe_assertion( r < size_E_cup_SN_);
  CGAL_qpe_assertion( c < size_BO_);
  return M[ r][ c];
}

// in-place update
template < class ET_, class Is_LP_ >  inline
void  QP_basis_inverse<ET_,Is_LP_>::
update_entry( ET& entry, const ET& d_new, const ET& y, const ET& d_old) const
{

  // TAG: DEBUG
//  ++CGAL::QP_solver_debug::timer.counter_integral_division;

  entry *= d_new;
  entry += y;
  entry = CGAL::integral_division(entry, d_old);
}

} //namespace CGAL

#include <CGAL/QP_solver/QP_basis_inverse_impl.h>

#endif // CGAL_QP_SOLVER_QP_BASIS_INVERSE_H

// ===== EOF ==================================================================
