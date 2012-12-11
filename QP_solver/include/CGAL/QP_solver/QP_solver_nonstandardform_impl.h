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
// $URL$
// $Id$
// 
//
// Author(s)     : Sven Schoenherr
//                 Bernd Gaertner <gaertner@inf.ethz.ch>
//                 Franz Wessendorp
//                 Kaspar Fischer
//                 Yves Brise

namespace CGAL {

// Looks in x_O_v_i which bound is present for variable i and returns
// the variable's value corresponding to this bound.
//
// Precondition: Is_nonnegative is Tag_false.
template < typename Q, typename ET, typename Tags >
ET QP_solver<Q, ET, Tags>::original_variable_value_under_bounds(int i) const
{
  CGAL_assertion(!Is_nonnegative::value && i<qp_n);
  switch (x_O_v_i[i]) {
  case UPPER:
    return *(qp_u+i);
  case ZERO:
    return et0;
  case LOWER:
  case FIXED:
    return *(qp_l+i);
  case BASIC:
    CGAL_qpe_assertion(false);
  }
  return et0; // dummy
}

template < typename Q, typename ET, typename Tags >
ET QP_solver<Q, ET, Tags>::variable_numerator_value(int i) const
{
  // Returns the current value of an *original* variable.
  CGAL_qpe_assertion( 0 <= i && i < qp_n );
  if (Is_nonnegative::value) {
    if (in_B[i] < 0) 
      return et0;
    else 
      return x_B_O[in_B[i]];
  }

  // now we have nonstandard form
  typedef QP_solver<Q, ET, Tags> QP;
  switch (x_O_v_i[i]) {
  case QP::UPPER:
    return static_cast<ET>(*(qp_u+i)) * denominator_;
  case QP::ZERO:
    return et0;
  case QP::LOWER:
  case QP::FIXED:
    return static_cast<ET>(*(qp_l+i)) * denominator_;
  case QP::BASIC:
    return x_B_O[in_B[i]];
  default: // never reached
    return et0;
  }
}

template < typename Q, typename ET, typename Tags >
ET QP_solver<Q, ET, Tags>::nonbasic_original_variable_value
(int i) const
{
  if (Is_nonnegative::value)
    return et0;

  CGAL_assertion(!is_basic(i));
  return original_variable_value_under_bounds(i);
}

// Computes r_i:= A_i x_init, for i=row, where x_init is the solution
// with which the solver starts the computation. I.e., computes the
// scalar product of the row-th row of A and the vector x_init which
// contains as its entries the values original_variable_value(i),
// 0<=i<qp_n.
template < typename Q, typename ET, typename Tags >
ET  QP_solver<Q, ET, Tags>::multiply__A_ixO(int row) const
{
  ET value = et0;
  ET temp = et0;
  
  A_sparse_column_iterator it_begin, it_end, it;
  
  for (int i = 0; i < qp_n; ++i) {
    it_begin = (*(qp_A_sparse+i)).begin();
    it_end = (*(qp_A_sparse+i)).end();
    while (it_begin != it_end && it_begin->first < row) {
      ++it_begin;
    }
       
    
    // Note: the following computes
    //
    //   value += original_variable_value(i) * qp_old_A[i][row];
    //
    // but for efficiency, we only add summands that are known to be
    // nonzero.
    switch (x_O_v_i[i]) {
      case UPPER:
        value += static_cast<ET>(*(qp_u+i)) * static_cast<ET>((it_begin != it_end && it_begin->first == row) ? it_begin->second : et0);
        break;
      case LOWER:
      case FIXED:
        value += static_cast<ET>(*(qp_l+i)) * static_cast<ET>((it_begin != it_end && it_begin->first == row) ? it_begin->second : et0);
        break;
      case BASIC:
        CGAL_qpe_assertion(false);
      default:
        break;
    }
  }
  
  return value;
}


// Computes r:= A x_init, where x_init is the solution
// with which the solver starts the computation. I.e., computes the
// scalar product of A and the vector x_init which
// contains as its entries the values.
// Note: Does the same job as calling multiply__A_ixO(i) for all
// 0 <= i < qp_m, but is way more efficient in the context of
// sparse iterator type input.
// PRE: out is random access with qp_m reserved entries, zero initialized
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::multiply__AxO(Value_iterator out) const
{  
  ET temp;
  
  A_sparse_column_iterator it, it_end;
  
  for (int i = 0; i < qp_n; ++i) {
  
    temp = et0;
 
    switch (x_O_v_i[i]) {
      case UPPER:
        temp = static_cast<ET>(*(qp_u+i));
        break;
      case LOWER:
      case FIXED:
        temp = static_cast<ET> (*(qp_l+i));
        break;
      case BASIC:
        CGAL_qpe_assertion(false);
      default:
        break;
    }

     
    it = (*(qp_A_sparse+i)).begin();
    it_end = (*(qp_A_sparse+i)).end();
    while (it != it_end) {
      *(out+it->first) += temp * static_cast<ET>(it->second);
      ++it;
    }
    
  }
  
}


// Computes r_{C}:= A_{C, N_O} x_{N_O}.
//
// Precondition: this routine should only be called for nonstandard form
// problems.
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
multiply__A_CxN_O(Value_iterator out) const
{
  CGAL_qpe_assertion(!Is_nonnegative::value);
  
  // initialize with zero vector:
  std::fill_n(out, C.size(), et0);
  
  A_sparse_column_iterator it;
  A_sparse_column_iterator it_end;
  
  if (no_ineq) { // in_C is not kept up to date, we have to set it up
    std::fill_n(in_C.begin(), qp_m, -1);
    for (int i = 0; i < qp_m; ++i) { // all equalities are in C
      in_C[ C[i] ] = i;
    }
  }  
  for (int i = 0; i < qp_n; ++i) {
    if (!is_basic(i)) {
      const ET x_i = nonbasic_original_variable_value(i);

      it = (*(qp_A_sparse+i)).begin();
      it_end = (*(qp_A_sparse+i)).end();
      while (it != it_end) {
        if (in_C[it->first] >= 0) {
          *(out + in_C[it->first]) += x_i * static_cast<ET>(it->second);
        }
        ++it;
      }
      
      /*
      Value_iterator out_it = out;
      for (Index_const_iterator row_it = C.begin(); row_it != C.end(); ++row_it, ++out_it)
        *out_it += x_i * static_cast<ET>(*(a_col+ *row_it));
        */
    }
  }
}

// Computes w:= 2D_{O, N_O} x_{N_O}.
//
// Precondition: this routine should only be called for nonstandard form
// problems.
//
// todo: In order to optimize this routine, we can if D is symmetric,
// multiply by two at the end of the computation instead of at each
// access to D. (Maybe its also faster to call
// nonbasic_original_variable_value() only O(n) times and not O(n^2)
// times.)
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
multiply__2D_OxN_O(Value_iterator out) const
{
  CGAL_qpe_assertion(!Is_nonnegative::value);

  // initialize with zero vector:
  std::fill_n(out, B_O.size(), et0);
  
  for (int row_it = 0; row_it < qp_n; ++row_it, ++out) {
    D_sparse_column_iterator it = (*(qp_D_sparse+row_it)).begin();
    D_sparse_column_iterator it_end = (*(qp_D_sparse+row_it)).end();
    while (it != it_end) {
      if (!is_basic(it->first)) {
        const ET value = nonbasic_original_variable_value(it->first);
        *out += static_cast<ET>(it->second) * value;
      }
      ++it;
    }
  }
}

// Computes r_{S_B}:= A_{S_B, N_O} x_{N_O}.
//
// Precondition: this routine should only be called for nonstandard form
// problems.
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
multiply__A_S_BxN_O(Value_iterator out) const
{
  // initialize with zero vector:
  std::fill_n(out, S_B.size(), et0);
  
  Indices in_S_B(qp_m, -1); // TAG: TODO maybe make this global
  int i = 0;
  for (Index_const_iterator S_B_it = S_B.begin(); S_B_it != S_B.end(); ++S_B_it) {
    in_S_B[*S_B_it] = i;
    ++i;
  }
  A_sparse_column_iterator it;
  A_sparse_column_iterator it_end;
  
  
  for (int i = 0; i < qp_n; ++i) {
    if (!is_basic(i)) {
      const ET x_i = nonbasic_original_variable_value(i);

      // reset A iterators
      it = (*(qp_A_sparse+i)).begin();
      it_end = (*(qp_A_sparse+i)).end();
      while (it != it_end) {
        if (in_S_B[it->first] >= 0) {
          *(out + in_S_B[it->first]) += x_i * static_cast<ET>(it->second);
        }
        ++it;
      }
    }
  }
}

// Initialize r_B_O.
//
// Note: this routine is called from transition() (and not during the
// initialization of the QP-solver).
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
init_r_B_O()
{
  CGAL_qpe_assertion(!Is_nonnegative::value &&
			!Is_linear::value);
  r_B_O.resize(B_O.size());
  multiply__2D_B_OxN_O(r_B_O.begin());
}

// Initialize w.
//
// Note: this routine is called from transition() (and not during the
// initialization of the QP-solver).
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
init_w()
{
  CGAL_qpe_assertion(!Is_nonnegative::value &&
			!Is_linear::value);
  w.resize(qp_n);
  multiply__2D_OxN_O(w.begin());
}

} //namespace CGAL

// ===== EOF ==================================================================
