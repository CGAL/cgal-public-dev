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

#include <CGAL/QP_solver/Initialization.h>
#include <CGAL/Timer.h>

namespace CGAL {

// =============================
// class implementation (cont'd)
// =============================

// transition (to phase II)
// ------------------------
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
transition( )
{
  CGAL_qpe_debug {
    if ( vout.verbose()) {
      vout2 << std::endl
      << "----------------------" << std::endl
      << 'T';
      vout1 << "[ t"; vout  << "ransition to phase II"; vout1 << " ]";
      vout  << std::endl;
      vout2 << "----------------------";
    }
  }
  
  
  // reset cycle bookkeping
  fingerprint_map_.clear();
  
  // update status
  m_phase    = 2;
  is_phaseI  = false;
  is_phaseII = true;
  
  // remove artificial variables
  in_B.erase( in_B.begin()+qp_n+slack_A.size(), in_B.end());
  //ensure_size(tmp_x_2, tmp_x.size());
  // update basis inverse
  CGAL_qpe_debug {
    vout4 << std::endl << "basis-inverse:" << std::endl;
  }
  
  transition( Is_linear());

  
  // initialize exact version of `-qp_c' (implicit conversion to ET)
  C_by_index_accessor  c_accessor( qp_c);
  std::transform( C_by_index_iterator( B_O.begin(), c_accessor),
                 C_by_index_iterator( B_O.end  (), c_accessor),
                 minus_c_B.begin(), std::negate<ET>());
  
  // compute initial solution of phase II
  compute_solution(Is_nonnegative());

  
  // diagnostic output
  CGAL_qpe_debug {
    if ( vout.verbose()) print_solution();
  }
  
  // notify pricing strategy
  strategyP->transition();
}
  
// access
// ------
// numerator of current solution; the denominator is 2*d*d, so we should
// compute here d*d*(x^T2Dx + x^T2c + 2c0)
// TODO: The following is not really sparse... temp vector is zero-constructed
// in every iteration of the loop. It's not so paramount though, because
// solution_numerator() is only called in debug output.
template < typename Q, typename ET, typename Tags >
ET QP_solver<Q, ET, Tags>::
solution_numerator( ) const
{
  ET   s, z = et0;
  int  i, j;
  
  if (Is_nonnegative::value || is_phaseI) {
    // standard form or phase I; it suffices to go
    // through the basic variables; all D- and c-entries
    // are obtained through the appropriate iterators 
    Index_const_iterator  i_it;
    Value_const_iterator  x_i_it, c_it;
    
    // foreach i
    x_i_it =       x_B_O.begin();
    c_it = minus_c_B  .begin();
    for ( i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++x_i_it, ++c_it){
      
      i = *i_it;
      
      // compute quadratic part: 2D_i x
      s = et0;
      if (is_QP && is_phaseII) {
        D_sparse_column_iterator it = (*(qp_D_sparse+i)).begin();
        D_sparse_column_iterator it_end = (*(qp_D_sparse+i)).end();
        Values temp(B_O.size(), et0);
        ET diag_temp = et0;
        while (it != it_end) {
          if (it->first == i) diag_temp = it->second;
          if (in_B[it->first] > -1) {
            temp[in_B[it->first]] = it->second;
          }
          ++it;
        }
        // half the off-diagonal contribution
        s += std::inner_product(x_B_O.begin(), x_i_it, temp.begin(), et0);
        // the other half
        s *= et2;
        // diagonal contribution
        s += diag_temp * *x_i_it;
      }
      s -= denominator_ * et2 * static_cast<ET>(*c_it );
      
      // add x_i(2D_i x + 2c_i)
      z += s * *x_i_it; // endowed with a factor of d*d now
    }
  } else {
    // nonstandard form and phase II, 
    // take all original variables into account; all D- and c-entries
    // are obtained from the input data directly
    // order in i_it and j_it matches original variable order
    if (is_QP) {
    
      // quadratic part
      i=0;
      for (Variable_numerator_iterator 
           i_it = this->original_variables_numerator_begin(); 
           i_it < this->original_variables_numerator_end(); ++i_it, ++i) {
        // do something only if *i_it != 0
        if (*i_it == et0) continue;
        s = et0; // contribution of i-th row
        
        D_sparse_column_iterator it = (*(qp_D_sparse+i)).begin();
        D_sparse_column_iterator it_end = (*(qp_D_sparse+i)).end();
        Variable_numerator_iterator j_it = this->original_variables_numerator_begin();

        // half the off-diagonal contribution        
        while (it != it_end && it->first < i) {
          s += static_cast<ET>(it->second) * *(j_it+it->first);
          ++it;
        }
        // the other half
        s *= et2;
        // the diagonal entry
        if (it != it_end && it->first == i) s += static_cast<ET>(it->second) * *(j_it+i);
        
        // accumulate
        z += s * *i_it;
      }
    }
    // linear part
    j=0; s = et0;
    for (Variable_numerator_iterator 
         j_it = this->original_variables_numerator_begin();
         j_it < this->original_variables_numerator_end(); ++j_it, ++j)
      s +=  et2 * static_cast<ET>(*(qp_c+j)) * *j_it;
    z += denominator_ * s;
  }
  
  // finally, add the constant term (phase II only)
  if (is_phaseII) z += et2 * static_cast<ET>(qp_c0) * denominator_ * denominator_;
  return z;
}

// pivot step
// ----------
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
pivot_step( )
{
  ++m_pivots;
  
  // diagnostic output
  CGAL_qpe_debug {
    vout2 << std::endl
    << "==========" << std::endl
    << "Pivot Step" << std::endl
    << "==========" << std::endl;
  }
  vout  << "[ phase " << ( is_phaseI ? "I" : "II")
  << ", iteration " << m_pivots << ", LU ]" << std::endl;
  
  
  // pricing
  // -------
  
  CGAL_qpe_debug {
    //CGAL::QP_solver_debug::timer.pricing.start();
  }
  
  // we need to make sure that lambda_sorted is properly set up.
  // we assume that it is filled with ET(0)'s at the beginning.
  for (int i = 0; i < static_cast<int>(C.size()); ++i) {
    lambda_sorted[ C[i] ] = lambda[i];
  }
  is_lambda_sorted = true;
  
  // actual pricing step
  pricing();
  
  // after the pricing we need to reset lambda_sorted to contain ET(0)'s
  for (int i = 0; i < static_cast<int>(C.size()); ++i) {
    lambda_sorted[ C[i] ] = et0;
  }
  is_lambda_sorted = false;
  
  
  CGAL_qpe_debug { 
    //CGAL::QP_solver_debug::timer.pricing.stop();
  }
  
  // check for optimality
  if ( index_entering < 0) {
    
    if ( is_phaseI) {                               // phase I
      // since we no longer assume full row rank and subsys assumption
      // we have to strengthen the precondition for infeasibility
      if (this->solution_numerator() > et0) {    
        // problem is infeasible
        m_phase  = 3;
        m_status = QP_INFEASIBLE;
        
        vout << "  ";
        vout << "problem is INFEASIBLE" << std::endl;
        
      } else {  // Drive/remove artificials out of basis
        expel_artificial_variables_from_basis();
        transition();
      }
    } else {                                        // phase II
      
      // optimal solution found
      m_phase  = 3;
      m_status = QP_OPTIMAL;
      
      vout << "  ";
      vout  << "solution is OPTIMAL" << std::endl;
      
    }
    return;
  }
  
  CGAL_qpe_debug {
    //CGAL::QP_solver_debug::timer.ratio_test.start();
  }
    
  
  
  // ratio test & update (step 1)
  // ----------------------------
  // initialize ratio test
  ratio_test_init();
  
  
  
  // diagnostic output
  CGAL_qpe_debug {
    if ( vout2.verbose() && is_QP && is_phaseII) {
      vout2.out() << std::endl
      << "----------------------------" << std::endl
      << "Ratio Test & Update (Step 1)" << std::endl
      << "----------------------------" << std::endl;
    }
  }
  
  
  CGAL_qpe_debug {
    //CGAL::QP_solver_debug::timer.ratio_test_1.start();
  }
  
  // loop (step 1)
  do {
    // ratio test
    ratio_test_1();
    
  
    // check for unboundedness
    if ( q_i == et0) {
      m_phase  = 3;
      m_status = QP_UNBOUNDED;
      
      vout << "  ";
      vout << "problem is UNBOUNDED" << std::endl;
      
      CGAL_qpe_debug {
        //nu should be zero in this case
        // note: (-1)/hat{\nu} is stored instead of \hat{\nu}
        // todo kf: as this is just used for an assertion check,
        // all the following lines should only be executed if
        // assertions are enabled...
        
        nu = lu_fact_.inner_product(     A_Cj.begin(), two_D_Bj.begin(),
                                    q_lambda.begin(),    q_x_O.begin());    
        
        if (is_QP) {
          // TAG: INEFFICIENT maybe we need the whole column earlier, then we would not have
          // to traverse it for each entry separately.
          if (index_entering < qp_n) { 
            D_sparse_column_iterator it = (*(qp_D_sparse+index_entering)).begin();
            D_sparse_column_iterator it_end = (*(qp_D_sparse+index_entering)).end();
            ET tmp(0);
            while (it != it_end && it->first < index_entering) ++it;
            if (it != it_end && it->first == index_entering) tmp = it->second;
            nu -= denominator_ * tmp;
          }
        }
        CGAL_qpe_assertion(nu == et0);
      }
      return;
    }
    // update
    update_1();
  } while ( index_entering >= 0);
  
  
  CGAL_qpe_debug {
    //CGAL::QP_solver_debug::timer.ratio_test_1.stop();
    //CGAL::QP_solver_debug::timer.ratio_test_2.start();
  }
  
  // ratio test & update (step 2)
  // ----------------------------
  /*    
   if ( i >= 0) {
   
   // diagnostic output
   CGAL_qpe_debug {
   vout2 << std::endl
   << "----------------------------" << std::endl
   << "Ratio Test & Update (Step 2)" << std::endl
   << "----------------------------" << std::endl;
   }
   
   // compute index of entering variable
   j += in_B.size();
   
   // loop (step 2)
   while ( ( i >= 0) && basis_matrix_stays_regular()) {
   
   // update
   update_2( Is_linear());
   
   // ratio test
   ratio_test_2( Is_linear());
   }
   }
   */
  // instead of the above piece of code we now have
  // diagnostic output
  if (is_RTS_transition) {
    is_RTS_transition = false;
    
    CGAL_qpe_debug {
      vout2 << std::endl
      << "----------------------------" << std::endl
      << "Ratio Test & Update (Step 2)" << std::endl
      << "----------------------------" << std::endl;
    }
    
    // compute index of entering variable
    index_entering += static_cast<int>(in_B.size());
    
    ratio_test_2( Is_linear());
    
    while ((index_leaving >= 0) && basis_matrix_stays_regular()) {
      
      update_2(Is_linear());
      
      ratio_test_2(Is_linear());
      
    }
  } 
  
  CGAL_qpe_debug { 
    //CGAL::QP_solver_debug::timer.ratio_test_2.stop();
    //CGAL::QP_solver_debug::timer.ratio_test_3.start();
  }
  
  // ratio test & update (step 3)
  // ----------------------------
  CGAL_qpe_assertion_msg( index_leaving < 0, "Step 3 should never be reached!");
  
  // diagnostic output
  
  if ( vout.verbose()) print_basis();
  if ( vout.verbose()) print_solution();
  
  // transition to phase II (if possible)
  // ------------------------------------
  if ( is_phaseI && ( art_basic == 0)) {
    CGAL_qpe_debug {
      if ( vout2.verbose()) {
        vout2.out() << std::endl
        << "all artificial variables are nonbasic"
        << std::endl;
      }
    }
    transition();
  }
  
  CGAL_qpe_debug { 
    //CGAL::QP_solver_debug::timer.ratio_test_3.stop();
    //CGAL::QP_solver_debug::timer.ratio_test.stop();
  }
  
  // Cycle detection bookkeeping
  if (!bland_flag) {
    if (!progress_flag) {
      QP_solver_impl::Iteration_fingerprint<Indices> current_fingerprint(B_O, B_S);
      current_fingerprint_hash_ = QP_solver_impl::hash_value<Indices>(current_fingerprint);
      fingerprint_map_[current_fingerprint_hash_].push_back(std::pair<int, QP_solver_impl::Iteration_fingerprint<Indices> >(m_pivots, current_fingerprint));
    } else {
      fingerprint_map_.clear();
    }
  }
  
  
}

// pricing
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
pricing( )
{
  // diagnostic output
  CGAL_qpe_debug {
	  if ( vout2.verbose() ) {
	    vout2 << std::endl
		  << "-------" << std::endl
		  << "Pricing" << std::endl
		  << "-------" << std::endl;
	  }
  }

  // call pricing strategy
  index_entering = strategyP->pricing(direction);

  // diagnostic output
  if ( vout.verbose() ) {
      if ( index_entering < 0) {
    	  CGAL_qpe_debug {
	        vout2 << "entering variable: none" << std::endl;
	      }
      } else {
        vout  << "  ";
	    vout  << "entering: ";
	    vout  << index_entering;
	    CGAL_qpe_debug {
	      vout2 << " (" << variable_type( index_entering) << ')' << std::endl;
	      vout2 << "direction: "
		    << ((direction == 1) ? "positive" : "negative") << std::endl;
	    }
    }
  }
}


// cycle detection
// ----------
template < typename Q, typename ET, typename Tags >
int
QP_solver<Q, ET, Tags>::
cycle_detection( )
{
  int ret(0);
  
  int n_items = fingerprint_map_[current_fingerprint_hash_].size();
  for (int i = n_items - 2; i >= 0; --i) {
    if (fingerprint_map_[current_fingerprint_hash_][i].second == fingerprint_map_[current_fingerprint_hash_][n_items - 1].second) {
      ret = fingerprint_map_[current_fingerprint_hash_][n_items - 1].first - fingerprint_map_[current_fingerprint_hash_][i].first;
      break;
    }
  }
  
  return ret;
}

// initialization of ratio-test
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
ratio_test_init( )
{
    // store exact version of `A_Cj' (implicit conversion)
    ratio_test_init__A_Cj( A_Cj.begin(), index_entering, no_ineq);

    // store exact version of `2 D_{B_O,j}'
    ratio_test_init__2_D_Bj( two_D_Bj.begin(), index_entering, Is_linear());
}

template < typename Q, typename ET, typename Tags >                                         // no ineq.
void  QP_solver<Q, ET, Tags>::
ratio_test_init__A_Cj( Value_iterator A_Cj_it, int j_, Tag_true)
{
  // store exact version of `A_Cj' (implicit conversion)
  if ( j_ < qp_n) {                                   // original variable
    
    std::fill_n( A_Cj_it, qp_m, et0); // TAG: INEFFICIENT... it's a pity that we have to fill with zeros first.
    A_sparse_column_iterator it = (*(qp_A_sparse+j_)).begin();
    A_sparse_column_iterator it_end = (*(qp_A_sparse+j_)).end();
    while (it != it_end) {
      *(A_Cj_it + it->first) = it->second;
      ++it;
    }
  } else {                                            // artificial variable
    
    unsigned int  k = j_;
    k -= qp_n;
    std::fill_n( A_Cj_it, qp_m, et0);
    A_Cj_it[ k] = ( art_A[ k].second ? -et1 : et1);
  }
}

template < typename Q, typename ET, typename Tags >     // has ineq.
void  QP_solver<Q, ET, Tags>::
ratio_test_init__A_Cj( Value_iterator A_Cj_it, int j_, Tag_false)
{
  // store exact version of `A_Cj' (implicit conversion)
  if ( j_ < qp_n) {                                   // original 
    
    // TAG: INEFFICIENT because it's not C-filtered
    std::fill_n( A_Cj_it, C.size(), et0);
    A_sparse_column_iterator it = (*(qp_A_sparse+j_)).begin();
    A_sparse_column_iterator it_end = (*(qp_A_sparse+j_)).end();
    while (it != it_end) {
      if (in_C[it->first] >= 0) {
        *(A_Cj_it + in_C[it->first]) = it->second;
      }
      ++it;
    }
  } else {
    unsigned int  k = j_;
    k -= qp_n;
    std::fill_n( A_Cj_it, C.size(), et0);
    
    if ( k < static_cast<unsigned int>(slack_A.size()) ) {                      // slack variable
      A_Cj_it[ in_C[ slack_A[ k].first]] = ( slack_A[ k].second ? -et1
                                            :  et1);
    } else {                                        // artificial variable
      k -= static_cast<unsigned int>(slack_A.size());
      
      if ( j_ != art_s_i) {                           // normal art.
        A_Cj_it[ in_C[ art_A[ k].first]] = ( art_A[ k].second ? -et1
                                            :  et1);
        
      } else {                                        // special art.
        S_by_index_accessor  s_accessor( art_s.begin());
        std::copy( S_by_index_iterator( C.begin(), s_accessor),
                  S_by_index_iterator( C.end  (), s_accessor),
                  A_Cj_it);
      }	
    }
  }
}

    
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
init__A_Ri( Value_iterator A_Ri_it, int i_, bool no_ineq)
{

  // store exact version of `A_Ri' (implicit conversion)
  std::fill_n( A_Ri_it, B_O.size(), et0);
  
  
  // original variables
  for (int j = 0; j < B_O.size(); ++j) {
    if (B_O[j] < qp_n) { // original variable
      A_sparse_column_iterator it = (*(qp_A_sparse+B_O[j])).begin();
      A_sparse_column_iterator it_end = (*(qp_A_sparse+B_O[j])).end();
      while (it != it_end && it->first < i_) {
        ++it;
      }
      if (it != it_end && it->first == i_) {
        *(A_Ri_it + j) = it->second;
      }
    } else {
      if (B_O[j] != art_s_i) { // regular artificial
        unsigned int k = B_O[j] - qp_n;
        if (!no_ineq) k -= slack_A.size();
        if (art_A[k].first == i_) {
          *(A_Ri_it + j) = ( art_A[ k].second ? -et1 : et1);
        }
      } else { // special artificial
        *(A_Ri_it + in_B[art_s_i]) = static_cast<ET>(art_s[i_]);
      }
    }
  }
}

template < typename Q, typename ET, typename Tags >                                         // no ineq.
void  QP_solver<Q, ET, Tags>::
init__A_ij( ET& val, int i, int j, Tag_true)
{
  // store exact version of `A_ij' (implicit conversion)
  val = et0;
  if ( j < qp_n) {                                   // original variable
    A_sparse_column_iterator it = (*(qp_A_sparse+j)).begin();
    A_sparse_column_iterator it_end = (*(qp_A_sparse+j)).end();
    while (it != it_end && it->first <= i) {
      if (it->first == i) {
        val = it->second;
      }
      ++it;
    }
  } else {                                            // artificial variable
    unsigned int  k = j;
    k -= qp_n;
    if (k == i) {
      val = ( art_A[ k].second ? -et1 : et1);
    }
  }
}

template < typename Q, typename ET, typename Tags >     // has ineq.
void  QP_solver<Q, ET, Tags>::
init__A_ij( ET& val, int i, int j, Tag_false)
{
  // store exact version of `A_Cj' (implicit conversion)
  val = et0;
  if ( j < qp_n) {                                   // original 
    A_sparse_column_iterator it = (*(qp_A_sparse+j)).begin();
    A_sparse_column_iterator it_end = (*(qp_A_sparse+j)).end();
    while (it != it_end) {
      if (in_C[it->first] == i) {
        val = it->second;
        return;
      }
      ++it;
    }
  } else {
    unsigned int  k = j;
    k -= qp_n;
    
    if ( k < static_cast<unsigned int>(slack_A.size()) ) {                      // slack variable
    
      if (in_C[ slack_A[ k].first] == i) {
        val = ( slack_A[ k].second ? -et1 : et1);
      }
    } else {                                        // artificial variable
      k -= static_cast<unsigned int>(slack_A.size());
      if ( j != art_s_i) {                           // normal art.
        if (in_C[ art_A[ k].first] == i) {
          val = ( art_A[ k].second ? -et1 : et1);
        }        
      } else {                                        // special art.
        val = art_s[C[i]];
      }	
    }
  }
}

/////////////////////////////////////////////


// ratio test (step 1)
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
ratio_test_1( )
{
  
  // diagnostic output
  CGAL_qpe_debug {
    if ( vout2.verbose()) {
      vout2.out() << std::endl;
      if ( is_LP || is_phaseI) {
        vout2.out() << "----------" << std::endl
        << "Ratio Test" << std::endl
        << "----------" << std::endl;
      } else {
        vout2.out() << "Ratio Test (Step 1)" << std::endl
        << "-------------------" << std::endl;
      }
      if ( vout3.verbose()) {
        vout3.out() << "    A_Cj: ";
        std::copy( A_Cj.begin(), A_Cj.begin()+C.size(),
                  std::ostream_iterator<ET>( vout3.out()," "));
        vout3.out() << std::endl;
        if ( is_QP && is_phaseII) {
          vout3.out() << "  2 D_Bj: ";
          std::copy( two_D_Bj.begin(), two_D_Bj.begin()+B_O.size(),
                    std::ostream_iterator<ET>( vout3.out()," "));
          vout3.out() << std::endl;
        }
        vout3.out() << std::endl;
      }
    }
  }
  
  
  // compute `q_lambda' and `q_x'
  ratio_test_1__q_x_O( Is_linear());
  ratio_test_1__q_x_S( no_ineq);
  
  
  // diagnostic output
  CGAL_qpe_debug {
    if ( vout3.verbose()) {
      if ( is_QP && is_phaseII) {
        vout3.out() << "q_lambda: ";
        std::copy( q_lambda.begin(), q_lambda.begin()+C.size(),
                  std::ostream_iterator<ET>( vout3.out()," "));
        vout3.out() << std::endl;
      }
      vout3.out() << "   q_x_O: ";
      std::copy( q_x_O.begin(), q_x_O.begin()+B_O.size(),
                std::ostream_iterator<ET>( vout3.out()," "));
      vout3.out() << std::endl;
      
      if ( has_ineq) {
        vout3.out() << "   q_x_S: ";
        std::copy( q_x_S.begin(), q_x_S.begin()+B_S.size(),
                  std::ostream_iterator<ET>( vout3.out()," "));
        vout3.out() << std::endl;
      }
      vout3.out() << std::endl;
    }
  }
  
  // check `t_i's
  x_i = et1;                                          // trick: initialize
  q_i = et0;                                          // minimum with +oo
  
  // computation of t_{min}^{j}
  ratio_test_1__t_min_j(Is_nonnegative());
  CGAL_qpe_debug { // todo kf: at first sight, this debug message should
    // only be output for problems in nonstandard form...
    if (vout2.verbose()) {
      vout2.out() << "t_min_j: " << x_i << '/' << q_i << std::endl;
      vout2.out() << std::endl;
    }
  }    
  
  
  // what happens, if all original variables are nonbasic?
  /*
   ratio_test_1__t_i(   B_O.begin(),   B_O.end(),
   x_B_O.begin(), q_x_O.begin(), Tag_false());
   ratio_test_1__t_i(   B_S.begin(),   B_S.end(),
   x_B_S.begin(), q_x_S.begin(), no_ineq);
   */		       
  ratio_test_1__t_min_B(no_ineq);    
  
  
  // check `t_j'
  ratio_test_1__t_j( Is_linear());
  
  // set progress flag for cycle detection
  progress_flag = ((q_i != et0 && x_i != 0) ? true : false);
  
  // diagnostic output
  CGAL_qpe_debug {
    if ( vout2.verbose()) {
      for ( unsigned int k = 0; k < static_cast<unsigned int>(B_O.size()); ++k) {
        print_ratio_1_original(k, x_B_O[k], q_x_O[k]);
      }     
      if ( has_ineq) {
        for ( unsigned int k = 0; k < static_cast<unsigned int>(B_S.size()); ++k) {
          /*
           vout2.out() << "t_S_" << k << ": "
           << x_B_S[ k] << '/' << q_x_S[ k]
           << ( ( q_i > et0) && ( i == B_S[ k]) ? " *":"")
           << std::endl;
           */
          print_ratio_1_slack(k, x_B_S[k], q_x_S[k]);
        }
      }
      if ( is_QP && is_phaseII) {
        vout2.out() << std::endl
        << "  t_j: " << mu << '/' << nu
        << ( ( q_i > et0) && ( index_leaving < 0) ? " *" : "")
        << std::endl;
      }
      vout2.out() << std::endl;
    }
  }
  if ( q_i > et0) {
    if ( index_leaving < 0) {
      vout2 << "leaving variable: none" << std::endl;
    } else {
      vout << ", ";
      vout  << "leaving: ";
      vout  << index_leaving;
      CGAL_qpe_debug {
        if ( vout2.verbose()) {
          if ( ( index_leaving < qp_n) || ( index_leaving >= static_cast<int>( qp_n+slack_A.size()))) {
            vout2.out() << " (= B_O[ " << in_B[ index_leaving] << "]: "
            << variable_type( index_leaving) << ')';
          } else {
            vout2.out() << " (= B_S[ " << in_B[ index_leaving] << "]: slack)";
          }
        }
        vout2 << std::endl;
      }
      
    }
  }
  
}

  
template < typename Q, typename ET, typename Tags >                         // Standard form
void  QP_solver<Q, ET, Tags>::
ratio_test_1__t_min_j(Tag_true /*is_nonnegative*/)
{
}

// By the pricing step we have the following precondition
// direction == +1 => x_O_v_i[j] == (LOWER v ZERO)
// direction == -1 => x_O_v_i[j] == (UPPER v ZERO) 
template < typename Q, typename ET, typename Tags >                         // Upper bounded
void  QP_solver<Q, ET, Tags>::
ratio_test_1__t_min_j(Tag_false /*is_nonnegative*/)
{
    if (index_entering < qp_n) {                                 // original variable
        if (direction == 1) {
            if (x_O_v_i[index_entering] == LOWER) {              // has lower bound value
                if (*(qp_fu+index_entering)) {                   // has finite upper bound
                    x_i = (*(qp_u+index_entering) - *(qp_l+index_entering));
                    q_i = et1;
                    index_leaving = index_entering;
                    ratio_test_bound_index = UPPER;
                } else {                            // has infinite upper bound
                    x_i = et1;
                    q_i = et0;
                }
            } else {                                // has value zero
                if (*(qp_fu+index_entering)) {                   // has finite upper bound
                    x_i = *(qp_u+index_entering);
                    q_i = et1;
                    index_leaving = index_entering;
                    ratio_test_bound_index = UPPER;
                } else {                            // has infinite upper bound
                    x_i = et1;
                    q_i = et0;                    
                }
            }
        } else {                                    // direction == -1
            if (x_O_v_i[index_entering] == UPPER) {              // has upper bound value
                if (*(qp_fl+index_entering)) {                   // has finite lower bound
                    x_i = (*(qp_u+index_entering) - *(qp_l+index_entering));
                    q_i = et1;
                    index_leaving = index_entering;
                    ratio_test_bound_index = LOWER;
                } else {                            // has infinite lower bound
                    x_i = et1;
                    q_i = et0;
                }
            } else {                                // has value zero
                if (*(qp_fl+index_entering)) {                   // has finite lower bound
                    x_i = -(*(qp_l+index_entering));
                    q_i = et1;
                    index_leaving = index_entering;
                    ratio_test_bound_index = LOWER;
                } else {                            // has infinite lower bound
                    x_i = et1;
                    q_i = et0;
                }
            }
        }
    } else {                                        // slack or artificial var
        x_i = et1;
        q_i = et0;
    }
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
ratio_test_1__t_min_B(Tag_true  /*has_equalities_only_and_full_rank*/)
{
    ratio_test_1_B_O__t_i(B_O.begin(), B_O.end(), x_B_O.begin(),
                        q_x_O.begin(), Is_nonnegative());
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
ratio_test_1__t_min_B(Tag_false /*has_equalities_only_and_full_rank*/)
{
    ratio_test_1_B_O__t_i(B_O.begin(), B_O.end(), x_B_O.begin(),
                        q_x_O.begin(), Is_nonnegative());
    ratio_test_1_B_S__t_i(B_S.begin(), B_S.end(), x_B_S.begin(),
                        q_x_S.begin(), Is_nonnegative());
}    

// ratio test for the basic original variables
template < typename Q, typename ET, typename Tags >                         // Standard form
void  QP_solver<Q, ET, Tags>::
ratio_test_1_B_O__t_i(Index_iterator i_it, Index_iterator end_it,
                    Value_iterator x_it, Value_iterator q_it,
                    Tag_true  /*is_nonnegative*/)
{    
    for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
        test_implicit_bounds_dir_pos(*i_it, *x_it, *q_it, index_leaving, x_i, q_i);
    }
}

// ratio test for the basic original variables                    
template < typename Q, typename ET, typename Tags >                         // Upper bounded
void  QP_solver<Q, ET, Tags>::
ratio_test_1_B_O__t_i(Index_iterator i_it, Index_iterator end_it,
                    Value_iterator x_it, Value_iterator q_it,
                    Tag_false /*is_nonnegative*/)
{
    if (is_phaseI) {
        if (direction == 1) {
            for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
                test_mixed_bounds_dir_pos(*i_it, *x_it, *q_it, index_leaving, x_i, q_i);
            }
        } else {
            for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
                test_mixed_bounds_dir_neg(*i_it, *x_it, *q_it, index_leaving, x_i, q_i);
            }
        }
    } else {
        if (direction == 1) {
            for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
                test_explicit_bounds_dir_pos(*i_it, *x_it, *q_it, index_leaving, x_i, q_i);
            }
        } else {
            for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
                test_explicit_bounds_dir_neg(*i_it, *x_it, *q_it, index_leaving, x_i, q_i);
            }
        }
    }
}

// ratio test for the basic slack variables
template < typename Q, typename ET, typename Tags >                         // Standard form
void  QP_solver<Q, ET, Tags>::
ratio_test_1_B_S__t_i(Index_iterator i_it, Index_iterator end_it,
                Value_iterator x_it, Value_iterator q_it,
                Tag_true  /*is_nonnegative*/)
{
    for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
        test_implicit_bounds_dir_pos(*i_it, *x_it, *q_it, index_leaving, x_i, q_i);
    }
}

// ratio test for the basic slack variables
template < typename Q, typename ET, typename Tags >                         // Upper bounded
void  QP_solver<Q, ET, Tags>::
ratio_test_1_B_S__t_i(Index_iterator i_it, Index_iterator end_it,
                Value_iterator x_it, Value_iterator q_it,
                Tag_false /*is_nonnegative*/)
{
    if (direction == 1) {
        for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
            test_implicit_bounds_dir_pos(*i_it, *x_it, *q_it, index_leaving, x_i, q_i);
        }
    } else {
        for ( ; i_it != end_it; ++i_it, ++x_it, ++q_it ) {
            test_implicit_bounds_dir_neg(*i_it, *x_it, *q_it, index_leaving, x_i, q_i);
        }    
    }
}

// test for one basic variable with implicit bounds only,
// note that this function writes the member variables i, x_i, q_i
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
test_implicit_bounds_dir_pos(int k, const ET& x_k, const ET& q_k, 
                                int& i_min, ET& d_min, ET& q_min)
{
  if (q_k > et0) {
    // BLAND rule: In case the ratios are the same, only update if the new index
    // is smaller. The special artificial variable is always made to leave first.
    if ((x_k * q_min < d_min * q_k) || ((k < i_min) && (i_min != art_s_i) && (x_k * q_min == d_min * q_k)) ) {
      i_min = k;
      d_min = x_k;
      q_min = q_k;
    }
  }
}

// test for one basic variable with implicit bounds only,
// note that this function writes the member variables i, x_i, q_i
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
test_implicit_bounds_dir_neg(int k, const ET& x_k, const ET& q_k, 
                                int& i_min, ET& d_min, ET& q_min)
{
  if (q_k < et0) {
    // BLAND rule: In case the ratios are the same, only update if the new index
    // is smaller. The special artificial variable is always made to leave first.
    if ((x_k * q_min < -(d_min * q_k)) || ((k < i_min) && (i_min != art_s_i) && (x_k * q_min == -(d_min * q_k))) ) {
      i_min = k;
      d_min = x_k;
      q_min = -q_k;
    }
  }
}

// test for one basic variable with explicit bounds only,
// note that this function writes the member variables i, x_i, q_i and
// ratio_test_bound_index, although the second and third variable name
// are in the context of upper bounding misnomers
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
test_explicit_bounds_dir_pos(int k, const ET& x_k, const ET& q_k, 
                                int& i_min, ET& d_min, ET& q_min)
{
  if (q_k > et0) {                                // check for lower bound
    if (*(qp_fl+k)) {
      ET  diff = x_k - (denominator_ * static_cast<ET>(*(qp_l+k)));
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ((diff * q_min < d_min * q_k) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == d_min * q_k)) ) {
        i_min = k;
        d_min = diff;
        q_min = q_k;
        ratio_test_bound_index = LOWER;
      }
    }
  } else {                                        // check for upper bound
    if ((q_k < et0) && (*(qp_fu+k))) {
      ET  diff = (denominator_ * static_cast<ET>(*(qp_u+k))) - x_k;
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ((diff * q_min < -(d_min * q_k)) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == -(d_min * q_k))) ) {
        i_min = k;
        d_min = diff;
        q_min = -q_k;
        ratio_test_bound_index = UPPER;
      }    
    }
  }
}

// test for one basic variable with explicit bounds only,
// note that this function writes the member variables i, x_i, q_i and
// ratio_test_bound_index, although the second and third variable name
// are in the context of upper bounding misnomers
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
test_explicit_bounds_dir_neg(int k, const ET& x_k, const ET& q_k, 
                                int& i_min, ET& d_min, ET& q_min)
{
  if (q_k < et0) {                                // check for lower bound
    if (*(qp_fl+k)) {
      ET  diff = x_k - (denominator_ * static_cast<ET>(*(qp_l+k)));
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ((diff * q_min < -(d_min * q_k)) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == -(d_min * q_k))) ) {
        i_min = k;
        d_min = diff;
        q_min = -q_k;
        ratio_test_bound_index = LOWER;
      }
    }
  } else {                                        // check for upper bound
    if ((q_k > et0) && (*(qp_fu+k))) {
      ET  diff = (denominator_ * static_cast<ET>(*(qp_u+k))) - x_k;
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ((diff * q_min < d_min * q_k) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == d_min * q_k)) ) {
        i_min = k;
        d_min = diff;
        q_min = q_k;
        ratio_test_bound_index = UPPER;
      }    
    }
  }
}

// test for one basic variable with mixed bounds,
// note that this function writes the member variables i, x_i, q_i and
// ratio_test_bound_index, although the second and third variable name
// are in the context of upper bounding misnomers
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
test_mixed_bounds_dir_pos(int k, const ET& x_k, const ET& q_k, 
                                int& i_min, ET& d_min, ET& q_min)
{
  if (q_k > et0) {                                // check for lower bound
    if (k < qp_n) {                             // original variable
      if (*(qp_fl+k)) {
        ET  diff = x_k - (denominator_ * static_cast<ET>(*(qp_l+k)));
        // BLAND rule: In case the ratios are the same, only update if the new index
        // is smaller. The special artificial variable is always made to leave first.
        if ((diff * q_min < d_min * q_k) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == d_min * q_k))) {
          i_min = k;
          d_min = diff;
          q_min = q_k;
          ratio_test_bound_index = LOWER;
        } // phase  I II switch
        
      }
    } else {                                    // artificial variable
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ((x_k * q_min < d_min * q_k) || ((k < i_min) && (i_min != art_s_i) && (x_k * q_min == d_min * q_k))) {
        i_min = k;
        d_min = x_k;
        q_min = q_k;
      }
    }
  } else {                                        // check for upper bound
    if ((q_k < et0) && (k < qp_n) && *(qp_fu+k)) {
      ET  diff = (denominator_ * static_cast<ET>(*(qp_u+k))) - x_k;
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ((diff * q_min < -(d_min * q_k)) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == -(d_min * q_k))) ) {
        i_min = k;
        d_min = diff;
        q_min = -q_k;
        ratio_test_bound_index = UPPER;
      }
    }
  }
}

// test for one basic variable with mixed bounds,
// note that this function writes the member variables i, x_i, q_i and
// ratio_test_bound_index, although the second and third variable name
// are in the context of upper bounding misnomers
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
test_mixed_bounds_dir_neg(int k, const ET& x_k, const ET& q_k, 
                                int& i_min, ET& d_min, ET& q_min)
{
  if (q_k < et0) {                                // check for lower bound
    if (k < qp_n) {                             // original variable
      if (*(qp_fl+k)) {
        ET  diff = x_k - (denominator_ * static_cast<ET>(*(qp_l+k)));
        // BLAND rule: In case the ratios are the same, only update if the new index
        // is smaller. The special artificial variable is always made to leave first.
        if ((diff * q_min < -(d_min * q_k)) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == -(d_min * q_k))) ) {
          i_min = k;
          d_min = diff;
          q_min = -q_k;
          ratio_test_bound_index = LOWER;
        }
      }
    } else {                                    // artificial variable
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ((x_k * q_min < -(d_min * q_k)) || ((k < i_min) && (i_min != art_s_i) && (x_k * q_min == -(d_min * q_k))) ) {
        i_min = k;
        d_min = x_k;
        q_min = -q_k;
      }
    }
  } else {                                        // check for upper bound
    if ((q_k > et0) && (k < qp_n) && *(qp_fu+k)) {
      ET  diff = (denominator_ * static_cast<ET>(*(qp_u+k))) - x_k;
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ((diff * q_min < d_min * q_k) || ((k < i_min) && (i_min != art_s_i) && (diff * q_min == d_min * q_k)) ) {
        i_min = k;
        d_min = diff;
        q_min = q_k;
        ratio_test_bound_index = UPPER;
      }
    }
  }
}    


template < typename Q, typename ET, typename Tags >                                         // QP case
void
QP_solver<Q, ET, Tags>::
ratio_test_2( Tag_false)
{
    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) {
	    vout2.out() << std::endl
			<< "Ratio Test (Step 2)" << std::endl
			<< "-------------------" << std::endl;
	}
    }

    // compute `p_lambda' and `p_x' (Note: `p_...' is stored in `q_...')
    ratio_test_2__p( no_ineq);
 
    // diagnostic output
    CGAL_qpe_debug {
	if ( vout3.verbose()) {
	    vout3.out() << "p_lambda: ";
	    std::copy( q_lambda.begin(), q_lambda.begin()+C.size(),
		       std::ostream_iterator<ET>( vout3.out()," "));
	    vout3.out() << std::endl;
	    vout3.out() << "   p_x_O: ";
	    std::copy( q_x_O.begin(), q_x_O.begin()+B_O.size(),
		       std::ostream_iterator<ET>( vout3.out()," "));
	    vout3.out() << std::endl;
	    if ( has_ineq) {
		vout3.out() << "   p_x_S: ";
		std::copy( q_x_S.begin(), q_x_S.begin()+B_S.size(),
			   std::ostream_iterator<ET>( vout3.out()," "));
		vout3.out() << std::endl;
	    }
	    vout3.out() << std::endl;
	}
    }
    
    // Idea here: At this point, the goal is to increase \mu_j until either we
    // become optimal (\mu_j=0), or one of the variables in x^*_\hat{B} drops
    // down to zero.
    //
    // Let us see first how this is done in the standard-form case (where
    // Sven's thesis applies).  Eq. (2.11) in Sven's thesis holds, and by
    // multlying it by $M_\hat{B}^{-1}$ we obtain an equation for \lambda and
    // x^*_\hat{B}.  The interesting equation (the one for x^*_\hat{B}) looks
    // more or less as follows:
    //
    //    x(mu_j)      = x(0) + mu_j      q_it                          (1)
    //
    // where q_it is the vector from (2.12).  In paritcular, for
    // mu_j=mu_j(t_1) (i.e., if we plug the value of mu_j at the beginning of
    // this ratio step 2 into (1)) we have
    //
    //    x(mu_j(t_1)) = x(0) + mu_j(t_1) q_it                          (2)
    //
    // where x(mu_j(t_1)) is the current solution of the solver at this point
    // (i.e., at the beginning of ratio step 2).
    //
    // By subtracting (2) from (1) we can thus eliminate the "unkown" x(0)
    // (which is cheaper than computing it):
    //
    //    x(mu_j) = x(mu_j(t_1)) + (mu_j-mu_j(t_1)) q_it
    //                             ----------------
    //                                  := delta
    //
    // In order to compute for each variable x_k in \hat{B} the value of mu_j
    // for which x_k(mu_j) = 0, we thus evaluate
    //
    //                x(mu_j(t_1))_k
    //    delta_k:= - --------------
    //                    q_it_k
    //
    // The first variable in \hat{B} that hits zero "in the future" is then
    // the one whose delta_k equals
    //
    //    delta_min:= min {delta_k | k in \hat{B} and (q_it)_k < 0 }
    //    
    // So in order to handle the standard-form case, we would compute this
    // minimum.  Once we have delta_min, we need to check whether we get
    // optimal BEFORE a variable drops to zero.  As delta = mu_j - mu_j(t_1),
    // the latter is precisely the case if delta_min >= -mu_j(t_1).
    //
    // (Note: please forget the crap identitiy between (2.11) and (2.12); the
    // notation is misleading.)
    //
    // Now to the nonstandard-form case.
    
    // fw: By definition delta_min >= 0, such that initializing
    // delta_min with -mu_j(t_1) has the desired effect that a basic variable
    // is leaving only if 0 <= delta_min < -mu_j(t_1).
    //
    // The only initialization of delta_min as fraction x_i/q_i that works is
    // x_i=mu_j(t_1); q_i=-1; (see below).
    //
    // Since mu_j(t_1) has been computed in ratio test step 1 we can
    // reuse it.
      
    x_i = mu;                                     // initialize minimum
    q_i = -et1;                                        // with -mu_j(t_1) 

    Value_iterator  x_it = x_B_O.begin();
    Value_iterator  q_it = q_x_O.begin();
    Index_iterator  i_it;
    for ( i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++x_it, ++q_it) {
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ( (*q_it < et0) && (
                             (( *x_it * q_i) < ( x_i * *q_it)) ||
                             ( (*i_it < index_leaving) && (index_leaving != art_s_i) && (( *x_it * q_i) == ( x_i * *q_it)) )
                             )
          ) {
        index_leaving = *i_it; x_i = *x_it; q_i = *q_it;
      }
    }
    x_it = x_B_S.begin();
    q_it = q_x_S.begin();
    for ( i_it = B_S.begin(); i_it != B_S.end(); ++i_it, ++x_it, ++q_it) {
      // BLAND rule: In case the ratios are the same, only update if the new index
      // is smaller. The special artificial variable is always made to leave first.
      if ( ( *q_it < et0) && (
                              (( *x_it * q_i) < ( x_i * *q_it)) ||
                              ( (*i_it < index_leaving) && (index_leaving != art_s_i) && (( *x_it * q_i) == ( x_i * *q_it)) )
                              )
          ){
        index_leaving = *i_it; x_i = *x_it; q_i = *q_it;
      }
    }

    CGAL_qpe_debug {
	if ( vout2.verbose()) {
	    for ( unsigned int k = 0; k < static_cast<unsigned int>(B_O.size()); ++k) {
		vout2.out() << "mu_j_O_" << k << ": - "
			    << x_B_O[ k] << '/' << q_x_O[ k]
			    << ( ( q_i < et0) && ( index_leaving == B_O[ k]) ? " *" : "")
			    << std::endl;
	    }
	    for ( unsigned int k = 0; k < static_cast<unsigned int>(B_S.size()); ++k) {
		vout2.out() << "mu_j_S_" << k << ": - "
			    << x_B_S[ k] << '/' << q_x_S[ k]
			    << ( ( q_i < et0) && ( index_leaving == B_S[ k]) ? " *" : "")
			    << std::endl;
	    }
	    vout2.out() << std::endl;
	}
    }
    if ( index_leaving < 0) {
      vout2 << "leaving variable: none" << std::endl;
    } else {
      vout1 << ", ";
      vout  << "leaving"; vout2 << " variable"; vout << ": ";
      vout  << index_leaving;
      if ( vout2.verbose()) {
	if ( index_leaving < qp_n) {
	  vout2.out() << " (= B_O[ " << in_B[ index_leaving] << "]: original)"
		      << std::endl;
	} else {
	  vout2.out() << " (= B_S[ " << in_B[ index_leaving] << "]: slack)"
		      << std::endl;
	}
      }
    }
    
}

// update (step 1)
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
update_1( )
{
  CGAL_qpe_debug {
    if ( vout2.verbose()) {
      vout2.out() << std::endl;
      if ( is_LP || is_phaseI) {
        vout2.out() << "------" << std::endl
        << "Update" << std::endl
        << "------" << std::endl;
      } else {
        vout2.out() << "Update (Step 1)" << std::endl
        << "---------------" << std::endl;
      }
    }
  }
  
  
  // update basis & basis inverse
  update_1( Is_linear());
  
  
  // check the updated vectors r_C and r_S_B
  CGAL_expensive_assertion(check_r_C(Is_nonnegative()));
  CGAL_expensive_assertion(check_r_S_B(Is_nonnegative()));
  
  // check the vectors r_B_O and w in phaseII for QPs
  CGAL_qpe_debug {
    if (is_phaseII && is_QP) {
      CGAL_expensive_assertion(check_r_B_O(Is_nonnegative()));
      CGAL_expensive_assertion(check_w(Is_nonnegative()));
    }
  }
  
  // compute current solution
  compute_solution(Is_nonnegative());
  
  
  // check feasibility 
  CGAL_qpe_debug {
    if (index_entering < 0 && !is_RTS_transition) // todo kf: is this too conservative?
      // Note: the above condition is necessary because of the
      // following.  In theory, it is true that the current
      // solution is at this point in the solver always
      // feasible. However, the solution has its x_j-entry equal to
      // the current t from the pricing, and is not zero (in the
      // standard-form case) or the current lower/upper bound
      // of the variable (in the non-standard-form case), resp., as
      // the routines is_solution_feasible() and
      // is_solution_feasible_for_auxiliary_problem() assume.
      if (is_phaseI) {
        CGAL_expensive_assertion(is_solution_feasible_for_auxiliary_problem());
      } else {
        CGAL_expensive_assertion(is_solution_feasible());
      }
      else
        vout2 << "(feasibility not checked in intermediate step)" << std::endl;
    CGAL_expensive_assertion(Is_nonnegative::value || r_C.size() == C.size());
  }
  
}

// update (step 2)
template < typename Q, typename ET, typename Tags >                                         // QP case
void
QP_solver<Q, ET, Tags>::
update_2( Tag_false)
{
    CGAL_qpe_debug {
	vout2 << std::endl
	      << "Update (Step 2)" << std::endl
	      << "---------------" << std::endl;
    }

    // leave variable from basis
    leave_variable();
    
    // check the updated vectors r_C, r_S_B, r_B_O and w
    CGAL_expensive_assertion(check_r_C(Is_nonnegative()));
    CGAL_expensive_assertion(check_r_S_B(Is_nonnegative()));
    CGAL_expensive_assertion(check_r_B_O(Is_nonnegative()));
    CGAL_expensive_assertion(check_w(Is_nonnegative()));

    // compute current solution
    compute_solution(Is_nonnegative());
}
 
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
expel_artificial_variables_from_basis( )
{
    int row_ind;
    ET r_A_Cj;
    
    CGAL_qpe_debug {
        vout2 << std::endl
	      << "---------------------------------------------" << std::endl
	      << "Expelling artificial variables from the basis" << std::endl
	      << "---------------------------------------------" << std::endl;
    }
    
    // try to pivot the artificials out of the basis
    // Note that we do not notify the pricing strategy about variables
    // leaving the basis, furthermore the pricing strategy does not
    // know about variables entering the basis.
    // The partial pricing strategies that keep the set of nonbasic vars
    // explicitly are synchronized during transition from phaseI to phaseII
    
    for (unsigned int i_ = static_cast<unsigned int>(qp_n + slack_A.size()); i_ < static_cast<unsigned int>(in_B.size()); ++i_) {
      if (is_basic(i_)) { 					// is basic
          row_ind = in_B[i_];
		  Values colvector(C.size());
		  Values einheit(C.size()+B_O.size(),et0);
	      Values dummy(B_O.size());
	      einheit[C.size() + row_ind] = et1;
	      lu_fact_.solve(einheit.begin(), einheit.begin() + C.size(), colvector.begin(), dummy.begin(), Is_linear::value, is_phaseI);
			
			  std::vector<int> nonzeros_of_colvector;  //for sparse representation for a faster inner_product
			  for (unsigned int k = 0; k < static_cast<unsigned int>(C.size()); ++k){
				  if (colvector[k]!=et0) nonzeros_of_colvector.push_back(k);
			  }
			
			  // determine first possible entering variable, if there is any
			  for (unsigned int j_ = 0; j_ < static_cast<unsigned int>(qp_n + slack_A.size()); ++j_) {
				  if (!is_basic(j_)) {  				// is nonbasic 
					  ratio_test_init__A_Cj( A_Cj.begin(), j_, no_ineq);
 					  r_A_Cj = et0;
					  r_A_Cj = lu_fact_.inner_product(colvector.begin(), A_Cj.begin(), C.size());  // bottle neck of this function so sparse inner_product and representation of the vectors
            
					  if (r_A_Cj != et0) {
						  ratio_test_1__q_x_O(Is_linear());
						  index_leaving = i_;
						  index_entering = j_;
						  update_1(Is_linear());
						  break;
					  } 
				  }
			  }
		  }
    }
    
        if ((art_basic != 0) && no_ineq) {
          // the vector in_C was not used in phase I, but now we remove redundant
          // constraints and switch to has_ineq treatment, hence we need it to
          // be correct at this stage
          for (int i=0; i<qp_m; ++i)
	          in_C.push_back(i);
        }
        diagnostics.redundant_equations = (art_basic != 0);

        // now reset the no_ineq and has_ineq flags to match the situation
        no_ineq = no_ineq && !diagnostics.redundant_equations;
        has_ineq = !no_ineq;
    
        // remove the remaining ones with their corresponding equality constraints
        // Note: the special artificial variable can always be driven out of the
        // basis
        for (unsigned int i_ = static_cast<unsigned int>(qp_n + slack_A.size()); i_ < static_cast<unsigned int>(in_B.size()); ++i_) {
          if (in_B[i_] >= 0) {
	          index_leaving = i_;
	          CGAL_qpe_debug {
	            vout2 << std::endl
		          << "~~> removing artificial variable " << index_leaving
		          << " and its equality constraint" << std::endl
		          << std::endl;
	          }
	          remove_artificial_variable_and_constraint();
	        }
        }
} // expel_artificial_variables_from_basis


// replace variable in basis
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
replace_variable( )
{
    CGAL_qpe_debug {
	vout2 <<   "<--> nonbasic (" << variable_type( index_entering) << ") variable " << index_entering
	      << " replaces basic (" << variable_type( index_leaving) << ") variable " << index_leaving
	      << std::endl << std::endl;
    }

    // replace variable
    replace_variable( no_ineq);

    // pivot step done
    index_leaving = index_entering = -1;
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
replace_variable_original_original( )
{
  // updates for the upper bounded case
  replace_variable_original_original_upd_r(Is_nonnegative());
  
  // store old column of basis matrix for update
  Values tmp(C.size());
  ratio_test_init__A_Cj(tmp.begin(), index_leaving, no_ineq); //update acj
  
  int  k = in_B[ index_leaving];
  
  // replace original variable [ in: index_entering | out: index_leaving]
  in_B  [ index_leaving] = -1;
  in_B  [ index_entering] = k;
  B_O[ k] = index_entering;
  
  minus_c_B[ k] = 
  ( is_phaseI ? 
   ( index_entering < qp_n ? et0 : -aux_c[index_entering-qp_n-slack_A.size()]) : -static_cast<ET>( *(qp_c+ index_entering)));
  
  if ( is_phaseI) {
    if ( index_entering >= qp_n) ++art_basic;
    if ( index_leaving >= qp_n) --art_basic;
  }
  
  // diagnostic output
  CGAL_qpe_debug {
    if ( vout2.verbose()) print_basis();
  }
  
  bool success(false);
  success = update_basis_matrix_O_O(k, tmp, Tag_true()); // linear case, type U5
  //lu_fact_.set_invalid();
  
  
  CGAL_qpe_debug {
      if (success) ++CGAL::QP_solver_debug::timer.counter_O_O;
      else ++CGAL::QP_solver_debug::timer.counter_O_O_fail;
  }
}

// update of the vector r for U_5 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
replace_variable_original_original_upd_r(Tag_true )
{
}

// update of the vector r for U_5 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded      
void  QP_solver<Q, ET, Tags>::
replace_variable_original_original_upd_r(Tag_false )
{
    ET      x_j, x_i;
    
    if (is_artificial(index_entering)) {
        if (!is_artificial(index_leaving)) {
            x_i = (ratio_test_bound_index == LOWER) ? *(qp_l+index_leaving) : *(qp_u+index_leaving);
            update_r_C_r_S_B__i(x_i);
            // update x_O_v_i
            x_O_v_i[index_leaving] = ratio_test_bound_index;
        }
    } else {
        x_j = nonbasic_original_variable_value(index_entering);
        if (is_artificial(index_leaving)) {
            update_r_C_r_S_B__j(x_j);
        } else {
            x_i = (ratio_test_bound_index == LOWER) ? *(qp_l+index_leaving) : *(qp_u+index_leaving);
            update_r_C_r_S_B__j_i(x_j, x_i);
            // update x_O_v_i
            x_O_v_i[index_leaving] = ratio_test_bound_index;
        }
        // update x_O_v_i
        x_O_v_i[index_entering] = BASIC;
    }
}


template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
replace_variable_slack_slack( )
{
  
  // updates for the upper bounded case
  replace_variable_slack_slack_upd_r(Is_nonnegative()); 
  
  // store old row of basis matrix for update
  Values tmp(B_O.size());
  init__A_Ri(tmp.begin(), slack_A[ index_entering-qp_n].first, no_ineq);
  
  
  int  k = in_B[ index_leaving];
  
  // replace slack variable [ in: index_entering | out: index_leaving]
  in_B  [ index_leaving] = -1;
  in_B  [ index_entering] = k;
  B_S[ k] = index_entering;
  S_B[ k] = slack_A[ index_entering-qp_n].first;
  
  // replace inequality constraint [ in: index_leaving | out: index_entering ]
  int old_row = S_B[ k];
  
  int new_row = slack_A[ index_leaving-qp_n].first;
  k = in_C[ old_row];
  
  in_C[ old_row] = -1;
  in_C[ new_row] = k;
  C[ k      ] = new_row;
  
  b_C[ k] = static_cast<ET>( *(qp_b+ new_row));
  
  // diagnostic output
  CGAL_qpe_debug {
    if ( vout2.verbose()) print_basis();
  }


  bool success(false);
  success = update_basis_matrix_LP_S_S(k, tmp); // type U6
  //lu_fact_.set_invalid(); // bail out variant (just refactor)
  
  
  CGAL_qpe_debug {
      if (success) ++CGAL::QP_solver_debug::timer.counter_S_S;
      else ++CGAL::QP_solver_debug::timer.counter_S_S_fail;
  }
  
}

// update of the vector r for U_6 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
replace_variable_slack_slack_upd_r(Tag_true )
{
}

// update of the vector r for U_6 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded      
void  QP_solver<Q, ET, Tags>::
replace_variable_slack_slack_upd_r(Tag_false )
{
    int     sigma_j = slack_A[ index_entering-qp_n].first;
    
    // swap r_gamma_C(sigma_j) in r_C with r_gamma_S_B(sigma_i) in r_S_B
    std::swap(r_C[in_C[sigma_j]], r_S_B[in_B[index_leaving]]); 
}


template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
replace_variable_slack_original( )
{
  
  Values old_row(B_O.size());
  Values old_col(C.size());
  Values last_row(B_O.size());
  Values last_col(C.size());

  // get old row and column
  int  old_row_index = slack_A[ index_entering-qp_n].first;
  
  init__A_Ri(old_row.begin(), old_row_index, no_ineq);
  ratio_test_init__A_Cj(old_col.begin(), index_leaving, no_ineq);
  init__A_Ri(last_row.begin(), C.back(), no_ineq);
  ratio_test_init__A_Cj(last_col.begin(), B_O.back(), no_ineq);

  // updates for the upper bounded case
  replace_variable_slack_original_upd_r(Is_nonnegative()); 
  
  int  k = in_B[ index_leaving];
  
  // leave original variable [ out: index_leaving]
  in_B  [ B_O.back()] = k;
  B_O[ k] = B_O.back();
  in_B  [ index_leaving] = -1;
  B_O.pop_back();
  
  minus_c_B[ k] = minus_c_B[ B_O.size()];
  
  if ( is_phaseI && ( index_leaving >= qp_n)) --art_basic;
  
  // enter slack variable [ in: index_entering ]
  
  in_B  [ index_entering] = static_cast<int>(B_S.size());
  B_S.push_back( index_entering);
  S_B.push_back( old_row_index);
  
  // leave inequality constraint [ out: index_entering ]
  int  l = in_C[ old_row_index];
  b_C[ l       ] = b_C[ C.size()-1];
  C[ l       ] = C.back();
  in_C[ C.back()] = l;
  in_C[ old_row_index ] = -1;
  C.pop_back();
  // diagnostic output
  CGAL_qpe_debug {
    if ( vout2.verbose()) print_basis();
  }
  
  bool success(false);
  success = update_basis_matrix_LP_S_O(l, old_row, k, old_col, last_row, last_col); // type U8
  //lu_fact_.set_invalid(); // bail out variant (just refactor)
  
  CGAL_qpe_debug {
      if (success) ++CGAL::QP_solver_debug::timer.counter_S_O;
      else ++CGAL::QP_solver_debug::timer.counter_S_O_fail;
  }
}

// update of the vector r for U_8 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
replace_variable_slack_original_upd_r(Tag_true )
{
}

// update of the vector r for U_8 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded      
void  QP_solver<Q, ET, Tags>::
replace_variable_slack_original_upd_r(Tag_false )
{
  if (!is_artificial(index_leaving)) {
    ET  x_i = (ratio_test_bound_index == LOWER) ? *(qp_l+index_leaving) : *(qp_u+index_leaving);
    update_r_C_r_S_B__i(x_i);
  }
  
  int     sigma_j = slack_A[ index_entering-qp_n].first;
  
  // append r_gamma_C(sigma_j) from r_C to r_S_B:
  r_S_B.push_back(r_C[in_C[sigma_j]]);
  
  // remove r_gamma_C(sigma_j) from r_C:
  r_C[in_C[sigma_j]] = r_C.back();
  r_C.pop_back();
  
  // update x_O_v_i
  if (!is_artificial(index_leaving)) // original and not artificial?
    x_O_v_i[index_leaving] = ratio_test_bound_index;
}


template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
replace_variable_original_slack( )
{
  // updates for the upper bounded case
  replace_variable_original_slack_upd_r(Is_nonnegative());
  
  int  k = in_B[ index_leaving];
  
  // enter original variable [ in: index_entering ]
  
  minus_c_B[ B_O.size()]
  = ( is_phaseI ? 
     ( index_entering < qp_n ? et0 : -aux_c[index_entering-qp_n-slack_A.size()]) 
     : -static_cast<ET>( *(qp_c+ index_entering)));
  
  
  in_B  [ index_entering] = static_cast<int>(B_O.size());
  B_O.push_back( index_entering);
  
  if ( is_phaseI && ( index_entering >= qp_n)) ++art_basic;
  
  // leave slack variable [ out: index_leaving]
  B_S[ k         ] = B_S.back();
  S_B[ k         ] = S_B.back();
  in_B  [ B_S.back()] = k;
  in_B  [ index_leaving] = -1; 
  B_S.pop_back();
  S_B.pop_back();
  
  // enter inequality constraint [ in: index_leaving ]
  int new_row = slack_A[ index_leaving-qp_n].first;
  
  b_C[ C.size()] = static_cast<ET>( *(qp_b+ new_row));
  in_C[ new_row ] = static_cast<int>(C.size());
  C.push_back( new_row);
  
  // diagnostic output
  CGAL_qpe_debug {
    if ( vout2.verbose()) print_basis();
  }

  CGAL_qpe_assertion(B_O.size() == C.size());
  
  bool success(false);
  success = update_basis_matrix_LP_O_S(); // type U7
  //lu_fact_.set_invalid(); // bail out variant (just refactor)
  
  CGAL_qpe_debug {
      if (success) ++CGAL::QP_solver_debug::timer.counter_O_S;
      else ++CGAL::QP_solver_debug::timer.counter_O_S_fail;
  }
}

// update of the vector r for U_7 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
replace_variable_original_slack_upd_r(Tag_true )
{
}

// update of the vector r for U_7 with upper bounding, note that we 
// need the headings C, and S_{B} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded      
void  QP_solver<Q, ET, Tags>::
replace_variable_original_slack_upd_r(Tag_false )
{
    if (!is_artificial(index_entering)) {
        ET x_j = nonbasic_original_variable_value(index_entering);
        update_r_C_r_S_B__j(x_j);
    }
    
    // append r_gamma_S_B(sigma_i) from r_S_B to r_C
    r_C.push_back(r_S_B[in_B[index_leaving]]);
    
    // remove r_gamma_S_B(sigma_i) from r_S_B
    r_S_B[in_B[index_leaving]] = r_S_B.back();
    r_S_B.pop_back();
    
    // update x_O_v_i
    if (!is_artificial(index_entering)) {
        x_O_v_i[index_entering] = BASIC;
    }
}


template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
remove_artificial_variable_and_constraint( )
{

  Values old_row(B_O.size());
  Values old_col(C.size());
  Values last_row(B_O.size());
  Values last_col(C.size());

  // get old row and column
  int old_row_index = art_A[index_leaving - qp_n - slack_A.size()].first;
  
  init__A_Ri(old_row.begin(), old_row_index, no_ineq);
  ratio_test_init__A_Cj(old_col.begin(), index_leaving, no_ineq);
  init__A_Ri(last_row.begin(), C.back(), no_ineq);
  ratio_test_init__A_Cj(last_col.begin(), B_O.back(), no_ineq);

  // updates for the upper bounded case
  remove_artificial_variable_and_constraint_upd_r(Is_nonnegative());
  
  int  k = in_B[ index_leaving];
  
  // leave artificial (original) variable [ out: index_leaving]
  in_B  [ B_O.back()] = k;
  B_O[ k] = B_O.back();
  in_B  [ index_leaving] = -1;
  B_O.pop_back();
  
  minus_c_B[ k] = minus_c_B[ B_O.size()];
  
  if ( is_phaseI && ( index_leaving >= qp_n)) --art_basic;
  
  
  
  // leave its equality constraint 
  int  l = in_C[ old_row_index];
  b_C[ l       ] = b_C[ C.size()-1];
  C[ l       ] = C.back();
  in_C[ C.back()] = l;
  in_C[ old_row_index ] = -1;
  C.pop_back();
  // diagnostic output
  CGAL_qpe_debug {
    if ( vout2.verbose()) print_basis();
  }
  
  bool success(false);
  success = update_basis_matrix_LP_S_O(l, old_row, k, old_col, last_row, last_col); // type U8
  //lu_fact_.set_invalid(); // bail out variant (just refactor)
  
  CGAL_qpe_debug {
    if (success) ++CGAL::QP_solver_debug::timer.counter_S_O;
    else ++CGAL::QP_solver_debug::timer.counter_S_O_fail;
  }
}

// update of the vector r with upper bounding for the removal of an
// artificial variable with its equality constraint, note that we 
// need the headings C before it is updated
template < typename Q, typename ET, typename Tags >                                 // Standard form
void  QP_solver<Q, ET, Tags>::
remove_artificial_variable_and_constraint_upd_r(Tag_true )
{
}

// update of the vector r with upper bounding for the removal of an
// artificial variable with its equality constraint, note that we 
// need the headings C before it is updated
template < typename Q, typename ET, typename Tags >                                 // Upper bounded
void  QP_solver<Q, ET, Tags>::
remove_artificial_variable_and_constraint_upd_r(Tag_false )
{
    int sigma_i = art_A[index_leaving - qp_n - slack_A.size()].first;
    
    // remove r_gamma_C(sigma_i) from r_C
    r_C[in_C[sigma_i]] = r_C.back();
    r_C.pop_back();
}

// update that occurs only with upper bounding in ratio test step 1
template < typename Q, typename ET, typename Tags >            
void  QP_solver<Q, ET, Tags>::
enter_and_leave_variable( )
{
    
    CGAL_qpe_assertion((index_leaving == index_entering) && (index_leaving >= 0));
    
    CGAL_qpe_debug {
	vout2 <<   "<--> nonbasic (" << variable_type( index_entering) << ") variable " << index_entering
	      << " enters and leaves basis" << std::endl << std::endl;
    }

    
    ET diff;
    ET x_j = nonbasic_original_variable_value(index_entering);
    
    if (ratio_test_bound_index == LOWER) {
        diff = x_j - static_cast<ET>(*(qp_l+index_entering));
    } else {
        diff = x_j - static_cast<ET>(*(qp_u+index_entering));
    }
    
    if (is_phaseI) {
        update_r_C_r_S_B__j(diff);
    } else {
        update_w_r_B_O__j(diff);
        update_r_C_r_S_B__j(diff);
    }
    
    x_O_v_i[index_entering] = ratio_test_bound_index;
    
    // notify pricing strategy (it has called enter_basis on i before)
    strategyP->leaving_basis (index_leaving);

    // variable entered and left basis
    index_leaving = -1; index_entering = -1;
}



// enter variable into basis
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
enter_variable( )
{
  CGAL_qpe_assertion (is_phaseII);
  CGAL_qpe_debug {
    vout2 << "--> nonbasic (" << variable_type( index_entering) << ") variable "
    << index_entering << " enters basis" << std::endl << std::endl;
  }
  
  // update basis & basis inverse:
  if (no_ineq || (index_entering < qp_n)) {              // original variable
    
    // updates for the upper bounded case:
    enter_variable_original_upd_w_r(Is_nonnegative());
    
    // enter original variable [ in: index_entering ]:
    if (minus_c_B.size() <= B_O.size()) { // Note: minus_c_B and the
      // containers resized in this
      // if-block are only enlarged
      // and never made smaller
      // (while B_O always has the
      // correct size). We check here
      // whether we need to enlarge
      // them.
      CGAL_qpe_assertion(minus_c_B.size() == B_O.size());
      minus_c_B.push_back(et0);
      q_x_O.push_back(et0);
      tmp_x  .push_back(et0);
      tmp_x_2.push_back(et0);
      two_D_Bj.push_back(et0);
      x_B_O.push_back(et0);
    }
    minus_c_B[B_O.size()] = -static_cast<ET>(*(qp_c+ index_entering)); // Note: B_O has always the
    // correct size.
    
    in_B[index_entering] = static_cast<int>(B_O.size());
    B_O.push_back(index_entering);
    
    // diagnostic output
    CGAL_qpe_debug {
      if (vout2.verbose())
        print_basis();
    }
    
    // update basis inverse
    // note: (-1)\hat{\nu} is stored instead of \hat{\nu}

    
    bool success(false);
    success = update_basis_matrix_QP_O_in(); // type U1
    //lu_fact_.set_invalid(); // bail out variant (just refactor)
    
    CGAL_qpe_debug {
      if (success) ++CGAL::QP_solver_debug::timer.counter_Oin;
      else ++CGAL::QP_solver_debug::timer.counter_Oin_fail;
    }
  
  } else {                                  // slack variable
    
    // store old column of basis matrix for update
    int old_row = slack_A[ index_entering-qp_n].first;
    int k = in_C[old_row];
    
    int last_row = C.back(); // equality constraint
    int k2 = in_C[last_row];
    Values tmp(C.size()+B_O.size(), et0);
    Values last_constraint(C.size()+B_O.size(), et0);
    init__A_Ri(tmp.begin()+C.size(), slack_A[ index_entering-qp_n].first, no_ineq); //update ari
    init__A_Ri(last_constraint.begin()+C.size(), last_row, no_ineq); //update ari
    
    // updates for the upper bounded case:
    enter_variable_slack_upd_w_r(Is_nonnegative());
    
    // enter slack variable [ in: index_entering ]:
    in_B  [ index_entering] = static_cast<int>(B_S.size());
    B_S.push_back( index_entering);
    S_B.push_back( slack_A[ index_entering-qp_n].first);    
    
    
    // reflect change of active constraints heading C in b_C:
    b_C[ k] = b_C[C.size()-1];
    
    C[ k] = C.back();
    in_C[ C.back()      ] = k;
    in_C[ old_row       ] = -1;
    C.pop_back();
    
    // diagnostic output:
    CGAL_qpe_debug {
      if (vout2.verbose())
        print_basis();
    }
    
    
    bool success(false);
    success = update_basis_matrix_QP_S_in(k, tmp, k2, last_constraint); // type U3
    //lu_fact_.set_invalid(); // bail out variant (just refactor)
    
    CGAL_qpe_debug {
      if (success) ++CGAL::QP_solver_debug::timer.counter_Sin;
      else ++CGAL::QP_solver_debug::timer.counter_Sin_fail;
    }
    
  }
  
  // variable entered:
  index_entering -= static_cast<int>(in_B.size());
}
  
// update of the vectors w and r for U_1 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
enter_variable_original_upd_w_r(Tag_true )
{
}

// update of the vectors w and r for U_1 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded     
void  QP_solver<Q, ET, Tags>::
enter_variable_original_upd_w_r(Tag_false )
{

    ET x_j = nonbasic_original_variable_value(index_entering);

    // Note: w needs to be updated before r_C, r_S_B
    update_w_r_B_O__j(x_j);
    update_r_C_r_S_B__j(x_j);
    
    // append w_j to r_B_O
    if (!Is_linear::value) // (kf.)
      r_B_O.push_back(w[index_entering]);
    
    // update x_O_v_i
    x_O_v_i[index_entering] = BASIC;
}

// update of the vectors w and r for U_3 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
enter_variable_slack_upd_w_r(Tag_true )
{
}

// update of the vectors w and r for U_3 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded     
void  QP_solver<Q, ET, Tags>::
enter_variable_slack_upd_w_r(Tag_false )
{
    
    int     sigma_j = slack_A[ index_entering-qp_n].first;       
    
    // append r_gamma_C(sigma_j) to r_S_B
    r_S_B.push_back(r_C[in_C[sigma_j]]);
    
    // remove r_gamma_C(sigma_j) from r_C   
    r_C[in_C[sigma_j]] = r_C.back();
    r_C.pop_back();
}

// leave variable from basis
template < typename Q, typename ET, typename Tags >
void
QP_solver<Q, ET, Tags>::
leave_variable( )
{
  CGAL_qpe_debug {
    vout2 << "<-- basic (" << variable_type( index_leaving) << ") variable "
    << index_leaving << " leaves basis" << std::endl << std::endl;
  }
  
  
  // update basis & basis inverse
  int  k = in_B[ index_leaving];
  if ( no_ineq || ( index_leaving < qp_n)) {                      // original variable
  
    // store old column of basis matrix for update
    Values tmp(C.size()+B_O.size(), et0);
    Values last_col(C.size()+B_O.size(), et0);
    ratio_test_init__A_Cj(tmp.begin(), index_leaving, no_ineq); //update acj
    ratio_test_init__2_D_Bj(tmp.begin()+C.size(), index_leaving, Tag_false() /*quadratic*/);
    ratio_test_init__A_Cj(last_col.begin(), B_O.back(), no_ineq); //update acj
    ratio_test_init__2_D_Bj(last_col.begin()+C.size(), B_O.back(), Tag_false() /*quadratic*/);
    
    
    // updates for the upper bounded case
    leave_variable_original_upd_w_r(Is_nonnegative());
    
    // leave original variable [ out: index_leaving]
    in_B  [ B_O.back()] = k;
    in_B  [ index_leaving] = -1; 
    //in_B  [ B_O.back()] = k;
    B_O[ k] = B_O.back(); B_O.pop_back();
    
    minus_c_B [ k] = minus_c_B [ B_O.size()];
    two_D_Bj[ k] =   two_D_Bj[ B_O.size()];
    
    
    // diagnostic output
    CGAL_qpe_debug {
      if ( vout2.verbose()) print_basis();
    }
    
    bool success(false);
    success = update_basis_matrix_QP_O_out(k, tmp, last_col); // type U2
    //lu_fact_.set_invalid(); // bail out variant (just refactor)
    
    CGAL_qpe_debug {
      if (success) ++CGAL::QP_solver_debug::timer.counter_Oout;
      else ++CGAL::QP_solver_debug::timer.counter_Oout_fail;
    }
    
  } else {                                            // slack variable
  
    // updates for the upper bounded case
    leave_variable_slack_upd_w_r(Is_nonnegative());
    
    // leave slack variable [ out: index_leaving]
    in_B  [ B_S.back()] = k;      // former last var moves to position k
    in_B  [ index_leaving] = -1;     // index_leaving gets deleted
    B_S[ k] = B_S.back(); B_S.pop_back();
    S_B[ k] = S_B.back(); S_B.pop_back();
    
    // enter inequality constraint [ in: index_leaving ]
    int new_row = slack_A[ index_leaving-qp_n].first;
    
    ET tmp = et0;
    if (index_entering < qp_n) {
      A_sparse_column_iterator it = (*(qp_A_sparse+index_entering)).begin();
      A_sparse_column_iterator it_end = (*(qp_A_sparse+index_entering)).end();
      while (it != it_end) {
        if (it->first >= new_row) {
          if (it->first == new_row) tmp = it->second;
          break;
        }
        ++it;
      }    
    }
    A_Cj[ C.size()] = tmp;
    
    b_C[ C.size()] = static_cast<ET>( *(qp_b+ new_row));
    in_C[ new_row ] = static_cast<int>(C.size());
    C.push_back( new_row);
    
    // diagnostic output
    CGAL_qpe_debug {
      if ( vout2.verbose()) print_basis();
    }
    
    bool success(false);
    success = update_basis_matrix_QP_S_out(); // type U4
    //lu_fact_.set_invalid(); // bail out variant (just refactor)
    
    CGAL_qpe_debug {
      if (success) ++CGAL::QP_solver_debug::timer.counter_Sout;
      else ++CGAL::QP_solver_debug::timer.counter_Sout_fail;
    }
    
  }
  
  // notify pricing strategy
  strategyP->leaving_basis( index_leaving);
  
  // variable left
  index_leaving = -1;
}

// update of the vectors w and r for U_2 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
leave_variable_original_upd_w_r(Tag_true )
{
}

// update of the vectors w and r for U_2 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded      
void  QP_solver<Q, ET, Tags>::
leave_variable_original_upd_w_r(Tag_false )
{

    ET      x_i = (ratio_test_bound_index == LOWER) ? *(qp_l+index_leaving) : *(qp_u+index_leaving);
    
    // Note: w needs to be updated before r_C, r_S_B
    update_w_r_B_O__i(x_i);
    update_r_C_r_S_B__i(x_i);    
    
    // remove r_beta_O(i) from r_B_O
    if (!Is_linear::value) { // (kf.)
      r_B_O[in_B[index_leaving]] = r_B_O.back();
      r_B_O.pop_back();
    }
    
    // update x_O_v_i
    x_O_v_i[index_leaving] = ratio_test_bound_index;
}

// update of the vectors w and r for U_4 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
leave_variable_slack_upd_w_r(Tag_true )
{
}

// update of the vectors w and r for U_4 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Upper bounded      
void  QP_solver<Q, ET, Tags>::
leave_variable_slack_upd_w_r(Tag_false )
{
    
    // append r_gamma_S_B(sigma_i) to r_C
    r_C.push_back(r_S_B[in_B[index_leaving]]);
    
    // remove r_gamma_S_B(sigma_i) from r_S_B
    r_S_B[in_B[index_leaving]] = r_S_B.back();
    r_S_B.pop_back();
}


// replace variable in basis QP-case, transition to Ratio Test Step 2
template < typename Q, typename ET, typename Tags >
void QP_solver<Q, ET, Tags>::
z_replace_variable( )
{
    CGAL_qpe_debug {
	vout2 <<   "<--> nonbasic (" << variable_type( index_entering) << ") variable " << index_entering
	      << " z_replaces basic (" << variable_type( index_leaving) << ") variable " << index_leaving
	      << std::endl << std::endl;
    }

    // replace variable
    z_replace_variable( no_ineq);

    // pivot step not yet completely done
    index_leaving = -1;
    index_entering -= static_cast<int>(in_B.size());
    is_RTS_transition = true;
}


template < typename Q, typename ET, typename Tags >  inline                           // no inequalities
void QP_solver<Q, ET, Tags>::
z_replace_variable( Tag_true)
{
  
  z_replace_variable_original_by_original();
  strategyP->leaving_basis(index_leaving);
  
}


template < typename Q, typename ET, typename Tags >  inline                          // has inequalities
void QP_solver<Q, ET, Tags>::
z_replace_variable( Tag_false)
{
  // determine type of variables
  bool  enter_original = ( (index_entering < qp_n) || (index_entering >= static_cast<int>( qp_n+slack_A.size())));
  bool  leave_original = ( (index_leaving < qp_n) || (index_leaving >= static_cast<int>( qp_n+slack_A.size())));
  
  // update basis and basis inverse
  if ( leave_original) {
    if ( enter_original) {               
      z_replace_variable_original_by_original();
    } else {                             
      z_replace_variable_original_by_slack();
    }
  } else {
    if ( enter_original) {
      z_replace_variable_slack_by_original();
    } else {
      z_replace_variable_slack_by_slack();
    }
  }
  strategyP->leaving_basis( index_leaving);
}


// replacement with precond det(M_{B \setminus \{i\}})=0
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
z_replace_variable_original_by_original( )
{
  // updates for the upper bounded case
  z_replace_variable_original_by_original_upd_w_r(Is_nonnegative());
  
  // store old column of basis matrix for update
  // TODO: replace variable names such as tmp
  Values tmp(C.size()+B_O.size(), et0);
  ratio_test_init__A_Cj(tmp.begin(), index_leaving, no_ineq); //update acj
  ratio_test_init__2_D_Bj(tmp.begin()+C.size(), index_leaving, Tag_false() /*quadratic*/);
  
  int  k = in_B[ index_leaving];
  
  // replace original variable [ in: index_entering | out: index_leaving]
  in_B  [ index_leaving] = -1;
  in_B  [ index_entering] = k;
  B_O[ k] = index_entering;
  
  minus_c_B[ k] = -static_cast<ET>( *(qp_c+ index_entering));
  
  // diagnostic output
  CGAL_qpe_debug {
    if ( vout2.verbose()) print_basis();
  }
  
  bool success(false);
  success = update_basis_matrix_O_O(k, tmp, Tag_false()); // quadratic case, type UZ1
  //lu_fact_.set_invalid();
  
  //CGAL_qpe_debug {
      if (success) ++CGAL::QP_solver_debug::timer.counter_UZ1;
      else ++CGAL::QP_solver_debug::timer.counter_UZ1_fail;
  //}
  
}

// update of the vectors w and r for U_Z_1 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form      
void  QP_solver<Q, ET, Tags>::
z_replace_variable_original_by_original_upd_w_r(Tag_true )
{
}

// update of the vectors w and r for U_Z_1 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                           
// Upper bounded      
void  QP_solver<Q, ET, Tags>::
z_replace_variable_original_by_original_upd_w_r(Tag_false )
{

    ET      x_j = nonbasic_original_variable_value(index_entering);
    ET      x_i = (ratio_test_bound_index == LOWER) ? *(qp_l+index_leaving) : *(qp_u+index_leaving);
    
    // Note: w needs to be updated before r_C, r_S_B
    update_w_r_B_O__j_i(x_j, x_i);
    update_r_C_r_S_B__j_i(x_j, x_i);
    
    // replace r_beta_O(i) with w_j    
    r_B_O[in_B[index_leaving]] = w[index_entering];
    
    // update x_O_v_i
    x_O_v_i[index_entering] = BASIC;
    x_O_v_i[index_leaving] = ratio_test_bound_index;    
}


// replacement with precond det(M_{B \setminus \{i\}})=0
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
z_replace_variable_original_by_slack( )
{

    Values old_constraint_row(C.size()+B_O.size(), et0);
    Values old_variable_row(C.size()+B_O.size(), et0);
    Values last_constraint_row(C.size()+B_O.size(), et0);
    Values last_variable_row(C.size()+B_O.size(), et0);

    // get old row and column
    int  k = in_B[ index_leaving];
    
    int  old_row_index = slack_A[ index_entering-qp_n].first;
    int last_row_index = C.back();  
    int k1 = in_C[old_row_index];
    int k2 = in_C[last_row_index];
    init__A_Ri(old_constraint_row.begin()+C.size(), old_row_index, no_ineq); //update ari
    init__A_Ri(last_constraint_row.begin()+C.size(), last_row_index, no_ineq); //update ari
    ratio_test_init__A_Cj(old_variable_row.begin(), index_leaving, no_ineq); //update acj
    ratio_test_init__2_D_Bj(old_variable_row.begin()+C.size(), index_leaving, Tag_false() /*quadratic*/);
    ratio_test_init__A_Cj(last_variable_row.begin(), B_O.back(), no_ineq); //update acj
    ratio_test_init__2_D_Bj(last_variable_row.begin()+C.size(), B_O.back(), Tag_false() /*quadratic*/);
    

    // updates for the upper bounded case
    z_replace_variable_original_by_slack_upd_w_r(Is_nonnegative());
    

    // leave original variable [ out: index_leaving]
    in_B  [ B_O.back()] = k;
       B_O[ k] = B_O.back();
       in_B  [ index_leaving] = -1;
       B_O.pop_back();

    minus_c_B[ k] = minus_c_B[ B_O.size()];

    // enter slack variable [ in: index_entering ]
    int  old_row = slack_A[ index_entering-qp_n].first;
    in_B  [ index_entering] = static_cast<int>(B_S.size());
       B_S.push_back( index_entering);
       S_B.push_back( old_row);

    // leave inequality constraint [ out: index_entering ]
    int  l = in_C[ old_row];
     b_C[ l       ] = b_C[ C.size()-1];
       C[ l       ] = C.back();
    in_C[ C.back()] = l;
    in_C[ old_row ] = -1;
       C.pop_back();
    
    // diagnostic output
    CGAL_qpe_debug {
	if ( vout2.verbose()) print_basis();
    }
     
  
  bool success(false);
  success = update_basis_matrix_QP_S_O(k1, old_constraint_row, k, old_variable_row, k2, last_constraint_row, last_variable_row); // type UZ2
  //lu_fact_.set_invalid(); // bail out variant (just refactor)
  
  CGAL_qpe_debug {
    if (success) ++CGAL::QP_solver_debug::timer.counter_UZ2;
    else ++CGAL::QP_solver_debug::timer.counter_UZ2_fail;
  }

}

// update of the vectors w and r for U_Z_2 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                         // Standard form
void  QP_solver<Q, ET, Tags>::
z_replace_variable_original_by_slack_upd_w_r(Tag_true )
{
}


// update of the vectors w and r for U_Z_2 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                         // Upper bounded
void  QP_solver<Q, ET, Tags>::
z_replace_variable_original_by_slack_upd_w_r(Tag_false )
{

    ET      x_i = (ratio_test_bound_index == LOWER) ? *(qp_l+index_leaving) : *(qp_u+index_leaving);
    
    // Note: w needs to be updated before r_C, r_S_B
    update_w_r_B_O__i(x_i);
    update_r_C_r_S_B__i(x_i);
    
    int     sigma_j = slack_A[ index_entering-qp_n].first;
    
    // append r_gamma_C(sigma_j) to r_S_B
    r_S_B.push_back(r_C[in_C[sigma_j]]);
    
    // remove r_gamma_C(sigma_j) from r_C
    r_C[in_C[sigma_j]] = r_C.back();
    r_C.pop_back();
    
    // remove r_beta_O(i) from r_B_O    
    if (!Is_linear::value) { // (kf.)
      r_B_O[in_B[index_leaving]] = r_B_O.back();
      r_B_O.pop_back();
    }
    
    // update x_O_v_i
    x_O_v_i[index_leaving] = ratio_test_bound_index;
}


// replacement with precond det(M_{B \setminus \{i\}})=0
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
z_replace_variable_slack_by_original( )
{
  // updates for the upper bounded case
  z_replace_variable_slack_by_original_upd_w_r(Is_nonnegative());
  
  int  k = in_B[ index_leaving];
  if (minus_c_B.size() <= B_O.size()) {  // Note: minus_c_B and the
    // containers resized in this
    // if-block are only enlarged
    // and never made smaller
    // (while B_O always has the
    // correct size). We check here
    // whether we need to enlarge
    // them.
    CGAL_qpe_assertion(minus_c_B.size() == B_O.size());
    minus_c_B.push_back(et0);
    q_x_O.push_back(et0);
    tmp_x  .push_back(et0);
    tmp_x_2.push_back(et0);
    two_D_Bj.push_back(et0);
    x_B_O.push_back(et0);
  }
  
  // enter original variable [ in: index_entering ]
  
  minus_c_B[ B_O.size()] = -static_cast<ET>( *(qp_c+ index_entering));
  
  
  in_B  [ index_entering] = static_cast<int>(B_O.size());
  B_O.push_back( index_entering);
  
  // leave slack variable [ out: index_leaving]
  B_S[ k         ] = B_S.back();
  S_B[ k         ] = S_B.back();
  in_B  [ B_S.back()] = k;
  in_B  [ index_leaving] = -1; 
  B_S.pop_back();
  S_B.pop_back();
  
  // enter inequality constraint [ in: index_leaving ]
  int new_row = slack_A[ index_leaving-qp_n].first;
  
  b_C[ C.size()] = static_cast<ET>( *(qp_b+ new_row));
  in_C[ new_row ] = static_cast<int>(C.size());
  C.push_back( new_row);
  
  // diagnostic output
  CGAL_qpe_debug {
    if ( vout2.verbose()) print_basis();
  }
  
  bool success(false);
  success = update_basis_matrix_QP_O_S(); // type UZ3
  //lu_fact_.set_invalid(); // bail out variant (just refactor)
  
  CGAL_qpe_debug {
    if (success) ++CGAL::QP_solver_debug::timer.counter_UZ3;
    else ++CGAL::QP_solver_debug::timer.counter_UZ3_fail;
  }
  
}

// update of the vectors w and r for U_Z_3 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                            // Standard form
void  QP_solver<Q, ET, Tags>::
z_replace_variable_slack_by_original_upd_w_r(Tag_true )
{
}

// update of the vectors w and r for U_Z_3 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                             // Upper bounded
void  QP_solver<Q, ET, Tags>::
z_replace_variable_slack_by_original_upd_w_r(Tag_false )
{

    ET      x_j = nonbasic_original_variable_value(index_entering);

    // Note: w needs to be updated before r_C, r_S_B    
    update_w_r_B_O__j(x_j);
    update_r_C_r_S_B__j(x_j);
        
    // append r_gamma_S_B(sigma_i) to r_C
    r_C.push_back(r_S_B[in_B[index_leaving]]);
    
    // remove r_gamma_S_B(sigma_i) from r_S_B
    r_S_B[in_B[index_leaving]] = r_S_B.back();
    r_S_B.pop_back();
    
    // append w_j to r_B_O    
    if (!Is_linear::value) // (kf.)
      r_B_O.push_back(w[index_entering]);
    
    // update x_O_v_i
    x_O_v_i[index_entering] = BASIC;
}


// replacement with precond det(M_{B \setminus \{i\}})=0
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
z_replace_variable_slack_by_slack( )
{

  // get old column
  Values old_column(B_O.size());
  init__A_Ri(old_column.begin(), slack_A[ index_entering-qp_n].first, no_ineq); //update ari
  
  // updates for the upper bounded case
  z_replace_variable_slack_by_slack_upd_w_r(Is_nonnegative());
  
  int  k = in_B[ index_leaving];
  
  // replace slack variable [ in: index_entering | out: index_leaving]
  in_B  [ index_leaving] = -1;
  in_B  [ index_entering] = k;
  B_S[ k] = index_entering;
  S_B[ k] = slack_A[ index_entering-qp_n].first;
  
  // replace inequality constraint [ in: index_leaving | out: index_entering ]
  int old_row = S_B[ k];
  int new_row = slack_A[ index_leaving-qp_n].first;
  k = in_C[ old_row];
  
  in_C[ old_row] = -1;
  in_C[ new_row] = k;
  C[ k      ] = new_row;
  
  b_C[ k] = static_cast<ET>( *(qp_b+ new_row));
  
  // diagnostic output
  CGAL_qpe_debug {
    if ( vout2.verbose()) print_basis();
  }
  
  
  bool success(false);
  success = update_basis_matrix_QP_S_S(k, old_column); // type UZ4
  //lu_fact_.set_invalid(); // bail out variant (just refactor)
  
  CGAL_qpe_debug {
    if (success) ++CGAL::QP_solver_debug::timer.counter_UZ4;
    else ++CGAL::QP_solver_debug::timer.counter_UZ4_fail;
  }
  
}

// update of the vectors w and r for U_Z_4 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                             // Standard form
void  QP_solver<Q, ET, Tags>::
z_replace_variable_slack_by_slack_upd_w_r(Tag_true )
{
}


// update of the vectors w and r for U_Z_4 with upper bounding, note that we 
// need the headings C, S_{B}, B_{O} before they are updated
template < typename Q, typename ET, typename Tags >                             // Upper bounded
void  QP_solver<Q, ET, Tags>::
z_replace_variable_slack_by_slack_upd_w_r(Tag_false )
{
    
    int     sigma_j = slack_A[ index_entering-qp_n].first;
    
    // swap r_sigma_j in r_C with r_sigma_i in r_S_B
    std::swap(r_C[in_C[sigma_j]], r_S_B[in_B[index_leaving]]);
}

// update of the vectors r_C and r_S_B with "x_j" column
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
update_r_C_r_S_B__j(ET& x_j)
{
  // TAG: INEFFICIENT because not C-filtered
  // TAG: TODO maybe make the no_ineq a template switch
  // update of vector r_{C}
  
  A_sparse_column_iterator it = (*(qp_A_sparse+index_entering)).begin();
  A_sparse_column_iterator it_end = (*(qp_A_sparse+index_entering)).end();
  
  if (no_ineq) { // in_C is not kept up to date, we have to set it up
    std::fill_n(in_C.begin(), qp_m, -1);
    for (int i = 0; i < qp_m; ++i) { // all equalities are in C
      in_C[ C[i] ] = i;
    }
  }
  
  while (it != it_end) {
    if (in_C[it->first] >= 0) {
      r_C[in_C[it->first]] -= x_j * static_cast<ET>(it->second);
    }
    ++it;
  }
  
  Indices in_S_B(qp_m, -1); // TAG: TODO maybe make this global
  int i = 0;
  for (Index_iterator S_B_it = S_B.begin(); S_B_it != S_B.end(); ++S_B_it) {
    in_S_B[*S_B_it] = i;
    ++i;
  }
  
  // reset A iterators
  it = (*(qp_A_sparse+index_entering)).begin();
  it_end = (*(qp_A_sparse+index_entering)).end();
  
  // update of r_{S_{B}}
  while (it != it_end) {
    if (in_S_B[it->first] >= 0) {
      r_S_B[in_S_B[it->first]] -= x_j * static_cast<ET>(it->second);
    }
    ++it;
  }
}

// update of the vectors r_C and r_S_B with "x_j" and "x_i" column
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
update_r_C_r_S_B__j_i(ET& x_j, ET& x_i)
{

  // TAG: INEFFICIENT because not C-filtered
  // TAG: TODO maybe make the no_ineq a template switch
  // update of vector r_{C}
  A_sparse_column_iterator it_entering = (*(qp_A_sparse+index_entering)).begin();
  A_sparse_column_iterator it_entering_end = (*(qp_A_sparse+index_entering)).end();
  A_sparse_column_iterator it_leaving = (*(qp_A_sparse+index_leaving)).begin();
  A_sparse_column_iterator it_leaving_end = (*(qp_A_sparse+index_leaving)).end();
  
  if (no_ineq) { // in_C is not kept up to date, we have to set it up
    std::fill_n(in_C.begin(), qp_m, -1);
    for (int i = 0; i < qp_m; ++i) { // all equalities are in C
      in_C[ C[i] ] = i;
    }
  }
  
  while (it_entering != it_entering_end) {
    if (in_C[it_entering->first] >= 0) {
      r_C[in_C[it_entering->first]] -= x_j * static_cast<ET>(it_entering->second);
    }
    ++it_entering;
  }
  while (it_leaving != it_leaving_end) {
    if (in_C[it_leaving->first] >= 0) {
      r_C[in_C[it_leaving->first]] += x_i * static_cast<ET>(it_leaving->second);
    }
    ++it_leaving;
  }
  
  
  Indices in_S_B(qp_m, -1); // TAG: TODO maybe make this global
  int i = 0;
  for (Index_iterator S_B_it = S_B.begin(); S_B_it != S_B.end(); ++S_B_it) {
    in_S_B[*S_B_it] = i;
    ++i;
  }
  
  // reset A iterators
  it_entering = (*(qp_A_sparse+index_entering)).begin();
  it_entering_end = (*(qp_A_sparse+index_entering)).end();
  it_leaving = (*(qp_A_sparse+index_leaving)).begin();
  it_leaving_end = (*(qp_A_sparse+index_leaving)).end();
  
  // update of r_{S_{B}}
  while (it_entering != it_entering_end) {
    if (in_S_B[it_entering->first] >= 0) {
      r_S_B[in_S_B[it_entering->first]] -= x_j * static_cast<ET>(it_entering->second);
    }
    ++it_entering;
  }
  while (it_leaving != it_leaving_end) {
    if (in_S_B[it_leaving->first] >= 0) {
      r_S_B[in_S_B[it_leaving->first]] += x_i * static_cast<ET>(it_leaving->second);
    }
    ++it_leaving;
  }
}

// update of the vectors r_C and r_S_B with "x_i'" column
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
update_r_C_r_S_B__i(ET& x_i)
{
  // TAG: INEFFICIENT because not C-filtered
  // TAG: TODO maybe make the no_ineq a template switch
  // update of vector r_{C}
  
  A_sparse_column_iterator it = (*(qp_A_sparse+index_leaving)).begin();
  A_sparse_column_iterator it_end = (*(qp_A_sparse+index_leaving)).end();
  
  if (no_ineq) { // in_C is not kept up to date, we have to set it up
    std::fill_n(in_C.begin(), qp_m, -1);
    for (int i = 0; i < qp_m; ++i) { // all equalities are in C
      in_C[ C[i] ] = i;
    }
  }
  
  while (it != it_end) {
    if (in_C[it->first] >= 0) {
      r_C[in_C[it->first]] += x_i * static_cast<ET>(it->second);
    }
    ++it;
  }

  Indices in_S_B(qp_m, -1); // TAG: TODO maybe make this global
  int i = 0;
  for (Index_iterator S_B_it = S_B.begin(); S_B_it != S_B.end(); ++S_B_it) {
    in_S_B[*S_B_it] = i;
    ++i;
  }
  
  // reset A iterators
  it = (*(qp_A_sparse+index_leaving)).begin();
  it_end = (*(qp_A_sparse+index_leaving)).end();
  
  // update of r_{S_{B}}
  while (it != it_end) {
    if (in_S_B[it->first] >= 0) {
      r_S_B[in_S_B[it->first]] += x_i * static_cast<ET>(it->second);
    }
    ++it;
  }
}
  

// Update of w and r_B_O with "x_j" column.
//
// todo: could be optimized slightly by factoring out the factor 2
// (which is implicitly contained in the pairwise accessor for D).
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
update_w_r_B_O__j(ET& x_j)
{
  // assertion checking:
  CGAL_expensive_assertion(!Is_nonnegative::value);

  // Note: we only do anything it we are dealing with a QP.
  if (!Is_linear::value) {

    D_sparse_column_iterator it = (*(qp_D_sparse+index_entering)).begin();
    D_sparse_column_iterator it_end = (*(qp_D_sparse+index_entering)).end();
    Values temp(B_O.size(), et0);
    while (it != it_end) {
      // update of vector w:
      w[it->first] -= static_cast<ET>(it->second) * x_j;
      if (in_B[it->first] > -1) {
        temp[in_B[it->first]] = it->second;
      }
      ++it;
    }
    for (int i = 0; i < r_B_O.size(); ++i) {
      // update of r_B_O:
      r_B_O[i] -= temp[i] * x_j;
    }
  }
}

// update of w and r_B_O with "x_j" and "x_i" column
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
update_w_r_B_O__j_i(ET& x_j, ET& x_i)
{
  // assertion checking:
  CGAL_expensive_assertion(!Is_nonnegative::value);

  // Note: we only do anything it we are dealing with a QP.
  if (!Is_linear::value) {
  
    D_sparse_column_iterator it_entering = (*(qp_D_sparse+index_entering)).begin();
    D_sparse_column_iterator it_entering_end = (*(qp_D_sparse+index_entering)).end();
    D_sparse_column_iterator it_leaving = (*(qp_D_sparse+index_leaving)).begin();
    D_sparse_column_iterator it_leaving_end = (*(qp_D_sparse+index_leaving)).end();
    Values temp_entering(B_O.size(), et0);
    Values temp_leaving(B_O.size(), et0);
    
    while (it_entering != it_entering_end && it_leaving != it_leaving_end) {
      if (it_entering->first == it_leaving->first) {
        // update of vector w
        w[it_leaving->first] += static_cast<ET>(it_leaving->second) * x_i - static_cast<ET>(it_entering->second) * x_j;
        if (in_B[it_entering->first] > -1) {
          temp_entering[in_B[it_entering->first]] = it_entering->second;
          temp_leaving[in_B[it_leaving->first]] = it_leaving->second;
        }        
        ++it_leaving;
        ++it_entering;
      } else {
        if (in_B[it_entering->first] > -1) {
          temp_entering[in_B[it_entering->first]] = it_entering->second;
        }
        if (in_B[it_leaving->first] > -1) {
          temp_leaving[in_B[it_leaving->first]] = it_leaving->second;
        }
        if (it_entering->first < it_leaving->first) {
          // update of vector w
          w[it_entering->first] -= static_cast<ET>(it_entering->second) * x_j;
          ++it_entering;
        } else {
          // update of vector w
          w[it_leaving->first] += static_cast<ET>(it_leaving->second) * x_i;
          ++it_leaving;
        }
      }
    }
    if (it_entering != it_entering_end) {
      while (it_entering != it_entering_end) {
        // update of vector w
        w[it_entering->first] -= static_cast<ET>(it_entering->second) * x_j;
        if (in_B[it_entering->first] > -1) {
          temp_entering[in_B[it_entering->first]] = it_entering->second;
        }
        ++it_entering;
      }
    }
    if (it_leaving != it_leaving_end) {
      while (it_leaving != it_leaving_end) {
        // update of vector w
        w[it_leaving->first] += static_cast<ET>(it_leaving->second) * x_i;
        if (in_B[it_leaving->first] > -1) {
          temp_leaving[in_B[it_leaving->first]] = it_leaving->second;
        }
        ++it_leaving;
      }
    }
    
    for (int i = 0; i < r_B_O.size(); ++i) {
      // update of r_B_O:
      r_B_O[i] += temp_leaving[i] * x_i - temp_entering[i] * x_j;
    }
  }
}

// update of w and r_B_O with "x_i" column
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
update_w_r_B_O__i(ET& x_i)
{
  CGAL_expensive_assertion(!Is_nonnegative::value);

  // Note: we only do anything it we are dealing with a QP.
  if (!Is_linear::value) {

    D_sparse_column_iterator it = (*(qp_D_sparse+index_leaving)).begin();
    D_sparse_column_iterator it_end = (*(qp_D_sparse+index_leaving)).end();
    Values temp(B_O.size(), et0);
    while (it != it_end) {
      // update of vector w:
      w[it->first] += static_cast<ET>(it->second) * x_i;
      if (in_B[it->first] > -1) {
        temp[in_B[it->first]] = it->second;
      }
      ++it;
    }
    for (int i = 0; i < r_B_O.size(); ++i) {
      // update of r_B_O:
      r_B_O[i] += temp[i] * x_i;
    }
  }
}

// Compute solution, meaning compute the solution vector x and the KKT
// coefficients lambda.
template < typename Q, typename ET, typename Tags >                                      // Standard form
void  QP_solver<Q, ET, Tags>::
compute_solution(Tag_true)
{
    
  // compute current solution, original variables and lambdas
  lu_fact_.solve( b_C.begin(), minus_c_B.begin(),
		    lambda.begin(), x_B_O.begin(), Is_linear::value, is_phaseI);
        
  
  // compute current solution, slack variables
  compute__x_B_S(no_ineq, Is_nonnegative());
}

// Compute solution, meaning compute the solution vector x and the KKT
// coefficients lambda.
template < typename Q, typename ET, typename Tags >                                      // Upper bounded
void  QP_solver<Q, ET, Tags>::
compute_solution(Tag_false)
{ 
  
  // compute the difference b_C - r_C
  std::transform(b_C.begin(), b_C.begin()+C.size(), // Note: r_C.size() ==
						    // C.size() always holds,
						    // whereas b_C.size() >=
						    // C.size() in general.
		 r_C.begin(),  tmp_l.begin(), std::minus<ET>());
  
  // compute the difference minus_c_B - r_B_O:
  if (is_phaseII && is_QP) {
    std::transform(minus_c_B.begin(), minus_c_B.begin() + B_O.size(), 
		   r_B_O.begin(), tmp_x.begin(), std::minus<ET>());
    
    // compute current solution, original variables and lambdas:
    lu_fact_.solve( tmp_l.begin(), tmp_x.begin(),
		      lambda.begin(), x_B_O.begin(), Is_linear::value, is_phaseI);
  } else {                                          // r_B_O == 0
    
    // compute current solution, original variables and lambdas        
    lu_fact_.solve( tmp_l.begin(), minus_c_B.begin(),
		      lambda.begin(), x_B_O.begin(), Is_linear::value, is_phaseI);
  }
   
  
  // compute current solution, slack variables
  compute__x_B_S( no_ineq, Is_nonnegative());
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
multiply__A_S_BxB_O(Value_iterator in, Value_iterator out) const
{
  // initialize with zero vector
  std::fill_n( out, B_S.size(), et0);
  
  // foreach original column of A in B_O (artificial columns are zero in S_B except special)
  Value_iterator        out_it;
  Index_const_iterator  row_it, col_it;
  A_sparse_column_iterator it;
  A_sparse_column_iterator it_end;
  
  Indices in_S_B(qp_m, -1); // TAG: TODO maybe make this global
  int i = 0;
  for (Index_const_iterator S_B_it = S_B.begin(); S_B_it != S_B.end(); ++S_B_it) {
    in_S_B[*S_B_it] = i;
    ++i;
  }
  
  for ( col_it = B_O.begin(); col_it != B_O.end(); ++col_it, ++in) {
    const ET in_value = *in;
    out_it   = out;
    
    if ( *col_it < qp_n) {	                        // original variable
      
      it = (*(qp_A_sparse+*col_it)).begin();
      it_end = (*(qp_A_sparse+*col_it)).end();
      
      while (it != it_end) {
        if (in_S_B[it->first] >= 0) {
          *(out_it + in_S_B[it->first]) += static_cast<ET>( it->second ) * in_value;
        }
        ++it;
      }
    } else {
      if ( *col_it == art_s_i) {                  // special artificial
        
        // foreach row of 'art_s'
        for ( row_it = S_B.begin(); row_it != S_B.end(); ++row_it, ++out_it) {
          *out_it += static_cast<ET>( art_s[ *row_it]) * in_value;
        }
      }
    }
  }
}

// compare the updated vector r_{C} with t_r_C=A_{C, N_O}x_{N_O}
template < typename Q, typename ET, typename Tags >                         // Standard form
bool  QP_solver<Q, ET, Tags>::
check_r_C(Tag_true) const
{
    return true;
}

// compare the updated vector r_{C} with t_r_C=A_{C, N_O}x_{N_O}
template < typename Q, typename ET, typename Tags >                         // Upper bounded
bool  QP_solver<Q, ET, Tags>::
check_r_C(Tag_false) const
{
    Values                  t_r_C;
    // compute t_r_C from scratch
    t_r_C.resize(C.size(), et0);
    multiply__A_CxN_O(t_r_C.begin());
    
    // compare r_C and the from scratch computed t_r_C
    bool failed = false;
    int count = 0;
    Value_const_iterator t_r_C_it = r_C.begin();
    for (Value_const_iterator r_C_it = r_C.begin(); r_C_it != r_C.end();
                                    ++r_C_it, ++t_r_C_it, ++count) {
        if (*r_C_it != *t_r_C_it) {
            failed = true;
        }
    }
    return (!failed);
}

// compare the updated vector r_{S_B} with t_r_S_B=A_{S_B, N_O}x_{N_O}
template < typename Q, typename ET, typename Tags >                             // Standard form
bool  QP_solver<Q, ET, Tags>::
check_r_S_B(Tag_true) const
{
    return true;
}

// compare the updated vector r_{S_B} with t_r_S_B=A_{S_B, N_O}x_{N_O}
template < typename Q, typename ET, typename Tags >                             // Upper bounded
bool  QP_solver<Q, ET, Tags>::
check_r_S_B(Tag_false) const
{
    Values                  t_r_S_B;
    // compute t_r_S_B from scratch
    t_r_S_B.resize(S_B.size(), et0);
    multiply__A_S_BxN_O(t_r_S_B.begin());
    
    // compare r_S_B and the from scratch computed t_r_S_B
    bool failed = false;
    int count = 0;
    Value_const_iterator    t_r_S_B_it = t_r_S_B.begin();
    for (Value_const_iterator r_S_B_it = r_S_B.begin(); r_S_B_it != r_S_B.end();
                                    ++r_S_B_it, ++t_r_S_B_it, ++count) {
        if (*r_S_B_it != *t_r_S_B_it) {
            failed = true;
        }
    }
    return (!failed);
}


// computes r_{B_{O}}:=2D_{B_O, N_O}x_{N_O} with upper bounding
// OPTIMIZATION: If D is symmetric we can multiply by two at the end of the
// computation of entry of r_B_O instead of each access to D
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
multiply__2D_B_OxN_O(Value_iterator out) const
{
    //initialize
    std::fill_n( out, B_O.size(), et0);

    Value_iterator          out_it;
    ET                      value;
    
    // foreach entry in r_B_O
    out_it = out;
    
    for (Index_const_iterator row_it = B_O.begin(); row_it != B_O.end(); ++row_it, ++out_it) {
      D_sparse_column_iterator it = (*(qp_D_sparse+*row_it)).begin();
      D_sparse_column_iterator it_end = (*(qp_D_sparse+*row_it)).end();
      while (it != it_end) {
        if (!is_basic(it->first)) {
          value = nonbasic_original_variable_value(it->first);
          *out_it += static_cast<ET>(it->second) * value;
        }
        ++it;
      } 
    }
}

// compares the updated vector r_{B_O} with t_r_B_O=2D_{B_O, N_O}x_{N_O}
template < typename Q, typename ET, typename Tags >                                // Standard form
bool  QP_solver<Q, ET, Tags>::
check_r_B_O(Tag_true) const
{
    return true;
}

// compares the updated vector r_{B_O} with t_r_B_O=2D_{B_O, N_O}x_{N_O}
template < typename Q, typename ET, typename Tags >                                 // Upper bounded
bool  QP_solver<Q, ET, Tags>::
check_r_B_O(Tag_false) const
{
    Values                  t_r_B_O;
    // compute t_r_B_O from scratch
    t_r_B_O.resize(B_O.size(), et0);
    multiply__2D_B_OxN_O(t_r_B_O.begin());
    
    // compare r_B_O and the from scratch computed t_r_B_O
    bool failed = false;
    int count = 0;
    Value_const_iterator    t_r_B_O_it = t_r_B_O.begin();
    for (Value_const_iterator r_B_O_it = r_B_O.begin(); r_B_O_it != r_B_O.end();
                                    ++r_B_O_it, ++t_r_B_O_it, ++count) {
        if (*r_B_O_it != *t_r_B_O_it) {
            failed = true;
        }
    }
    return (!failed);   
}

// compares the updated vector w with t_w=2D_{O,N_O}*x_N_O
template < typename Q, typename ET, typename Tags >                             // Standard form
bool  QP_solver<Q, ET, Tags>::
check_w(Tag_true) const
{
    return true;
}

// compares the updated vector w with t_w=2D_O_N_O*x_N_O
template < typename Q, typename ET, typename Tags >                             // Upper bounded
bool  QP_solver<Q, ET, Tags>::
check_w(Tag_false) const
{
    Values              t_w;
    // compute t_w from scratch
    t_w.resize( qp_n, et0);
    multiply__2D_OxN_O(t_w.begin());
    
    // compare w and the from scratch computed t_w
    bool  failed = false;
    Value_const_iterator    t_w_it = t_w.begin();
    for (int i = 0; i < qp_n; ++i, ++t_w_it) {
        if (w[i] != *t_w_it) {
            failed = true;
        }
    }
    return (!failed);
}

// check basis inverse
template < typename Q, typename ET, typename Tags >
bool
QP_solver<Q, ET, Tags>::
check_basis_inverse()
{
    // diagnostic output
    CGAL_qpe_debug {
	vout4 << "check: " << std::flush;
    }
    bool ok;
    if (is_phaseI) {
    	ok = check_basis_inverse(Tag_true());
    } else {
        ok = check_basis_inverse( Is_linear());
    }

    // diagnostic output
    CGAL_qpe_debug {
	if ( ok) {
	    vout4 << "check ok";
	} else {
	    vout4 << "check failed";
	}
	vout4 << std::endl;
    }

    return ok;
}

template < typename Q, typename ET, typename Tags >                                         // LP case
bool  QP_solver<Q, ET, Tags>::
check_basis_inverse( Tag_true)
{
    CGAL_qpe_debug {
	vout4 << std::endl;
    }
    bool res = true;
    unsigned int    row, rows =   static_cast<unsigned int>(C.size());
    unsigned int    col, cols = static_cast<unsigned int>(B_O.size());
    Index_iterator  i_it = B_O.begin();
    Value_iterator  q_it;

    
    // BG: is this a real check?? How does the special artifical
    // variable come in, e.g.? OK: it comes in through
    // ratio_test_init__A_Cj
    for ( col = 0; col < cols; ++col, ++i_it) {
	ratio_test_init__A_Cj( tmp_l.begin(), *i_it,
			       no_ineq);
	lu_fact_.solve_x( tmp_l.begin(), q_x_O.begin());

	CGAL_qpe_debug {
       	    if ( vout4.verbose()) {
		std::copy( tmp_l.begin(), tmp_l.begin()+rows,
			   std::ostream_iterator<ET>( vout4.out(), " "));
		vout4.out() << " ||  ";
		std::copy( q_x_O.begin(), q_x_O.begin()+cols,
			   std::ostream_iterator<ET>( vout4.out(), " "));
		vout4.out() << std::endl;
	    }
	}

	q_it = q_x_O.begin();
	for ( row = 0; row < rows; ++row, ++q_it) {
	    if ( *q_it != ( row == col ? denominator_ : et0)) {
		if ( ! vout4.verbose()) {
		    std::cerr << std::endl << "basis-inverse check: ";
		}
		std::cerr << "failed ( row=" << row << " | col=" << col << " )"
		          << std::endl;
		res = false;
	    }
	}
    }

    return res;
}

template < typename Q, typename ET, typename Tags >                                         // QP case
bool  QP_solver<Q, ET, Tags>::
check_basis_inverse( Tag_false)
{
  bool res = true;
  unsigned int    row, rows =   static_cast<unsigned int>(C.size());
  unsigned int    col, cols = static_cast<unsigned int>(B_O.size());
  Value_iterator  v_it;
  Index_iterator  i_it;
  
  CGAL_qpe_debug {
    vout4 << std::endl;
  }
  A_sparse_column_iterator it, it_end;
  
  // left part of M_B
  std::fill_n( tmp_l.begin(), rows, et0);
  for ( col = 0; col < rows; ++col) {

    // get column of A_B^T (i.e. row of A_B)
    row = ( has_ineq ? C[ col] : col);
    v_it = tmp_x.begin();
    for ( i_it = B_O.begin(); i_it != B_O.end(); ++i_it, ++v_it) {
      if (*i_it < qp_n) { // original
        // TAG: INEFFICIENT, binary search...?
        it = (*(qp_A_sparse+*i_it)).begin();
        it_end = (*(qp_A_sparse+*i_it)).end();
        while (it != it_end && it->first < row) {
          ++it;
        }
        if (it != it_end && it->first == row) {
          *v_it = it->second;
        }
      } else { //artificial
        *v_it = (art_A[ *i_it - qp_n].first != static_cast<int>(row) ? et0 :// artific.
                (art_A[ *i_it - qp_n].second ? -et1 : et1));
      }
    }
    
    CGAL_qpe_assertion (art_s_i < 0);
    
    lu_fact_.solve( tmp_l.begin(), tmp_x.begin(),
                   q_lambda.begin(), q_x_O.begin(), Is_linear::value, is_phaseI);
    
    CGAL_qpe_debug {
      if ( vout4.verbose()) {
        std::copy( tmp_l.begin(), tmp_l.begin()+rows,
                  std::ostream_iterator<ET>( vout4.out(), " "));
        vout4.out() << "| ";
        std::copy( tmp_x.begin(), tmp_x.begin()+cols,
                  std::ostream_iterator<ET>( vout4.out(), " "));
        vout4.out() << " ||  ";
        std::copy( q_lambda.begin(), q_lambda.begin()+rows,
                  std::ostream_iterator<ET>( vout4.out(), " "));
        vout4.out() << " |  ";
        std::copy( q_x_O.begin(), q_x_O.begin()+cols,
                  std::ostream_iterator<ET>( vout4.out(), " "));
        vout4.out() << std::endl;
      }
    }
    
    v_it = q_lambda.begin();
    for ( row = 0; row < rows; ++row, ++v_it) {
      if ( *v_it != ( row == col ? denominator_ : et0)) {
        if ( ! vout4.verbose()) {
          std::cerr << std::endl << "basis-inverse check: ";
        }
        std::cerr << "failed ( row=" << row << " | col=" << col << " )"
        << std::endl;
        //		return false;
        res = false;
      }
    }
    v_it = std::find_if( q_x_O.begin(), q_x_O.begin()+cols,
                        std::bind2nd( std::not_equal_to<ET>(), et0));
    if ( v_it != q_x_O.begin()+cols) {
      if ( ! vout4.verbose()) {
        std::cerr << std::endl << "basis-inverse check: ";
      }
      std::cerr << "failed ( row=" << rows+(v_it-q_x_O.begin())
      << " | col=" << col << " )" << std::endl;
      // ToDo: return false;
      res = false;
    }
  }
  vout4 << "= = = = = = = = = =" << std::endl;
  
  // right part of M_B
  if ( is_phaseI) std::fill_n( tmp_x.begin(), B_O.size(), et0);
  i_it = B_O.begin();
  for ( col = 0; col < cols; ++col, ++i_it) {
    ratio_test_init__A_Cj  ( tmp_l.begin(), *i_it, 
                            no_ineq);
    ratio_test_init__2_D_Bj( tmp_x.begin(), *i_it, Tag_false());
    
    lu_fact_.solve( tmp_l.begin(), tmp_x.begin(),
                   q_lambda.begin(), q_x_O.begin(), Is_linear::value, is_phaseI);
    
    
    
    CGAL_qpe_debug {
      if ( vout4.verbose()) {
        std::copy( tmp_l.begin(), tmp_l.begin()+rows,
                  std::ostream_iterator<ET>( vout4.out(), " "));
        vout4.out() << "| ";
        std::copy( tmp_x.begin(), tmp_x.begin()+cols,
                  std::ostream_iterator<ET>( vout4.out(), " "));
        vout4.out() << " ||  ";
        std::copy( q_lambda.begin(), q_lambda.begin()+rows,
                  std::ostream_iterator<ET>( vout4.out(), " "));
        vout4.out() << " |  ";
        std::copy( q_x_O.begin(), q_x_O.begin()+cols,
                  std::ostream_iterator<ET>( vout4.out(), " "));
        vout4.out() << std::endl;
      }
    }
    
    v_it = std::find_if( q_lambda.begin(), q_lambda.begin()+rows,
                        std::bind2nd( std::not_equal_to<ET>(), et0));
    if ( v_it != q_lambda.begin()+rows) {
      if ( ! vout4.verbose()) {
        std::cerr << std::endl << "basis-inverse check: ";
      }
      std::cerr << "failed ( row=" << v_it-q_lambda.begin()
      << " | col=" << col << " )" << std::endl;
      //	    return false;
      res = false;
    }
    v_it = q_x_O.begin();
    for ( row = 0; row < cols; ++row, ++v_it) {
      if ( *v_it != ( row == col ? denominator_ : et0)) {
        if ( ! vout4.verbose()) {
          std::cerr << std::endl << "basis-inverse check: ";
        }
        std::cerr << "failed ( row=" << row+rows << " | col="
        << col << " )" << std::endl;
        //		return false;
        res = false;
      }
    }
  }
  return res;
  }

// filtered strategies are only allowed if the input type as
// indicated by C_entry is double; this may still fail if the
// types of some other program entries are not double, but 
// since the filtered stuff is internal anyway, we can probably
// live with this simple check
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
set_pricing_strategy
( Quadratic_program_pricing_strategy strategy)
{
    CGAL_qpe_assertion( phase() != 1);
    CGAL_qpe_assertion( phase() != 2);

    if (strategy == QP_DANTZIG)
      strategyP = new QP_full_exact_pricing<Q, ET, Tags>;
    else if (strategy == QP_FILTERED_DANTZIG)    
      // choose between FF (double) and FE (anything else)
      strategyP = 
	new typename QP_solver_impl::Filtered_pricing_strategy_selector
	<Q, ET, Tags, C_entry>::FF;
    else if (strategy == QP_PARTIAL_DANTZIG)
      strategyP = new QP_partial_exact_pricing<Q, ET, Tags>;
    else if (strategy == QP_PARTIAL_FILTERED_DANTZIG
	     || strategy == QP_CHOOSE_DEFAULT)   
      // choose between PF (double) and PE (anything else)
      strategyP = 
	new typename QP_solver_impl::Filtered_pricing_strategy_selector
	<Q, ET, Tags, C_entry>::PF;
    else if (strategy == QP_BLAND) {
      strategyP = new QP_exact_bland_pricing<Q, ET, Tags>;
      bland_flag = true;
    }
    CGAL_qpe_assertion(strategyP != static_cast<Pricing_strategy*>(0));
   
    if ( phase() != -1) strategyP->set( *this, vout2);
}

// diagnostic output
// -----------------
template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
set_verbosity( int verbose, std::ostream& stream)
{
    vout  = Verbose_ostream( verbose >  0, stream);
    vout1 = Verbose_ostream( verbose == 1, stream);
    vout2 = Verbose_ostream( verbose >= 2, stream);
    vout3 = Verbose_ostream( verbose >= 3, stream);
    vout4 = Verbose_ostream( verbose == 4, stream);
    vout5 = Verbose_ostream( verbose == 5, stream);
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
print_program( ) const
{
  int  row, i;
  
  // objective function
  vout4.out() << std::endl << "objective function:" << std::endl;
  if ( is_QP) {
    vout4.out() << "--> output of MATRIX must go here <--";
    vout4.out() << std::endl;
  }
  std::copy( qp_c, qp_c+qp_n,
            std::ostream_iterator<C_entry>( vout4.out(), " "));
  vout4.out() << "(+ " << qp_c0 << ") ";
  vout4.out() << std::endl;
  vout4.out() << std::endl;
  
  std::vector<A_sparse_column_iterator> it(qp_n);
  std::vector<A_sparse_column_iterator> it_end(qp_n);
  for (int i = 0; i < qp_n; ++i) {
    it[i] = (*(qp_A_sparse+i)).begin();
    it_end[i] = (*(qp_A_sparse+i)).end();
  }
  
  // constraints
  vout4.out() << "constraints:" << std::endl;
  for ( row = 0; row < qp_m; ++row) {
    // original variables
    for (i = 0; i < qp_n; ++i) {
      if (it[i] != it_end[i]) {
        if (it[i]->first == row) {
          vout4.out() << it[i]->second << ' ';
        } else {
          vout4.out() << et0 << ' ';
        }
        ++(it[i]);
      } else {
        vout4.out() << et0 << ' ';
      }
    }
    
    // slack variables
    if ( ! slack_A.empty()) {
      vout4.out() << " |  ";
      for ( i = 0; i < static_cast<int>(slack_A.size()); ++i) {
        vout4.out() << ( slack_A[ i].first != row ? " 0" :
                        ( slack_A[ i].second ? "-1" : "+1")) << ' ';
      }
    }
    
    // artificial variables
    if ( ! art_A.empty()) {
      vout4.out() << " |  ";
      for ( i = 0; i < static_cast<int>(art_A.size()); ++i) {
        if (art_s_i == i+qp_n+static_cast<int>(slack_A.size()))
          vout4.out() << " * ";          // for special artificial column
        vout4.out() << ( art_A[ i].first != row ? " 0" :
                        ( art_A[ i].second ? "-1" : "+1")) << ' ';
      }
    }
    if ( ! art_s.empty()) vout4.out() << " |  " << art_s[ row] << ' ';
    
    // rhs
    vout4.out() << " |  "
    << ( *(qp_r+ row) == CGAL::EQUAL      ? ' ' :
        ( *(qp_r+ row) == CGAL::SMALLER ? '<' : '>')) << "=  "
    << *(qp_b+ row);
    if (!is_nonnegative) {
      vout4.out() << " - " << multiply__A_ixO(row);
    }
    vout4.out() << std::endl;
  }
  vout4.out() << std::endl;
  
  // explicit bounds
  if (!is_nonnegative) {
    vout4.out() << "explicit bounds:" << std::endl; 
    for (int i = 0; i < qp_n; ++i) {
      if (*(qp_fl+i)) {                   // finite lower bound
        vout4.out() << *(qp_l+i);
      } else {                            // infinite lower bound
        vout4.out() << "-inf";
      }
      vout4.out() << " <= x_" << i << " <= ";
      if (*(qp_fu+i)) {
        vout4.out() << *(qp_u+i);
      } else {
        vout4.out() << "inf";
      }
      vout4.out() << std::endl;
    }
  }
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
print_basis( ) const
{
  char label;
  vout << "  basis: ";
  CGAL_qpe_debug {
    vout2 << "basic variables" << ( has_ineq ? "  " : "") << ":  ";
  }
  std::copy( B_O.begin(), B_O.end  (),
            std::ostream_iterator<int>( vout.out(), " "));
  CGAL_qpe_debug {
    if ( vout2.verbose()) {
      if ( has_ineq && ( ! slack_A.empty())) {
        vout2.out() << " |  ";
        std::copy( B_S.begin(), B_S.end(),
                  std::ostream_iterator<int>( vout2.out(), " "));
      }
      if ( is_phaseI) {
        vout2.out() << " (# of artificials: " << art_basic << ')';
      }
      if ( has_ineq) {
        vout2.out() << std::endl
        << "basic constraints:  ";
        for (Index_const_iterator i_it = 
             C.begin(); i_it != C.end(); ++i_it) {
          label = (*(qp_r+ *i_it) == CGAL::EQUAL) ? 'e' : 'i';
          vout2.out() << *i_it << ":" << label << " ";
        }
        /*
         std::copy( C.begin(), C.begin()+(C.size()-slack_A.size()),
         std::ostream_iterator<int>( vout2.out(), " "));
         if ( ! slack_A.empty()) {
         vout2.out() << " |  ";
         std::copy( C.end() - slack_A.size(), C.end(),
         std::ostream_iterator<int>( vout2.out(), " "));
         }
         */
      }
      if ( vout3.verbose()) {
        vout3.out() << std::endl
        << std::endl
        << "    in_B: ";
        std::copy( in_B.begin(), in_B.end(),
                  std::ostream_iterator<int>( vout3.out(), " "));
        if ( has_ineq) {
          vout3.out() << std::endl
          << "    in_C: ";
          std::copy( in_C.begin(), in_C.end(),
                    std::ostream_iterator<int>( vout3.out(), " "));
        }
      }
    }
  }
  vout.out() << std::endl;
  CGAL_qpe_debug {
    vout4 << std::endl << "basis-inverse:" << std::endl;
  }
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
print_solution( ) const
{
  CGAL_qpe_debug {
    if ( vout3.verbose()) {
	   vout3.out() << std::endl
            << "     b_C: ";
	   std::copy( b_C.begin(), b_C.begin()+C.size(),
		   std::ostream_iterator<ET>( vout3.out()," "));
	   vout3.out() << std::endl
            << "  -c_B_O: ";
	   std::copy( minus_c_B.begin(), minus_c_B.begin()+B_O.size(),
		   std::ostream_iterator<ET>( vout3.out()," "));
	   if (!is_nonnegative) {
	       vout3.out() << std::endl
                << "     r_C: ";
	       std::copy( r_C.begin(), r_C.begin()+r_C.size(),
	           std::ostream_iterator<ET>( vout3.out(), " "));
	       vout3.out() << std::endl
	           << "   r_B_O: ";
	       std::copy( r_B_O.begin(), r_B_O.begin()+r_B_O.size(),
	           std::ostream_iterator<ET>( vout3.out(), " "));
	       if (r_B_O.size() == 0)
		 vout3.out() << "< will only be allocated in phase II>";

	   }
	   vout3.out() << std::endl;
    }
    if ( vout2.verbose()) {
        vout2.out() << std::endl << "  lambda: ";
        std::copy( lambda.begin(), lambda.begin()+C.size(),
            std::ostream_iterator<ET>( vout2.out(), " "));
        vout2.out() << std::endl << "   x_B_O: ";
        std::copy( x_B_O.begin(), x_B_O.begin()+B_O.size(),
            std::ostream_iterator<ET>( vout2.out(), " "));
        vout2.out() << std::endl;
        if (!is_nonnegative) {
            vout2.out() << "   x_N_O: ";
            for (int i = 0; i < qp_n; ++i) {
                if (!is_basic(i)) {
                    vout2.out() << nonbasic_original_variable_value(i);
                    vout2.out() << " ";    
                }
            }
            vout2.out() << std::endl;
        }
        if ( has_ineq) {
            vout2.out() << "   x_B_S: ";
            std::copy( x_B_S.begin(), x_B_S.begin()+B_S.size(),
                std::ostream_iterator<ET>( vout2.out()," "));
            vout2.out() << std::endl;
        }
    }
  }
  vout << "  ";
  vout.out() 
    << "solution: " 
    << solution_numerator() << " / " << solution_denominator() 
    << "  ~= " 
    << to_double
    (CGAL::Quotient<ET>(solution_numerator(), solution_denominator()))
    << std::endl;
  CGAL_qpe_debug {
      vout2 << std::endl;
  }
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
print_ratio_1_original(int k, const ET& x_k, const ET& q_k)
{
    if (is_nonnegative) {                      // => direction== 1
        if (q_k > et0) {                            // check for lower bound
            vout2.out() << "t_O_" << k << ": "
            << x_k << '/' << q_k
            << ( ( q_i != et0) && ( index_leaving == B_O[ k]) ? " *" : "")
            << std::endl;
        } else if (q_k < et0) {                     // check for upper bound
            vout2.out() << "t_O_" << k << ": "
            << "inf" << '/' << q_k
            << ( ( q_i != et0) && ( index_leaving == B_O[ k]) ? " *" : "")
            << std::endl;
        } else {                                    // q_k == 0
            vout2.out() << "t_O_" << k << ": "
            << "inf"
            << ( ( q_i != et0) && ( index_leaving == B_O[ k]) ? " *" : "")
            << std::endl;
        }
    } else {                                        // upper bounded
        if (q_k * static_cast<ET>(direction) > et0) {            // check for lower bound
            if (B_O[k] < qp_n) {                         // original variable
                if (*(qp_fl+B_O[k])) {                   // finite lower bound
                    vout2.out() << "t_O_" << k << ": "
                    << x_k - (denominator_ * static_cast<ET>(*(qp_l+B_O[k]))) << '/' << q_k
                    << ( ( q_i != et0) && ( index_leaving == B_O[ k]) ? " *" : "")
                    << std::endl;
                } else {                            // lower bound -infinity
                    vout2.out() << "t_O_" << k << ": "
                    << "-inf" << '/' << q_k
                    << ( ( q_i != et0) && ( index_leaving == B_O[ k]) ? " *" : "")
                    << std::endl;                
                }
            } else {                                // artificial variable
                vout2.out() << "t_O_" << k << ": "
                << x_k << '/' << q_k
                << ( ( q_i != et0) && ( index_leaving == B_O[ k]) ? " *" : "")
                << std::endl;
            }
        } else if (q_k * static_cast<ET>(direction) < et0) {     // check for upper bound
            if (B_O[k] < qp_n) {                         // original variable
                if (*(qp_fu+B_O[k])) {                   // finite upper bound
                    vout2.out() << "t_O_" << k << ": "
                    << (denominator_ * static_cast<ET>(*(qp_l+B_O[k]))) - x_k << '/' << q_k
                    << ( ( q_i != et0) && ( index_leaving == B_O[ k]) ? " *" : "")
                    << std::endl;                    
                } else {                            // upper bound infinity
                    vout2.out() << "t_O_" << k << ": "
                    << "inf" << '/' << q_k
                    << ( ( q_i != et0) && ( index_leaving == B_O[ k]) ? " *" : "")
                    << std::endl;
                }
            } else {                                // artificial variable
                vout2.out() << "t_O_" << k << ": "
                << "inf" << '/' << q_k
                << ( ( q_i != et0) && ( index_leaving == B_O[ k]) ? " *" : "")
                << std::endl;
            }
        } else {                                    // q_k == 0
            vout2.out() << "t_O_" << k << ": "
            <<  "inf"
            << ( ( q_i != et0) && ( index_leaving == B_O[ k]) ? " *" : "")
            << std::endl;
        }
    }
}

template < typename Q, typename ET, typename Tags >
void  QP_solver<Q, ET, Tags>::
print_ratio_1_slack(int k, const ET& x_k, const ET& q_k)
{
    if (is_nonnegative) {                      // => direction == 1
        if (q_k > et0) {                            // check for lower bound
            vout2.out() << "t_S_" << k << ": "
            << x_k << '/' << q_k
            << ( ( q_i != et0) && ( index_leaving == B_S[ k]) ? " *" : "")
            << std::endl;
        } else {                                    // check for upper bound
            vout2.out() << "t_S_" << k << ": "
            << "inf" << '/' << q_k
            << ( ( q_i != et0) && ( index_leaving == B_S[ k]) ? " *" : "")
            << std::endl;        
        }
    } else {                                        // upper bounded
        if (q_k * static_cast<ET>(direction) > et0) {            // check for lower bound
            vout2.out() << "t_S_" << k << ": "
            << x_k << '/' << q_k
            << ( ( q_i != et0) && ( index_leaving == B_S[ k]) ? " *" : "")
            << std::endl;            
        } else if (q_k * static_cast<ET>(direction) < et0) {     // check for upper bound
            vout2.out() << "t_S_" << k << ": "
            << "inf" << '/' << q_k
            << ( ( q_i != et0) && ( index_leaving == B_S[ k]) ? " *" : "")
            << std::endl;
        } else {                                    // q_k == 0
            vout2.out() << "t_S_" << k << ": "
            << "inf"
            << ( ( q_i != et0) && ( index_leaving == B_S[ k]) ? " *" : "")
            << std::endl;
        }
    }
}


template < typename Q, typename ET, typename Tags >
const char*  QP_solver<Q, ET, Tags>::
variable_type( int k) const
{
    return ( k <        qp_n                 ? "original"  :
	   ( k < static_cast<int>( qp_n+slack_A.size()) ? "slack"     :
	                                       "artificial"));
}

template < typename Q, typename ET, typename Tags > 
bool QP_solver<Q, ET, Tags>::
is_artificial(int k) const
{
    return (k >= static_cast<int>(qp_n+slack_A.size())); 
}

template < typename Q, typename ET, typename Tags > 
int QP_solver<Q, ET, Tags>::
get_l() const
{
    return min_N_M_;
}

template <typename Q, typename ET, typename Tags>
boost::shared_ptr<QP_sparse_matrix<ET> >
QP_solver<Q,ET,Tags>::get_basis_matrix( int& csize, int& bosize,
                                       bool is_linear) {
  
  CGAL_qpe_debug {
    //CGAL::QP_solver_debug::timer.get_basis_matrix.start();
  }
  
  if (!is_linear) {
    csize = C.size();
    bosize = B_O.size();
    
    basis_matrix_ = boost::shared_ptr<QP_sparse_matrix<ET> >(
        new QP_sparse_matrix<ET>(csize + bosize, csize + bosize, et0)
    );
    
    // Get A_B part of M_B
    A_sparse_column_iterator it, it_end;
    for (int col = 0; col < bosize; ++col) { // process columns of A_B
      if (B_O[col] < qp_n ) {
        // TAG: INEFFICIENT, binary search...?
        it = (*(qp_A_sparse+B_O[col])).begin();
        it_end = (*(qp_A_sparse+B_O[col])).end();
        while (it != it_end) {
          if (has_ineq) {
            if (in_C[it->first] >= 0) {
              basis_matrix_->set_entry(csize + col, in_C[it->first], it->second);
              basis_matrix_->set_entry(in_C[it->first], csize + col, it->second);
            }
          } else {
            basis_matrix_->set_entry(csize + col, it->first, it->second);
            basis_matrix_->set_entry(it->first, csize + col, it->second);
          }
          ++it;
        }
      } else { // articicial columns of A
        
        for (int row = 0; row < csize; ++row) {
          int ineqrow = (has_ineq ? C[row] : row);
          
          if (art_A[B_O[col] - qp_n].first == ineqrow){
            if( art_A[B_O[col] - qp_n].second){
              basis_matrix_->set_entry(row, col + csize, -et1);
              basis_matrix_->set_entry(col + csize, row, -et1);
            } else {
              basis_matrix_->set_entry(row, col + csize, et1);
              basis_matrix_->set_entry(col + csize, row, et1);
            }
          }
        }
      }
    }
    
    // Get D_B part of M_B
    D_sparse_column_iterator d_it, d_it_end;
    for (int row = 0; row < bosize; ++row) {
      d_it = (*(qp_D_sparse+B_O[row])).begin();
      d_it_end = (*(qp_D_sparse+B_O[row])).end();
      while (d_it != d_it_end) {
        if (in_B[d_it->first] > -1 && B_O[in_B[d_it->first]] <= B_O[row]) {
          basis_matrix_->set_entry(csize + row, csize + in_B[d_it->first], d_it->second);
          basis_matrix_->set_entry(csize + in_B[d_it->first], csize + row, d_it->second);          
        }
        ++d_it;
      }
    }
  } else { // is_linear
    csize = C.size();
    bosize = B_O.size();
    
    basis_matrix_ = boost::shared_ptr<QP_sparse_matrix<ET> >(new QP_sparse_matrix<ET>(bosize, bosize, et0));
    
    Values tmp(csize);  
    
    for (int col = 0; col < B_O.size(); ++col) {
      ratio_test_init__A_Cj(tmp.begin(), B_O[col], no_ineq); //update acj
      for (int row = 0; row < C.size(); ++row){
        basis_matrix_->set_entry(row, col, tmp[row]);
      }
    }  
  }

  CGAL_qpe_debug {
    //CGAL::QP_solver_debug::timer.get_basis_matrix.stop();
  }

  return basis_matrix_;
}



// PRE: index_entering is original
template <typename Q, typename ET, typename Tags>
bool
QP_solver<Q,ET,Tags>::update_basis_matrix_QP_O_in() {

  CGAL_qpe_assertion(!is_phaseI && !Is_linear::value);
  
  bool ret(false);
  
  int csize(C.size()), bosize(B_O.size());
  QP_sparse_vector<ET> y(csize+bosize), z(csize+bosize);
  Values tmp(csize+bosize);
  
  // TODO: try different insertion strategies... for a and k in enlarge
  ret = lu_fact_.enlarge(csize+bosize-1, csize+bosize-1, true /*enlarge basic set*/);
  if (!ret) {
    lu_fact_.set_invalid();
    return ret; // bail out
  }
  
  // get new row/col for basis matrix
  ratio_test_init__A_Cj(tmp.begin(), index_entering, no_ineq); //update acj
  ratio_test_init__2_D_Bj(tmp.begin()+csize, index_entering, Tag_false()/*is linear?*/);
  
  // do row update
  y.set_entry(csize+bosize-1, et1);
  
  // construct update row z
  for (int col = 0; col < csize+bosize; ++col){
    z.set_entry(col, tmp[col]);
  }
  // TODO: check what happens if we subtract et1 in the following update not here
  z.set_entry(csize+bosize-1, tmp[csize+bosize-1]-et1);
  
  
  ret = lu_fact_.rank_1_update(et1, y, z);
  
  if (!ret) {
    lu_fact_.set_invalid();
    return ret; // bail out
  }
  
  // do column update
  z.set_entry(csize+bosize-1, et0);
  
  ret = lu_fact_.rank_1_update(et1, z, y);
  if (!ret) {
    lu_fact_.set_invalid();
  }
  
  return ret;
}

// PRE: index_leaving is original
template <typename Q, typename ET, typename Tags>
bool
QP_solver<Q,ET,Tags>::update_basis_matrix_QP_O_out(unsigned int k, Values old_column, Values last_column) {

  CGAL_qpe_assertion(!is_phaseI && !Is_linear::value);
  
  bool ret(false);
  
  int csize(C.size());
  int bosize(B_O.size()+1); // correct for reduced size
  
  
  QP_sparse_vector<ET> y(csize+bosize),  z(csize+bosize);
  int index;
  
  
  // TODO: better handling of symmetric/unsymmetric permutation
  if ( (index = lu_fact_.perm_col_inv(bosize+csize-1)) == lu_fact_.perm_row_inv(bosize+csize-1) ) { // pivot permuation symmetric
  
    // reset last column/row to unit col/row
    // do row update
    y.set_entry(bosize+csize-1, et1);
    
    // construct update row z
    for (int col = 0; col < bosize+csize; ++col){
      z.set_entry(col, -last_column[col]);
    }
    z.set_entry(bosize+csize-1, et1 - last_column[bosize+csize-1]);

    ret = lu_fact_.rank_1_update(et1, y, z);
    
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
     
    
    // do column update
    z.set_entry(bosize+csize-1, et0);

    ret = lu_fact_.rank_1_update(et1, z, y);
    
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
    
    if (k != bosize+csize-1) {
      
      // install last column/row in place of column/row that will be dropped
      y.clear();
      z.clear();
      
      
      
      // reset last column/row to unit col/row
      // do row update
      y.set_entry(csize+k, et1);
      
      // construct update row z
      for (int col = 0; col < bosize+csize; ++col){
        z.set_entry(col, -old_column[col]+last_column[col]);
      }
      z.set_entry(bosize+csize-1, et0);
      z.set_entry(csize+k, last_column[bosize+csize-1]-old_column[csize+k]);

      ret = lu_fact_.rank_1_update(et1, y, z);

      if (!ret) {
        lu_fact_.set_invalid();
        return ret; // bail out
      }
      
      
      
      // do column update
      z.set_entry(csize+k, et0);

      ret = lu_fact_.rank_1_update(et1, z, y);

      if (!ret) {
        lu_fact_.set_invalid();
        return ret; // bail out
      }
      
    }
    
    ret = lu_fact_.shrink(index, true /*is_basic*/);
    
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
    
        
  } else { // pivot permutation not symmetric
    // no swap strategy implemented yet
  }
    
  
  if (!ret) {
    lu_fact_.set_invalid();
  }

  return ret;
}


// PRE: index_leaving is slack
template <typename Q, typename ET, typename Tags>
bool
QP_solver<Q,ET,Tags>::update_basis_matrix_QP_S_in(unsigned int k, Values old_column, unsigned int last_index, Values last_constraint) {

  CGAL_qpe_assertion(!is_phaseI && !Is_linear::value);
  
  bool ret(false);
  
  int csize(C.size()+1); // correct for reduced size
  int bosize(B_O.size()); 
  
  
  QP_sparse_vector<ET> y(csize+bosize), z(csize+bosize);
  int index;
 
  // TODO: better handling of symmetric/unsymmetric permutation
  if ( (index = lu_fact_.perm_col_inv(last_index)) == lu_fact_.perm_row_inv(last_index) ) { // pivot permuation symmetric

    // reset last column/row to unit col/row
    // do row update
    y.set_entry(last_index, et1);
    
    // construct update row z
    for (int col = 0; col < bosize+csize; ++col){
      z.set_entry(col, -last_constraint[col]);
    }
    z.set_entry(last_index, et1 - last_constraint[last_index]);
    ret = lu_fact_.rank_1_update(et1, y, z);
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
    // do column update
    z.set_entry(last_index, et0);
    ret = lu_fact_.rank_1_update(et1, z, y);
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
    // install last column/row in place of column/row that will be dropped
    y.clear();
    z.clear();
  
    // reset last column/row to unit col/row
    // do row update
    y.set_entry(k, et1);
    
    // construct update row z
    for (int col = 0; col < bosize+csize; ++col){
      z.set_entry(col, -old_column[col]+last_constraint[col]);
    }
    z.set_entry(last_index, et0);
    z.set_entry(k, last_constraint[last_index]-old_column[k]);
    ret = lu_fact_.rank_1_update(et1, y, z);
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
    // do column update
    z.set_entry(k, et0);
    ret = lu_fact_.rank_1_update(et1, z, y);
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
    ret = lu_fact_.shrink(index, false /*is slack*/);
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }     
  } else { // pivot permutation not symmetric
    // no swap strategy implemented yet
  }
    
  
  if (!ret) {
    lu_fact_.set_invalid();
  }
  return ret;
}


// PRE: index_leaving is slack
template <typename Q, typename ET, typename Tags>
bool
QP_solver<Q,ET,Tags>::update_basis_matrix_QP_S_out() {

  CGAL_qpe_assertion(!is_phaseI && !Is_linear::value);

  bool ret(false);
  int csize(C.size()), bosize(B_O.size());
  QP_sparse_vector<ET> y(csize+bosize), z(csize+bosize);
  Values tmp(bosize, et0);
  
  ret = lu_fact_.enlarge(csize-1, csize+bosize-1, false /*enlarge constraint set*/);
  if (!ret) {
    lu_fact_.set_invalid();
    return ret; // bail out
  }
  
  // get new row for basis matrix
  init__A_Ri(tmp.begin(), slack_A[ index_leaving-qp_n].first, no_ineq); //update ari

  // do row update
  y.set_entry(csize-1, et1);
  
  // construct update row z
  for (int col = csize; col < csize+bosize; ++col){
    z.set_entry(col, tmp[col-csize]);
  }
  // TODO: check what happens if we subtract et1 in the following update not here
  ret = lu_fact_.rank_1_update(et1, y, z);
  
  if (!ret) {
    lu_fact_.set_invalid();
    return ret; // bail out
  }
  
  // do column update
  z.set_entry(csize-1, -et1);
  ret = lu_fact_.rank_1_update(et1, z, y);
    if (!ret) {
    lu_fact_.set_invalid();
  }
  return ret;
}


// PRE: index_entering and index_leaving are original variables
template <typename Q, typename ET, typename Tags>
bool
QP_solver<Q,ET,Tags>::update_basis_matrix_O_O(unsigned int k, Values old_column, Tag_true /*linear case*/) {
  
  CGAL_qpe_assertion(old_column.size() == C.size());
  CGAL_qpe_assertion(is_phaseI || Is_linear::value);
  
  bool ret(false);
  int csize(C.size()), bosize(B_O.size());
  int swap_index_col(-1);
  
  QP_sparse_vector<ET> y(csize), z(bosize);
  
  
  Values tmp(csize);

  ET test_val(et0);  
  
  // prepare update
  ratio_test_init__A_Cj(tmp.begin(), index_entering, no_ineq); //update acj
    
  if (tmp[ lu_fact_.perm_row(lu_fact_.perm_col_inv(k)) ] == et0) { // pivot is bad
  //if (false) { // don't try swap
    
    Values swap_col(csize, et0);
    
    QP_sparse_vector<ET> y1(csize), z1(csize);
    QP_sparse_vector<ET> y2(csize), z2(csize);
    
    // TAG: TODO more efficient finding of suitable column, i.e., there can't
    // be no extraction of A_Ri in the following loop, because this makes it n^3
    // We have to extract the corresponding element individually.
    // TODO: replace perm_row[perm_column_inv[csize+k]] by stored value
    
    // find suitable swap column
    for (int j = 0; j < bosize; ++j) {
      // check if pivots are alright
      //if (tmp[ perm_row[perm_column_inv[j]] ] != et0) {
      if (tmp[ lu_fact_.perm_row(lu_fact_.perm_col_inv(j)) ] != et0) {
        test_val = et0;
        //init__A_ij(test_val, perm_row[perm_column_inv[k]], B_O[j], no_ineq);
        init__A_ij(test_val, lu_fact_.perm_row(lu_fact_.perm_col_inv(k)), B_O[j], no_ineq);
        if (test_val != et0) {
          swap_index_col = j;
          break;
        }
      }
    }
    
    
    if (swap_index_col < 0) {
      lu_fact_.set_invalid();
      return ret; // bail out
    } else { // get swap column
      ratio_test_init__A_Cj(swap_col.begin(), B_O[swap_index_col], no_ineq); //update acj
    }
    // set up new column in place of the swap column
    z1.set_entry(swap_index_col, et1);
    for (int row = 0; row < bosize; ++row){
      y1.set_entry(row, tmp[row] - swap_col[row]);
    }
    
    ret = lu_fact_.rank_1_update(et1, y1, z1);
    
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }

    // set up swap column in new place
    z2.set_entry(k, et1);
    for (int row = 0; row < bosize; ++row) {
      y2.set_entry(row, swap_col[row] - old_column[row]);
    }
    ret = lu_fact_.rank_1_update(et1, y2, z2);
    
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    // adjust internal permuations
    lu_fact_.swap_columns_physically(k, swap_index_col);
  } else {

    z.set_entry(k, et1);
    // construct update column y
    for (int row = 0; row < csize; ++row){
      y.set_entry(row, tmp[row] - old_column[row]);
    }
    ret = lu_fact_.rank_1_update(et1, y, z);
    
  }
  
  return ret;
}
  
// PRE: index_entering and index_leaving are original variables
template <typename Q, typename ET, typename Tags>
bool
QP_solver<Q,ET,Tags>::update_basis_matrix_O_O(unsigned int k, Values old_column, Tag_false /*quadratic case*/) {

  CGAL_qpe_assertion(old_column.size() == C.size()+B_O.size());
  CGAL_qpe_assertion(!is_phaseI && !Is_linear::value);
  
  bool ret(false);
  int swap_index_col(-1);
  int swap_index_row(-1);
  int csize(C.size()), bosize(B_O.size());
  
  QP_sparse_vector<ET> y(csize+bosize), z(csize+bosize);
  Values tmp(csize+bosize);

  ET test_val(et0);
  
  // get new column/row
  ratio_test_init__A_Cj(tmp.begin(), index_entering, no_ineq); //update acj
  ratio_test_init__2_D_Bj(tmp.begin()+csize, index_entering, Tag_false()/*is linear?*/);
  
  
  // prepare update
  z.set_entry(csize+k, et1);
  for (int row = 0; row < csize+bosize; ++row){
    y.set_entry(row, tmp[row] - old_column[row]);
  }
  
  if (tmp[ lu_fact_.perm_col(lu_fact_.perm_row_inv(csize+k)) ] == et0 || tmp[ lu_fact_.perm_col(lu_fact_.perm_row_inv(csize+k)) ] == et0) { // pivot is bad
  //if (false) { // don't try swap

      Values swap_col(csize+bosize, et0);
      Values swap_row(csize+bosize, et0);
  
      QP_sparse_vector<ET> y1(csize+bosize), z1(csize+bosize);
      QP_sparse_vector<ET> y2(csize+bosize), z2(csize+bosize);
  
      // TAG: TODO more efficient finding of suitable column, i.e., there can't
      // be no extraction of A_Ri in the following loop, because this makes it n^3
      // We have to extract the corresponding element individually.
      // TODO: replace perm_row[perm_column_inv[csize+k]] by stored value
      
      
      // find suitable swap column
      for (int j = bosize+csize-1; j >= 0; --j) {
        // check if pivots are alright
           if (tmp[ lu_fact_.perm_row(lu_fact_.perm_col_inv(j)) ] != et0) {
          test_val = et0;
          if (j < csize){
            if (lu_fact_.perm_row(lu_fact_.perm_col_inv(csize+k)) >= csize) {
              init__A_ij( test_val, j, B_O[ lu_fact_.perm_row(lu_fact_.perm_col_inv(csize+k)) - csize ], no_ineq);
            }
          } else {
            if (lu_fact_.perm_row(lu_fact_.perm_col_inv(csize+k)) < csize) {
              init__A_ij(test_val, lu_fact_.perm_row(lu_fact_.perm_col_inv(csize+k)), B_O[j-csize], no_ineq);
            } else {
              init__2_D_ij(test_val, B_O[lu_fact_.perm_row(lu_fact_.perm_col_inv(csize+k))-csize], B_O[j-csize], Tag_false()/*quadratic case*/);
            }
          }
          
          //CGAL_qpe_assertion(test_val == test_col[ perm_row[perm_column_inv[csize+k]] ]);
          
          if (test_val != et0) {
            swap_index_col = j;
            break;
          }
        }
      }
      
      if (swap_index_col < 0) {
        lu_fact_.set_invalid();
        return ret; // bail out
      } else if (swap_index_col < csize) { // get swap column from left part
        init__A_Ri(swap_col.begin()+csize, C[swap_index_col], no_ineq);
      } else { // get swap column from right part
        ratio_test_init__A_Cj(swap_col.begin(), B_O[swap_index_col-csize], no_ineq); //update acj
        ratio_test_init__2_D_Bj(swap_col.begin()+csize, B_O[swap_index_col-csize], Tag_false()/*is linear?*/);
      }
      swap_col[csize+k] = old_column[swap_index_col]; // we retrieved the new element already... so change back to old value

      // set up new column in place of the swap column
      z1.set_entry(swap_index_col, et1);
      for (int row = 0; row < csize+bosize; ++row){
        y1.set_entry(row, tmp[row] - swap_col[row]);
      }
      
      ret = lu_fact_.rank_1_update(et1, y1, z1);
      
      if (!ret) {
        lu_fact_.set_invalid();
        return ret; // bail out
      }
    
      // set up swap column in new place
      z2.set_entry(csize+k, et1);
      for (int row = 0; row < csize+bosize; ++row) {
        y2.set_entry(row, swap_col[row] - old_column[row]);
      }
      
      ret = lu_fact_.rank_1_update(et1, y2, z2);
      
      if (!ret) {
        lu_fact_.set_invalid();
        return ret; // bail out
      }

      // adjust internal permuations
      lu_fact_.swap_columns_physically(csize+k, swap_index_col);

      // do row update
      //if (tmp[ perm_column[perm_row_inv[csize+k]] ] == et0) {
      if (tmp[ lu_fact_.perm_col(lu_fact_.perm_row_inv(csize+k)) ] == et0) {
      
        
        QP_sparse_vector<ET> y3(csize+bosize), z3(csize+bosize);
        QP_sparse_vector<ET> y4(csize+bosize), z4(csize+bosize);
        
        
        // find suitable swap row
        for (int i = bosize+csize-1; i >= 0; --i) {
          test_val = et0;
          // check if pivots are alright
          if (tmp[ lu_fact_.perm_col(lu_fact_.perm_row_inv(i)) ] != et0) {
            if (i < csize){
              if (lu_fact_.perm_col(lu_fact_.perm_row_inv(csize+k)) >= csize) {
                init__A_ij(test_val, i, B_O[ lu_fact_.perm_col(lu_fact_.perm_row_inv(i))-csize ] , no_ineq);
              }
            } else {
              if (lu_fact_.perm_col(lu_fact_.perm_row_inv(csize+k)) < csize) {
                init__A_ij(test_val, lu_fact_.perm_col(lu_fact_.perm_row_inv(csize+k)), B_O[i-csize] , no_ineq);
              } else {
                init__2_D_ij(test_val, B_O[i-csize], B_O[lu_fact_.perm_col(lu_fact_.perm_row_inv(csize+k))-csize] , Tag_false());
              }
            }
            if (test_val != et0) {
              swap_index_row = i;
              break;
            }
          }
        }
        
        
        if (swap_index_row < 0) {
          lu_fact_.set_invalid();
          return ret; // bail out
        } else if (swap_index_row < csize) { // get swap column from left part
          init__A_Ri(swap_row.begin()+csize, swap_index_row, no_ineq);
        } else { // get swap column from right part
          ratio_test_init__A_Cj(swap_row.begin(), B_O[swap_index_row-csize], no_ineq); //update acj
          ratio_test_init__2_D_Bj(swap_row.begin()+csize, B_O[swap_index_row-csize], Tag_false()/*is linear?*/);
        }
        
        // set up new row in place of the swap row
        y3.set_entry(swap_index_row, et1);
        for (int col = 0; col < csize+bosize; ++col){
          z3.set_entry(col, tmp[col] - swap_row[col]);
        }
        z3.set_entry(csize+k, tmp[csize+k]-tmp[swap_index_row]);
        ret = lu_fact_.rank_1_update(et1, y3, z3);
        
        if (!ret) {
          lu_fact_.set_invalid();
          return ret; // bail out
        }
        
        // set up swap row in new place
        y4.set_entry(csize+k, et1);
        for (int col = 0; col < csize+bosize; ++col) {
          z4.set_entry(col, swap_row[col] - old_column[col]);
        }
        z4.set_entry(csize+k, swap_row[csize+k]-tmp[csize+k]);
        
        ret = lu_fact_.rank_1_update(et1, y4, z4);
        
        if (!ret) {
          lu_fact_.set_invalid();
          return ret; // bail out
        }
        
        // adjust internal permuations
        lu_fact_.swap_rows_physically(csize+k, swap_index_row);
      } else { // row update pivot alright
    
        y.set_entry(csize+k, et0);
        
        ret = lu_fact_.rank_1_update(et1, z, y);
        
        if (!ret) {
          lu_fact_.set_invalid();
          return ret; // bail out
        }
      }
  } else { // pivot is ok
    
    
    ret = lu_fact_.rank_1_update(et1, y, z);
    
    if (!ret) {
      lu_fact_.set_invalid();
      return ret;
    }
    
    // do second upate
    y.set_entry(csize+k, et0);
    ret = lu_fact_.rank_1_update(et1, z, y);
    
    if (!ret) {
      lu_fact_.set_invalid();
      return ret;
    }
  }
  return ret;
}

// PRE: index_entering and index_leaving are slack variables
template <typename Q, typename ET, typename Tags>
bool
QP_solver<Q,ET,Tags>::update_basis_matrix_LP_S_S(unsigned int k, Values old_row) {

  CGAL_qpe_assertion(old_row.size() == B_O.size());
  CGAL_qpe_assertion(is_phaseI || Is_linear::value);
  
  int csize(C.size()), bosize(B_O.size());
  
  QP_sparse_vector<ET> y(csize), z(bosize);
  Values tmp(bosize);
  
  y.set_entry(k, et1);
  
  // construct update column z
  init__A_Ri(tmp.begin(), slack_A[ index_leaving-qp_n].first, no_ineq); //update ari
  for (int col = 0; col < B_O.size(); ++col){
        z.set_entry(col, tmp[col] - old_row[col]);
  }
  return lu_fact_.rank_1_update(et1, y, z);
}

// PRE: index_entering is original and index_leaving is slack
template <typename Q, typename ET, typename Tags>
bool
QP_solver<Q,ET,Tags>::update_basis_matrix_LP_O_S() {

  CGAL_qpe_assertion(is_phaseI || Is_linear::value);
  
  bool ret = false;
  int swap_index = -1;
  int csize(C.size()), bosize(B_O.size());
  Permutation perm_row, perm_column_inv;
  
  QP_sparse_vector<ET> y1(csize), y2(csize), z1(bosize), z2(bosize);
  Values tmp_col(csize), tmp_row(bosize);
  Values swap_col(csize);
  
  // enlarge lu factorization (in rare cases this might not work, e.g. in the first
  // iteration of phase II, when the basis matrix is not valid initially)
  // TAG: TODO try different enlargement strategies for a and k in enlarge
  ret = lu_fact_.enlarge(bosize-1, bosize-1, true /*dummy*/);
  if (!ret) {
    lu_fact_.set_invalid();
    return ret; // bail out
  }
    

  // get new row for basis matrix
  init__A_Ri(tmp_row.begin(), slack_A[ index_leaving-qp_n].first, no_ineq); //update ari

  // get new column for basis matrix
  ratio_test_init__A_Cj(tmp_col.begin(), index_entering, no_ineq); //update acj
  
  CGAL_qpe_assertion(tmp_col[csize-1] == tmp_row[bosize-1]);
  
  if (tmp_col[csize-1] == et0) { // pivot is 0
  
      QP_sparse_vector<ET> y3(csize), z3(bosize);
      // get permutations
      lu_fact_.get_perm_row(perm_row);
      lu_fact_.get_perm_col_inv(perm_column_inv);
      
      // find suitable swap column
      for (int i = bosize-2; i >= 0; --i) {
        // check if pivots are alright
        if (tmp_row[i] != et0 && tmp_col[ perm_row[perm_column_inv[i]] ] != et0) {
          swap_index = i;
          break;
        }
      }
      if (swap_index < 0) {
        lu_fact_.set_invalid();
        return ret; // bail out
      }
      
      // get swap column
      ratio_test_init__A_Cj(swap_col.begin(), B_O[swap_index], no_ineq);

      // set up new column in place of suitable column
      z1.set_entry(swap_index, et1);
      for (int row = 0; row < csize-1; ++row) {
        y1.set_entry(row, -swap_col[row]+tmp_col[row]);
      }
      y1.set_entry(csize-1, tmp_col[csize-1]);
      ret = lu_fact_.rank_1_update(et1, y1, z1);
      
      if (!ret) {
        lu_fact_.set_invalid();
        return ret; // bail out
      }
      
      // set up last column
      z2.set_entry(bosize-1, et1);
      for (int row = 0; row < csize; ++row) {
        y2.set_entry(row, swap_col[row]);
      }
      y2.set_entry(csize-1, y2.get_entry(csize-1)-et1);
      ret = lu_fact_.rank_1_update(et1, y2, z2);
      if (!ret) {
        lu_fact_.set_invalid();
        return ret; // bail out
      }
      // TODO: put et0s in z3 instead of tmp_row...
      // do row update
      tmp_row[csize-1] = et0;
      tmp_row[swap_index] = et0;
      y3.set_entry(csize-1, et1);
      for (int col = 0; col < bosize; ++col){
        z3.set_entry(col, tmp_row[col]);
      }
      ret = lu_fact_.rank_1_update(et1, y3, z3);
      
      if (!ret) {
        lu_fact_.set_invalid();
        return ret; // bail out
      }
      
      // adjust internal permuations
      lu_fact_.swap_columns_physically(bosize-1, swap_index);
  } else { // pivot is ok
    
    // do row update
    y1.set_entry(bosize-1, et1);
    
    // construct update row z1
    for (int col = 0; col < bosize; ++col){
      z1.set_entry(col, tmp_row[col]);
    }
    z1.set_entry(bosize-1, z1.get_entry(bosize-1)-et1);
    ret = lu_fact_.rank_1_update(et1, y1, z1);
    
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
    // do column update
    z2.set_entry(bosize-1, et1);
    
    // construct update column y2
    for (int row = 0; row < csize; ++row) {
      y2.set_entry(row, tmp_col[row]);
    }
    y2.set_entry(bosize-1, et0);
    ret = lu_fact_.rank_1_update(et1, y2, z2);
  }
  
  if (!ret) lu_fact_.set_invalid();
  return ret;
}


// PRE: index_entering is original and index_leaving is slack.
//      B_O.size() and C.size() have already been reduced by 1, i.e., old_row
//      and old_col have size B_O.size()+1.
template <typename Q, typename ET, typename Tags>
bool
QP_solver<Q,ET,Tags>::update_basis_matrix_LP_S_O(unsigned int old_row_index, Values old_row, unsigned int old_col_index, Values old_col, Values last_row, Values last_col) {
  
  CGAL_qpe_assertion(is_phaseI || Is_linear::value);
  bool ret(false);
  int csize(C.size()+1), bosize(B_O.size()+1); // correct for reduced size
  QP_sparse_vector<ET> y1(csize), y2(csize), z1(bosize), z2(bosize);
  int index;
  
  if ( (index = lu_fact_.perm_col_inv(bosize-1)) == lu_fact_.perm_row_inv(csize-1) ) { // pivot permuation symmetric
    // reset last column/row to unit col/row
    // do row update
    y1.set_entry(csize-1, et1);
    
    // construct update row z1
    for (int col = 0; col < bosize; ++col){
      z1.set_entry(col, -last_row[col]);
    }
    z1.set_entry(bosize-1, et1 - last_row[bosize-1]);
    ret = lu_fact_.rank_1_update(et1, y1, z1);
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
     
    
    // do column update
    z2.set_entry(bosize-1, et1);
    
    // construct update column y2
    for (int row = 0; row < csize; ++row) {
      y2.set_entry(row, -last_col[row]);
    }
    y2.set_entry(csize-1, et0);
    ret = lu_fact_.rank_1_update(et1, y2, z2);
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
     
    
    // install last column/row in place of column/row that will be dropped
    y1.clear();
    y2.clear();
    z1.clear();
    z2.clear();
    
    // do row update
    y1.set_entry(old_row_index, et1);
    
    // construct update row z1
    for (int col = 0; col < bosize; ++col){
      z1.set_entry(col, last_row[col]-old_row[col]);
    }
    z1.set_entry(bosize-1, et0);
    z1.set_entry(old_col_index, last_row[bosize-1]-old_row[old_col_index]);
    ret = lu_fact_.rank_1_update(et1, y1, z1);
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
    
    
    // do column update
    z2.set_entry(old_col_index, et1);
    
    // construct update column y2
    for (int row = 0; row < csize; ++row) {
      y2.set_entry(row, last_col[row]-old_col[row]);
    }
    y2.set_entry(old_row_index, et0);
    y2.set_entry(csize-1, et0);
    ret = lu_fact_.rank_1_update(et1, y2, z2);
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
    ret = lu_fact_.shrink(index, true/*dummy*/);
  } else { // pivot permutation not symmetric
    // apply swap trick
  }
    
  if (!ret) {
    lu_fact_.set_invalid();
  }
  
  return ret;
  
}

template <typename Q, typename ET, typename Tags>
bool
QP_solver<Q,ET,Tags>::update_basis_matrix_QP_S_O(unsigned int old_constraint_row_index, Values old_constraint_row, unsigned int old_variable_row_index, Values old_variable_row, unsigned int last_constraint_row_index, Values last_constraint_row, Values last_variable_row) {
  
  CGAL_qpe_assertion(!is_phaseI && !Is_linear::value);
  bool ret(false);
  int csize(C.size()+1), bosize(B_O.size()+1); // correct for reduced size
  QP_sparse_vector<ET> y1(csize+bosize), y2(csize+bosize), z1(csize+bosize), z2(csize+bosize);
  
  int index1, index2;
  
  
  if ( (index1 = lu_fact_.perm_col_inv(last_constraint_row_index)) == lu_fact_.perm_row_inv(last_constraint_row_index) &&
       (index2 = lu_fact_.perm_col_inv(bosize+csize-1)) == lu_fact_.perm_row_inv(bosize+csize-1)) { // pivot permuation symmetric
  
  
    /////////////////////////////
    // do constraint update first
  
    // reset last column/row to unit col/row
    // do row update
    y1.set_entry(last_constraint_row_index, et1);
    
    // construct update row z1
    for (int col = 0; col < csize+bosize; ++col){
      z1.set_entry(col, -last_constraint_row[col]);
    }
    z1.set_entry(last_constraint_row_index, et1 - last_constraint_row[last_constraint_row_index]);
    ret = lu_fact_.rank_1_update(et1, y1, z1);
    
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }

    // do column update
    z1.set_entry(last_constraint_row_index, et0);
    ret = lu_fact_.rank_1_update(et1, z1, y1);
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
    // install last constraint column/row in place of column/row that will be dropped
    if (old_constraint_row_index != last_constraint_row_index) {
      y1.clear();
      z1.clear();
      
      // do row update
      y1.set_entry(old_constraint_row_index, et1);
      
      // construct update row z1
      for (int col = 0; col < csize+bosize; ++col){
        z1.set_entry(col, last_constraint_row[col]-old_constraint_row[col]);
      }
      z1.set_entry(last_constraint_row_index, et0);
      z1.set_entry(old_constraint_row_index, last_constraint_row[last_constraint_row_index]-old_constraint_row[old_constraint_row_index]);
      ret = lu_fact_.rank_1_update(et1, y1, z1);
      if (!ret) {
        lu_fact_.set_invalid();
        return ret; // bail out
      }

      // do column update
      z1.set_entry(old_constraint_row_index, et0);
      ret = lu_fact_.rank_1_update(et1, z1, y1);
      if (!ret) {
        lu_fact_.set_invalid();
        return ret; // bail out
      }
    }
    
    
    
    ////////////////////////////
    // do variable update second

    // reflect changes of constraint swap
    old_variable_row[old_constraint_row_index] = last_constraint_row[csize+old_variable_row_index];
    old_variable_row[last_constraint_row_index] = et0;
    last_variable_row[old_constraint_row_index] = last_constraint_row[csize+bosize-1];
    last_variable_row[last_constraint_row_index] = et0;

    // reset last column/row to unit col/row
    // do row update
    y2.set_entry(csize+bosize-1, et1);
    
    // construct update row z2
    for (int col = 0; col < csize+bosize; ++col){
      z2.set_entry(col, -last_variable_row[col]);
    }
    z2.set_entry(csize+bosize-1, et1 - last_variable_row[csize+bosize-1]);
    ret = lu_fact_.rank_1_update(et1, y2, z2);
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
    // do column update
    z2.set_entry(csize+bosize-1, et0);
    ret = lu_fact_.rank_1_update(et1, z2, y2);
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    
    // install last constraint column/row in place of column/row that will be dropped
    if (old_variable_row_index != bosize-1) {
      y2.clear();
      z2.clear();
      
      // do row update
      y2.set_entry(csize+old_variable_row_index, et1);
      
      // construct update row z2
      for (int col = 0; col < csize+bosize; ++col){
        z2.set_entry(col, last_variable_row[col]-old_variable_row[col]);
      }
      z2.set_entry(csize+bosize-1, et0);
      z2.set_entry(csize+old_variable_row_index, last_variable_row[csize+bosize-1]-old_variable_row[csize+old_variable_row_index]);
      ret = lu_fact_.rank_1_update(et1, y2, z2);
      
      
      if (!ret) {
        lu_fact_.set_invalid();
        return ret; // bail out
      }
      
      // do column update
      z2.set_entry(csize+old_variable_row_index, et0);
      ret = lu_fact_.rank_1_update(et1, z2, y2);
      if (!ret) {
        lu_fact_.set_invalid();
        return ret; // bail out
      }
    }
    ret = lu_fact_.shrink((index2>index1 ? index2 : index1), true);
    if (!ret) {
      lu_fact_.set_invalid();
      return ret; // bail out
    }
    ret = lu_fact_.shrink((index2>index1 ? index1 : index2), false);  
  } else { // last constraint row permutation not symmetric
   // apply swap trick
  }
  if (!ret) {
    lu_fact_.set_invalid();
  }
  
  return ret;
  
}
  


// PRE: index_entering is original and index_leaving is slack
template <typename Q, typename ET, typename Tags>
bool
QP_solver<Q,ET,Tags>::update_basis_matrix_QP_O_S() {

  CGAL_qpe_assertion(!is_phaseI && !Is_linear::value);
  
  bool ret(false);
  int csize(C.size()), bosize(B_O.size());
  
  //Permutation perm_row, perm_column_inv;
  QP_sparse_vector<ET> y1(csize+bosize), y2(csize+bosize), z1(csize+bosize), z2(csize+bosize);
    
  // enlarge lu factorization (in rare cases this might not work, e.g. in the first
  // iteration of phase II, when the basis matrix is not valid initially)
  // TAG: TODO try different enlargement strategies for a and k in enlarge
  ret = lu_fact_.enlarge(csize+bosize-2, csize+bosize-2, true);
  if (!ret) {
    lu_fact_.set_invalid();
    return ret; // bail out
  }
  ret = lu_fact_.enlarge(csize-1, csize+bosize-1, false);
  if (!ret) {
    lu_fact_.set_invalid();
    return ret; // bail out
  }
  Values new_variable_row(csize+bosize, et0), new_constraint_row(csize+bosize, et0);
  
  //TODO: think about sparsifying the fetching procedures for those rows and columns

  // get new column for basis matrix
  ratio_test_init__A_Cj(new_variable_row.begin(), index_entering, no_ineq); //update acj
  ratio_test_init__2_D_Bj( new_variable_row.begin()+csize, index_entering, Tag_false());
    
  // do variable row update
  y1.set_entry(csize+bosize-1, et1);
  
  // construct update row z
  for (int col = 0; col < csize+bosize; ++col){
    z1.set_entry(col, new_variable_row[col]);
  }
  // TODO: check what happens if we subtract et1 in the following update not here
  z1.set_entry(csize+bosize-1, new_variable_row[csize+bosize-1]-et1);
  
  ret = lu_fact_.rank_1_update(et1, y1, z1);
  if (!ret) {
    lu_fact_.set_invalid();
    return ret; // bail out
  }

  // do column update
  z1.set_entry(csize+bosize-1, et0);
  ret = lu_fact_.rank_1_update(et1, z1, y1);
  if (!ret) {
    lu_fact_.set_invalid();
    return ret; // bail out
  }
  // get new row for basis matrix
  init__A_Ri(new_constraint_row.begin()+csize, slack_A[ index_leaving-qp_n].first, no_ineq); //update ari
  
  // do constraint row update
  y2.set_entry(csize-1, et1);
  
  // construct update row z
  for (int col = 0; col < csize+bosize; ++col){
    z2.set_entry(col, new_constraint_row[col]);
  }
  // TODO: check what happens if we subtract et1 in the following update not here
  z2.set_entry(csize-1, new_constraint_row[csize-1]);
  z2.set_entry(csize+bosize-1, et0);
  ret = lu_fact_.rank_1_update(et1, y2, z2);
  if (!ret) {
    lu_fact_.set_invalid();
    return ret; // bail out
  }
  // do column update
  z2.set_entry(csize-1, -et1);
  ret = lu_fact_.rank_1_update(et1, z2, y2);
  if (!ret) lu_fact_.set_invalid();
  return ret;
}


// PRE: index_entering is original and index_leaving is slack
template <typename Q, typename ET, typename Tags>
bool
QP_solver<Q,ET,Tags>::update_basis_matrix_QP_S_S(unsigned int k, Values old_column) {

  CGAL_qpe_assertion(k < C.size());
  CGAL_qpe_assertion(!is_phaseI && !Is_linear::value);
  
  bool ret(false);
  int csize(C.size()), bosize(B_O.size());
  
  QP_sparse_vector<ET> y(csize+bosize), z(csize+bosize);
  Values new_column(bosize);
  
  init__A_Ri(new_column.begin(), slack_A[ index_leaving-qp_n].first, no_ineq); //update ari
  
  // do row update
  y.set_entry(k, et1);
  
  // construct update column z
  for (int col = csize; col < csize+bosize; ++col){
        z.set_entry(col, new_column[col-csize]-old_column[col-csize]);
  }
  ret = lu_fact_.rank_1_update(et1, y, z);
  if (!ret) {
    lu_fact_.set_invalid();
    return ret;
  }
  ret = lu_fact_.rank_1_update(et1, z, y);
  if (!ret) {
    lu_fact_.set_invalid();
  }
  
  return ret;
}

} //namespace CGAL

// ===== EOF ==================================================================
