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

namespace CGAL {


template <typename ET, typename Is_LP, typename Matrix_Provider>
QP_LU_factorization<ET, Is_LP, Matrix_Provider>::QP_LU_factorization (
  CGAL::Verbose_ostream& verbose_ostream,
  Matrix_Provider* matrix_provider):
      valid_(false),
      et0_(0),
      et1_(1),
      et2_(2),
      n_(0),
      matrix_u_(),
      matrix_l_(),
      d_row_(),
      d_(ET(1)),
      d_neg_(false),
      csize_(0),
      bosize_(0),
      is_linear_( Is_LP::value ),
      vout_(verbose_ostream),
      perm_row_(),
      perm_column_(),
      perm_row_inv_(),
      perm_column_inv_(),
      matrix_provider_(matrix_provider)
{

  CGAL_qpe_assertion(matrix_provider_ != 0);
  
    
  matrix_u_ = matrix_provider->get_basis_matrix(csize_, bosize_, is_linear_ || is_phase_I_);
  CGAL_qpe_assertion( matrix_u_->get_size().first == matrix_u_->get_size().second);  // We require square matrices
  n_ = matrix_u_->get_size().first;
  
  perm_row_ = permutation_t(n_);
  perm_column_ = permutation_t(n_);
  perm_row_inv_ = permutation_t(n_);
  perm_column_inv_ = permutation_t(n_);
  
  d_row_ = dense_vector_t(n_);
  std::fill(d_row_.begin(), d_row_.end(), et1_);
                        
  // initialize the permutations with the identity
  std::generate(perm_row_.begin(), perm_row_.end(), Init_Natural_Numbers(0));
  std::generate(perm_column_.begin(), perm_column_.end(),
                Init_Natural_Numbers(0));
  std::generate(perm_row_inv_.begin(), perm_row_inv_.end(),
                Init_Natural_Numbers(0));
  std::generate(perm_column_inv_.begin(), perm_column_inv_.end(),
                Init_Natural_Numbers(0));            
}

// =========================================================================
// Getter & setter methods
// =========================================================================
    

// =========================================================================
// Debug methods
// =========================================================================

template <typename ET, typename Is_LP, typename Matrix_Provider>
void QP_LU_factorization<ET, Is_LP, Matrix_Provider>::test() {

  CGAL_qpe_debug {
    if (vout_.verbose()) {
      vout_.out() << "Testing..." << std::endl;
    }
  }
  
  compute_factorization(is_linear_ || is_phase_I_);
  
  CGAL_qpe_debug {
    if (vout_.verbose()) {
      vout_.out() << " Factorization computed.\n";
    }
  }
  
}

// =========================================================================
// Operations
// =========================================================================
// public

template <typename ET, typename Is_LP, typename Matrix_Provider>
void QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::set_invalid() {
  valid_ = false;
}

template <typename ET, typename Is_LP, typename Matrix_Provider>
std::ostream&
QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::print(std::ostream& out) const {
  if (valid_) {
    out << "Matrix U:\n" << *matrix_u_;
    out << "Matrix L:\n" << *matrix_l_;
    out << "Denominator: " << d_ << std::endl;
    
    out << "n_: " << n_ << std::endl;
    out << "bosize_: " << bosize_ << std::endl;
    out << "csize_: " << csize_ << std::endl;
    
    out << "Row permutation:\n(";
    std::copy(perm_row_.begin(), perm_row_.end()-1, std::ostream_iterator<int>(out, ", "));
    if (n_ > 0) out << perm_row_.at(n_-1) << ")\n";
    
    out << "Column permutation:\n(";
    std::copy(perm_column_.begin(), perm_column_.end()-1, std::ostream_iterator<int>(out, ", "));
    
    if (n_ > 0) out << perm_column_.at(n_-1) << ")\n";
    out << "Inverse row permutation:\n(";
    std::copy(perm_row_inv_.begin(), perm_row_inv_.end()-1, std::ostream_iterator<int>(out, ", "));
    if (n_ > 0) out << perm_row_inv_.at(n_-1) << ")\n";
    
    out << "Inverse column permutation:\n(";
    std::copy(perm_column_inv_.begin(), perm_column_inv_.end()-1, std::ostream_iterator<int>(out, ", "));
    if (n_ > 0) out << perm_column_inv_.at(n_-1) << ")\n";
    
  } else { // !valid_
    out << "invalid matrix!\n"; 
  }
  return out;
}


// private

//  getPivotElement(const int k) const
//  Wrapper function to call the actual pivot finder
template <typename ET, typename Is_LP, typename Matrix_Provider>
typename QP_LU_factorization<ET, Is_LP, Matrix_Provider>::index_pair_t
QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::getPivotElement(const int k)
const {
  //return getPivotElement_no_pivoting(k);
  //return getPivotElement_markowitz_extended(k);
  return getPivotElement_markowitz(k);
  //return getPivotElement_partial_pivoting(k);
}

// getPivotElement_no_pivoting(const int k) const 
// Simply chooses the diagonal element (k,k) as pivot.
// PRE: If element (k,k) is zero, then the subdiagonal
//      elements below it must be zero too.
template <typename ET, typename Is_LP, typename Matrix_Provider>
typename QP_LU_factorization<ET, Is_LP, Matrix_Provider>::index_pair_t
QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::getPivotElement_no_pivoting(const int k) const {
  return std::make_pair(k, k);
}

  
// getPivotElement_partial_pivoting(const int k) const
// Chooses the largest subdiagonal (or diagonal) element in column k as pivot.
template <typename ET, typename Is_LP, typename Matrix_Provider>
typename QP_LU_factorization<ET, Is_LP, Matrix_Provider>::index_pair_t
QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::getPivotElement_partial_pivoting(const int k) const {
  int index = k;
  
  typename QP_sparse_matrix<ET>::value_const_iterator it_begin = matrix_u_->begin_column_const_value(perm_column_[k]);
  typename QP_sparse_matrix<ET>::value_const_iterator it_end = matrix_u_->end_column_const_value(perm_column_[k]);
  dense_vector_t temp(n_, et0_);
  // TAG: TODO handle permutation nicer O(log n) vs. O(n)
  while (it_begin != it_end) {
    temp[perm_row_inv_[it_begin->first]] = *it_begin->second;
    ++it_begin;
  }
  
  ET min_abs = et0_; 
  for (int i = k; i < n_; ++i) {
    // TAG: TODO not just the first one (e.g. largest one)
    if (temp[i] != et0_) {
      if ((temp[i]>et0_ ? temp[i] : -temp[i]) > min_abs) {
        index = perm_row_[i];
        min_abs = temp[i];
      }
    }
  }
  return std::make_pair(index, perm_column_[k]);
}

// getPivotElement_markowitz(const int k) const
// Chooses the the element that minimizes the Markowitz count from the
// remaining matrix.
template <typename ET, typename Is_LP, typename Matrix_Provider>
typename QP_LU_factorization<ET, Is_LP, Matrix_Provider>::index_pair_t
QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::getPivotElement_markowitz(const int k) const {

  int index = k;
     
  std::vector<unsigned int> row_count(n_-k, -1);
  std::vector<unsigned int> column_count(n_-k, -1);

  
  typename QP_sparse_matrix<ET>::value_const_iterator it_begin;
  typename QP_sparse_matrix<ET>::value_const_iterator it_end;
  for (int i = k; i < n_; ++i) {
    it_begin = matrix_u_->begin_column_const_value(perm_column_[i]);
    it_end  = matrix_u_->end_column_const_value(perm_column_[i]);
    while (it_begin != it_end) {
      if ( (index = perm_row_inv_[it_begin->first]) >= k && *it_begin->second != 0) {
        ++row_count[index-k];
        ++column_count[i-k];
      }
      ++it_begin;
    }
  }
  
  int index_row = k;
  int index_column = k;
  unsigned int temp = 0;
  unsigned int min = static_cast<unsigned int>(-1);
  for (int i = 0; i < n_-k; ++i) {
    for (int j = 0; j < n_-k; ++j) {
        if ((temp = row_count[i]*column_count[j]) < min && matrix_u_->get_entry(perm_row_[i+k], perm_column_[j+k]) != 0) {
          min = temp;
          index_row = i;
          index_column = j;
        }
    }
  }
  return std::make_pair(perm_row_[index_row+k], perm_column_[index_column+k]);
}


// getPivotElement_markowitz(const int k) const
// Chooses the the element that minimizes the "real" Markowitz count from the
// remaining matrix.
template <typename ET, typename Is_LP, typename Matrix_Provider>
typename QP_LU_factorization<ET, Is_LP, Matrix_Provider>::index_pair_t
QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::getPivotElement_markowitz_extended(const int k) const {
  
  int index = k;
     
  std::vector<unsigned int> row_count(n_-k, 0);
  std::vector<unsigned int> column_count(n_-k, 0);
  
  QP_sparse_matrix<unsigned int> pattern(n_-k, n_-k);
  
  typename QP_sparse_matrix<ET>::value_const_iterator it_begin;
  typename QP_sparse_matrix<ET>::value_const_iterator it_end;
  for (int i = k; i < n_; ++i) {
    it_begin = matrix_u_->begin_column_const_value(perm_column_[i]);
    it_end  = matrix_u_->end_column_const_value(perm_column_[i]);
    while (it_begin != it_end) {
      if ( (index = perm_row_inv_[it_begin->first]) >= k && *it_begin->second != 0) {
        pattern.set_entry(index-k, i-k, 1);
      }
      ++it_begin;
    }
  }
  
  QP_sparse_matrix<unsigned int> pattern_2(pattern);  
  pattern_2 = *pattern.matrix_matrix_multiplication(pattern_2);

  typename QP_sparse_matrix<unsigned int>::value_const_iterator it_pattern_begin;
  typename QP_sparse_matrix<unsigned int>::value_const_iterator it_pattern_end;
  int index_row = k;
  int index_column = k;
  unsigned int temp = 0;
  unsigned int min = static_cast<unsigned int>(-1);
  for (int i = 0; i < n_-k; ++i) {
    it_pattern_begin = pattern.begin_column_const_value(i);
    it_pattern_end = pattern.end_column_const_value(i);
    while (it_pattern_begin != it_pattern_end) {
      if ((temp = pattern_2.get_entry(it_pattern_begin->first, i)) < min) {
        min = temp;
        index_row = it_pattern_begin->first;
        index_column = i;
      }
      ++it_pattern_begin;
    }
  }
  return std::make_pair(perm_row_[index_row+k], perm_column_[index_column+k]);
}
 

// compute_factorization()
// Wrapper function to call the actual factorization routine.
template <typename ET, typename Is_LP, typename Matrix_Provider>
void QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::compute_factorization(bool is_linear) {
  
  //CGAL_qpe_assertion(matrix_provider_ != 0);
  
  if (matrix_provider_ != 0) { // TAG: ARTIFACT
    
    matrix_u_ = matrix_provider_->get_basis_matrix(csize_, bosize_, is_linear);
    
    // We require square matrices
    CGAL_qpe_assertion( matrix_u_->get_size().first == matrix_u_->get_size().second);
    n_ = matrix_u_->get_size().first;
    
    perm_row_ = permutation_t(n_);
    perm_column_ = permutation_t(n_);
    perm_row_inv_ = permutation_t(n_);
    perm_column_inv_ = permutation_t(n_);
    d_row_ = dense_vector_t(n_, et1_);
    
    if (n_ <= 0) {
      valid_ = false;
      return;
    }
    
    // initialize the permutations with the identity
    std::generate(perm_row_.begin(), perm_row_.end(), Init_Natural_Numbers(0));
    std::generate(perm_column_.begin(), perm_column_.end(),
                  Init_Natural_Numbers(0));
    std::generate(perm_row_inv_.begin(), perm_row_inv_.end(),
                  Init_Natural_Numbers(0));
    std::generate(perm_column_inv_.begin(), perm_column_inv_.end(),
                  Init_Natural_Numbers(0));    
    
    CGAL_qpe_debug {
      if (vout_.verbose()) {
        vout_.out() << "==========================================" << std::endl;
        vout_.out() << "We are in LU compute_factorization" << std::endl;
        vout_.out() << "csize_: " << csize_ << std::endl;
        vout_.out() << "bosize_: " << bosize_ << std::endl;
        vout_.out() << "n_: " << n_ << std::endl;
        vout_.out() << "Original matrix:\n" << *matrix_u_ << std::endl;
       }
     }
    
    compute_factorization_gauss();
    valid_ = true;
    
    CGAL_qpe_debug {
      if (vout_.verbose()) {
        vout_.out() << "Factorization output:\n";
        print(vout_.out());
        vout_.out() << "==========================================" << std::endl;
      }
    }
    
  }
}

// TODO: improved algorithm according to write-up
template <typename ET, typename Is_LP, typename Matrix_Provider>
void QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::compute_factorization_gauss() {


  CGAL_qpe_debug {
    if (vout_.verbose()) {
      vout_.out() << "||========================================" << std::endl;
      vout_.out() << "BEGINNING of compute_factorization_gauss" << std::endl;
      print(vout_.out());
      vout_.out() << "\nn: " << n_;
      vout_.out() << "\nDensity BEFORE factorization: " << matrix_u_->get_density() << std::endl;
      CGAL::QP_solver_debug::timer.compute_factorization.start();
    }
  }

  // Compute FTRAN and BTRAN factorization
  matrix_l_ = QP_sparse_matrix<ET>::get_identity(n_);
  
  typename QP_sparse_matrix<ET>::value_iterator it_column, it_column_end;
  typename QP_sparse_matrix<ET>::value_iterator it_row, it_row_end;
  typename QP_sparse_matrix<ET>::index_iterator it_column_index, it_column_index_end;
  typename QP_sparse_matrix<ET>::index_iterator it_row_index, it_row_index_end;
  
  std::vector<std::pair<typename vector_t::key_type, typename vector_t::value_type> > row, column;
  typedef typename std::vector<std::pair<typename vector_t::key_type, typename vector_t::value_type> >::iterator Iterator;
 
 
  // Fill the left hand part of matrix_u_
  for (int k = 1; k < n_; ++k) {
      
    index_pair_t pivot_index = getPivotElement(k-1);
  
    d_row_[k] = matrix_u_->get_entry(pivot_index);
    
    if (perm_row_inv_[pivot_index.first] > k-1) swap_rows(perm_row_inv_[pivot_index.first], k-1);
    if (perm_column_inv_[pivot_index.second] > k-1) swap_columns(perm_column_inv_[pivot_index.second], k-1);
       
    
    // premultiply with d_row_[k]
    for (int i = k; i < n_; ++i) {
      it_column = matrix_u_->begin_column_value(perm_column_[i]);
      it_column_end = matrix_u_->end_column_value(perm_column_[i]);
      
      while (it_column != it_column_end) {     
        if ( perm_row_inv_[it_column->first] >= k) {
          *it_column->second *= d_row_[k];
        }
        ++it_column;
      }
    }
    
    // get pivot column
    column.clear();
    it_column = matrix_u_->begin_column_value(perm_column_[k-1]);
    it_column_end = matrix_u_->end_column_value(perm_column_[k-1]);
    while (it_column != it_column_end) {
      if ( perm_row_inv_[it_column->first] >= k) {
        column.push_back(std::make_pair(it_column->first, *it_column->second));
      }
      ++it_column;
    }
    
    // get pivot row
    row.clear();
    it_row = matrix_u_->begin_row_value(perm_row_[k-1]);
    it_row_end = matrix_u_->end_row_value(perm_row_[k-1]);
    while (it_row != it_row_end && it_row->first < n_) { // better stopping criterion? This relies on at least one other entry following in the right part.
      if (perm_column_inv_[it_row->first] >= k) {
        row.push_back(std::make_pair(it_row->first, *it_row->second));
      }
      ++it_row;
    }
        
    for (Iterator it_r = column.begin(); it_r != column.end(); ++it_r) {
      for (Iterator it_c = row.begin(); it_c != row.end(); ++it_c) {
        matrix_u_->set_entry(it_r->first, it_c->first, matrix_u_->get_entry(it_r->first, it_c->first) - it_r->second * it_c->second);
      }
    }

    CGAL_qpe_assertion(d_row_[k-1] != 0);
    
    // divide by d_row_[k-1]
    for (int i = k; i < n_; ++i) {
      it_column = matrix_u_->begin_column_value(perm_column_[i]);
      it_column_end = matrix_u_->end_column_value(perm_column_[i]);
      while (it_column != it_column_end) {
        if ( perm_row_inv_[it_column->first] >= k) {
          
          CGAL_qpe_debug{
            ++CGAL::QP_solver_debug::timer.counter_integral_division;
          }
          
          *it_column->second = integral_division(*it_column->second, d_row_[k-1]);
        }
        ++it_column;
      }
    }
  }
    
    
  // Fill the right hand part of matrix_u_. It's important that this
  // is done in a separate loop after the previous one, bacause the
  // permutation has to be fully known (possibly this may be done
  // more efficiently TODO YB)
  for (int k = 1; k < n_; ++k) {
    
    // premultiply with d_row_[k]
    for (int i = k; i < n_; ++i) {
      it_column = matrix_l_->begin_column_value(i);
      it_column_end = matrix_l_->end_column_value(i);
      while (it_column != it_column_end && it_column->first < k) {
        *it_column->second *= d_row_[k];
        ++it_column;
      }
      it_row = matrix_l_->begin_row_value(i);
      it_row_end = matrix_l_->end_row_value(i);
      while (it_row != it_row_end && it_row->first < k) {
        *it_row->second *= d_row_[k];
        ++it_row;
      }
    }
    
    
    // process diagonal
    for (int i = k; i < n_; ++i) {
      matrix_l_->set_entry(i, i, matrix_l_->get_entry(i, i) * d_row_[k]);
      
      CGAL_qpe_debug{
        ++CGAL::QP_solver_debug::timer.counter_integral_division;
      }
      
      matrix_l_->set_entry(i, i, integral_division(matrix_l_->get_entry(i,i), d_row_[k-1]) );
    }
    
    
    // get pivot column
    column.clear();
    it_column = matrix_u_->begin_column_value(perm_column_[k-1]);
    it_column_end = matrix_u_->end_column_value(perm_column_[k-1]);
    while (it_column != it_column_end) {
      if ( perm_row_inv_[it_column->first] >= k) {
        column.push_back(std::make_pair(it_column->first, *it_column->second));
      }
      ++it_column;
    }
    
    it_row_end = matrix_l_->end_row_value(k-1);
    for (Iterator it = column.begin(); it != column.end(); ++it) {
      it_row = matrix_l_->begin_row_value(k-1);
      while (it_row != it_row_end && it_row->first < k) { // TODO: better stopping criterion? This relies on at least one other entry following in the right part.
        matrix_l_->set_entry(perm_row_inv_[it->first], it_row->first, matrix_l_->get_entry(perm_row_inv_[it->first], it_row->first) - it->second * *it_row->second );
        ++it_row;
      }
    }
    
    // get pivot row
    row.clear();
    it_row = matrix_u_->begin_row_value(perm_row_[k-1]);
    it_row_end = matrix_u_->end_row_value(perm_row_[k-1]);
    while (it_row != it_row_end && it_row->first < n_) { // TODO: get rid of < n_
      if (perm_column_inv_[it_row->first] >= k) {
        row.push_back(std::make_pair(it_row->first, *it_row->second) );
      }
      ++it_row;
    }
    
    // fill upper part
    it_column_end = matrix_l_->end_column_value(k-1);
    for (Iterator it = row.begin(); it != row.end(); ++it) {
      it_column = matrix_l_->begin_column_value(k-1);
      while (it_column != it_column_end && it_column->first < k) { // better stopping criterion? This relies on the diagonal of the right part being non-zero (which we may assume in this implementation)
        matrix_l_->set_entry(it_column->first, perm_column_inv_[it->first], matrix_l_->get_entry(it_column->first, perm_column_inv_[it->first]) - it->second * *it_column->second);
        ++it_column;
      }
    }
    
    // divide by d_row_[k-1]
    
    for (int i = 0; i < k; ++i) {
      it_column_index = matrix_l_->get_column_indices(i).lower_bound(k);
      it_column_index_end = matrix_l_->get_column_indices(i).end();
      while (it_column_index != it_column_index_end) {
              
        CGAL_qpe_debug{
          ++CGAL::QP_solver_debug::timer.counter_integral_division;
        }
        
        matrix_l_->set_element_by_index(it_column_index->second,
          integral_division(matrix_l_->get_element_by_index(it_column_index->second), d_row_[k-1]) );
        ++it_column_index;
      }
    }
    
    for (int i = k; i < n_; ++i) {
      it_column = matrix_l_->begin_column_value(i);
      it_column_end = matrix_l_->end_column_value(i);
      while (it_column != it_column_end && it_column->first < k) {
      
        CGAL_qpe_debug{
          ++CGAL::QP_solver_debug::timer.counter_integral_division;
        }
        
        *it_column->second = integral_division(*it_column->second, d_row_[k-1]);
        ++it_column;
      }
    }
  }  // loop over k
  
  d_row_.push_back(matrix_u_->get_entry(perm_row_[n_-1], perm_column_[n_-1]));

  if (d_row_[n_] > 0) {
    d_ = d_row_[n_];
    d_neg_ = false;
  } else {
    d_ = -d_row_[n_];
    d_neg_ = true; 
  }
  
  CGAL_qpe_debug {
    if (vout_.verbose()) {
      vout_.out() << "END of compute_factorization_gauss" << std::endl;
      vout_.out() << "Matrix U:\n" << *matrix_u_ << std::endl;
      vout_.out() << "Denominator: " << d_ << std::endl;
      vout_.out() << "Density AFTER factorization: " << matrix_u_->get_density() << std::endl;
      CGAL::QP_solver_debug::timer.compute_factorization.stop();
      vout_.out() << "||========================================" << std::endl;
    }
  }
}

// TAG: TODO make template switch
template <typename ET, typename Is_LP, typename Matrix_Provider>
void QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::swap_rows(int i, int j) {
  CGAL_qpe_assertion(i < n_  && i >= 0 && j < n_ && j >= 0);
  std::swap(perm_row_[i], perm_row_[j]);
  std::swap(perm_row_inv_[perm_row_[i]], perm_row_inv_[perm_row_[j]]);
}


// TAG: TODO make template switch
template <typename ET, typename Is_LP, typename Matrix_Provider>
void QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::swap_columns(int i, int j) {
  CGAL_qpe_assertion(i < n_  && i >= 0 && j < n_ && j >= 0);
  std::swap(perm_column_[i], perm_column_[j]);
  std::swap(perm_column_inv_[perm_column_[i]],
            perm_column_inv_[perm_column_[j]]);
}

// TAG: TODO make template switch
template <typename ET, typename Is_LP, typename Matrix_Provider>
void QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::swap_columns_physically(int i, int j) {
  CGAL_qpe_assertion(i < n_  && i >= 0 && j < n_ && j >= 0);
  std::swap(perm_column_inv_[i], perm_column_inv_[j]);
  std::swap(perm_column_[perm_column_inv_[i]],
            perm_column_[perm_column_inv_[j]]);
  matrix_u_->swap_columns(i, j);
}

// TAG: TODO make template switch
template <typename ET, typename Is_LP, typename Matrix_Provider>
void QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::swap_rows_physically(int i, int j) {
  CGAL_qpe_assertion(i < n_  && i >= 0 && j < n_ && j >= 0);
  std::swap(perm_row_inv_[i], perm_row_inv_[j]);
  std::swap(perm_row_[perm_row_inv_[i]],
            perm_row_[perm_row_inv_[j]]);
  matrix_u_->swap_rows(i, j);
}

// ATTENTION: This is for debugging purposes only. It is highly inefficient,
// because the involved numbers grow hugely.
template <typename ET, typename Is_LP, typename Matrix_Provider>
boost::shared_ptr<QP_sparse_matrix<ET> >
QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::recover_original_matrix() {

  boost::shared_ptr<QP_sparse_matrix<ET> > ret;

  if (!valid_) {
    compute_factorization(is_linear_ || is_phase_I_);
  }
    
  int dim = ((is_linear_ || is_phase_I_) ? bosize_ : bosize_ + csize_);
  ret = boost::shared_ptr<QP_sparse_matrix<ET> >(new QP_sparse_matrix<ET>(dim, dim, et0_));

  QP_sparse_matrix<ET> triu(dim, dim, et0_);
  QP_sparse_matrix<ET> tril(dim, dim, et0_);

  ET divisor(et1_);
  for (int i = 0; i <= dim; ++i) divisor *= d_row_[i];
  
  for (int i=0; i < dim; ++i) {
    for (int j=0; j < dim; ++j) {
      if (i <= j) {
        triu.set_entry( i, j, matrix_u_->get_entry(perm_row_[i], perm_column_[j]) );
      }
    }
  }

  
  for (int k = 0; k < dim; ++k) {
    for (int i = k; i < dim; ++i) {
      ET tmp = (i == k ? divisor : et0_);
      for (int j = k; j < i; ++j) { // k instead of 0
        tmp -= tril.get_entry(j, k) * ( matrix_l_->get_entry(i, j) );
        tmp *= d_row_[j+1];
      }
      tril.set_entry(i, k, integral_division(tmp, matrix_l_->get_entry(i, i)) );
    }
    // I think there's something unnecessary about the following double loop. I use
    // it to get rid of surplus factors, which could probably be handled in the previous loop
    // or by not starting with such large numbers to begin with. However, it does not change the
    // fact that there is no way we can avoid going to divisor in number size, because
    // we have to invert the Linv... too bad... that's why it's for debugging purposes.
    for (int i = k+1; i < dim; ++i) {
      for (int j = i; j < dim; ++j) { // k instead of 0
          tril.set_entry(j, k, integral_division(tril.get_entry(j, k), matrix_l_->get_entry(i, i)) );
      }
    }
  }
  
  
  QP_sparse_matrix<ET> temp_matrix(*(tril*triu));
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      ret->set_entry(perm_row_[i], perm_column_[j], integral_division( temp_matrix.get_entry(i,j), divisor) );
    }
  }

  return ret;
}

// public solve FTRAN ( LU y = v )
template <typename ET, typename Is_LP, typename Matrix_Provider>
template < typename ForIt, typename OutIt >
void
QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::solve( ForIt v_l_it, ForIt v_x_it,
         OutIt y_l_it,  OutIt y_x_it,
         bool is_linear, bool is_phase_I) {
         
  // TODO: CODE DESIGN move tags to other place (not as parameters)
  // Also: There's a potential problem about the consistency of is_linear_
  // and phase tags. They are only set here at this point, when solve is called.
  // However, in some rare cases the recomputation of the factorization is triggered
  // in routines such as enlarge, rank_1_update and such. This happens if the
  // initial matrix is not valid for example. Then, the current state of the
  // tags is used to recompute the factorization... that might not be correct.
  // Need to fix this / make this consistent.
  is_phase_I_ = is_phase_I;
  is_phase_II_ = !is_phase_I;
  is_linear_ = is_linear;
  
  if (!valid_) {
    compute_factorization(is_linear_ || is_phase_I_);
  }

   CGAL_qpe_debug {
     if (vout_.verbose()) {
       vout_.out() << "==========================================" << std::endl;
       vout_.out() << "BEGINNING of LU solve" << std::endl;
       vout_.out() << "is_phase_I_: " << is_phase_I_ << std::endl;
       vout_.out() << "is_phase_II_: " << is_phase_II_ << std::endl;
       vout_.out() << "is_linear_: " << is_linear_ << std::endl;
       vout_.out() << "valid_: " << valid_ << std::endl;
       vout_.out() << "csize_: " << csize_ << std::endl;
       vout_.out() << "bosize_: " << bosize_ << std::endl;
       vout_.out() << "n_: " << n_ << std::endl;
       //print(vout_.out());
     }
   }
  
  if (n_ <= 0) return;
  
  // clear output iterators, i.e. initialize with et0_
  for (int i = 0; i < bosize_; ++i) {
    *(y_x_it+i) = et0_;
  }
  for (int i = 0; i < csize_; ++i) {
    *(y_l_it+i) = et0_;
  }
  
  
  // TAG: CODE DESIGN are both necessary (is_linear_ and is_phase_I)?
  if (is_linear_ || is_phase_I_) {
    solve_LP(v_l_it, v_x_it, y_l_it, y_x_it);
  } else {
    solve_QP(v_l_it, v_x_it, y_l_it, y_x_it);
  }
  
  
   CGAL_qpe_debug {
     if (vout_.verbose()) {
       vout_.out() << "END of LU solve" << std::endl;
       vout_.out() << "is_phase_I_: " << is_phase_I_ << std::endl;
       vout_.out() << "is_phase_II_: " << is_phase_II_ << std::endl;
       vout_.out() << "is_linear_: " << is_linear_ << std::endl;
       vout_.out() << "valid_: " << valid_ << std::endl;
       vout_.out() << "csize_: " << csize_ << std::endl;
       vout_.out() << "bosize_: " << bosize_ << std::endl;
       vout_.out() << "n_: " << n_ << std::endl;
       //print(vout_.out());
       vout_.out() << "==========================================" << std::endl;
     }
   }
}


// private FTRAN QP case (for integral types)
// PRE: symmetric matrix (except in phase_I_)
template <typename ET, typename Is_LP, typename Matrix_Provider>
template < class ForIt, class OutIt>      // QP case
void QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::solve_QP( ForIt v_l_it, ForIt v_x_it,
           OutIt y_l_it, OutIt y_x_it) {
  
  CGAL_qpe_assertion(csize_+bosize_ == n_);

  
  // use 'LP' functions in phase I
  if (is_phase_I_) {
    solve__l(v_x_it, y_l_it);
    solve__x(v_l_it, y_x_it);
    return;
  }
  
  // copy input vector
  QP_sparse_vector<ET> b(n_);
  for (int i = 0; i < csize_; ++i) {
    b.set_entry(i, *(v_l_it + i));
  }
  for (int i = 0; i < bosize_; ++i) {
    b.set_entry(csize_+i, *(v_x_it + i));
  }
  
  typename QP_sparse_vector<ET>::sparse_iterator_t it_b, it_b_end;
  
  b.permute(perm_row_inv_);
  
  // First compute y2 = L^-1 v
  QP_sparse_vector<ET> y2(csize_ + bosize_);
  
  typename QP_sparse_matrix<ET>::value_const_iterator it_row_const_value, it_row_const_value_end;
  it_b_end = b.end();
  for (int i = 0; i < n_; ++i) {
    ET temp = et0_;
    it_b = b.begin();
    it_row_const_value = matrix_l_->begin_row_const_value(i);
    it_row_const_value_end = matrix_l_->end_row_const_value(i);
    
    while (it_row_const_value->first < i+1 && it_b->first < i+1 && it_row_const_value != it_row_const_value_end && it_b != it_b_end) {
      if (it_row_const_value->first == it_b->first) {
        temp += it_b->second * (*it_row_const_value->second);
        ++it_row_const_value;
        ++it_b;
      } else if (it_row_const_value->first < it_b->first) {
        ++it_row_const_value;
      } else {
        ++it_b;
      }
    }
    
    
    y2.set_entry(i, temp);
  }

  // Compute y = U^-1 y2
  QP_sparse_vector<ET> y(csize_ + bosize_);
  
  /* // Inefficient, but straight forward method
  for (int i = n_ - 1; i >= 0; --i) {
    ET temp = y2.get_entry(i) * d_row_[n_];
    for (int j = i + 1; j < n_; ++j) {
      temp -= y.get_entry(j) * matrix_u_->get_entry(perm_row_[i],j);
    }
    y.set_entry(i, integral_division(temp, d_row_[i+1]));
  }*/

  Sparse_const_vector_iterator_t it_vec_y, it_vec_y_end;
  
  for (int i = n_-1; i >= 0; --i) {
    it_vec_y = y.begin();
    it_vec_y_end = y.end();
    it_row_const_value = matrix_u_->begin_row_const_value(perm_row_[i]);
    it_row_const_value_end = matrix_u_->end_row_const_value(perm_row_[i]);
    ET temp = y2.get_entry(i) * d_row_[n_];
    while (it_vec_y != it_vec_y_end &&  it_row_const_value != it_row_const_value_end && it_row_const_value->first < n_) {
      if (it_vec_y->first < it_row_const_value->first) {
        ++it_vec_y;
      } else if (it_vec_y->first > it_row_const_value->first) {
        ++it_row_const_value;
      } else {
        temp -= it_vec_y->second * *it_row_const_value->second; 
        ++it_row_const_value;
        ++it_vec_y;
      }
    }
    
    CGAL_qpe_debug{
      ++CGAL::QP_solver_debug::timer.counter_integral_division;
    }
    
    y.set_entry(perm_column_[i], integral_division(temp, d_row_[i+1]));
  }
  
  // copy result to the output iterator
  Sparse_const_vector_iterator_t it = y.begin();
  while (it != y.end() && it->first < csize_) {
    *(y_l_it + it->first) = d_neg_ ? -it->second : it->second;
    ++it;
  }
  it = y.lower_bound(csize_); // iterator to first nonzero emelent in x part
  while (it != y.end()) {
    *(y_x_it + (it->first - csize_)) = d_neg_ ? -it->second : it->second;
    ++it;
  }
}


// private FTRAN LP case
template <typename ET, typename Is_LP, typename Matrix_Provider>
template < class ForIt, class OutIt>      // LP case
void QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::solve_LP( ForIt v_l_it, ForIt v_x_it,
		     OutIt y_l_it, OutIt y_x_it) {
  solve__l(v_x_it, y_l_it);
  solve__x(v_l_it, y_x_it);
}

    

// special matrix-vector multiplication function for LPs (for integral types)
// private BTRAN (y^T LU = v^T), x part, (y^T LU = v^T)
template <typename ET, typename Is_LP, typename Matrix_Provider>
template < class ForIt, class OutIt >
void
QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::solve__l( ForIt v_x_it, OutIt y_l_it) {
  
  CGAL_qpe_assertion(n_ == csize_ && n_ == bosize_);
  
  
  // TAG: INEFFICIENT (do update instead of recomputing every time)
  if (!valid_) {
    compute_factorization(is_linear_ || is_phase_I_);
  }

  
  if (n_ <= 0) {
    return;
  }
  
  // First compute y2^T = v_x^T U^-1
  QP_sparse_vector<ET> y2(n_);
  
  typename QP_sparse_matrix<ET>::value_const_iterator it_col_const_value, it_col_const_value_end;
  for (int j = 0; j < n_; ++j) {
    it_col_const_value = matrix_l_->begin_column_const_value(j);
    it_col_const_value_end = matrix_l_->end_column_const_value(j);
    ET temp = et0_;
    while (it_col_const_value != it_col_const_value_end && it_col_const_value->first <= j) {
      temp += *(v_x_it + perm_column_[it_col_const_value->first]) * (*it_col_const_value->second);
      ++it_col_const_value;
    }
    y2.set_entry(j, temp);
  }

  // Compute y^T = y2^T L^-1
  QP_sparse_vector<ET> y(n_);
  
  Sparse_const_vector_iterator_t it_vec_y, it_vec_y_end;
  for (int j = n_-1; j >= 0; --j) {
    it_vec_y = y.begin();
    it_vec_y_end = y.end();
    it_col_const_value = matrix_u_->begin_column_const_value(perm_column_[j]);
    it_col_const_value_end = matrix_u_->end_column_const_value(perm_column_[j]);
    ET temp = y2.get_entry(j) * d_row_[n_];
    while (it_vec_y != it_vec_y_end && it_col_const_value != it_col_const_value_end) {
      if (it_vec_y->first < it_col_const_value->first) {
        ++it_vec_y;
      } else if (it_vec_y->first > it_col_const_value->first) {
        ++it_col_const_value;
      } else {
        temp -= it_vec_y->second * *it_col_const_value->second; 
        ++it_col_const_value;
        ++it_vec_y;
      }
    }
    
    CGAL_qpe_debug{
      ++CGAL::QP_solver_debug::timer.counter_integral_division;
    }
    
    y.set_entry(perm_row_[j], integral_division(temp, d_row_[j+1]));
  }
  
  
  // TAG: INEFFICIENT2, because it's done several times...
  // clear y_l
  for (int i = 0; i < n_; ++i){
    *(y_l_it + i) = et0_;
  }
  
  // copy result to the output iterator
  Sparse_const_vector_iterator_t it = y.begin();
  while (it != y.end()) {
    *(y_l_it + it->first) = d_neg_ ? -it->second : it->second;
    ++it;
  }
}


// special matrix-vector multiplication function for LPs (for integral types)
// private FTRAN LP case, lambda part, (LU y = v)
template <typename ET, typename Is_LP, typename Matrix_Provider>
template < class ForIt, class OutIt >
void
QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::solve__x( ForIt v_l_it, OutIt y_x_it) {
  
  CGAL_qpe_assertion(n_ == csize_ && n_ == bosize_);
  
  // TAG: INEFFICIENT (do update instead of recomputing every time)
  if (!valid_) {
    compute_factorization(is_linear_ || is_phase_I_);
  }
  
  if (n_ <= 0) {
    return;
  }

  // First compute y2 = L^-1 v
  QP_sparse_vector<ET> y2(n_);

  typename QP_sparse_matrix<ET>::value_const_iterator it_row_const, it_row_const_end;
  for (int i = 0; i < n_; ++i) {    
    it_row_const = matrix_l_->begin_row_const_value(i);
    it_row_const_end = matrix_l_->end_row_const_value(i);

    ET temp = et0_;

    while (it_row_const != it_row_const_end && it_row_const->first <= i) {
      temp += *(v_l_it + perm_row_[it_row_const->first]) * *it_row_const->second;
      ++it_row_const;
    }
    y2.set_entry(i, temp);
  }
  
  // Compute y = U^-1 y2
  QP_sparse_vector<ET> y(n_);
  
  
  Sparse_const_vector_iterator_t it_vec_y, it_vec_y_end;
  typename QP_sparse_matrix<ET>::value_const_iterator it_row_const_value, it_row_const_value_end;
  for (int i = n_-1; i >= 0; --i) {
    it_vec_y = y.begin();
    it_vec_y_end = y.end();
    it_row_const_value = matrix_u_->begin_row_const_value(perm_row_[i]);
    it_row_const_value_end = matrix_u_->end_row_const_value(perm_row_[i]);
    ET temp = y2.get_entry(i) * d_row_[n_];
    while (it_vec_y != it_vec_y_end &&  it_row_const_value != it_row_const_value_end && it_row_const_value->first < n_) {
      if (it_vec_y->first < it_row_const_value->first) {
        ++it_vec_y;
      } else if (it_vec_y->first > it_row_const_value->first) {
        ++it_row_const_value;
      } else {
        temp -= it_vec_y->second * *it_row_const_value->second; 
        ++it_row_const_value;
        ++it_vec_y;
      }
    }
    
    CGAL_qpe_debug{
      ++CGAL::QP_solver_debug::timer.counter_integral_division;
    }
    
    y.set_entry(perm_column_[i], integral_division(temp, d_row_[i+1]));
  }
  
  
  // TAG: INEFFICIENT2, because it's done several times...
  // clear y_x
  for (int i = 0; i < n_; ++i) {
    *(y_x_it + i) = et0_;
  }
  
  
  // copy result to the output iterator
  Sparse_const_vector_iterator_t it = y.begin();
  while (it != y.end()) {
    *(y_x_it + it->first) = d_neg_ ? -it->second : it->second;
    ++it;
  }
}

// Rank 1 Update
// A' = A + a y z^t
// where a is a scalar, y and z are vectors of appropriate dimensions
template < typename ET, typename Is_LP, typename Matrix_Provider >
bool
QP_LU_factorization<ET, Is_LP, Matrix_Provider>::rank_1_update(ET alpha, QP_sparse_vector<ET> y, QP_sparse_vector<ET> z) {

   CGAL_qpe_debug {
     if (vout_.verbose()) {
       vout_.out() << "==========================================" << std::endl;
       vout_.out() << "BEGINNING of LU rank_1_update" << std::endl;
       //print(vout_.out());
     }
   }

  if (!valid_) {
    CGAL_qpe_debug {
      if (vout_.verbose()) {
        vout_.out() << "Matrix not valid in rank_1_update" << std::endl;
        vout_.out() << "==========================================" << std::endl;
      }
    }
    compute_factorization(is_linear_ || is_phase_I_);
    return false;
  }
  
  
    /*
    % Matlab code, no pivoting, L, Linv, U, Uinv are provided
    P = [tril(L,-1)+U, tril(Linv,-1)+Uinv];
    P2 = [tril(L,-1)+U, tril(Linv,-1)+Uinv];
    Q2 = ones(1,n+1);
    
    yl1 = alpha * Linv * y;
    zl1 = [z,zeros(1,n)];
    yl2 = alpha * Uinv' * z';
    zl2 = [y',zeros(1,n)];
    
    % division-free
    % compute Linv and U AND L and Uinv at the same time
    for i=1:n,
        u = P(i,i);
        urow = P(i,i+1:n+i); % store old values of row
        ucol = [P(i+1:n,i); P(1:i,n+i)];
        
        P(i,i:n+i) = idivide_enforce_int( Q2(i)*P(i,i:n+i) + yl1(i)*zl1(i:n+i), Q(i) ); % update U & Linv
        
        P(i+1:n,i) = idivide_enforce_int(Q2(i)*P(i+1:n,i) + yl2(i)*zl2(i+1:n)', Q(i));
        P(1:i-1,n+i) = idivide_enforce_int( Q2(i)*P(1:i-1,n+i) + yl2(i)*zl2(n+1:n+i-1)', Q(i));
        
        Q2(i+1) = P(i,i);
        zl1(i+1:n+i) = idivide_enforce_int(u*zl1(i+1:n+i) - zl1(i)*urow, Q(i)); % update z
        zl2(i+1:n+i) = idivide_enforce_int(u*zl2(i+1:n+i) - zl2(i)*ucol', Q(i));
        
        
    end % for
    */
      
  
  CGAL_qpe_assertion(y.get_size() == n_ && z.get_size() == n_);
  
  
  QP_sparse_vector<ET> y1(n_), z1(z), y2(n_), z2(y);
  dense_vector_t d_row_tmp(n_+1);
  std::fill(d_row_tmp.begin(), d_row_tmp.end(), et1_);
  
  QP_sparse_vector<ET> zeros(n_);
  z1.append(zeros);
  z2.append(zeros);
  
  
  typename QP_sparse_matrix<ET>::value_const_iterator it_row_const_value, it_row_const_value_end;
  typename QP_sparse_matrix<ET>::value_const_iterator it_column_const_value, it_column_const_value_end;
  
  typedef typename std::vector<std::pair<typename vector_t::key_type, typename vector_t::value_type> >::iterator Iterator;
  std::vector<std::pair<typename vector_t::key_type, typename vector_t::value_type> > urow_u, urow_l, ucol_u, ucol_l;
  
  
  
  // prepare y vectors, i.e., compute L^-1 y and y' U^-1
  typename QP_sparse_vector<ET>::sparse_iterator_t it_y, it_y_end;
  y.permute(perm_row_inv_);
  it_y_end = y.end();
  for (int i = 0; i < n_; ++i) {
    it_row_const_value = matrix_l_->begin_row_const_value(i);
    it_row_const_value_end = matrix_l_->end_row_const_value(i);
    it_y = y.begin();
    ET temp = et0_;
    while (it_y != it_y_end && it_row_const_value != it_row_const_value_end && it_row_const_value->first <= i) {
      if (it_row_const_value->first == it_y->first) {
        temp += alpha * it_y->second * *it_row_const_value->second;
        ++it_row_const_value;
        ++it_y;
      } else if (it_row_const_value->first < it_y->first) {
        ++it_row_const_value;
      } else {
        ++it_y;
      }
    }
    y1.set_entry(perm_row_[i], alpha*temp);
  }
  
  typename QP_sparse_matrix<ET>::value_const_iterator it_col_const_value, it_col_const_value_end;
  typename QP_sparse_vector<ET>::sparse_iterator_t it_z, it_z_end;
  z.permute(perm_column_inv_);
  it_z_end = z.end();
  for (int j = 0; j < n_; ++j) {
    it_col_const_value = matrix_l_->begin_column_const_value(j);
    it_col_const_value_end = matrix_l_->end_column_const_value(j);
    it_z = z.begin();
    ET temp = et0_;
    while (it_z != it_z_end && it_col_const_value != it_col_const_value_end && it_col_const_value->first <= j) {
      if (it_col_const_value->first == it_z->first) {
        temp += it_z->second * *it_col_const_value->second;
        ++it_col_const_value;
        ++it_z;
      } else if (it_col_const_value->first < it_z->first) {
        ++it_col_const_value;
      } else {
        ++it_z;
      }
    }
    y2.set_entry(perm_column_[j], alpha*temp);
  }
  
  for (int i = 0; i < n_; ++i) {
  
    ET u = matrix_u_->get_entry(perm_row_[i], perm_column_[i]);
    
    if (u == 0) {
      set_invalid();
      CGAL_qpe_debug {
        if (vout_.verbose()) {
          vout_.out() << "Pivot failure rank_1_update" << std::endl;
          vout_.out() << "==========================================" << std::endl;
        }
      }
      return false;
    }
    
    // update U and L
    typename QP_sparse_matrix<ET>::value_iterator it_row_value, it_row_value_end;
    typename QP_sparse_matrix<ET>::value_iterator it_column_value, it_column_value_end;
    
    // get old values of pivot rows and columns AND
    // premultiply entries with d_row_tmp[i]
    urow_u.clear();
    urow_l.clear();
    ET tmp = z1.get_entry(perm_column_[i]);
  
    it_row_value = matrix_u_->begin_row_value(perm_row_[i]);
    it_row_value_end = matrix_u_->end_row_value(perm_row_[i]);
    
    
    if (tmp != et0_) {
      while (it_row_value != it_row_value_end) {
        if (perm_column_inv_[it_row_value->first] > i) {
          urow_u.push_back(std::make_pair(it_row_value->first, tmp * *it_row_value->second));
          *it_row_value->second *= d_row_tmp[i];
        }
        ++it_row_value;
      }
      it_row_value = matrix_l_->begin_row_value(i);
      it_row_value_end = matrix_l_->end_row_value(i);
      while (it_row_value != it_row_value_end && it_row_value->first < i) {
        urow_l.push_back(std::make_pair(it_row_value->first, tmp * *it_row_value->second));
        *it_row_value->second *= d_row_tmp[i];
        ++it_row_value;
      }
      if (it_row_value != it_row_value_end && it_row_value->first == i) {
        urow_l.push_back(std::make_pair(it_row_value->first, tmp * *it_row_value->second));
      }
    } else {
      while (it_row_value != it_row_value_end) {
        if (perm_column_inv_[it_row_value->first] > i) {
          *it_row_value->second *= d_row_tmp[i];
        }
        ++it_row_value;
      }
      it_row_value = matrix_l_->begin_row_value(i);
      it_row_value_end = matrix_l_->end_row_value(i);
      while (it_row_value != it_row_value_end &&  it_row_value->first < i) {
        *it_row_value->second *= d_row_tmp[i];
        ++it_row_value;
      }
    }
      
      
    ucol_u.clear();
    ucol_l.clear();
    tmp = z2.get_entry(perm_row_[i]);
    
    it_column_value = matrix_u_->begin_column_value(perm_column_[i]);
    it_column_value_end = matrix_u_->end_column_value(perm_column_[i]);
    
    if (tmp != et0_) {
      while (it_column_value != it_column_value_end) {
        if (perm_row_inv_[it_column_value->first] > i) {
          ucol_u.push_back(std::make_pair(it_column_value->first, tmp * *it_column_value->second));
          *it_column_value->second *= d_row_tmp[i];
        } else if (perm_row_inv_[it_column_value->first] == i) {
          *it_column_value->second *= d_row_tmp[i];
        }
        ++it_column_value;
      }
      it_column_value = matrix_l_->begin_column_value(i);
      it_column_value_end = matrix_l_->end_column_value(i);
      while (it_column_value != it_column_value_end && it_column_value->first <= i) {
        ucol_l.push_back(std::make_pair(it_column_value->first, tmp * *it_column_value->second));
        *it_column_value->second *= d_row_tmp[i];
        ++it_column_value;
      }
    } else {
      while (it_column_value != it_column_value_end) {
        if (perm_row_inv_[it_column_value->first] >= i) {
          *it_column_value->second *= d_row_tmp[i];
        }
        ++it_column_value;
      }
      it_column_value = matrix_l_->begin_column_value(i);
      it_column_value_end = matrix_l_->end_column_value(i);
      while (it_column_value != it_column_value_end && it_column_value->first <= i) {
        *it_column_value->second *= d_row_tmp[i];
        ++it_column_value;
      }
    }

    
    // add y*z term
    typename QP_sparse_vector<ET>::sparse_iterator_t it_z, it_z_end;
    it_z = z1.begin();
    it_z_end = z1.end();
    ET y_tmp = y1.get_entry(perm_row_[i]);
    if (y_tmp != et0_) {
      while (it_z != it_z_end && it_z->first < n_) {
        if (perm_column_inv_[it_z->first] >= i) {
          matrix_u_->set_entry(perm_row_[i], it_z->first, matrix_u_->get_entry(perm_row_[i], it_z->first) + y_tmp * it_z->second );
        }
        ++it_z;
      }
      while (it_z != it_z_end && it_z->first <= n_ + i) {
        matrix_l_->set_entry(i, it_z->first - n_, matrix_l_->get_entry(i, it_z->first - n_) + y_tmp * it_z->second );
        ++it_z;
      }
    }

    it_z = z2.begin();
    it_z_end = z2.end();
    y_tmp = y2.get_entry(perm_column_[i]);
    if (y_tmp != et0_) {
      while (it_z != it_z_end && it_z->first < n_) {
        if (perm_row_inv_[it_z->first] > i) {
          matrix_u_->set_entry(it_z->first, perm_column_[i], matrix_u_->get_entry(it_z->first, perm_column_[i]) + y_tmp * it_z->second );
        }
        ++it_z;
      }
      while (it_z != it_z_end && it_z->first < n_ + i) {
        matrix_l_->set_entry(it_z->first - n_, i, matrix_l_->get_entry(it_z->first - n_, i) + y_tmp * it_z->second );
        ++it_z;
      }
    }
    
    if (d_row_[i] == 0) {
      set_invalid();
      CGAL_qpe_debug {
        if (vout_.verbose()) {
          vout_.out() << "Pivot failure rank_1_update" << std::endl;
          vout_.out() << "==========================================" << std::endl;
        }
      }
      return false;
    }
    
    // divide by d_row_[i]
    it_row_value = matrix_u_->begin_row_value(perm_row_[i]);
    it_row_value_end = matrix_u_->end_row_value(perm_row_[i]);
    while (it_row_value != it_row_value_end) {
      if (perm_column_inv_[it_row_value->first] >= i) {
        *it_row_value->second = integral_division(*it_row_value->second, d_row_[i]);
      }
      ++it_row_value;
    }
    it_row_value = matrix_l_->begin_row_value(i);
    it_row_value_end = matrix_l_->end_row_value(i);
    while (it_row_value != it_row_value_end) {
      if (it_row_value->first <= i) {
        *it_row_value->second = integral_division(*it_row_value->second, d_row_[i]);
      }
      ++it_row_value;
    }
    
    it_column_value = matrix_u_->begin_column_value(perm_column_[i]);
    it_column_value_end = matrix_u_->end_column_value(perm_column_[i]);
    while (it_column_value != it_column_value_end) {
      if (perm_row_inv_[it_column_value->first] > i) {
        *it_column_value->second = integral_division(*it_column_value->second, d_row_[i]);
      }
      ++it_column_value;
    }
    it_column_value = matrix_l_->begin_column_value(i);
    it_column_value_end = matrix_l_->end_column_value(i);
    while (it_column_value != it_column_value_end) {
      if (it_column_value->first < i) {
        *it_column_value->second = integral_division(*it_column_value->second, d_row_[i]);
      }
      ++it_column_value;
    }
  
    // update z1
    Iterator it_u, it_u_end;
    it_z = z1.begin();
    it_z_end = z1.end();
    ET zi = z1.get_entry(perm_column_[i]);
    while (it_z != it_z_end && it_z->first < n_) {
      if (perm_column_inv_[it_z->first] > i) {
        it_z->second *= u;
      }
      ++it_z;
    }
    while (it_z != it_z_end && it_z->first < n_+i+1) {
      it_z->second *= u;
      ++it_z;
    }
    it_u = urow_u.begin();
    it_u_end = urow_u.end();
    while (it_u != it_u_end) {
      z1.set_entry(it_u->first, z1.get_entry(it_u->first) - it_u->second);
      ++it_u;
    }
    it_u = urow_l.begin();
    it_u_end = urow_l.end();
    while (it_u != it_u_end) {
      z1.set_entry(it_u->first + n_, z1.get_entry(it_u->first + n_) - it_u->second);
      ++it_u;
    }
    it_z = z1.begin();
    while (it_z != it_z_end && it_z->first < n_) {
      if (perm_column_inv_[it_z->first] > i) {
        it_z->second = integral_division(it_z->second, d_row_[i]);
      }
      ++it_z;
    }
    while (it_z != it_z_end && it_z->first < n_+i+1) {
      it_z->second = integral_division(it_z->second, d_row_[i]);
      ++it_z;
    }
    
    // update z2
    it_z = z2.begin();
    it_z_end = z2.end();
    zi = z2.get_entry(perm_row_[i]);
    while (it_z != it_z_end && it_z->first < n_) {
      if (perm_row_inv_[it_z->first] > i) {
        it_z->second *= u;
      }
      ++it_z;
    }
    while (it_z != it_z_end && it_z->first < n_+i+1) {
      it_z->second *= u;
      ++it_z;
    }
    it_u = ucol_u.begin();
    it_u_end = ucol_u.end();
    while (it_u != it_u_end) {
      z2.set_entry(it_u->first, z2.get_entry(it_u->first) - it_u->second);
      ++it_u;
    }
    it_u = ucol_l.begin();
    it_u_end = ucol_l.end();
    while (it_u != it_u_end) {
      z2.set_entry(it_u->first + n_, z2.get_entry(it_u->first + n_) - it_u->second);
      ++it_u;
    }
    it_z = z2.begin();
    while (it_z != it_z_end && it_z->first < n_) {
      if (perm_row_inv_[it_z->first] > i) {
        it_z->second = integral_division(it_z->second, d_row_[i]);
      }
      ++it_z;
    }
    while (it_z != it_z_end && it_z->first < n_+i+1) {
      it_z->second = integral_division(it_z->second, d_row_[i]);
      ++it_z;
    }
      
    // get new denominator for next pivot row
    d_row_tmp[i+1] = matrix_u_->get_entry(perm_row_[i], perm_column_[i]);
    
    // check if pivot is bad...
    // either stop and refactor or introduce pivoting in rank 1 update
    // as of now: stop and refactor
    if (d_row_tmp[i+1] == 0) {
      set_invalid();
      CGAL_qpe_debug {
        if (vout_.verbose()) {
          vout_.out() << "Pivot failure rank_1_update" << std::endl;
          vout_.out() << "==========================================" << std::endl;
        }
      }
      return false;
    }
  }
  
  // copy denominators
  d_row_.clear();
  std::copy(d_row_tmp.begin(), d_row_tmp.end(), std::back_inserter(d_row_));
  
  if (d_row_[n_] > 0) {
    d_ = d_row_[n_];
    d_neg_ = false;
  } else {
    d_ = -d_row_[n_];
    d_neg_ = true; 
  }
  
  CGAL_qpe_debug {
     if (vout_.verbose()) {
       vout_.out() << "END of LU rank_1_update" << std::endl;
       //print(vout_.out());
       vout_.out() << "==========================================" << std::endl;
     }
   }
  
  return true;
}

// Enlarge matrix by unit row/column
// PRE: We are in linear case!
// POST: zero row and column are added at position k
// PARAM: k is the row/col at which the matrix should be enlarged. It pertains to
//        the input ordering of the matrix. On the other hand, a is the index
//        where the row/col should be added in the factorization ordering.
//        is_basic indicates whether the set of basic variables or the set of
//        constraints should be enlarged (only applies in QP case).
// TODO YB: add functionality for the quadratic case as well
// TODO YB: Make all is_linear distinction template parameters and not bools
template < typename ET, typename Is_LP, typename Matrix_Provider >
bool
QP_LU_factorization<ET, Is_LP, Matrix_Provider>::enlarge(unsigned int k, unsigned int a, bool is_basic) {
  
  CGAL_qpe_assertion(k <= n_);
  CGAL_qpe_assertion(a <= n_);
  
  if (!valid_) {
    compute_factorization(is_linear_ || is_phase_I_);
    return false;
  }
  
  CGAL_qpe_debug {
    if (vout_.verbose()) {
      vout_.out() << "==========================================" << std::endl;
      vout_.out() << "BEGINNING of LU enlarge matrix" << std::endl;
      //print(vout_.out());
      vout_.out() << "a: " << a << std::endl;
      vout_.out() << "k: " << k << std::endl;
      vout_.out() << "Original matrix:\n" << *recover_original_matrix();
    }
  }
  
  if (a < n_) {
    matrix_l_->insert_row(a, vector_t(n_, et0_));
    matrix_l_->insert_column(a, vector_t(n_+1, et0_));
  } else {
    matrix_l_->append_row(vector_t(n_, et0_));
    matrix_l_->append_column(vector_t(n_+1, et0_));
  }
  
  if (k < n_) {
    matrix_u_->insert_row(k, vector_t(n_, et0_));
    matrix_u_->insert_column(k, vector_t(n_+1, et0_));
  } else {
    matrix_u_->append_row(vector_t(n_, et0_));
    matrix_u_->append_column(vector_t(n_+1, et0_));
  
  }
  
  ET new_denominator = d_row_[a];
  d_row_.insert(d_row_.begin()+a, new_denominator);
  
  
  matrix_u_->set_entry(k, k, new_denominator);
  matrix_l_->set_entry(a, a, new_denominator);

  // TODO YB: do more efficiently by traversing counterpart structure...
  std::transform(perm_row_.begin(), perm_row_.end(), perm_row_.begin(), std::bind1st(Conditional_inc(), k));
  std::transform(perm_column_.begin(), perm_column_.end(), perm_column_.begin(), std::bind1st(Conditional_inc(), k));
  std::transform(perm_row_inv_.begin(), perm_row_inv_.end(), perm_row_inv_.begin(), std::bind1st(Conditional_inc(), a));
  std::transform(perm_column_inv_.begin(), perm_column_inv_.end(), perm_column_inv_.begin(), std::bind1st(Conditional_inc(), a));

  perm_row_.insert(perm_row_.begin()+a, k);
  perm_row_inv_.insert(perm_row_inv_.begin()+k, a);
  perm_column_.insert(perm_column_.begin()+a, k);
  perm_column_inv_.insert(perm_column_inv_.begin()+k, a);
  
  if (is_linear_ || is_phase_I_) {
    ++bosize_;
    ++csize_;
  } else {
    if (is_basic) {
      ++bosize_;
    } else {
      ++csize_;
    }
  }
  ++n_;
  
  if (d_row_[n_] > 0) {
    d_ = d_row_[n_];
    d_neg_ = false;
  } else {
    d_ = -d_row_[n_];
    d_neg_ = true; 
  }

  CGAL_qpe_debug {
    if (vout_.verbose()) {
      vout_.out() << "END of LU enlarge matrix" << std::endl;
      //print(vout_.out());
      vout_.out() << "Modified matrix:\n" << *recover_original_matrix();
      vout_.out() << "==========================================" << std::endl;
    }
  }
  
  return true;
}

// Shrink matrix by row/column
// PRE: We are in linear case and the columne/row to be deleted is a unit vector
// POST: zero row and column are deleted at position k
// TODO: add functionality for the quadratic case as well
// TODO: Make all is_linear distinction template parameters and not bools
template < typename ET, typename Is_LP, typename Matrix_Provider >
bool
QP_LU_factorization<ET, Is_LP, Matrix_Provider>::shrink(unsigned int k, bool is_basic) {
  
  CGAL_qpe_assertion(k < n_);
  CGAL_qpe_assertion(valid_); // should never be invalid. Not like enlarge, where
                              // we might start with empty matrix
  
  /*
  // remove/replace the following. Refactoring at this stage is not a clever idea...
  if (!valid_) {
    compute_factorization(is_linear_ || is_phase_I_);
    return false;
  }*/

  
  CGAL_qpe_debug {
    if (vout_.verbose()) {
      vout_.out() << "==========================================" << std::endl;
      vout_.out() << "BEGINNING of LU shrink matrix" << std::endl;
      //print(vout_.out());
      vout_.out() << "k: " << k << std::endl;
      vout_.out() << "Original matrix:\n" << *recover_original_matrix();
    }
  }
   
  
  int actual_row_index = perm_row_[k];
  int actual_col_index = perm_column_[k];
   
  
  if (actual_row_index == n_-1) {
    matrix_u_->strip_row();
  } else {
    matrix_u_->delete_row(actual_row_index);
  }
  
  if  (actual_col_index == n_-1) {
    matrix_u_->strip_column();
  } else {
    matrix_u_->delete_column(actual_col_index);  
  }
  
  if (k == n_-1) {
    matrix_l_->strip_row();
    matrix_l_->strip_column();
    d_row_.pop_back();
  } else {
    matrix_l_->delete_row(k);
    matrix_l_->delete_column(k);
    d_row_.erase(d_row_.begin()+k);
  }
  
  perm_row_.erase(perm_row_.begin()+k);
  perm_column_.erase(perm_column_.begin()+k);
  std::transform(perm_row_.begin(), perm_row_.end(), perm_row_.begin(), std::bind1st(Conditional_dec(), actual_row_index));
  std::transform(perm_column_.begin(), perm_column_.end(), perm_column_.begin(), std::bind1st(Conditional_dec(), actual_col_index));
  
  
  perm_row_inv_.erase(perm_row_inv_.begin()+actual_row_index);
  perm_column_inv_.erase(perm_column_inv_.begin()+actual_col_index);
  std::transform(perm_row_inv_.begin(), perm_row_inv_.end(), perm_row_inv_.begin(), std::bind1st(Conditional_dec(), k));
  std::transform(perm_column_inv_.begin(), perm_column_inv_.end(), perm_column_inv_.begin(), std::bind1st(Conditional_dec(), k));  
  

  if (is_linear_ || is_phase_I_) {
    --bosize_;
    --csize_;
  } else {
    if (is_basic) {
      --bosize_;
    } else {
      --csize_;
    }
  }
  --n_;  
  
  if (d_row_[n_] > 0) {
    d_ = d_row_[n_];
    d_neg_ = false;
  } else {
    d_ = -d_row_[n_];
    d_neg_ = true; 
  }
  


  CGAL_qpe_debug {
    if (vout_.verbose()) {
      vout_.out() << "END of LU shrink matrix" << std::endl;
      //print(vout_.out());
      vout_.out() << "Modified matrix:\n" << *recover_original_matrix();
      vout_.out() << "==========================================" << std::endl;
    }
  }
   
  return true;
}


// To retreive the internal permutation with regards to LU pivoting
template < typename ET, typename Is_LP, typename Matrix_Provider >
void
QP_LU_factorization<ET, Is_LP, Matrix_Provider>::get_perm_row(permutation_t& perm) const {
  perm.clear();
  int n = perm_row_.size();
  perm.resize(n);
  std::copy(perm_row_.begin(), perm_row_.end(), perm.begin());
}

template < typename ET, typename Is_LP, typename Matrix_Provider >
void
QP_LU_factorization<ET, Is_LP, Matrix_Provider>::get_perm_col(permutation_t& perm) const {
  perm.clear();
  int n = perm_column_.size();
  perm.resize(n);
  std::copy(perm_column_.begin(), perm_column_.end(), perm.begin());
}

template < typename ET, typename Is_LP, typename Matrix_Provider >
void
QP_LU_factorization<ET, Is_LP, Matrix_Provider>::get_perm_row_inv(permutation_t& perm) const {
  perm.clear();
  int n = perm_row_inv_.size();
  perm.resize(n);
  std::copy(perm_row_inv_.begin(), perm_row_inv_.end(), perm.begin());
}

template < typename ET, typename Is_LP, typename Matrix_Provider >
void
QP_LU_factorization<ET, Is_LP, Matrix_Provider>::get_perm_col_inv(permutation_t& perm) const {
  perm.clear();
  int n = perm_column_inv_.size();
  perm.resize(n);
  std::copy(perm_column_inv_.begin(), perm_column_inv_.end(), perm.begin());
}


template < typename ET, typename Is_LP, typename Matrix_Provider >
typename Algebraic_structure_traits< typename Coercion_traits<ET,ET>::Type >
::Gcd::result_type
QP_LU_factorization<ET, Is_LP, Matrix_Provider>::lcm(ET a, ET b) {
    return lcm( a, b, Category());
}

template < typename ET, typename Is_LP, typename Matrix_Provider >
typename Algebraic_structure_traits< typename Coercion_traits<ET,ET>::Type >
::Gcd::result_type
QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::lcm(ET a, ET b, Field_tag ) {
  std::cout << "field" << std::endl;
  std::cout << Coercion_traits<ET,ET>::Are_explicit_interoperable::value << std::endl;
  std::cout << typeid(typename Algebraic_structure_traits< typename Coercion_traits<ET,ET>::Type >
                      ::Gcd::result_type).name() << std::endl;
  std::cout << typeid(typename Coercion_traits<ET,ET>::Type).name() << std::endl;
//  return a > b ? a : b;
  Null_tag ret;
  return ret;
}

template < typename ET, typename Is_LP, typename Matrix_Provider >
typename Algebraic_structure_traits< typename Coercion_traits<ET,ET>::Type >
::Gcd::result_type
QP_LU_factorization<ET, Is_LP, Matrix_Provider>
::lcm(ET a, ET b, Unique_factorization_domain_tag ) {
  std::cout << "unique\n";  
  return CGAL::gcd<ET, ET>( a, b );
}

} //namespace CGAL

// ===== EOF ==================================================================
