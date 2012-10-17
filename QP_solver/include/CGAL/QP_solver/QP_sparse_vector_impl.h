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


#include <utility>

namespace CGAL {


// =============================================================================
// Construction
// =============================================================================
template <typename NT>
QP_sparse_vector<NT>::QP_sparse_vector(): n_(0), nt0_(NT(0)), entries_(map_t())
{}

template <typename NT>
QP_sparse_vector<NT>::QP_sparse_vector(int n, NT default_value):
  n_(n),
  nt0_(default_value),
  entries_(map_t())
{}
  
template <typename NT>
QP_sparse_vector<NT>::QP_sparse_vector(const QP_sparse_vector<NT>& v):
  n_(v.n_),
  nt0_(v.nt0_),
  entries_(map_t())
{
  this->entries_.insert(v.entries_.begin(), v.entries_.end());
  /*std::copy(v.entries_.begin(), v.entries_.end(),
              std::inserter(this->entries_, this->entries_.begin()));*/
}


// =============================================================================
// Info functions
// =============================================================================

template <typename NT>
int QP_sparse_vector<NT>::get_size() const {
  return n_;
}

template <typename NT>
int QP_sparse_vector<NT>::get_number_of_nonzeros() const {
  return entries_.size();
}


// =============================================================================
// Getter & setter methods
// =============================================================================

template <typename NT>  
void QP_sparse_vector<NT>::set_entry(key_type n, NT val) {
  if (val != nt0_) {
    entries_[n] = val;
  } else {
    entries_.erase(n);
  }
}

template <typename NT>
template <typename Iterator>
void QP_sparse_vector<NT>::insert(Iterator begin, Iterator end) {
  entries_.insert(begin, end);
}

template <typename NT>  
void QP_sparse_vector<NT>::clear_entry(key_type n) {
  sparse_iterator_t it = entries_.find(n);
  if (it != entries_.end()) {
    entries_.erase(it);
  }
}


template <typename NT>
const NT& QP_sparse_vector<NT>::get_entry(key_type n) const {
  sparse_const_iterator_t it = entries_.find(n);
  return (it == entries_.end() ? nt0_ : (*it).second);
}

template <typename NT>  
void QP_sparse_vector<NT>::swap_entries(key_type i, key_type j) {
  typename map_t::iterator it_i = entries_.find(i);
  typename map_t::iterator it_j = entries_.find(j);
  typename map_t::iterator it_end = entries_.end();
  
  if (i != j) {
    if (it_i != it_end) {
      if (it_j != it_end) {
        std::swap(it_i->second, it_j->second);
      } else {
        NT tmp = it_i->second;
        entries_.erase(it_i);
        entries_[j] = tmp;
      }
    } else {
      if (it_j != it_end) {
        NT tmp = it_j->second;
        entries_.erase(it_j);
        entries_[i] = tmp;
      }
    }
  }
}


// =============================================================================
// Operations
// =============================================================================
template <typename NT>
NT QP_sparse_vector<NT>::inner_product(const QP_sparse_vector& v) const {

  CGAL_qpe_assertion(n_ == v.n_);
  
  NT ret = nt0_;
  sparse_const_iterator_t itA(begin()), itAend(end()),
                        itV(v.begin()), itVend(v.end());
  // TAG: TODO: remove special treatment of nt0_!=NT(0 )
  if (nt0_ == NT(0) && v.nt0_ == NT(0)) {
     while (itA != itAend && itV != itVend) {
      // DEBUG output
      // std::cout << (*itA).first << " " << (*itV).first << std::endl;
      if ((*itA).first < (*itV).first) {
        ++itA;
      } else if ((*itA).first > (*itV).first) {
        ++itV;
      } else {
        ret += (*itA).second * (*itV).second;
        ++itA;
        ++itV;
      }
    }
  } else {
    for (int i = 0; i < n_; ++i) {
      ret += get_entry(i) * v.get_entry(i);
    }
  }
  return ret;
}



// TAG: TODO: ATTENTION this applies inverse permutation...
template <typename NT>
QP_sparse_vector<NT>&
QP_sparse_vector<NT>::permute(const std::vector<int>& permutation) {
  
  CGAL_qpe_assertion(permutation.size() == n_);
  
  map_t temp;
  sparse_const_iterator_t it = this->begin();
  sparse_const_iterator_t it_end = this->end();
  while (it != it_end) {
    temp.insert(std::pair<int, NT>(permutation[(*it).first], (*it).second));
    ++it;
  }
  entries_ = temp;
  return *this;
}

// PRE: this != &v
template <typename NT>
void
QP_sparse_vector<NT>::append(const QP_sparse_vector& v) {

  CGAL_qpe_assertion(this != &v);
  
  sparse_const_iterator_t it;
  for (it = v.begin(); it != v.end(); ++it) {
    entries_.insert(std::pair<int, NT>((*it).first + n_, (*it).second));
  }
  n_ += v.n_;
}

template <typename NT>
void
QP_sparse_vector<NT>::clear() {
  entries_ = map_t();
}



template <typename NT>
std::ostream& QP_sparse_vector<NT>::print(std::ostream& out) const {
  out << "(";
  for (int i = 0; i < n_ - 1; ++i) {
    out << get_entry(i) << ", ";
  }
  out << get_entry(n_-1) << ")";
  return out;
}


// =============================================================================
// Operators
// =============================================================================
template <typename NT>
const NT& QP_sparse_vector<NT>::operator[] (int n) const {
  return get_entry(n);
}

template <typename NT>
NT QP_sparse_vector<NT>::operator* (const QP_sparse_vector<NT>& v) const {
  return this->inner_product(v);
}

template <typename NT>
QP_sparse_vector<NT>&
QP_sparse_vector<NT>::operator= (const QP_sparse_vector<NT>& rhs) {
  if (this != &rhs) {
    this->n_ = rhs.n_;
    this->nt0_ = rhs.nt0_;
    this->entries_.clear();
    std::copy(rhs.entries_.begin(), rhs.entries_.end(),
              std::inserter(this->entries_, this->entries_.begin()));
  }
  return *this;
}

// =============================================================================
// Iterators
// =============================================================================

template <typename NT>
typename QP_sparse_vector<NT>::dense_iterator_t
QP_sparse_vector<NT>::begin_dense() {
  return dense_iterator_t(boost::counting_iterator<int>(0),
                  CGAL::Map_with_default<map_t>(&entries_, NT(0)));
}

template <typename NT>
typename QP_sparse_vector<NT>::sparse_iterator_t QP_sparse_vector<NT>::begin() {
  return entries_.begin();
}

template <typename NT>
typename QP_sparse_vector<NT>::sparse_iterator_t QP_sparse_vector<NT>::end() {
  return entries_.end();
}

template <typename NT>
typename QP_sparse_vector<NT>::sparse_iterator_t
QP_sparse_vector<NT>::lower_bound(const key_type i) {
  return entries_.lower_bound(i);
}

template <typename NT>
typename QP_sparse_vector<NT>::sparse_const_iterator_t
QP_sparse_vector<NT>::begin() const {
  return entries_.begin();
}

template <typename NT>
typename QP_sparse_vector<NT>::sparse_const_iterator_t
QP_sparse_vector<NT>::end() const {
  return entries_.end();
}

template <typename NT>
typename QP_sparse_vector<NT>::sparse_const_iterator_t
QP_sparse_vector<NT>::lower_bound(const key_type i) const {
  return entries_.lower_bound(i);
}



} //namespace CGAL

// ===== EOF ==================================================================