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
QP_sparse_matrix<NT>::QP_sparse_matrix():
  m_(0),
  n_(0),
  nt0_(NT(0)),
  rows_(std::vector<QP_sparse_vector<int> >()),
  columns_(std::vector<QP_sparse_vector<int> >()),
  data_(std::vector<NT>()),
  next_index_(0)
{
}

// TODO try initializer list contruction of rows_ and columns_
template <typename NT>
QP_sparse_matrix<NT>::QP_sparse_matrix(int m, int n, NT default_value):
  m_(m),
  n_(n),
  nt0_(default_value),
  rows_(std::vector<QP_sparse_vector<int> >()),
  columns_(std::vector<QP_sparse_vector<int> >()),
  data_(std::vector<NT>()),
  next_index_(0)
{
  std::fill_n(std::back_inserter(rows_), m_,
              QP_sparse_vector<int>(n_, -1));
  std::fill_n(std::back_inserter(columns_), n_,
              QP_sparse_vector<int>(m_, -1));              
}

template <typename NT>
QP_sparse_matrix<NT>::QP_sparse_matrix(const QP_sparse_matrix& m) {
  this->operator=(m);
}



// =============================================================================
// Info functions
// =============================================================================
template <typename NT>
typename QP_sparse_matrix<NT>::index_pair_t
QP_sparse_matrix<NT>::get_size() const {
  return std::make_pair(m_, n_);
}



// =============================================================================
// Getter & setter methods
// =============================================================================

template <typename NT>
const NT& QP_sparse_matrix<NT>::get_entry(int m, int n) const {
  int index = rows_.at(m).get_entry(n);
  if (index >= 0) {
    return data_[index];
  } else {
    return nt0_;
  }
}

template <typename NT>
const NT& QP_sparse_matrix<NT>::get_entry(index_pair_t p) const {
  int index = rows_.at(p.first).get_entry(p.second);
  if (index >= 0) {
    return data_[index];
  } else {
    return nt0_;
  }
}

template <typename NT>
bool QP_sparse_matrix<NT>::set_entry(int m, int n, const NT& val) {
  bool ret = true;
  int index = rows_.at(m).get_entry(n);
  
  if (index < 0) {
    if (val != nt0_) {
      rows_.at(m).set_entry(n, next_index_);
      columns_.at(n).set_entry(m, next_index_);
      if (data_.size() >= data_.capacity()) { // resize data_
        data_.reserve(data_.capacity() * 2);
      }
      ++next_index_;
      data_.push_back(val);
    }
  } else {
    data_[index] = val;
    /*
    if (val != nt0_) {
      data_[index] = val;
    } else if (val == nt0_) { 
      // TAG: TODO Not clear whether this should be adopted. Test this version
      // against the version in which we don't do anything at all in this case,
      // and just leave the unused entries in the middle of data_. This version
      // could be advantageous if matrices are kept for a long time.
      rows_.at(m).clear_entry(n);
      columns_.at(n).clear_entry(m);
      data_.erase(data_.begin() + index);
      if (next_index_ > 0) --next_index_;
      erase_index(index);
    }*/
  }
  
  return ret;
}

template <typename NT>
const QP_sparse_vector<int>&
QP_sparse_matrix<NT>::get_row_indices(int m) const {
  CGAL_qpe_assertion(m < m_ && m > -1);
  return rows_.at(m);
}

template <typename NT>
const QP_sparse_vector<int>&
QP_sparse_matrix<NT>::get_column_indices(int n) const {
  CGAL_qpe_assertion(n < n_ && n > -1);
  return columns_.at(n);
}

template <typename NT>
QP_sparse_vector<int>&
QP_sparse_matrix<NT>::get_row_indices(int m) {
  CGAL_qpe_assertion(m < m_ && m > -1);
  return rows_.at(m);
}

template <typename NT>
QP_sparse_vector<int>&
QP_sparse_matrix<NT>::get_column_indices(int n) {
  CGAL_qpe_assertion(n < n_ && n > -1);
  return columns_.at(n);
}


template <typename NT>
const NT&
QP_sparse_matrix<NT>::get_element_by_index(int index) const {
  CGAL_qpe_assertion(index < data_.size() && index > -1);
  return data_[index];
}

template <typename NT>
NT*
QP_sparse_matrix<NT>::get_element_pointer_by_index(int index) {
  CGAL_qpe_assertion(index < data_.size() && index > -1);
  return &data_[index];
}

template <typename NT>
void
QP_sparse_matrix<NT>::set_element_by_index(int index, const NT& val) {
  CGAL_qpe_assertion(index < data_.size() && index > -1);
  data_[index] = val;
}


// POST: This computes the programmatic density, i.e., it counts the number of
//       sparse entries in the matrix. Note that this may differ from the non-zero
//       density, because some of the entries may be 0.
template <typename NT>
double
QP_sparse_matrix<NT>::get_density() const {

  double count = 0.0;
  
  for (int i = 0; i < m_; ++i) {
    value_const_iterator it_row_begin = begin_row_const_value(i);
    value_const_iterator it_row_end = end_row_const_value(i);
    while (it_row_begin != it_row_end) {
      count += 1.0;
      ++it_row_begin;
    }
  }
  
  return count / (m_ * n_);
}


// =============================================================================
// Operations
// =============================================================================

template <typename NT>
QP_sparse_vector<NT>
QP_sparse_matrix<NT>::matrix_vector_multiplication
(const QP_sparse_vector<NT>& v) const {
  CGAL_qpe_assertion(n_ == v.get_size());
  QP_sparse_vector<NT> ret(m_, nt0_);
  for (int i = 0; i < m_; ++i) {
    NT temp = inner_product_row(i, v);
    ret.set_entry(i, temp);
  }
  return ret;
}

template <typename NT>
boost::shared_ptr<QP_sparse_matrix<NT> >
QP_sparse_matrix<NT>::matrix_matrix_multiplication
(const QP_sparse_matrix<NT>& m) const {
  CGAL_qpe_assertion( n_ == m.m_ );
  boost::shared_ptr<QP_sparse_matrix<NT> >
    ret( new QP_sparse_matrix<NT>(m_, m.n_, nt0_) );
  for (int i = 0; i < m_; ++i) {
      value_const_iterator it_row_end = end_row_const_value(i);
    for (int j = 0; j < m.n_; ++j) {
      value_const_iterator it_row_begin = begin_row_const_value(i);
      value_const_iterator it_column_begin = m.begin_column_const_value(j);
      value_const_iterator it_column_end = m.end_column_const_value(j);
      NT tmp = nt0_;
      while (it_row_begin != it_row_end && it_column_begin != it_column_end) {
        if (it_row_begin->first < it_column_begin->first) {
          ++it_row_begin;
        } else if (it_row_begin->first > it_column_begin->first) {
          ++it_column_begin;
        } else { // it_row_begin->first == it_column_begin->first
          tmp += (*it_row_begin->second) * (*it_column_begin->second);
          ++it_row_begin;
          ++it_column_begin;
        }
      }
      ret->set_entry(i, j, tmp);
    } // j
  } // i
  return ret;
}

template <typename NT>
NT
QP_sparse_matrix<NT>::inner_product_row(int i, const QP_sparse_vector<NT>& v) const {
  value_const_iterator it_row_begin = begin_row_const_value(i);
  value_const_iterator it_row_end = end_row_const_value(i);
  return inner_product_by_iterators(it_row_begin, it_row_end, v);
}

template <typename NT>
NT
QP_sparse_matrix<NT>::inner_product_column(int j, const QP_sparse_vector<NT>& v) const {
  value_const_iterator it_column_begin = begin_column_const_value(j);
  value_const_iterator it_column_end = end_column_const_value(j);
  return inner_product_by_iterators(it_column_begin, it_column_end, v);
}


template <typename NT>
void
QP_sparse_matrix<NT>::append_column(const QP_sparse_vector<NT>& v) {
  CGAL_qpe_assertion(m_ == v.n_);
  // set up data structures
  columns_.push_back(QP_sparse_vector<int>(m_,-1));
  for (int i = 0; i < m_; ++i) {
    ++rows_[i].n_;
  }
  ++n_;  
  
  // set entries
  typename QP_sparse_vector<NT>::sparse_const_iterator_t it;
  for (it = v.begin(); it != v.end(); ++it) {
    //rows_[(*it).first].entries_.insert(std::pair<int, NT>(n_, (*it).second));
    set_entry(it->first, n_-1, it->second);
  }

}

template <typename NT>
void
QP_sparse_matrix<NT>::strip_column() {
  // set up data structures
  columns_.pop_back();
  for (int i = 0; i < m_; ++i) {
    rows_[i].clear_entry(n_-1);
    --rows_[i].n_;
  }
  --n_;  
}

template <typename NT>
void
QP_sparse_matrix<NT>::strip_row() {
  // set up data structures
  rows_.pop_back();
  for (int i = 0; i < n_; ++i) {
    columns_[i].clear_entry(m_-1);
    --columns_[i].n_;
  }
  --m_;  
}

template <typename NT>
void
QP_sparse_matrix<NT>::insert_column(int index, const QP_sparse_vector<NT>& v) {

  CGAL_qpe_assertion(index < n_);
  CGAL_qpe_assertion(m_ == v.n_);
  
  std::vector<QP_sparse_vector<int> >::iterator col_it = columns_.begin() + index;
  
  // set up data structures
  columns_.insert(col_it, QP_sparse_vector<int>(m_,-1));
  ++n_;
  
  //QP_sparse_vector<int>::sparse_const_iterator_t index_it, index_it_end;
  for (int i = 0; i < m_; ++i) {
    QP_sparse_vector<int> new_row(n_, -1);
    /*
    index_it = rows_[i].begin();
    index_it_end = rows_[i].end();
    // TODO: efficiency... inserting only instead of set_entry...
    while (index_it != index_it_end) {
      new_row.set_entry((index_it->first < index ? index_it->first : index_it->first+1), index_it->second);
      ++index_it;
    }
    */
    std::transform(rows_[i].begin(), rows_[i].end(), std::inserter(new_row.entries_, new_row.entries_.end()), std::bind1st(Conditional_pair_inc(), index));
    rows_[i] = new_row;
  }
  
  // set entries
  typename QP_sparse_vector<NT>::sparse_const_iterator_t it;
  for (it = v.begin(); it != v.end(); ++it) {
    set_entry(it->first, index, it->second);
  }
   
}

template <typename NT>
void
QP_sparse_matrix<NT>::append_row(const QP_sparse_vector<NT>& v) {
  CGAL_qpe_assertion(n_ == v.n_);
  // set up data structures
  rows_.push_back(QP_sparse_vector<int>(n_,-1));
  for (int i = 0; i < n_; ++i) {
    ++columns_[i].n_;
  }
  ++m_;
  
  // set entries
  typename QP_sparse_vector<NT>::sparse_const_iterator_t it;
  for (it = v.begin(); it != v.end(); ++it) {
    //columns_[(*it).first].entries_.insert(std::pair<int, NT>(m_, (*it).second));
    set_entry(m_-1, it->first, it->second);
  }

}

template <typename NT>
void
QP_sparse_matrix<NT>::insert_row(int index, const QP_sparse_vector<NT>& v) {

  CGAL_qpe_assertion(index < m_);
  CGAL_qpe_assertion(n_ == v.n_);
  
  std::vector<QP_sparse_vector<int> >::iterator row_it = rows_.begin() + index;
  
  // set up data structures  
  rows_.insert(row_it, QP_sparse_vector<int>(n_,-1));
  ++m_;
  
  //QP_sparse_vector<int>::sparse_iterator_t index_it, index_it_end;
  //std::vector<std::pair<int, int> > temp_elements;
  for (int i = 0; i < n_; ++i) {
    QP_sparse_vector<int> new_col(m_, -1);
    //index_it = columns_[i].begin();
    //index_it_end = columns_[i].end();
    //temp_elements.clear();
    // TODO: efficiency... inserting only instead of set_entry...

    std::transform(columns_[i].begin(), columns_[i].end(), std::inserter(new_col.entries_, new_col.entries_.end()), std::bind1st(Conditional_pair_inc(), index));
    columns_[i] = new_col;
  }
  
  // set entries
  typename QP_sparse_vector<NT>::sparse_const_iterator_t it;
  for (it = v.begin(); it != v.end(); ++it) {
    set_entry(index, it->first, it->second);
  }

}


template <typename NT>
void
QP_sparse_matrix<NT>::delete_column(int index) {

  CGAL_qpe_assertion(index < n_);
  
  std::vector<QP_sparse_vector<int> >::iterator col_it = columns_.begin() + index;
  
  columns_.erase(col_it);
  
  /*
  for (int i = 0; i < m_; ++i) {
    //rows_[i].set_entry(index, -1);
    rows_[i].clear_entry(index);
  }*/
  
  //QP_sparse_vector<int>::sparse_const_iterator_t index_it, index_it_end;
  for (int i = 0; i < m_; ++i) {
    QP_sparse_vector<int> new_row(n_-1, -1);
    rows_[i].clear_entry(index);
    /*
    index_it = rows_[i].begin();
    index_it_end = rows_[i].end();
    // TODO: efficiency... inserting only instead of set_entry... (with hints for the map)
    
    while (index_it != index_it_end) {
      new_row.set_entry((index_it->first < index ? index_it->first : index_it->first-1), index_it->second);
      ++index_it;
    }
    */
    std::transform(rows_[i].begin(), rows_[i].end(), std::inserter(new_row.entries_, new_row.entries_.end()), std::bind1st(Conditional_pair_dec(), index));
    rows_[i] = new_row;
  }
  --n_;
}

template <typename NT>
void
QP_sparse_matrix<NT>::delete_row(int index) {

  CGAL_qpe_assertion(index < m_);
  
  std::vector<QP_sparse_vector<int> >::iterator row_it = rows_.begin() + index;
  
  rows_.erase(row_it);
  
  //QP_sparse_vector<int>::sparse_const_iterator_t index_it, index_it_end;
  for (int i = 0; i < n_; ++i) {
    QP_sparse_vector<int> new_col(m_-1, -1);
    columns_[i].clear_entry(index);
    /*
    index_it = columns_[i].begin();
    index_it_end = columns_[i].end();
    // TODO: efficiency... inserting only instead of set_entry... (with hints for the map)
    while (index_it != index_it_end) {
      new_col.set_entry((index_it->first < index ? index_it->first : index_it->first-1), index_it->second);
      ++index_it;
    }
    */
    std::transform(columns_[i].begin(), columns_[i].end(), std::inserter(new_col.entries_, new_col.entries_.end()), std::bind1st(Conditional_pair_dec(), index));
    columns_[i] = new_col;
  }
  --m_; 
}


template <typename NT>
QP_sparse_matrix<NT>& QP_sparse_matrix<NT>::transpose() {
  std::swap(m_, n_);
  std::swap(rows_, columns_);
  return *this;
}

template <typename NT>
boost::shared_ptr<QP_sparse_matrix<rational_t<NT> > >
QP_sparse_matrix<NT>::transform_to_rational() {
  boost::shared_ptr<QP_sparse_matrix<rational_t<NT> > > ret
    (new QP_sparse_matrix<rational_t<NT> >(m_, n_, rational_t<NT>(nt0_)));
  for (int i = 0; i < m_; ++i) {
    typename QP_sparse_vector<int>::sparse_const_iterator_t it;
    for (it = rows_[i].begin(); it != rows_[i].end(); ++it) {
      ret->set_entry(i, (*it).first, rational_t<NT>(data_[(*it).second]));
    }
  }
  return ret;
}

template <typename NT>
std::ostream& QP_sparse_matrix<NT>::print(std::ostream& out) const {
  // TAG: TODO why does bind not work? && replace by bind
  /*
  std::for_each(m.rows_.begin(), m.rows_.end(),
          boost::bind(operator<< <NT>(), _1, out));
  std::for_each(m.rows_.begin(), m.rows_.end(),
          boost::bind1st(operator<< <NT>(), out));
  */      
          
  for (int i = 0; i < rows_.size(); ++i) {
    out << "(";
    for (int j = 0; j < columns_.size(); ++j) {
     out << get_entry(i,j) << ","; 
    }
    out << ")" << std::endl;
  }
  
  return out;
}


template <typename NT>
void QP_sparse_matrix<NT>::swap_rows(int i, int j) {
  std::swap(rows_[i], rows_[j]);
  for (int k = 0; k < n_; ++k) {
    columns_[k].swap_entries(i,j);
  }
}


template <typename NT>
void QP_sparse_matrix<NT>::swap_columns(int i, int j) {
  std::swap(columns_[i], columns_[j]);
  for (int k = 0; k < m_; ++k) {
    rows_[k].swap_entries(i,j);
  }
}

// PRE: This assumes that nt0_ is actually 0
template <typename NT>
NT
QP_sparse_matrix<NT>::inner_product_by_iterators(value_const_iterator begin, value_const_iterator end, const QP_sparse_vector<NT>& v) const {
  NT ret = nt0_;
  typename QP_sparse_vector<NT>::sparse_const_iterator_t itV_begin(v.begin()), itV_end(v.end());
  while (begin != end && itV_begin != itV_end) {
    if (begin->first < itV_begin->first) {
        ++begin;
      } else if (begin->first > itV_begin->first) {
        ++itV_begin;
      } else { // itV->first == begin->first
        ret += *begin->second * itV_begin->second;
        ++begin;
        ++itV_begin;
      }
  }
  return ret;
}

template <typename NT>
void
QP_sparse_matrix<NT>::erase_index(int index) {
    
  for (int i = 0; i < m_; ++i) {
    typename QP_sparse_vector<int>::sparse_iterator_t it = rows_[i].begin();
    typename QP_sparse_vector<int>::sparse_iterator_t it_end = rows_[i].end();
    while (it != it_end) {
      if (it->second > index) --it->second;
      ++it;
    }
  }
  
  for (int i = 0; i < n_; ++i) {
    typename QP_sparse_vector<int>::sparse_iterator_t it = columns_[i].begin();
    typename QP_sparse_vector<int>::sparse_iterator_t it_end = columns_[i].end();
    while (it != it_end) {
      if (it->second > index) --it->second;
      ++it;
    }
  }
}


// =============================================================================
// Operators
// =============================================================================
template <typename NT>
QP_sparse_matrix<NT>&
QP_sparse_matrix<NT>::operator= (const QP_sparse_matrix<NT>& rhs) {
  if (this != &rhs) {
    this->m_ = rhs.m_;
    this->n_ = rhs.n_;
    this->nt0_ = rhs.nt0_;
    this->rows_.resize(this->m_);
    this->columns_.resize(this->n_);
    this->data_.resize(rhs.data_.size());
    std::copy(rhs.rows_.begin(), rhs.rows_.end(), this->rows_.begin());
    std::copy(rhs.columns_.begin(), rhs.columns_.end(), this->columns_.begin());
    std::copy(rhs.data_.begin(), rhs.data_.end(), this->data_.begin());
    this->next_index_ = rhs.next_index_;
  }
  return *this;
}


template <typename NT>
QP_sparse_vector<NT>
QP_sparse_matrix<NT>::operator* (const QP_sparse_vector<NT>& v) const {
  return this->matrix_vector_multiplication(v);
}


template <typename NT>
boost::shared_ptr<QP_sparse_matrix<NT> >
QP_sparse_matrix<NT>::operator* (const QP_sparse_matrix<NT>& m) const {
  return this->matrix_matrix_multiplication(m);
}

template <typename NT>
QP_sparse_matrix<NT>&
QP_sparse_matrix<NT>::operator*= (const QP_sparse_matrix<NT>& m) {
  return *this = *(this->matrix_matrix_multiplication(m));
}

// =============================================================================
// Iterators
// =============================================================================

// value iterators
template <typename NT>
typename QP_sparse_matrix<NT>::value_iterator QP_sparse_matrix<NT>::begin_column_value(int i) {
  return value_iterator(columns_[i].begin(), Index_to_NT(&data_));
}

template <typename NT>
typename QP_sparse_matrix<NT>::value_iterator QP_sparse_matrix<NT>::lower_bound_column_value(int i, int k) {
  return value_iterator(columns_[i].lower_bound(k), Index_to_NT(&data_));
}

template <typename NT>
typename QP_sparse_matrix<NT>::value_iterator QP_sparse_matrix<NT>::end_column_value(int i) {
  return value_iterator(columns_[i].end(), Index_to_NT(&data_));
}

template <typename NT>
typename QP_sparse_matrix<NT>::value_const_iterator QP_sparse_matrix<NT>::begin_column_const_value(int i) const {
  return value_const_iterator(columns_[i].begin(), Index_to_const_NT(&data_));
}

template <typename NT>
typename QP_sparse_matrix<NT>::value_const_iterator QP_sparse_matrix<NT>::lower_bound_column_const_value(int i, int k) const {
  return value_const_iterator(columns_[i].lower_bound(k), Index_to_const_NT(&data_));
}

template <typename NT>
typename QP_sparse_matrix<NT>::value_const_iterator QP_sparse_matrix<NT>::end_column_const_value(int i) const{
  return value_const_iterator(columns_[i].end(), Index_to_const_NT(&data_));
}

template <typename NT>
typename QP_sparse_matrix<NT>::value_iterator QP_sparse_matrix<NT>::begin_row_value(int i) {
  return value_iterator(rows_[i].begin(), Index_to_NT(&data_));
}

template <typename NT>
typename QP_sparse_matrix<NT>::value_iterator QP_sparse_matrix<NT>::lower_bound_row_value(int i, int k) {
  return value_iterator(rows_[i].lower_bound(k), Index_to_NT(&data_));
}

template <typename NT>
typename QP_sparse_matrix<NT>::value_iterator QP_sparse_matrix<NT>::end_row_value(int i) {
  return value_iterator(rows_[i].end(), Index_to_NT(&data_));
}

template <typename NT>
typename QP_sparse_matrix<NT>::value_const_iterator QP_sparse_matrix<NT>::begin_row_const_value(int i) const {
  return value_const_iterator(rows_[i].begin(), Index_to_const_NT(&data_));
}

template <typename NT>
typename QP_sparse_matrix<NT>::value_const_iterator QP_sparse_matrix<NT>::lower_bound_row_const_value(int i, int k) const {
  return value_const_iterator(rows_[i].lower_bound(k), Index_to_const_NT(&data_));
}

template <typename NT>
typename QP_sparse_matrix<NT>::value_const_iterator QP_sparse_matrix<NT>::end_row_const_value(int i) const{
  return value_const_iterator(rows_[i].end(), Index_to_const_NT(&data_));
}



// index iterators
template <typename NT>
typename QP_sparse_matrix<NT>::row_iterator QP_sparse_matrix<NT>::begin_rows_index() {
  return rows_.begin();
}

template <typename NT>
typename QP_sparse_matrix<NT>::row_iterator QP_sparse_matrix<NT>::end_rows_index() {
  return rows_.end();
}

template <typename NT>
typename QP_sparse_matrix<NT>::row_const_iterator
QP_sparse_matrix<NT>::begin_rows_const_index() const {
  return rows_.begin();
}

template <typename NT>
typename QP_sparse_matrix<NT>::row_const_iterator
QP_sparse_matrix<NT>::end_rows_const_index() const {
  return rows_.end();
}

template <typename NT>
typename QP_sparse_matrix<NT>::column_iterator QP_sparse_matrix<NT>::begin_columns_index() {
  return columns_.begin();
}

template <typename NT>
typename QP_sparse_matrix<NT>::column_iterator QP_sparse_matrix<NT>::end_columns_index() {
  return columns_.end();
}

template <typename NT>
typename QP_sparse_matrix<NT>::column_const_iterator
QP_sparse_matrix<NT>::begin_columns_const_index() const {
  return columns_.begin();
}

template <typename NT>
typename QP_sparse_matrix<NT>::column_const_iterator
QP_sparse_matrix<NT>::end_columns_const_index() const {
  return columns_.end();
}


} //namespace CGAL

// ===== EOF ==================================================================