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

#ifndef CGAL_QP_SPARSE_MATRIX_H
#define CGAL_QP_SPARSE_MATRIX_H


#include <CGAL/QP_solver/basic.h>
#include <CGAL/QP_solver/functors.h>
#include <CGAL/QP_solver/QP_sparse_vector.h>

#include <ostream>
#include <map>
#include <utility>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>

namespace CGAL {

// =============================================================================
// helper objects
// =============================================================================
  
// Used to manipulate the index sets
struct Conditional_pair_inc : public std::binary_function<int,
                                                     std::pair<int, int>,
                                                     std::pair<int, int> > {
  std::pair<int,int> operator() (int k, std::pair<int,int> a) const  {
    return (a.first < k ? a : std::pair<int,int>(a.first+1, a.second));
  }
};

struct Conditional_pair_dec : public std::binary_function<int,
                                                     std::pair<int, int>,
                                                     std::pair<int, int> > {
  std::pair<int,int> operator() (int k, std::pair<int,int> a) const  {
    return (a.first < k ? a : std::pair<int,int>(a.first-1, a.second));
  }
};


// =============================================================================
// declarations
// =============================================================================

template <typename NT>
class QP_sparse_vector;
template <typename NT>
class QP_sparse_matrix;

template <typename NT>
std::ostream& operator<<(std::ostream& out, const QP_sparse_matrix<NT>& m);

template <typename NT>
struct rational_t {
  rational_t(): n_(NT(0)), d_(NT(1)) {};
  rational_t(NT val): n_(val), d_(NT(1)) {};
  
  bool operator== (const rational_t& b) const {
    return b.d_ * n_ == b.n_ * d_;
  };
  
  bool operator!= (const rational_t& b) const {
    return !(*this == b);
  };
  
  NT n_;
  NT d_;
};



// =============================================================================
// class interface
// =============================================================================
template <typename NT>
class QP_sparse_matrix {

  
  
  // ===========================================================================
  // friend functions
  // ===========================================================================
  friend
  std::ostream& operator<< <>(std::ostream& out, const QP_sparse_matrix<NT>& m);

  private:
  // =============================================================================
  // Accessor
  // =============================================================================
  struct Index_to_const_NT {
    typedef typename std::pair<int,int> input_type;
    typedef typename std::pair<int, const NT*> result_type;
    
    Index_to_const_NT(): data_pointer_(0) {}
    Index_to_const_NT(const std::vector<NT>* p): data_pointer_(p) {}
    
    result_type operator() (const input_type& p) const {
      return result_type(p.first, &data_pointer_->at(p.second));
    }
    
    private:
      const std::vector<NT>* data_pointer_;
  };
  
  // TODO: return reference instead of pointers?
  struct Index_to_NT {
    typedef typename std::pair<int,int> input_type;
    typedef typename std::pair<int, NT*> result_type;
    
    Index_to_NT(): data_pointer_(0) {}
    Index_to_NT(std::vector<NT>* p): data_pointer_(p) {}
    
    result_type operator() (const input_type& p) const {
      return result_type(p.first, &(data_pointer_->at(p.second)));
    }
    
    private:
      std::vector<NT>* data_pointer_;
  };

  private:
  
    // =========================================================================
    // typedefs
    // =========================================================================
    typedef typename std::map<int, NT> map_t;
  

  public:
    
    // =========================================================================
    // typedefs
    // =========================================================================
    typedef typename std::pair<int, int> index_pair_t;
    typedef typename std::pair<index_pair_t, NT> index_value_pair_t;
  
    // =========================================================================
    // Construction
    // =========================================================================
    QP_sparse_matrix ();
    QP_sparse_matrix (int m, int n, NT default_value = NT(0));
    QP_sparse_matrix (const QP_sparse_matrix<NT>& m);

    // =========================================================================
    // Info functions
    // =========================================================================
    index_pair_t get_size() const;
    
    // =========================================================================
    // Getter & setter methods
    // =========================================================================
    
    const NT& get_entry(int m, int n) const;
    const NT& get_entry(index_pair_t p) const;
    bool set_entry(int m, int n, const NT& val);
    
    const QP_sparse_vector<int>& get_row_indices(int m) const;
    const QP_sparse_vector<int>& get_column_indices(int n) const;
    QP_sparse_vector<int>& get_row_indices(int m);
    QP_sparse_vector<int>& get_column_indices(int n);
    
    const NT& get_element_by_index(int index) const;
    NT* get_element_pointer_by_index(int index);
    void set_element_by_index(int index, const NT& val);
    
    double get_density() const;
    
    // =========================================================================
    // Operations
    // =========================================================================
    // default values are inherited from the this pointer
    QP_sparse_vector<NT>
      matrix_vector_multiplication (const QP_sparse_vector<NT>& v) const;
    boost::shared_ptr<QP_sparse_matrix<NT> >
      matrix_matrix_multiplication (const QP_sparse_matrix<NT>& m) const;
    
    NT inner_product_row(int i, const QP_sparse_vector<NT>& v) const;
    NT inner_product_column(int j, const QP_sparse_vector<NT>& v) const;
     
    
    void append_column(const QP_sparse_vector<NT>& v);
    void append_row(const QP_sparse_vector<NT>& v);
    
    void strip_column();
    void strip_row();
  
    //void append_columns(const QP_sparse_matrix<NT>& m);
    //void append_rows(const QP_sparse_matrix<NT>& m);
            
    void insert_column(int index, const QP_sparse_vector<NT>& v);
    void insert_row(int index, const QP_sparse_vector<NT>& v);
    
    void delete_column(int index);
    void delete_row(int index);
    

    QP_sparse_matrix<NT>& transpose ();
    
    boost::shared_ptr<QP_sparse_matrix<rational_t<NT> > >
    transform_to_rational();
  
    std::ostream& print(std::ostream& out) const;
    
    
    void swap_rows(int i, int j);
    void swap_columns(int i, int j);
    
    QP_sparse_matrix<NT>& apply_column_permutation(const std::vector<int>& perm);

    
    
    // =========================================================================
    // Operators
    // =========================================================================
  public:
    QP_sparse_vector<NT> operator* (const QP_sparse_vector<NT>& v) const;
    boost::shared_ptr<QP_sparse_matrix<NT> >
      operator* (const QP_sparse_matrix<NT>& m) const;
    QP_sparse_matrix<NT>& operator*= (const QP_sparse_matrix<NT>& m);
    QP_sparse_matrix<NT>& operator= (const QP_sparse_matrix<NT>& rhs);
    
    // =========================================================================
    // Iterators
    // =========================================================================
    typedef typename std::vector<QP_sparse_vector<int> >::iterator   row_iterator;
    typedef typename std::vector<QP_sparse_vector<int> >::iterator   column_iterator;
    typedef typename std::vector<QP_sparse_vector<int> >::const_iterator   row_const_iterator;
    typedef typename std::vector<QP_sparse_vector<int> >::const_iterator   column_const_iterator;
    typedef typename QP_sparse_vector<int>::sparse_iterator_t        index_iterator;
    typedef typename QP_sparse_vector<int>::sparse_const_iterator_t  index_const_iterator;
    typedef typename boost::transform_iterator<Index_to_NT, index_iterator> value_iterator;
    typedef typename boost::transform_iterator<Index_to_const_NT, index_const_iterator> value_const_iterator;
    
    
    // TODO: can we make const vs non-const selection implicit?
    value_iterator begin_column_value(int i);
    value_iterator lower_bound_column_value(int i, int k);
    value_iterator end_column_value(int i);
    value_const_iterator begin_column_const_value(int i) const;
    value_const_iterator lower_bound_column_const_value(int i, int k) const;
    value_const_iterator end_column_const_value(int i) const;
    
    value_iterator begin_row_value(int i);
    value_iterator lower_bound_row_value(int i, int k);
    value_iterator end_row_value(int i);
    value_const_iterator begin_row_const_value(int i) const;
    value_const_iterator lower_bound_row_const_value(int i, int k) const;
    value_const_iterator end_row_const_value(int i) const;
    
    
    
	  row_iterator begin_rows_index();
	  row_iterator end_rows_index();
	  row_const_iterator begin_rows_const_index() const;
	  row_const_iterator end_rows_const_index() const;
		column_iterator begin_columns_index();
	  column_iterator end_columns_index();
	  column_const_iterator begin_columns_const_index() const;
	  column_const_iterator end_columns_const_index() const;
    
    
    // =========================================================================
    // Static funtions
    // =========================================================================
    static boost::shared_ptr<QP_sparse_matrix<NT> > get_identity(const int n) {
      boost::shared_ptr<QP_sparse_matrix<NT> > m( new QP_sparse_matrix<NT>(n, n, NT(0)) );
      for (int i = 0; i < n; ++i) {
        m->set_entry(i, i, NT(1));
      }
      return m;
    }
  
    
    // =========================================================================
    // Operations
    // =========================================================================
  private:
  
    NT inner_product_by_iterators(value_const_iterator begin, value_const_iterator end, const QP_sparse_vector<NT>& v) const;
    
    void erase_index(int index);
    
  private:
    // =========================================================================
    // data members
    // =========================================================================
    int m_; // row dimension
    int n_; // column dimension
    NT nt0_;
    std::vector<QP_sparse_vector<int> > rows_, columns_;
    
    // new data layout
    std::vector<NT> data_;
    int next_index_;
};



// =============================================================================
// Non-member functions
// =============================================================================

template <typename NT>
std::ostream& operator<<(std::ostream& out, const QP_sparse_matrix<NT>& m) {
  return m.print(out);
}


template <typename NT>
std::ostream& operator<<( std::ostream& out,
                         const rational_t<NT>& r) {
  return (out << "(" << r.n_ << "/" << r.d_ << ")");
}



// =============================================================================
// class implementation (inline)
// =============================================================================

} //namespace CGAL

#include <CGAL/QP_solver/QP_sparse_matrix_impl.h>

#endif //CGAL_QP_SPARSE_MATRIX_H

// ===== EOF ==================================================================