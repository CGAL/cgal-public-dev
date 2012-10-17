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



// TAG TODO: worry about mem management and resize stuff
// TAG COMMENT: if default value is not NT(0) then inner_product is not sparse
// any more

#ifndef CGAL_QP_SPARSE_VECTOR_H
#define CGAL_QP_SPARSE_VECTOR_H



#include <CGAL/QP_solver/basic.h>
#include <CGAL/QP_solver/functors.h>

#include <ostream>
#include <map>
#include <utility>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

namespace CGAL {

// =============================================================================
// declarations
// =============================================================================

template <typename NT>
class QP_sparse_vector;
template <typename NT>
class QP_sparse_matrix;


template <typename NT>
std::ostream& operator<<(std::ostream& out, const QP_sparse_vector<NT>& v);



// =============================================================================
// class interface
// =============================================================================
template <typename NT>
class QP_sparse_vector {

  // ===========================================================================
  // friends
  // ===========================================================================
  template <typename S> friend class QP_sparse_matrix;
  
  friend
  std::ostream& operator<< <>(std::ostream& out, const QP_sparse_vector<NT>& v);

  public:
  
    // =========================================================================
    // typedefs
    // =========================================================================
    typedef int key_type;
    typedef NT value_type;
    

  private:
    // =========================================================================
    // typedefs
    // =========================================================================
    typedef typename std::map<key_type, NT> map_t;
    typedef typename std::pair<key_type, NT> index_value_pair_t;
  
  

  public:

    // =========================================================================
    // Construction
    // =========================================================================
    QP_sparse_vector (); 
    QP_sparse_vector (key_type n, NT default_value = NT(0));
    QP_sparse_vector (const QP_sparse_vector& v);
    
    // =========================================================================
    // Info functions
    // =========================================================================
    int get_size() const;
    int get_number_of_nonzeros() const;
    
    // =========================================================================
    // Getter & setter methods
    // =========================================================================
    const NT& get_entry(key_type n) const;
    void set_entry(key_type n, NT val);
    template <typename Iterator> void insert(Iterator begin, Iterator end);
    void clear_entry(key_type n);
    void swap_entries(key_type i, key_type j);
    
    
    // =========================================================================
    // Operations
    // =========================================================================
    NT inner_product(const QP_sparse_vector& v) const;
    QP_sparse_vector<NT>& permute(const std::vector<int>& permutaion);
    
    // PRE: this != &v
    void append(const QP_sparse_vector& v);
    
    void clear();
  
    std::ostream& print(std::ostream& out) const;
    
    // =========================================================================
    // Operators
    // =========================================================================
    const NT& operator[] (key_type n) const;
    NT operator* (const QP_sparse_vector<NT>& v) const;
    QP_sparse_vector<NT>& operator= (const QP_sparse_vector<NT>& rhs);
    
    // =========================================================================
    // Iterators
    // =========================================================================
    typedef typename boost::transform_iterator<CGAL::Map_with_default<map_t>,
				    boost::counting_iterator<int> >
      dense_iterator_t;
	  typedef typename map_t::iterator sparse_iterator_t;
	  typedef typename map_t::const_iterator sparse_const_iterator_t;
	  
	  dense_iterator_t begin_dense();
	  sparse_iterator_t begin();
	  sparse_iterator_t end();
	  sparse_iterator_t lower_bound(const key_type i);
	  sparse_const_iterator_t begin() const;
	  sparse_const_iterator_t end() const;
	  sparse_const_iterator_t lower_bound(const key_type i) const;
	  
				    
    
    
  private:
    // =========================================================================
    // data members
    // =========================================================================
    int n_; // the size of the vector
    NT nt0_; // default value
    map_t entries_; // used to store the entries of the vector
};


// =============================================================================
// Non-member functions
// =============================================================================

template <typename NT>
std::ostream& operator<<(std::ostream& out, const QP_sparse_vector<NT>& v) {
  return v.print(out);
}



// =============================================================================
// class implementation (inline)
// =============================================================================

} //namespace CGAL

#include <CGAL/QP_solver/QP_sparse_vector_impl.h>

#endif //CGAL_QP_SPARSE_VECTOR_H

// ===== EOF ==================================================================