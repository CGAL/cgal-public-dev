// Copyright (c) 2007 GeometryFactory (France), Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/STL_Extension/include/CGAL/Iterator_transform.h $
// $Id: Iterator_transform.h 46206 2008-10-11 20:21:08Z spion $
// 
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion
//                 Fernando Cacciola <fernando.cacciola@geometryfactory.com> 
//                 Michael Hemmer

#ifndef CGAL_TRANSFORMED_OUTPUT_ITERATOR_H
#define CGAL_TRANSFORMED_OUTPUT_ITERATOR_H 1


namespace CGAL {

template < class OutputIterator, class Creator>
class Transformed_output_iterator {
protected:
  OutputIterator oit;    // The internal iterator.
  Creator        creator;          // The internal iterator.
public:
  typedef Transformed_output_iterator<OutputIterator,Creator>  Self;
  typedef OutputIterator                                       OIT;
  typedef std::iterator_traits<OIT>                            OITT;
  typedef typename OITT::iterator_category  iterator_category;
  typedef typename OITT::pointer            pointer;   // void
  typedef typename OITT::reference          reference; // void 
  typedef typename OITT::difference_type    difference_type; // void 

  typedef typename Creator::argument_type         value_type;
  typedef typename Creator::argument_type         argument_type;
  typedef typename Creator::result_type           result_type;

  // CREATION
  // --------
  Transformed_output_iterator(){};
  Transformed_output_iterator(OutputIterator oit_, Creator creator_) : oit(oit_),creator(creator_) {};

  // OPERATIONS Forward Category
  // ---------------------------

  OutputIterator  current_iterator() const { return oit;}
  bool      operator==( const Self& i) const { return ( oit == i.oit); }
  bool      operator!=( const Self& i) const { return !(*this == i); }
  
  struct Proxy
  {
    Self it; // should it be & ?
    Proxy(Self it_) : it(it_){}
    Proxy operator = (const argument_type& x ){
      *(it.oit) = it.creator(x);
    }  
  };

//  Proxy operator->() const
//  {
//      return Proxy(creator(*nt));
//  }
  
  Proxy operator* () const {
    return Proxy(*this);
  }
  
  Self&     operator++() {
    ++oit;
    return *this;
  }
  Self      operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }
};

template<class OutputIterator, class Creator>
Transformed_output_iterator<OutputIterator, Creator> 
make_transformed_output_iterator(OutputIterator oi, Creator creator){
  return Transformed_output_iterator<OutputIterator, Creator>(oi,creator);
}



} //namespace CGAL
#endif // CGAL_TRANSFORMED_OUTPUT_ITERATOR_H //
// EOF //
