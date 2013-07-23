#ifndef SV_FILTER_ITERATOR_H
#define SV_FILTER_ITERATOR_H

namespace SV {
namespace internal {


template <class Container>
struct Contains{
  Container const*  m_container; 
  // Contains(){};
  Contains(Container const* container):m_container(container){};
  template< class Arg>
  bool operator()(const Arg& arg) const {
    return m_container->contains(arg);
  }
};

template < class OI, class P > struct Filter_iterator;

template < class OI, class P >
Filter_iterator< OI, P >
make_filter_iterator(OI e, const P& p)
{ return Filter_iterator<OI,P>(e,p);}


// template < class OI, class P >
// bool operator==(const Filter_iterator<OI,P>&, const Filter_iterator<OI,P>&);

template < class OI, class P >
struct Filter_iterator {
  typedef OI                                OIterator;
  typedef P                                 Predicate;
  typedef Filter_iterator<OI,P>             Self;
  typedef std::iterator_traits<OI>          OITI;
  typedef Self&                             reference;
  typedef Self*                             pointer;
//typedef value_type                        value_type;
//typedef typename OITI::difference_type    difference_type;
//  typedef typename OITI::iterator_category  iterator_category;
  // Special for circulators.
//  typedef OI_Circulator_size_traits<iterator_category,OI> C_S_Traits;
//  typedef typename  C_S_Traits::size_type               size_type;

protected:
  OIterator c_;      // current position.
  Predicate p_;      // Leave out x <==> p_(x).
public:

  Filter_iterator() {}

  Filter_iterator(OIterator c, const Predicate& p)
  : c_(c), p_(p) {}

  

  template <class ARG>
  Self& operator=(const ARG& x){
    if(!this->p_(x)) *(this->c_)++=x;
    return *this; 
  }
  
  Self& operator=(const Self& s)
  {
    this->p_=s.p_;
    this->c_=s.c_;
    return *this;
  } 
  
  Self& operator++() { return *this; }
  Self& operator++(int) { return *this; }
  Self& operator*() { return *this; }
  
  // friend bool operator== <>(const Self&, const Self&);
};

// template < class OI, class P >
// inline
// bool operator==(const Filter_iterator<OI,P>& it1,
//                 const Filter_iterator<OI,P>& it2)
// {
//   CGAL_precondition(it1.e_ == it2.e_);
//   return it1.base() == it2.base();
// }

// template < class OI, class P >
// inline
// bool operator!=(const Filter_iterator<OI,P>& it1,
//                 const Filter_iterator<OI,P>& it2)
// { return !(it1 == it2); }
















// template < class OutputIterator, class BooleanPredicate> 
// struct Filter_iterator{
//   typedef Filter_iterator Self; 
//   OutputIterator m_oit; 
//   const BooleanPredicate& m_pred;
//   Filter_iterator(OutputIterator oit, const BooleanPredicate& pred):m_oit(oit),m_pred(pred){};
 
//   Filter_iterator& operator++(){return *this;}
//   Filter_iterator& operator++(int){return *this;}
//   Filter_iterator& operator*(){return *this;} 
  
//   template< class Arg> 
//   Filter_iterator& operator=(const Arg& x){
//     if(!m_pred(x)) *m_oit++=x;
//     return *this; 
//   }
// };  

// template <class OutputIterator, class BooleanPredicate> 
// Filter_iterator<OutputIterator,BooleanPredicate>
// make_filter_iterator(OutputIterator oit, const BooleanPredicate& pred){
//   return  Filter_iterator<OutputIterator,BooleanPredicate>(oit,pred);
// }


} // namespace internal 
} // namespace SV
#endif // SV_FILTER_ITERATOR_H
