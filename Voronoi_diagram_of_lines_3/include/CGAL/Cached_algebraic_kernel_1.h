// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//                 Sebastian Limbach <slimbach@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_CHACHED_ALGEBRAIC_KERNEL_1_H
#define CGAL_CHACHED_ALGEBRAIC_KERNEL_1_H

#include <CGAL/basic.h>
#include <map>

namespace CGAL {

// sofar this is a pure wrapper 
template< class AlgebraicKernel_1 > 
class Cached_algebraic_kernel_1{
  typedef AlgebraicKernel_1 AK_1;
public:
  typedef AlgebraicKernel_1 Algebraic_kernel_d_1;
  typedef Cached_algebraic_kernel_1<Algebraic_kernel_d_1> Self;
  typedef typename Algebraic_kernel_d_1::Coefficient      Coefficient; 
  typedef typename Algebraic_kernel_d_1::Polynomial_1     Polynomial_1; 
  typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real_1;
  typedef typename Algebraic_kernel_d_1::Boundary         Boundary;
  
#define CGAL_ALGEBRAIC_KERNEL_1_PRED(Y,Z)       \
  typedef typename Algebraic_kernel_d_1::Y Y;     \
    Y Z() const { return Y(); }
  
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Is_square_free_1,
      is_square_free_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Make_square_free_1,
      make_square_free_1_object);
//  CGAL_ALGEBRAIC_KERNEL_1_PRED(Square_free_factorize_1,
//      square_free_factorize_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Is_coprime_1,
      is_coprime_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Make_coprime_1,
      make_coprime_1_object);
//  CGAL_ALGEBRAIC_KERNEL_1_PRED(Solve_1,
//      solve_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Sign_at_1,
      sign_at_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Compare_1,
      compare_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Refine_1,
      refine_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Lower_boundary_1,
      lower_boundary_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Upper_boundary_1,
      upper_boundary_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Boundary_between_1,
      boundary_between_1_object);
      
#undef CGAL_ALGEBRAIC_KERNEL_1_PRED  

private:
  Algebraic_kernel_d_1 m_algebraic_kernel_1;
public:
  const Algebraic_kernel_d_1& algebraic_kernel_1() const {
    return m_algebraic_kernel_1; 
  };
  Algebraic_kernel_d_1&       algebraic_kernel_1(){
    return m_algebraic_kernel_1; 
  }

public:
  Cached_algebraic_kernel_1(){};
  Cached_algebraic_kernel_1(const Algebraic_kernel_d_1& ak_1)
    :m_algebraic_kernel_1(ak_1){};
  virtual ~Cached_algebraic_kernel_1(){};


private:
  typedef std::list<std::pair<Polynomial_1,int> > SQFF; 
  typedef boost::optional<SQFF> Sqff_cache_data;
  typedef std::map<Polynomial_1, Sqff_cache_data > Sqff_cache;
  mutable Sqff_cache m_sqff_cache;
  
  SQFF& get_sqff(const Polynomial_1& p) const {
    Sqff_cache_data& data = m_sqff_cache[p];
    if(!data.is_initialized()){
      data = Sqff_cache_data(std::list<std::pair<Polynomial_1,int> >());
      m_algebraic_kernel_1.square_free_factorize_1_object()
        (p,std::back_inserter(*data));
    }
    return *data;
  }
  
public:
  struct Square_free_factorize_1{
    typedef Polynomial_1 first_argument_type;
    const Cached_algebraic_kernel_1* m_cak_1;
    Square_free_factorize_1(const Cached_algebraic_kernel_1 * cak_1):m_cak_1(cak_1){}

    template <typename OutputIterator>
    OutputIterator operator()(const Polynomial_1& p, OutputIterator oit){
      SQFF& sqff = m_cak_1->get_sqff(p);
      std::copy(sqff.begin(),sqff.end(),oit);
    }
  };
  
  Square_free_factorize_1 square_free_factorize_1_object() const { 
    return Square_free_factorize_1(this); 
  }
  
private:
  typedef std::list<Algebraic_real_1 > ROOTS; 
  typedef boost::optional<ROOTS> Roots_cache_data;
  typedef std::map<Polynomial_1, Roots_cache_data > Roots_cache;
  mutable Roots_cache m_roots_cache;
  
  ROOTS& get_roots(const Polynomial_1& p) const {
    CGAL_assertion(CGAL::is_square_free(p));
    Roots_cache_data& data = this->m_roots_cache[p];
    if(!data.is_initialized()){
      data = Roots_cache_data(ROOTS());
      algebraic_kernel_1().solve_1_object()
        (p,std::back_inserter(*data),true);
    }
    return *data;
  }
  
public:
   struct Solve_1{
    typedef Polynomial_1 argument_type;

    const Cached_algebraic_kernel_1 *m_cak_1;
    Solve_1(const Cached_algebraic_kernel_1 *cak_1):m_cak_1(cak_1){}

    template <typename OutputIterator>
    OutputIterator operator()(const Polynomial_1& p, OutputIterator oit, 
        bool known_to_be_square_free = false){
      
      if(known_to_be_square_free){
        ROOTS& roots = m_cak_1->get_roots(p);
        return std::copy(roots.begin(),roots.end(),oit);
      }
      
      SQFF& sqff = m_cak_1->get_sqff(p);
      for(typename SQFF::iterator it = sqff.begin(); it != sqff.end(); it++){
        ROOTS roots = m_cak_1->get_roots(it->first);
        oit = std::copy(roots.begin(),roots.end(),oit);
      }
      return oit; 
    } 

     template< class OutputIteratorRoots, class OutputIteratorMults >
     std::pair< OutputIteratorRoots, OutputIteratorMults >
     operator()( 
         const Polynomial_1& p, 
         OutputIteratorRoots roit,  
         OutputIteratorMults moit ) const {
       SQFF& sqff = m_cak_1->get_sqff(p);
       for(typename SQFF::iterator it = sqff.begin(); it != sqff.end(); it++){
         ROOTS roots = m_cak_1->get_roots(it->first);
         for(typename ROOTS::iterator rit = roots.begin(); rit != roots.end(); rit++){
           *roit++ = *rit;
           *moit++ = it->second; 
         }
       }
       return std::make_pair(roit,moit);
     }
  };
  
  Solve_1 solve_1_object() const { 
    return Solve_1(this); 
  }

};


} //namespace CGAL



#endif // CGAL_CHACHED_ALGEBRAIC_KERNEL_1_H
