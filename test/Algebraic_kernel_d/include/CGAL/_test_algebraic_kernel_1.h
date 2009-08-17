// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Sebastian Limbach <slimbach@mpi-inf.mpg.de>
//                 Michael Hemmer    <hemmer@mpi-inf.mpg.de> 
//
// ============================================================================

// DONE: test solve_1_object(), etc.
// TODO: Use proper construction of Polynomials 

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Algebraic_kernel_1.h>
#include <CGAL/Test/_test_real_embeddable.h>
#include <CGAL/Test/_test_algebraic_structure.h>
#include <CGAL/Test/_test_coercion_traits.h>
#include <CGAL/Test/_test_polynomial_traits_d.h>


// Test for the Algebraic_kernel syntax
#ifndef CGAL_TEST_ALGEBRAIC_KERNEL_1_H
#define CGAL_TEST_ALGEBRAIC_KERNEL_1_H

CGAL_BEGIN_NAMESPACE


// Test for an exact AlgebraicKernel1
template <class AlgebraicKernel_d_1>
void test_algebraic_kernel_1(const AlgebraicKernel_d_1& ak_1){
  typedef AlgebraicKernel_d_1 Algebraic_kernel_d_1;

  typedef typename AlgebraicKernel_d_1::Coefficient Coefficient;
  typedef typename AlgebraicKernel_d_1::Polynomial_1 Polynomial_1; 
  typedef typename AlgebraicKernel_d_1::Algebraic_real_1 Algebraic_real_1;
  typedef typename AlgebraicKernel_d_1::Bound Bound;
  typedef std::pair<Bound,Bound> BInterval;
  typedef CGAL::Polynomial_traits_d<Polynomial_1> PT;
  
  { 
    // check Coefficient 
    typedef Algebraic_structure_traits<Coefficient> AST;
    typedef typename AST::Algebraic_category Algebraic_category; 
    test_algebraic_structure< Coefficient,Algebraic_category,Tag_true>();
    test_real_embeddable<Coefficient>();
  }
  {
    // check Polynomial_1
    typedef Polynomial_traits_d<Polynomial_1> PT;
    test_polynomial_traits_d(PT());

    // test not possible due to bug in test_algebraic_structure
    // div(3,2)=3/2 != 0 in case of Polynomial<Rational> 
    //  typedef Algebraic_structure_traits<Polynomial_1> AST;
    //  typedef typename AST::Algebraic_category Algebraic_category; 
    //  test_algebraic_structure< Polynomial_1,Algebraic_category,Tag_true>();
  }
  {
    // check Algebraic_real_1
    test_real_embeddable<Algebraic_real_1>();
  }
  { 
    typedef Algebraic_structure_traits<Bound> AST;
    typedef typename AST::Algebraic_category Algebraic_category; 
    test_algebraic_structure< Bound,Algebraic_category,Tag_true>();
    test_real_embeddable<Bound>();
  }
  
  test_explicit_interoperable_from_to<int, Coefficient>();
  test_explicit_interoperable_from_to<int, Bound>();

  test_explicit_interoperable_from_to<int        , Algebraic_real_1>();
  test_explicit_interoperable_from_to<Bound      , Algebraic_real_1>();
  test_explicit_interoperable_from_to<Coefficient, Algebraic_real_1>();
  
#define CGAL_GET_FTOR(Name,name,get_object)             \
  typedef typename AlgebraicKernel_d_1::Name Name;      \
  const Name name = ak_1.get_object();                    
  
  CGAL_GET_FTOR(Is_square_free_1, is_square_free_1, is_square_free_1_object);
  CGAL_GET_FTOR(Make_square_free_1, make_square_free_1, make_square_free_1_object);
  CGAL_GET_FTOR(Square_free_factorize_1, square_free_factorize_1,square_free_factorize_1_object);
  CGAL_GET_FTOR(Is_coprime_1, is_coprime_1,is_coprime_1_object);
  CGAL_GET_FTOR(Make_coprime_1, make_coprime_1, make_coprime_1_object);
  CGAL_GET_FTOR(Solve_1, solve_1, solve_1_object);
  CGAL_GET_FTOR(Sign_at_1, sign_at_1, sign_at_1_object);
  CGAL_GET_FTOR(Compare_1, compare_1, compare_1_object);
  CGAL_GET_FTOR(Bound_between_1, bound_between_1, bound_between_1_object);
  CGAL_GET_FTOR(Approximate_absolute_1, approximate_absolute_1, approximate_absolute_1_object);
  CGAL_GET_FTOR(Approximate_relative_1, approximate_relative_1, approximate_relative_1_object);
#undef CGAL_GET_FTOR

#define CGAL_CHECK_UFUNCTION(Name,AT,RT)                        \
  {                                                             \
    typedef typename Name::argument_type AT_;                   \
    typedef typename Name::result_type   RT_;                   \
    {BOOST_STATIC_ASSERT(( ::boost::is_same<AT,AT_>::value));}  \
    {BOOST_STATIC_ASSERT(( ::boost::is_same<RT,RT_>::value));}  \
  }                
#define CGAL_CHECK_BFUNCTION(Name,AT1,AT2,RT)                           \
  {                                                                     \
    typedef typename Name::first_argument_type AT1_;                    \
    typedef typename Name::second_argument_type AT2_;                   \
    typedef typename Name::result_type   RT_;                           \
    {BOOST_STATIC_ASSERT(( ::boost::is_same<AT1,AT1_>::value));}        \
    {BOOST_STATIC_ASSERT(( ::boost::is_same<AT2,AT2_>::value));}        \
    {BOOST_STATIC_ASSERT(( ::boost::is_same<RT,RT_>::value));}          \
  }                                                             

  CGAL_CHECK_UFUNCTION(Is_square_free_1,Polynomial_1,bool);
  CGAL_CHECK_UFUNCTION(Make_square_free_1,Polynomial_1,Polynomial_1);
  // Square_free_factorize_1 
  CGAL_CHECK_BFUNCTION(Is_coprime_1,Polynomial_1,Polynomial_1,bool);
  // Make_coprime_1
  // Solve_1
  CGAL_CHECK_BFUNCTION(Sign_at_1,Polynomial_1,Algebraic_real_1,Sign);
  CGAL_CHECK_BFUNCTION(Compare_1,Algebraic_real_1,Algebraic_real_1,Sign);
  CGAL_CHECK_BFUNCTION(Bound_between_1,Algebraic_real_1,Algebraic_real_1,Bound);
  CGAL_CHECK_BFUNCTION(Approximate_absolute_1,Algebraic_real_1,int,BInterval);
  CGAL_CHECK_BFUNCTION(Approximate_relative_1,Algebraic_real_1,int,BInterval);
#undef CGAL_CHECK_BFUNCTION
#undef CGAL_CHECK_UFUNCTION
  
  Polynomial_1 x = typename PT::Shift()(Polynomial_1(1),1);
  
  assert( is_square_free_1(ipower((x-1),1)));
  assert(!is_square_free_1(ipower((x-1),2)));

  assert( make_square_free_1(ipower((x-1),2))==ipower((x-1),1));

  {
    std::list< std::pair<Polynomial_1,int> > factors; 
    square_free_factorize_1((x-1)*(x-2)*(x-2),std::back_inserter(factors));
    assert(factors.size()==2);
    assert(factors.front() != factors.back());
    assert(
        factors.front() == std::make_pair((x-1),1) || 
        factors.front() == std::make_pair((x-2),2) );
    assert(
        factors.back()  == std::make_pair((x-1),1) || 
        factors.back()  == std::make_pair((x-2),2) );
  }
  
  assert( is_coprime_1((x-1),(x-2)));
  assert(!is_coprime_1((x-1)*(x-2),(x-1)*(x-3)));
  
  {
    Polynomial_1 a,b,c,d,e;
    a = (x-1)*(x-2);
    b = (x-1)*(x-3);
    assert( make_coprime_1(a,b,c,d,e));
    assert( c == (x-1) ); // gcd 
    assert( d == (x-2) );
    assert( e == (x-3) ); 

    a = (x-1);
    b = (x-2);
    assert(!make_coprime_1(a,b,c,d,e) );
    assert( c == (1) ); // gcd 
    assert( d == (x-1) );
    assert( e == (x-2) );     
  }
  assert(false); // TODO 
//   {
//     std::list<std::pair<Algebraic_real_1,int> > roots; 
//     solve_1((x-1)*(x-2)*(x-2),std::back_inserter(roots));
//     assert(roots.size()==2);
//     assert(roots.front() != roots.back());
//     assert(
//         roots.front() == std::make_pair(Algebraic_real_1(1),1) || 
//         roots.front() == std::make_pair(Algebraic_real_1(2),2) );
//     assert(
//         roots.back()  == std::make_pair(Algebraic_real_1(1),1) || 
//         roots.back()  == std::make_pair(Algebraic_real_1(2),2) );
//   }


//   // Solve_1
//   CGAL_CHECK_BFUNCTION(Sign_at_1,Polynomial_1,Algebraic_real_1,Sign);
//   CGAL_CHECK_BFUNCTION(Compare_1,Algebraic_real_1,Algebraic_real_1,Sign);
//   CGAL_CHECK_BFUNCTION(Bound_between_1,Algebraic_real_1,Algebraic_real_1,Bound);
//   CGAL_CHECK_BFUNCTION(Approximate_absolute_1,Algebraic_real_1,int,BInterval);
//   CGAL_CHECK_BFUNCTION(Approximate_relative_1,Algebraic_real_1,int,BInterval);

}


namespace CGALi {

template< 
class AK_, 
class AlgebraicReal1, 
class Isolator_, 
class Coefficient_, 
class Polynomial1, 
class Bound_  >
void old_test_algebraic_kernel_1() {
    typedef AK_            AK;
    typedef AlgebraicReal1 Algebraic_real_1;
    typedef Isolator_      Isolator;
    typedef Coefficient_   Coefficient;
    typedef Polynomial1    Polynomial_1;
    typedef Bound_      Bound;
        
    BOOST_STATIC_ASSERT( (::boost::is_same< 
            Algebraic_real_1, typename AK::Algebraic_real_1 >::value) );

    BOOST_STATIC_ASSERT((::boost::is_same<
            Isolator,
            typename AK::Isolator >::value) );
            
    BOOST_STATIC_ASSERT((::boost::is_same< 
            Coefficient, 
            typename AK::Coefficient >::value));
            
    BOOST_STATIC_ASSERT((::boost::is_same<
            Polynomial_1,
            typename AK::Polynomial_1 >::value));
    
    // Test of functors
    // Test AK::Solve_1...
    
    typename AK::Solve_1 solve_1 = AK().solve_1_object();
    Polynomial_1 poly1( -4,0,1 );
    Polynomial_1 poly2( 0, 0, 1 );
    std::vector< Algebraic_real_1 > roots_vec;
    std::vector< int > mults_vec;
    
    solve_1( poly1, std::back_inserter( roots_vec ) );        
    assert( roots_vec.size() == 2 );
    assert( CGAL::abs( roots_vec[0] ) == CGAL::abs( Algebraic_real_1(2) ) );
    roots_vec.clear();
    
    solve_1( poly1, std::back_inserter( roots_vec ), std::back_inserter( mults_vec ) );
    assert( roots_vec.size() == 2 );
    assert( mults_vec.size() == 2 );
    assert( CGAL::abs( roots_vec[0] ) == CGAL::abs( Algebraic_real_1(2) ) );
    assert( CGAL::abs( roots_vec[1] ) == CGAL::abs( Algebraic_real_1(2) ) );
    assert( mults_vec[0] == 1 );
    assert( mults_vec[1] == 1 );
    roots_vec.clear();
    mults_vec.clear();

    solve_1( poly2, std::back_inserter( roots_vec ), std::back_inserter( mults_vec ) );
    assert( roots_vec.size() == 1 );
    assert( mults_vec.size() == 1 );
    assert( CGAL::abs( roots_vec[0] ) == CGAL::abs( Algebraic_real_1(0) ) );
    assert( mults_vec[0] == 2 );        
    roots_vec.clear();
    mults_vec.clear();
    
    // Test AK::Sign_at_1
    typename AK::Sign_at_1 sign_at_1;
    typename AK::Polynomial_1 poly4( -2,0,1 );
    solve_1( poly4, std::back_inserter( roots_vec ) );
    typename AK::Polynomial_1 poly3( 0,0,0,1 ); 
    assert( sign_at_1( poly3, roots_vec[0] ) == CGAL::sign( roots_vec[0] ) );
    assert( sign_at_1( poly3, roots_vec[1] ) == CGAL::sign( roots_vec[1] ) );
    assert( sign_at_1( poly3, Algebraic_real_1(0) ) == CGAL::ZERO );  
    roots_vec.clear();
    
    solve_1( poly1, std::back_inserter( roots_vec ) );
    assert( sign_at_1( poly3, roots_vec[0] ) == CGAL::sign( roots_vec[0] ) );
    assert( sign_at_1( poly3, roots_vec[1] ) == CGAL::sign( roots_vec[1] ) );
    assert( sign_at_1( poly3, Algebraic_real_1(0) ) == CGAL::ZERO );  
    roots_vec.clear();
    
    typename AK::Polynomial_1 poly5( 0,0,-1,0,1 );
    typename AK::Algebraic_real_1 algreal1( poly1, Bound(-3), Bound(1) );
    typename AK::Algebraic_real_1 algreal2( poly1, Bound(-1), Bound(3) );
    assert( sign_at_1( poly5, algreal2 ) == CGAL::POSITIVE );
    
    
    // Just syntax tests... (TODO)
    // Test AK::Is_square_free_1...
    typename AK::Is_square_free_1 is_square_free_1 
      = AK().is_square_free_1_object();
    is_square_free_1( poly1 );
    
    // Test AK::Is_coprime_1...
    typename AK::Is_coprime_1 is_coprime_1
      = AK().is_coprime_1_object();
    is_coprime_1( poly1, poly2 );
        
    // Test AK::Make_square_free_1...
    typename AK::Make_square_free_1 make_square_free_1
      = AK().make_square_free_1_object();
    make_square_free_1( poly1 );
    
    // Test AK::Make_coprime_1...
    typename AK::Make_coprime_1 make_coprime_1 
      = AK().make_coprime_1_object();
    Polynomial_1 g, q1, q2;
    make_coprime_1( poly1, poly2, g, q1, q2 );
    
    // Test AK::Square_free_factorize_1...
    typename AK::Square_free_factorize_1 square_free_factorize_1 
      = AK().square_free_factorize_1_object();
    
    std::vector<std::pair<Polynomial_1,int> > fac_mult_pairs;
    square_free_factorize_1( poly1, std::back_inserter(fac_mult_pairs) );
        
    ////////////////////////////////////////////////////////////////////////////
    
    // (Not only) syntax tests for Algebraic_real_traits
            
    // Create test polynomial
    Polynomial_1 p2( -2,0,1 );
    std::vector< Algebraic_real_1 > roots_vec2;
    
    solve_1( p2, std::back_inserter( roots_vec2 ) );
    
    // Test AK::Bound_between...
    typename AK::Bound_between_1 bound_between 
      = AK().bound_between_1_object();
    assert( typename AK::Bound( -2 ) < 
        bound_between( roots_vec2[0], roots_vec2[1] ) );
    assert( typename AK::Bound(  2 ) > 
        bound_between( roots_vec2[0], roots_vec2[1] ) );
    
    // Test AK::Lower_bound
    typename AK::Lower_bound_1 lower_bound 
      = AK().lower_bound_1_object();
    assert( lower_bound( roots_vec2[0] ) < typename AK::Bound(-1) );
    assert( lower_bound( roots_vec2[1] ) < typename AK::Bound( 2) );

    // Test AK::Upper_bound
    typename AK::Upper_bound_1 upper_bound
      = AK().upper_bound_1_object();
    assert( upper_bound( roots_vec2[0] ) > typename AK::Bound(-1) );
    assert( upper_bound( roots_vec2[1] ) > typename AK::Bound( 1) );
    
    // Test AK::Refine
    typename AK::Refine_1 refine
      = AK().refine_1_object();
    Algebraic_real_1 ar = roots_vec2[1];
    typename AK::Bound old_lower_bound = ar.low();
    typename AK::Bound old_upper_bound = ar.high(); 

    refine( ar );
    
    assert( old_lower_bound <= lower_bound( ar ) );
    assert( old_upper_bound >= upper_bound( ar ) );
    typename AK::Bound interval_size_old 
      = CGAL::abs( old_upper_bound - old_lower_bound );
    typename AK::Bound interval_size_new 
      = CGAL::abs( upper_bound( ar ) - lower_bound( ar ) );
    assert( interval_size_new * typename AK::Bound(2) <= interval_size_old );
    
    refine( ar, 100 );
    assert( CGAL::abs( upper_bound( ar ) - lower_bound( ar ) ) < 
        (typename AK::Bound(1) / CGAL::ipower(typename AK::Bound(2), 99 )) );

}

} //namespace CGALi




CGAL_END_NAMESPACE

#endif //CGAL_TEST_ALGEBRAIC_KERNEL_1_H
