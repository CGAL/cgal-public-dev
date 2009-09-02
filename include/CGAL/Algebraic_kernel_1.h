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

#ifndef CGAL_ALGEBRAIC_KERNEL_1_H
#define CGAL_ALGEBRAIC_KERNEL_1_H

#ifndef CGAL_AK_ENABLE_DEPRECATED_INTERFACE 
#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 0
#endif 

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_pure.h>
#include <CGAL/Algebraic_kernel_d/Descartes.h>
#include <CGAL/Algebraic_kernel_d/Real_roots.h>
#include <CGAL/Algebraic_kernel_d/refine_zero_against.h>
#include <CGAL/ipower.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {
template< class AlgebraicReal1, class Isolator_ >
class Algebraic_kernel_1_base {
  
public:
  typedef AlgebraicReal1                              Algebraic_real_1;
  typedef Isolator_                                   Isolator;
            
  typedef typename Algebraic_real_1::Coefficient      Coefficient;
  typedef typename Algebraic_real_1::Rational         Bound;
  typedef typename 
      CGAL::Polynomial_type_generator< Coefficient,1 >::Type Polynomial_1;

private:
  typedef CGAL::Polynomial_traits_d< Polynomial_1 >   PT_1;

public:    
  class Algebraic_real_traits {
  public:
    typedef Algebraic_real_1                      Type;
    typedef typename Algebraic_real_1::Rational   Bound;
                
    struct Bound_between 
      : public std::binary_function< Type, Type, Bound > {
      Bound operator()( const Type& t1, 
          const Type& t2 ) const {
        return t1.rational_between( t2 );
      }
    };
                                
    struct Lower_bound
      : public std::unary_function< Type, Bound > {
      Bound operator()( const Type& t ) const {
        return t.low();
      }
    };
                
    struct Upper_bound
      : public std::unary_function< Type, Bound > {
      Bound operator()( const Type& t ) const {
        return t.high();
      }
    };
                
    struct Refine
      : public std::unary_function< Type, void > {
      void operator()( const Type& t ) const {
        t.refine();
      }
                    
      void operator()( Type& t, int rel_prec ) const {
        // If t is zero, we can refine the interval to
        //  infinite precission
        if( CGAL::is_zero( t ) ) {
          t = Type(0);                            
        } else {
          // Refine until both boundaries have the same sign
          while( CGAL::sign( t.high() ) != 
              CGAL::sign( t.low() ) )
            t.refine();
                            
          CGAL_assertion( CGAL::sign( t.high() ) != CGAL::ZERO &&
              CGAL::sign( t.low() ) != CGAL::ZERO );
                            
          // Calculate the needed precision
          Bound prec = Bound(1) / 
            CGAL::ipower( Bound(2), rel_prec );
                            
          // Refine until precision is reached
          while( CGAL::abs( t.high() - t.low() ) /
              CGAL::max( CGAL::abs( t.high() ),
                  CGAL::abs( t.low() ) ) > prec ) {
            t.refine();

            CGAL_assertion( CGAL::sign( t.high() ) != CGAL::ZERO &&
                CGAL::sign( t.low() ) != CGAL::ZERO );
                                
          } 
        }                        
      }
    };              
    
    struct Approximate_absolute_1:
      public std::binary_function<Algebraic_real_1,int,std::pair<Bound,Bound> >{
      std::pair<Bound,Bound> 
      operator()(const Algebraic_real_1& x, int prec) const {
        Lower_bound lower; 
        Upper_bound upper; 
        Refine refine; 
        Bound l = lower(x);  
        Bound u = upper(x);
        Bound error = CGAL::ipower(Bound(2),CGAL::abs(prec));
        while((prec>0)?((u-l)*error>Bound(1)):((u-l)>error)){
          refine(x);
          u = upper(x);
          l = lower(x);
        }
        return std::make_pair(l,u);
      }
    }; 
    
    struct Approximate_relative_1:
      public std::binary_function<Algebraic_real_1,int,std::pair<Bound,Bound> >{
      std::pair<Bound,Bound> 
      operator()(const Algebraic_real_1& x, int prec) const {
        
        if(CGAL::is_zero(x)) return std::make_pair(Bound(0),Bound(0));
        
        Lower_bound lower; 
        Upper_bound upper; 
        Refine refine; 
        Bound l = lower(x);  
        Bound u = upper(x);
        Bound error = CGAL::ipower(Bound(2),CGAL::abs(prec));
        Bound max_b = (CGAL::max)(CGAL::abs(u),CGAL::abs(l));
        while((prec>0)?((u-l)*error>max_b):((u-l)>error*max_b)){
          refine(x);
          u = upper(x);
          l = lower(x);
          max_b = (CGAL::max)(CGAL::abs(u),CGAL::abs(l));
        }
        return std::make_pair(l,u);
      }
    };

#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
    typedef Lower_bound Lower_boundary;
    typedef Upper_bound Upper_boundary;
    typedef Bound_between Boundary_between; 
#endif    


  }; // class Algebraic_real_traits
            
  // Functors of Algebraic_kernel_1
  struct Solve_1 {
  public:
    template <class OutputIterator>
    OutputIterator 
    operator()(const Polynomial_1& p, OutputIterator oi) const {  
#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
#else
      CGAL_precondition(!CGAL::is_zero(p));
#endif
      CGALi::Real_roots< Algebraic_real_1, Isolator > real_roots; 
      std::list< int > mults;
      std::list< Algebraic_real_1 > roots; 
      real_roots( p, std::back_inserter(roots), std::back_inserter( mults ) );
      CGAL_assertion(roots.size()==mults.size());
      std::list<int>::iterator mit =mults.begin();
      typename std::list< Algebraic_real_1 >::iterator rit = roots.begin();
      while(rit != roots.end()){
        //*oi++ = std::make_pair(*rit, (unsigned int)(*mit));
        *oi++ = std::make_pair(*rit, *mit);
        rit++;
        mit++;
      }
      return oi;
    }
    
    template< class OutputIterator >
    OutputIterator operator()( 
        const Polynomial_1& p, 
        OutputIterator oi , 
        bool known_to_be_square_free) const {

      CGALi::Real_roots< Algebraic_real_1, Isolator > real_roots; 
#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
#else
      CGAL_precondition(!CGAL::is_zero(p));
#endif
      std::list<Algebraic_real_1> roots; 
      if( known_to_be_square_free ){
        real_roots(p,std::back_inserter(roots));
      }else{
        std::list<int> dummy;
        real_roots(p,std::back_inserter(roots),std::back_inserter(dummy));
      }
      return std::copy(roots.begin(),roots.end(),oi);
    }

#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
    template< class OutputIteratorRoots , class OutputIteratorMults >
    std::pair<OutputIteratorRoots,OutputIteratorMults> 
    operator()( 
        const Polynomial_1& p, 
        OutputIteratorRoots roi,
        OutputIteratorMults moi) const {

      CGALi::Real_roots< Algebraic_real_1, Isolator > real_roots; 
      real_roots(p,roi,moi);
      return std::make_pair(roi,moi);
    }
#endif
  };
            
  struct Sign_at_1 
    : public std::binary_function< Polynomial_1, Algebraic_real_1, CGAL::Sign > {
    CGAL::Sign operator()( const Polynomial_1& p, const Algebraic_real_1& ar ) const {
      if(CGAL::is_zero(p)) return ZERO; 
      if(CGAL::degree(p)==0) return p.sign_at(0);
      if( ar.low() == ar.high() ) return p.sign_at( ar.low() );

      Polynomial_1 g = gcd_utcf(p,ar.polynomial());
      if (g.sign_at(ar.low()) != g.sign_at(ar.high())) return ZERO; 
     
      while(CGALi::descartes(p,ar.low(),ar.high()) > 0) ar.refine();
      while( p.sign_at(ar.low())  == ZERO )  ar.refine();
      while( p.sign_at(ar.high()) == ZERO )  ar.refine();

      CGAL::Sign result = p.sign_at(ar.low());
      CGAL_assertion(result == p.sign_at(ar.high()));
      return result; 
    }
  };                
  struct Is_zero_at_1 
    : public std::binary_function< Polynomial_1, Algebraic_real_1, bool > {
    bool operator()( const Polynomial_1& p, const Algebraic_real_1& ar ) const {
      if(CGAL::is_zero(p)) return true; 
      if( ar.low() == ar.high() ) return p.sign_at( ar.low() ) == ZERO;
      Polynomial_1 g = gcd_utcf(p,ar.polynomial());
      return g.sign_at(ar.low()) != g.sign_at(ar.high());
    }
  };                
                   
  struct Is_square_free_1 
    : public std::unary_function< Polynomial_1, bool > {
    bool operator()( const Polynomial_1& p ) const {
      typename CGAL::Polynomial_traits_d< Polynomial_1 >::Is_square_free isf;
      return isf(p);
    }
  };
            
  struct Is_coprime_1
    : public std::binary_function< Polynomial_1, Polynomial_1, bool > {
    bool operator()( const Polynomial_1& p1, const Polynomial_1& p2 ) const {
      typename CGAL::Polynomial_traits_d< Polynomial_1 >::Total_degree total_degree;
                        
      // TODO: Is GCD already filtered? 
      return( total_degree( gcd_utcf( p1, p2 ) ) == 0 );                        
    } 
  }; 
            
  struct Make_square_free_1
    : public std::unary_function< Polynomial_1, Polynomial_1 > {
    Polynomial_1 operator()( const Polynomial_1& p ) const {
      return typename CGAL::Polynomial_traits_d< Polynomial_1 >::Make_square_free()( p );
    }
  };
            
  struct Make_coprime_1 {
    typedef bool         result_type;
    typedef Polynomial_1 first_argument_type;
    typedef Polynomial_1 second_argument_type;
    typedef Polynomial_1 third_argument_type;
    typedef Polynomial_1 fourth_argument_type;
    typedef Polynomial_1 fifth_argument_type;
                
    bool operator()( const Polynomial_1& p1,
        const Polynomial_1& p2,
        Polynomial_1& g, // ggT utcf 
        Polynomial_1& q1, // Rest utcf
        Polynomial_1& q2 ) const {
      g = typename CGAL::Polynomial_traits_d< Polynomial_1 >::Gcd_up_to_constant_factor()( p1, p2 );
      q1 = p1 / g;
      q2 = p2 / g;
      return CGAL::is_one(g);
    }                 
  };
            

  struct Square_free_factorize_1 {
    template< class OutputIterator>
    OutputIterator operator()( const Polynomial_1& p, OutputIterator it) const {
      typename PT_1::Square_free_factorize_up_to_constant_factor sqff;
      return sqff(p,it);
    } 
  };
                
  typedef typename Real_embeddable_traits<Algebraic_real_1>::Compare Compare_1;
  typedef typename Algebraic_real_traits::Bound_between Bound_between_1;
  typedef typename Algebraic_real_traits::Approximate_absolute_1 Approximate_absolute_1;
  typedef typename Algebraic_real_traits::Approximate_relative_1 Approximate_relative_1;


  
      
#define CGAL_ALGEBRAIC_KERNEL_1_PRED(Y,Z) Y Z() const { return Y(); }

  CGAL_ALGEBRAIC_KERNEL_1_PRED(Is_square_free_1,
      is_square_free_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Make_square_free_1,
      make_square_free_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Square_free_factorize_1,
      square_free_factorize_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Is_coprime_1,
      is_coprime_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Make_coprime_1,
      make_coprime_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Solve_1,
      solve_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Sign_at_1,
      sign_at_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Is_zero_at_1,
      is_zero_at_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Compare_1,
      compare_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Bound_between_1,
      bound_between_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Approximate_absolute_1,
      approximate_absolute_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Approximate_relative_1,
      approximate_relative_1_object);

  // Deprecated 
#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
  typedef Bound Boundary;  
  typedef typename Algebraic_real_traits::Refine Refine_1;
  typedef typename Algebraic_real_traits::Lower_bound Lower_bound_1;
  typedef typename Algebraic_real_traits::Upper_bound Upper_bound_1;
  typedef typename Algebraic_real_traits::Lower_bound Lower_boundary_1;
  typedef typename Algebraic_real_traits::Upper_bound Upper_boundary_1;
  typedef Bound_between_1 Boundary_between_1;

  CGAL_ALGEBRAIC_KERNEL_1_PRED(Refine_1,      refine_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Lower_bound_1, lower_bound_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Upper_bound_1, upper_bound_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Lower_boundary_1, lower_boundary_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Upper_boundary_1, upper_boundary_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Boundary_between_1, boundary_between_1_object);
#endif
      
#undef CGAL_ALGEBRAIC_KERNEL_1_PRED  
          
};
} // namespace CGALi


template< class Coefficient,
          class Bound = typename CGAL::Get_arithmetic_kernel< Coefficient >::Arithmetic_kernel::Rational,
          class RepClass = CGALi::Algebraic_real_rep< Coefficient, Bound >,
          class Isolator = CGALi::Descartes< typename CGAL::Polynomial_type_generator<Coefficient,1>::Type, Bound > >
class Algebraic_kernel_1    
  : public CGALi::Algebraic_kernel_1_base< 

    // Template argument #1 (AlgebraicReal1)        
        CGALi::Algebraic_real_pure< 
            Coefficient, 
            Bound,
            ::CGAL::Handle_policy_no_union,     
            RepClass >,
        
    // Template argument #2 (Isolator_)
        Isolator >
            
{};


CGAL_END_NAMESPACE



#endif // CGAL_ALGEBRAIC_KERNEL_1_H
