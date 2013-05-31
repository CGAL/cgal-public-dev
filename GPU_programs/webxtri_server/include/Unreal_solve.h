
/*! \file Unreal_solve.h
  \brief Defines class CGAL::Real_solve
  
  Isolate UNreal roots of polynomials with Fabrice Roullier's Rs.

  The polynomial has to be a univariat polynomial over any number
  type which is contained in the real numbers.

*/

#ifndef CGAL_UNREAL_SOLVE_H
#define CGAL_UNREAL_SOLVE_H

#include <CGAL/config.h>
#include <CGAL/Algebraic_kernel_rs_gmpz_d_1.h>


namespace CGAL {

namespace internal {

/*! \brief A model of concept UnrealRootIsolator.
 */
template <class Polynomial_, class Rational_> 
class Unreal_solve {

public:
    //! First template parameter
    typedef Polynomial_ Polynomial;

    //! Second template parameter
    typedef Rational_ Rational;

    //! Bound type of the isolating intervals
    typedef Rational_ Bound;
private:
    
    //! Coefficient type of polynomial
    typedef typename Polynomial::NT Coefficient;

    //! internal algebraic kernel
    typedef ::Algebraic_kernel_rs_1< CGAL::Gmpz > AK_1;
//                 CGAL::NTLgcd_1

    //! type of algebraic real
    typedef typename AK_1::Algebraic_real_1 Algebraic_real_1;
    
public:
    /*! \brief Constructor from univariate square free polynomial.
    
    The UnrealRootIsolator provides isolating intervals for the real
    roots of the polynomial.
    \pre the polynomial is square free      
    */
   Unreal_solve(const Polynomial& p = Polynomial(Coefficient(0))) :
      _m_polynomial(p)
      //, _m_interval_given(false)
    {

       std::cout << "\nstarting RealSolve\n";

      typename AK_1::Solve_1()(_m_polynomial, true, std::back_inserter(_m_real_roots));

      if (_m_real_roots.size() > 1) {
        for (std::vector< Algebraic_real_1 >::size_type i = 0;
             i < _m_real_roots.size() - 1;
             ++i) {
          // Isolate them against each other
          typename AK_1::Compare_1()(_m_real_roots[i],_m_real_roots[i+1]);
        }
      }
      std::cout << "\nRealSolve done\n";
    }
    
    // TODO interval constructor?

public: // functions
    
    /*! \brief returns the defining polynomial*/ 
    Polynomial polynomial() const { 
      return _m_polynomial;
    }
    
    //! returns the number of UNreal roots
    int number_of_real_roots() const { 
      return _m_real_roots.size();
    }

    /*! \brief returns true if the isolating interval is degenerated to a 
      single point.
      
      If is_exact_root(i) is true, 
      then left_bound(int i) equals  \f$root_i\f$. \n
      If is_exact_root(i) is true, 
      then right_bound(int i) equals  \f$root_i\f$. \n 
    */
    bool is_exact_root(int i) const {
      return _m_real_roots[i].is_point();
    }
      
public:   
  
    /*! \brief returns  \f${l_i}\f$ the left bound of the isolating interval 
      for root  \f$root_{i}\f$.
      
      In case is_exact_root(i) is true,  \f$l_i = root_{i}\f$,\n
      otherwise:  \f$l_i < root_{i}\f$. 
         
      If  \f$i-1>=0\f$, then  \f$l_i > root_{i-1}\f$. \n
      If  \f$i-1>=0\f$, then  \f$l_i >= r_{i-1}\f$, 
      the right bound of  \f$root_{i-1}\f$\n

      \pre 0 <= i < number_of_real_roots()
    */
    Bound left_bound(int i) const { 
      CGAL_assertion(i >= 0);
      CGAL_assertion(i < this->number_of_real_roots());
      return Bound(_m_real_roots[i].left());
    }
    
    /*! \brief returns  \f${r_i}\f$ the right bound of the isolating interval 
      for root  \f$root_{i}\f$.
      
      In case is_exact_root(i) is true,  \f$r_i = root_{i}\f$,\n
      otherwise:  \f$r_i > root_{i}\f$. 
         
      If  \f$i+1< n \f$, then  \f$r_i < root_{i+1}\f$,
      where \f$n\f$ is number of real roots.\n
      If  \f$i+1< n \f$, then  \f$r_i <= l_{i+1}\f$, 
      the left bound of  \f$root_{i+1}\f$\n
        
      \pre 0 <= i < number_of_real_roots()
    */
    Bound right_bound(int i) const { 
      CGAL_assertion(i >= 0);
      CGAL_assertion(i < this->number_of_real_roots());
      return Bound(_m_real_roots[i].right());
    }  
    
private:
    
    //! the squarefree polynomial
    Polynomial _m_polynomial;
    
    //! the solutions
    std::vector<  Algebraic_real_1 > _m_real_roots;
   
};

} // namespace internal

} //namespace CGAL

#endif // CGAL_UNREAL_SOLVE_H
