// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     :  Sarah Schaeffer
//
// ============================================================================


/*! \file include/CGAL/Algebraic_kernel_d/Descartes_star.h
 * \brief defines class \c Descartes_star
 *
 * Isolates real roots of a bitstream polynomial. It is garanteed that the computed isolating
 * intervals contain the real roots of the given polynomial.
 */

#ifndef CGAL_BITSTREAM_ISOLATORS_DESCARTES_STAR_H
#define CGAL_BITSTREAM_ISOLATORS_DESCARTES_STAR_H

#if CGAL_DESCARTES_VERBOSE
#define CGAL_real_roots_log(x) std::clog << x;
#else
#define CGAL_real_roots_log(x) static_cast< void >(0);
#endif

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h> 
#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Algebraic_kernel_d_1_generator.h>
#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>

#include <CGAL/Algebraic_kernel_d/Star_neighborhood_descartes.h>
#include <CGAL/Bitstream_coefficient.h>

#include <CGAL/Timer.h>

namespace CGAL {

namespace internal {

template < class Coefficient_, class Bound_ > 
class Descartes_star {

public:

  //! First template parameter
  typedef Coefficient_ Coefficient;

  //! Second template parameter
  typedef Bound_ Bound;
  
private:

  //! Vector of coefficients
  typedef std::vector< Coefficient > Coefficient_vector;
  
  //! size type
  typedef typename Coefficient_vector::size_type size_type;

  //! Fraction traits of bound
  typedef CGAL::Fraction_traits< Bound > FT;

  //! Numerator type of bound.
  typedef typename FT::Numerator_type Numerator;

  //! Denominator type of bound.
  typedef typename FT::Numerator_type Denominator;

  //! Integer or Numerator/Denominator type of bound.
  typedef typename Bitstream_coefficient<Coefficient>::Rational_approx 
  Rational_approx;
  
  //! Is zero functor
  typedef typename Bitstream_coefficient<Coefficient>::Is_zero Is_zero;
    
private:
    //! Univariate rational polynomial
  typedef typename CGAL::Polynomial_type_generator< Bound, 1 >::Type 
  Poly_bound_1;

  //! Integer type 
  // TODO to be taken from Bitstream_traits model
  typedef Numerator Integer;

  typedef typename CGAL::Polynomial_type_generator< Integer, 1 >::Type 
  Poly_int_1;

  //! traits class for approximation instance
  typedef ::CGAL::internal::Star_neighborhood_descartes< Poly_int_1, Bound >
  Star_neighborhood_descartes;
  
public:

  /*!
   * \brief Constructor from univariate square free polynomial.
   * 
   * The RealRootIsolator provides isolating intervals for the real 
   * roots of the polynomial.
   * \pre the polynomial is square free      
   */
  template < class InputIterator >
  Descartes_star(InputIterator begin, InputIterator end) {

    CGAL_real_roots_log("\n Start Descartes_star with new polynomial:"
                        << std::endl)

		CGAL::Timer timer_descartes_star;
		CGAL::Timer timer_star_neighborhood_descartes;
		timer_descartes_star.start();

    _m_coefficients.reserve(std::distance(begin, end) - 1);
    std::copy(begin, end, std::back_inserter(_m_coefficients));
    _m_degree = _m_coefficients.size() - 1;
    
    // Check for leading coefficient
    Is_zero is_zero;
//    if (is_zero(_m_coefficients[_m_degree])) {
     //    TODO: throw exception
//    }
    
    Poly_bound_1 p_approx, r, q;
    int number_of_rounds = 0;
    int approximation_exponent = 1;
    Bound precision;
    bool precision_sufficient = false;

		p_approx = approximated_polynomial(approximation_exponent);
		CGAL_real_roots_log("\t Polynomial: " << p_approx << std::endl);

		//TODO: root bound is rounded up to next integer to prevent
    //      huge coefficients.
    //      would be better to approximate it with log(n) bits
    _m_root_bound = round_up(root_bound());
    CGAL_real_roots_log("\t Root bound: " << _m_root_bound << std::endl);

    while (!precision_sufficient && number_of_rounds < 30) {

		  timer_star_neighborhood_descartes.reset();
      timer_star_neighborhood_descartes.start();

      CGAL_real_roots_log("\t Start " << number_of_rounds << ". round:"
                          << std::endl)

      //TODO: more criteria to prevent endless loop
	    //compute approximation and precision for next round
      approximation_exponent = CGAL::ipower(2, number_of_rounds);     
      precision = compute_precision(number_of_rounds);

      CGAL_real_roots_log("\t \t Precision: " << precision << std::endl)
	  
      p_approx = approximated_polynomial(approximation_exponent);
	  
	    //scale polynomial so that roots lie in (-1/4,1/4)
      r = ::CGAL::scale_up(p_approx, Bound(8) * Bound(_m_root_bound));
	  
      typedef Fraction_traits< Poly_bound_1 > FT_poly;
      typename FT_poly::Decompose decompose;
      typedef typename FT_poly::Numerator_type Numerator_poly;
      typedef typename Numerator_poly::NT Coeff;
      typename FT_poly::Numerator_type r_num;
      typename FT_poly::Denominator_type common_denominator;
      
	    //multiplying r with common divisor common_denominator results in r_num
      decompose(r, r_num, common_denominator);

	    //start isolating real roots
      _m_snd = Star_neighborhood_descartes(r_num, precision,
																					 common_denominator);
	  
      precision_sufficient = _m_snd.precision_sufficient();
      number_of_rounds++;

			timer_star_neighborhood_descartes.stop();
      CGAL_real_roots_log("\t Runtime Star_neighborhood_descartes: "
                      << timer_star_neighborhood_descartes.time() << std::endl)
    }

		//TODO: leave out?
    if (number_of_rounds >= 30) {
      std::cout << "Can't guarantee that intervals are correct!" << std::endl;
    }

		timer_descartes_star.stop();
    CGAL_real_roots_log("Runtime Descartes_star: "
                        << timer_descartes_star.time() << std::endl)
  }
 
  //!returns the number of real roots
  int number_of_real_roots() const {
    return _m_snd.number_of_real_roots();
  }
  
private:

  /*!
   * computes the left bound of the isolating interval for root i
   * \param i number of root
   * \param numerator variable in which numerator of left bound is written
   * \param denominator variable in which denominator of left bound is written 
   */
  void left_bound(int i, Numerator& numerator, Denominator& denominator) {
    _m_snd.left_bound(i, numerator, denominator);
    numerator = Numerator(8) * numerator * _m_root_bound;
  }
    
  /*!
   * computes the right bound of the isolating interval for root i
   * \param i number of root
   * \param numerator variable in which numerator of right bound is written
   * \param denominator variable in which denominator of right bound is written 
   */
  void right_bound(int i,Numerator& numerator, Denominator& denominator) {
    _m_snd.right_bound(i, numerator, denominator);
    numerator = Numerator(8) * numerator * _m_root_bound;
  }

public: 
    
  /*!
   * returns the left bound of the isolating interval for root i
   * \param i number of root
   * \return left bound of isolating interval
   */
  Bound left_bound(int i) { 
    Numerator numerator;
    Denominator denominator;
    left_bound(i,numerator,denominator);
    return Bound(numerator) / Bound(denominator);
  }

  /*!
   * returns the right bound of the isolating interval for root i
   * \param i number of root
   * \return right bound of isolating interval
   */
  Bound right_bound(int i) { 
    Numerator numerator;
    Denominator denominator;
    right_bound(i,numerator,denominator);
    return Bound(numerator) / Bound(denominator);
  }

  //TODO: method is needed because of requirements for "isolator", but makes
  //      no sense in this algorithm
  //!returns true, if root i lies on boundary of isolating interval
  bool is_exact_root(int i) const { 
    return false;
  }


private:

  /*!
   * computes and returns an upper root bound for p as follows:
   * bound = 2 * max[i=1..n-1]{(n-i)th_root(|a_i / a_n|)} 
   * \return upper root bound of _m_poly
   */
  Bound root_bound () {
    Bound result;
    Bound nth_coefficient;
    Bound ith_coefficient;
    Bound mth_root;
    int j = 0;
    std::pair<Bound,Bound> approximation_bounds;
    Rational_approx rational_approx;

    approximation_bounds = rational_approx(_m_coefficients[_m_degree], 3);
    nth_coefficient = CGAL::max(
      CGAL::abs(approximation_bounds.first),
      CGAL::abs(approximation_bounds.second));
    for (typename Coefficient_vector::const_iterator it = 
           _m_coefficients.begin(); it != _m_coefficients.end() - 1; it++) { 
      approximation_bounds = rational_approx(*it, 3);
      ith_coefficient = CGAL::min(
          CGAL::abs(approximation_bounds.first),
          CGAL::abs(approximation_bounds.second));
      mth_root = m_th_root(CGAL::abs(ith_coefficient / nth_coefficient),
                           _m_degree - j);
      result = CGAL::max(result, mth_root);
      j++;
    }
    return Bound(2) * result;
  }

  Numerator floor_log_of_root_bound () {
    CGAL::internal::floor_log2_abs(_m_root_bound);
  }

  /*!
   * returns an approximation of the input polynom with precision
   * 2^approximation_exponent
   * \param approximation_exponent exponent of approximation
   * \return approximated polynomial with coefficients of type Bound
   */
  Poly_bound_1 approximated_polynomial (int approximation_exponent) {

    typedef CGAL::Polynomial_traits_d<Poly_bound_1> PT_bound_1;

    typename PT_bound_1::Construct_polynomial construct_polynomial;
    std::cout << "approximated_polynomial" << std::endl;
    std::pair<Bound,Bound> approximation_bounds;
    std::list<Bound> coefficients_approx;
    Rational_approx rational_approx;
    for (typename Coefficient_vector::const_iterator it = 
           _m_coefficients.begin(); it != _m_coefficients.end(); it++) { 
      approximation_bounds = rational_approx(*it,approximation_exponent);
      coefficients_approx.push_back(
          (approximation_bounds.first + approximation_bounds.second) / Bound(2)
      );
    }

    Poly_bound_1 p_approx = construct_polynomial(coefficients_approx.begin(),
                                                 coefficients_approx.end());
    return p_approx;
  }

  /*!
   * computes precision according to current round
   * \param number_of_rounds number of current round
   * \return precision
   */
  Bound compute_precision (int number_of_rounds) {
    return (Bound(power<Bound>(Bound(8 * _m_root_bound), int(_m_degree))) * 
            power<Bound>(Bound(2), Numerator(3))) / 
            power<Bound>(Bound(2),
            power<Numerator>(Numerator(2), number_of_rounds));
  }

  /*!
   * computes an upper bound of the m-th root of _m_precision as follows:
   * decompose radicand into numerator and denominator
   * compute root1 of y1(x) = x^m - numerator
   * compute root2 of y2(x) = x^m - denominator
   * compute root1(rounded up) / root2(rounded down)
   * \return m-th root of radicand
   */
  Bound m_th_root (Bound radicand, int m) {

    typedef typename CGAL::Algebraic_kernel_d_1_generator< Numerator, Bound >::
    Algebraic_kernel_with_qir_and_descartes_1  Algebraic_kernel_1;
    typedef typename Algebraic_kernel_1::Algebraic_real_1 Algebraic_real_1;
    typedef typename Algebraic_kernel_1::Approximate_relative_1 
        Approximate_relative_1;
    typedef std::list< std::pair< Algebraic_real_1, unsigned int > > Roots;

    std::list<Algebraic_real_1> coeff_real;
    Roots roots_num;
    Roots roots_den;
    typename Algebraic_kernel_1::Solve_1 solve;
    
    //TODO: use polynomial instead of polynomial_1
    typedef typename Algebraic_kernel_1::Polynomial_1 Polynomial_1;
    typedef CGAL::Polynomial_traits_d<Polynomial_1> PT_1;

    Numerator radicand_num;
    Denominator radicand_den;

    typename FT::Decompose decompose;
    decompose(radicand,radicand_num,radicand_den);

    Polynomial_1 x = typename PT_1::Shift()(Polynomial_1(1),1);
    Polynomial_1 m_th_root_polynomial_num = x;
    for (int j = 0; j < m-1; j++) {
      m_th_root_polynomial_num *= x;
    }
    m_th_root_polynomial_num = m_th_root_polynomial_num - radicand_num;

    Polynomial_1 m_th_root_polynomial_den = x;
    for (int j = 0; j < m-1; j++) {
      m_th_root_polynomial_den *= x;
    }
    m_th_root_polynomial_den = m_th_root_polynomial_den - radicand_den;
    
    solve(m_th_root_polynomial_num, std::back_inserter(roots_num));
    Algebraic_real_1 m_th_root_num = roots_num.back().first;  

    solve(m_th_root_polynomial_den, std::back_inserter(roots_den));
    Algebraic_real_1 m_th_root_den = roots_den.back().first;

    Approximate_relative_1 approximation;  
    std::pair<Bound,Bound> isolating_interval_num =
              approximation(m_th_root_num,3);
    std::pair<Bound,Bound> isolating_interval_den =
              approximation(m_th_root_den,3);

    return isolating_interval_num.second / isolating_interval_den.first;
  }

  /*!
   * rounds a value up to next integer value
   * \param value value which will be rounded
   * \return next upper integer of value
   */
  Numerator round_up (Bound value) {
    typename FT::Decompose decompose;
    Numerator quotient, rest;
    Numerator value_num;
    Denominator value_den;
    decompose(value,value_num,value_den);
    CGAL::div_mod(value_num,value_den,quotient,rest);
    if (CGAL::compare(rest,Numerator(0)) != CGAL::EQUAL) {
      quotient++;
    }
    return quotient;
  }

  
  template <class T>
  /*!
   *  computes base^e
   * \param base base
   * \param e exponent
   * \return base^e
   */
  T power (const T base, Numerator e) {
    T result = base;
    if (CGAL::compare(e, Numerator(0)) == CGAL::EQUAL) {
      return T(1);
    }
    for (Numerator i = Numerator(0); i < e - Numerator(1); i++) {
      result *= base;
    }
    return result;
  }

private:
  
  //! stores coefficients
  Coefficient_vector _m_coefficients;
  
  //! degree of 'polynomial'
  size_type _m_degree;

  //! approx instance
  Star_neighborhood_descartes _m_snd;

  //! rational root bound
  Numerator _m_root_bound;
  
};


} // namespace internal

} //namespace CGAL

#endif // CGAL_BITSTREAM_ISOLATORS_DESCARTES_STAR_H

