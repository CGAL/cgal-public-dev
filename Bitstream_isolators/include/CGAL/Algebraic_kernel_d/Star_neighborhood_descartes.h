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


// TODO: Instead of power try using ipower somehow


/*!
 * \file include/CGAL/Algebraic_kernel_d/Star_neighborhood_descartes.h
 * \brief defines class \c Star_neighborhood_descartes
 *
 * Isolates real roots of an approximated polynomial. Whenever the isolating
 * intervals get too small with respect to a given precision, it is reported.
 */

#ifndef CGAL_BITSTREAM_ISOLATORS_STAR_NEIGHBORHOOD_DESCARTES_H
#define CGAL_BITSTREAM_ISOLATORS_STAR_NEIGHBORHOOD_DESCARTES_H

#if CGAL_DESCARTES_VERBOSE
#define CGAL_real_roots_log(x) std::clog << x;
#else
#define CGAL_real_roots_log(x) static_cast< void >(0);
#endif

#if CGAL_DESCARTES_EXTENDED_VERBOSE
#define CGAL_real_roots_extended_log(x) std::clog << x;
#else
#define CGAL_real_roots_extended_log(x) static_cast< void >(0);
#endif


#include <CGAL/config.h>

#include <CGAL/Polynomial.h> 
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Coercion_traits.h>

#include <CGAL/Algebraic_kernel_d/univariate_polynomial_utils.h>
#include <CGAL/Algebraic_kernel_d/construct_binary.h>

#include <CGAL/Timer.h>

namespace CGAL {

namespace internal {


template < class Polynomial_, class Bound_ > 
class Star_neighborhood_descartes {

public:

  //! First template parameter
  typedef Polynomial_ Polynomial;
  
  //! Second template parameter
  typedef Bound_ Bound;

private:

  // TODO replace Fractions with Mantissa/Exponent types
  //! Fraction traits of Bound
  typedef Fraction_traits< Bound > FT_bound;

  //! Numerator of Bound
  typedef typename FT_bound::Numerator_type Numerator;
  
  //! Denominator of Bound
  typedef typename FT_bound::Denominator_type Denominator;
  
  //! Coefficient of polynomial
  typedef typename Polynomial::NT Coefficient;

  //! Polynomial traits
  typedef CGAL::Polynomial_traits_d<Polynomial> PT;
  
  //!\name Constructors
  //!@{
  
public:

  /*!
   * \brief Constructor from univariate square free polynomial.
    
   *  The RealRootIsolator provides isolating intervals for the real 
   *  roots of the polynomial.
   * \param p the given polynomial
   * \param precision precision 
   * \param common_denominator value with which the polynomial was multiplied previously
   * \pre the polynomial is square free      
   */
  Star_neighborhood_descartes(const Polynomial& p = Polynomial(Coefficient(0)),
                              Bound precision = Bound(2),
                              Denominator common_denominator = Denominator(1)):
    _m_poly(p),
    _m_precision(precision), 
    _m_common_denom(common_denominator),
    _m_number_of_real_roots(0),
    _m_precision_sufficient(true) {
        
    _m_degree = CGAL::degree(p);
    _m_numerator = new Numerator[_m_degree]; 
    _m_denominator_exponent = new Denominator[_m_degree];
    _m_radius_numerator = new Numerator[_m_degree];
    _m_radius_denominator = new Denominator[_m_degree];

    if (_m_degree == 0) { 
      if (p.is_zero()) {
        _m_number_of_real_roots = -1;
      }
      return;
    }
   
    isolate();
  }
  
  /*!
  * copy constructor
  * \param cnd copy instance
  */
  Star_neighborhood_descartes(const Star_neighborhood_descartes& t):
    _m_poly(t._m_poly), 
    _m_number_of_real_roots(t._m_number_of_real_roots),
    _m_precision(t._m_precision),
    _m_common_denom(t._m_common_denom),
    _m_precision_sufficient(t._m_precision_sufficient),
    _m_degree(t._m_degree) {
        
    _m_numerator = new Numerator[_m_degree]; 
    _m_denominator_exponent = new Denominator[_m_degree];
    _m_radius_numerator = new Numerator[_m_degree];
    _m_radius_denominator = new Denominator[_m_degree];

    for (int i = 0; i < _m_degree; i++) {
      _m_numerator[i] = t._m_numerator[i];
      _m_denominator_exponent[i] = t._m_denominator_exponent[i];
    }
  }

  //!@}

  //!\name Desctructors
  //!@{
  
  //! destructor
  ~Star_neighborhood_descartes() {
    delete[] _m_numerator;
    delete[] _m_denominator_exponent;
    delete[] _m_radius_numerator;
    delete[] _m_radius_denominator;
  }
  
  //!@}

  //! = operator
  void operator=(const Star_neighborhood_descartes& snd) {
    delete[] _m_numerator;
    delete[] _m_denominator_exponent;
    delete[] _m_radius_numerator;
    delete[] _m_radius_denominator;

    _m_degree = snd._m_degree;
    _m_poly = snd._m_poly;
    _m_common_denom = snd._m_common_denom;
    _m_number_of_real_roots = snd._m_number_of_real_roots;   

    _m_numerator = new Numerator[_m_degree]; 
    _m_denominator_exponent = new Denominator[_m_degree];
    _m_radius_numerator = new Numerator[_m_degree];
    _m_radius_denominator = new Denominator[_m_degree];
    
    for (int i = 0; i < _m_degree; i++) {
      _m_numerator[i] = snd._m_numerator[i];
      _m_denominator_exponent[i] = snd._m_denominator_exponent[i];
      _m_radius_numerator[i] = snd._m_radius_numerator[i];
      _m_radius_denominator[i] = snd._m_radius_denominator[i];
    }

    _m_precision = snd._m_precision;
    _m_precision_sufficient = snd._m_precision_sufficient;
  }   

public:
    
 //! returns whether used precision sufficed while executing/constructing
  bool precision_sufficient() const { 
    return _m_precision_sufficient; 
  }
  
  //! returns the number of real roots
  int number_of_real_roots() const { 
    return _m_number_of_real_roots; 
  }
  
  /*!
   * computes the left bound of the isolating interval for root i
   * \param i number of root
   * \param numerator variable in which numerator of left bound is written
   * \param denominator variable in which denominator of left bound is written 
   */
  void left_bound(int i, Numerator& numerator, Denominator& denominator) const {
    
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < _m_number_of_real_roots);
    
    construct_binary(_m_denominator_exponent[i], denominator);

    numerator = Numerator(2) * _m_numerator[i] - denominator;
    denominator = Denominator(4) * denominator;
    numerator = numerator * _m_radius_denominator[i] -
                _m_radius_numerator[i] * denominator;
    denominator = denominator * _m_radius_denominator[i];
  }
   
  /*!
   * computes the right bound of the isolating interval for root i
   * \param i number of root
   * \param numerator variable in which numerator of right bound is written
   * \param denominator variable in which denominator of right bound is written 
   */
  void right_bound(int i, Numerator& numerator, Denominator& denominator)const {
    
    CGAL_assertion(i >= 0);
    CGAL_assertion(i < _m_number_of_real_roots);
    
    construct_binary(_m_denominator_exponent[i],denominator);

    numerator = Numerator(2) * (_m_numerator[i] + Numerator(1)) - denominator;
    denominator = Denominator(4) * denominator;
    numerator = numerator * _m_radius_denominator[i] +
                _m_radius_numerator[i] * denominator;
    denominator = denominator * _m_radius_denominator[i];
  }

  /*! 
   * \brief returns  \f${l_i}\f$ the left bound of the isolating interval 
   * for root  \f$root_{i}\f$.
    
   * If  \f$i-1>=0\f$, then  \f$l_i > root_{i-1}\f$. \n
   * If  \f$i-1>=0\f$, then  \f$l_i >= r_{i-1}\f$, 
   * the right bound of  \f$root_{i-1}\f$\n

   * \pre 0 <= i < number_of_real_roots()
   */
  Bound left_bound(int i) const { 
    Numerator numerator;
    Denominator denominator;
    left_bound(i,numerator,denominator);
    return Bound(numerator) / Bound(denominator);
  }
    
  /*!
   * \brief returns  \f${r_i}\f$ the right bound of the isolating interval 
   * for root  \f$root_{i}\f$.

   * If  \f$i+1< n \f$, then  \f$r_i < root_{i+1}\f$,
   * where \f$n\f$ is number of real roots.\n
   * If  \f$i+1< n \f$, then  \f$r_i <= l_{i+1}\f$, 
   * the left bound of  \f$root_{i+1}\f$\n
     
   * \pre 0 <= i < number_of_real_roots()
   */
  Bound right_bound(int i) const { 
    Numerator numerator;
    Denominator denominator;
    right_bound(i,numerator,denominator);
    return Bound(numerator) / Bound(denominator);
  } 
    
private:

  /*!scales and translates P onto interval (0,1)
   * starts descartes algorithm with tests
   * starts test, which checks whether all real roots are covered
   */
  void isolate() {

    Polynomial r = ::CGAL::scale_down(_m_poly, Coefficient(4));
    Polynomial q = ::CGAL::translate(r, Coefficient(-1));
    Polynomial s = ::CGAL::scale_up(q, Coefficient(2));
		
	  _m_common_denom *= power<Denominator,int>(Denominator(4), _m_degree);

		CGAL::Timer timer_methods;
		timer_methods.start();

    zero_one_descartes(s, Numerator(0), Denominator(0),
                       Bound(0), Bound(1) / Bound(4), Denominator(1));

		timer_methods.stop();
		CGAL_real_roots_log("\t \t Runtime iso-apx: " << timer_methods.time()
                        << std::endl)
		timer_methods.reset();
    
    CGAL_real_roots_extended_log("\t \t \t Runtime apx transform polys: " <<
                        timer_polys.time() << std::endl)
    CGAL_real_roots_extended_log("\t \t \t Runtime apx detect 'no root': " <<
                        timer_no_root.time() << std::endl)
    CGAL_real_roots_extended_log("\t \t \t Runtime apx detect 'one root': " <<
                        timer_one_root.time() << std::endl)
    CGAL_real_roots_extended_log("\t \t \t Runtime apx (mu,k)-test: " <<
                        timer_mu_k_test.time() << std::endl)
    CGAL_real_roots_log("\t \t Runtime certify: " << timer_certify.time()
                        << std::endl)

    if (_m_precision_sufficient) {
			timer_methods.start();
			_m_common_denom /= power<Denominator,int>(Denominator(4), _m_degree);
      all_real_roots_covered();
			timer_methods.stop();
      CGAL_real_roots_log("\t \t Runtime iso-ex: " << timer_methods.time()
                          << std::endl)
    }

    CGAL_real_roots_extended_log("\t \t Runtime circle_test: " <<
                        timer_circle_test.time() << std::endl)

		CGAL_real_roots_extended_log("\t \t \t Runtime circle_test 1: " <<
                        timer_circle_test1.time() << std::endl)
    CGAL_real_roots_extended_log("\t \t \t Runtime circle_test 2: " <<
                        timer_circle_test2.time() << std::endl)
    CGAL_real_roots_extended_log("\t \t \t Runtime circle_test 3: " <<
                        timer_circle_test3.time() << std::endl)
    CGAL_real_roots_extended_log("\t \t \t Runtime circle_test 4: " <<
                        timer_circle_test4.time() << std::endl)
    CGAL_real_roots_extended_log("\t \t \t Runtime circle_test 5: " <<
                        timer_circle_test5.time() << std::endl)
    CGAL_real_roots_extended_log("\t \t \t Runtime circle_test 6: " <<
                        timer_circle_test6.time() << std::endl)

    CGAL_real_roots_extended_log("\t \t Runtime (mu,k)-test: " <<
                        timer_mu_k_test.time() << std::endl)

    CGAL_real_roots_extended_log("\t \t \t Runtime (mu,k)-test 1: " <<
                        timer_mu_k_test1.time() << std::endl)
    CGAL_real_roots_extended_log("\t \t \t Runtime (mu,k)-test 2: " <<
                        timer_mu_k_test2.time() << std::endl)
    CGAL_real_roots_extended_log("\t \t \t Runtime (mu,k)-test 3: " <<
                        timer_mu_k_test3.time() << std::endl)
    CGAL_real_roots_extended_log("\t \t \t Runtime (mu,k)-test 4: " <<
                        timer_mu_k_test4.time() << std::endl)
    CGAL_real_roots_extended_log("\t \t \t Runtime (mu,k)-test 5: " <<
                        timer_mu_k_test5.time() << std::endl)
  }
  
  //! returns the polynomial $(1 + x)^n P(1/(1 + x))$.
  Polynomial variation_transformation(const Polynomial& p) { 
    Polynomial r = reversal(p);
    return translate_by_one(r); 
  }
  
  /*! 
   * circle test can determine whether in a given interval exist no roots or exactly one,
   * moreover used in mu_k_test
   * \param k scale factor
   * \param p polynomial (already scaled and shifted)
   * \param midpoint midpoint of current interval
   * \param radius radius of current interval
   * \return true if circle test returns true
   */
  bool circle_test (Bound k, const Polynomial& p, const Coefficient& radius) {

		timer_circle_test.start();

    Coefficient result_sum = Coefficient(0);
    Coefficient new_radius = Coefficient(1);
		
    for (int i = 1; i <= CGAL::degree(p); i++) {

			timer_circle_test1.start();
			new_radius *= radius;
			timer_circle_test1.stop();

			timer_circle_test2.start();
			const Coefficient& coeff = get_coefficient(p,i);
			timer_circle_test2.stop();

			timer_circle_test3.start();
			Coefficient coeff_abs = CGAL::abs(coeff);
			timer_circle_test3.stop();

			timer_circle_test4.start();
      Coefficient mul =  coeff_abs * new_radius;
			timer_circle_test4.stop();

			timer_circle_test5.start();
      result_sum = result_sum + mul;
			timer_circle_test5.stop();
			
    }

		bool result = CGAL::compare(Bound(CGAL::abs(get_coefficient(p,0))),
                         k * Bound(result_sum)) == CGAL::LARGER;

		timer_circle_test.stop();

    return result;
  }
  
  /*!
   * mu_k_test can determine whether the current interval is too small according to the current precision
   * and the given polynomial  
   * \param midpoint midpoint of current interval
   * \param radius radius of current interval
   * \param p polynomial (already scaled and shifted)
   * \param k scale factor
   * \param denominator_correction value with wihich the polynomial was multiplied before
   * \return true if (mu,k)-test returns true
   */
  bool mu_k_test (const Coefficient& radius, const Bound& midpoint,
                  const Polynomial& p, const Bound& k,
                  const Denominator& denominator_correction) {

    typedef typename PT::Coefficient_const_iterator_range 
       Coefficient_const_iterator_range;
    
    typename PT::Construct_coefficient_const_iterator_range itrange;
    typename PT::Construct_polynomial construct_polynomial;
  
    typename FT_bound::Decompose decompose;

    bool test_1_succeeds, test_2_succeeds;
    Polynomial r = CGAL::scale_down(p,Numerator(2));
    
		timer_mu_k_test1.start();
    Polynomial poly_test_1, poly_test_2;
    poly_test_1 = r;
    poly_test_2 = CGAL::differentiate(r);
    Coefficient new_value;
    Bound precision_test_1 = _m_precision *
                      Bound(_m_common_denom * denominator_correction);
		Bound precision_test_2 = _m_precision *
                      Bound(_m_common_denom * denominator_correction * radius);
    timer_mu_k_test1.stop();
    
    timer_mu_k_test2.start();
		Numerator precision_num_test_1;
    Denominator precision_den_test_1;
    decompose(precision_test_1,precision_num_test_1,precision_den_test_1);

		Numerator precision_num_test_2;
    Denominator precision_den_test_2;
    decompose(precision_test_2,precision_num_test_2,precision_den_test_2);

    std::vector< Coefficient > coeffs_test_1;
    std::vector< Coefficient > coeffs_test_2;
    coeffs_test_1.reserve(_m_degree + 1);
    coeffs_test_2.reserve(_m_degree);
    timer_mu_k_test2.stop();

    timer_mu_k_test3.start();
    Coefficient_const_iterator_range range_test_1 = itrange(poly_test_1);
    Coefficient_const_iterator_range range_test_2 = itrange(poly_test_2);
    std::copy(range_test_1.first, range_test_1.second,
              std::back_inserter(coeffs_test_1));
    std::copy(range_test_2.first, range_test_2.second,
              std::back_inserter(coeffs_test_2));

    coeffs_test_1[0] = CGAL::abs(coeffs_test_1[0]) * precision_den_test_1 +
                       Numerator(2) * precision_num_test_1;
    coeffs_test_1[_m_degree] *= precision_den_test_1;
  
    for (int i = 1; i < _m_degree; i++) {
      new_value = Coefficient(CGAL::abs(coeffs_test_1[i]) * precision_den_test_1
									- Numerator(2) * precision_num_test_1);
      if (CGAL::compare(new_value,Coefficient(0)) == CGAL::SMALLER) {
        coeffs_test_1[i] = Coefficient(0);
      } else {
        coeffs_test_1[i] = new_value;        
      }
    }
    timer_mu_k_test3.stop();

    timer_mu_k_test4.start();
    coeffs_test_2[0] = CGAL::abs(coeffs_test_2[0]) * Coefficient(2) *
                       precision_den_test_2 + Numerator(2) * _m_degree *
											 precision_num_test_2;
    coeffs_test_2[_m_degree-1] *= precision_den_test_2;
    for (int i = 1; i < _m_degree-1; i++) {
      new_value = CGAL::abs(coeffs_test_2[i]) * Denominator(2) *
									precision_den_test_2 - Numerator(2) * _m_degree *
									precision_num_test_2;
      if (CGAL::compare(new_value,Coefficient(0)) == CGAL::SMALLER) {
        coeffs_test_2[i] = Coefficient(0);
      } else {
        coeffs_test_2[i] = new_value;
      }
    }
    timer_mu_k_test4.stop();

    timer_mu_k_test5.start();
    poly_test_1 = construct_polynomial(coeffs_test_1.begin(),
                                       coeffs_test_1.end());
    poly_test_2 = construct_polynomial(coeffs_test_2.begin(),
                                       coeffs_test_2.end());
    timer_mu_k_test5.stop();

    test_1_succeeds = circle_test(k, poly_test_1, Coefficient(2) * radius);
    test_2_succeeds = circle_test(k, poly_test_2, Coefficient(2) * radius);
    
    return test_1_succeeds || test_2_succeeds;
  }
  
  /*!
   * tests whether the current interval is adjacent to an already saved interval
   * \param i numerator of interval: [i/2^exp, (i+1)/2^exp]
   * \param exp denominator of interval: [i/2^exp, (i+1)/2^exp]
   * \return true if the current interval is adjacent to a saved interval
   */
  bool adjacency_test(const Numerator& i, const Denominator& exp) {
    Bound left_bound, right_bound;
    Bound left_bound_root, right_bound_root;
    left_bound = Bound(i) / power<Bound, Denominator>(Bound(2), exp);
    right_bound = Bound(i + 1) / power<Bound, Denominator>(Bound(2), exp);

    for (int j = 0; j < _m_number_of_real_roots; j++){
      left_bound_root = Bound(_m_numerator[j]) /
        power<Bound, Denominator>(Bound(2), _m_denominator_exponent[j]);
      right_bound_root = Bound(_m_numerator[j] + 1) /
        power<Bound, Denominator>(Bound(2), _m_denominator_exponent[j]);
      if (CGAL::compare(left_bound, left_bound_root) == CGAL::EQUAL) {
        return true;
      }
      if (CGAL::compare(left_bound, right_bound_root) == CGAL::EQUAL) {
        return true;
      }
      if (CGAL::compare(right_bound, left_bound_root) == CGAL::EQUAL) {
        return true;
      }
      if (CGAL::compare(right_bound, right_bound_root) == CGAL::EQUAL) {
        return true;
      }
    }
    return false;
  }
  
  /*! Algoritm to determine isolating intervals for the roots 
   * lying in the interval (0,1).
   * The parameters $(i,D)$ describe the interval $(i/2^D, (i+1)/2^D)$.
   * Here $0\leq i < 2^D$.
   */
  void zero_one_descartes(const Polynomial& p, 
                          Numerator i, Denominator exp, 
                          const Bound& midpoint, const Bound& radius,
                          const Denominator& denominator_correction) { 

    typename FT_bound::Decompose decompose;
    typename PT::Evaluate evaluate;
    
    if (!_m_precision_sufficient) {
      return;
    }
  
    Polynomial r = variation_transformation(p);

    int descarte = sign_variations(r);
    Coefficient new_radius = Coefficient(7 * _m_degree * _m_degree * (5 * _m_degree + 1));
    
    Numerator radius_num;
    Denominator radius_den;
    decompose(radius, radius_num, radius_den);
    
		timer_polys.start();
    Polynomial q = scale_down(p, Coefficient(2));
    Polynomial s = translate_by_one(q);
    
		Polynomial poly_circle_test = s;
		timer_polys.stop();
    
		//no root
		timer_no_root.start();
    if (descarte == 0 ||
				circle_test(Bound(1), poly_circle_test, Coefficient(1))) {
      if (CGAL::compare(get_coefficient(p,0) * evaluate(p,1), Coefficient(0))
		      != CGAL::EQUAL ) {
				timer_no_root.stop();
	      return;
      }
    }
		timer_no_root.stop();

    // exactly one root
		timer_one_root.start();
    if ( circle_test(Bound(3) / Bound(2), CGAL::differentiate(poly_circle_test),
										 Coefficient(5*_m_degree + 1)) ) { 
      if (CGAL::compare(get_coefficient(p,0) * evaluate(p,1), Coefficient(0))
                != CGAL::LARGER )  {
        if (_m_number_of_real_roots == 0 || !adjacency_test(i,exp)) {
          if (isolating_intervals(poly_circle_test, radius, midpoint,
                                  denominator_correction)) {
            _m_numerator[_m_number_of_real_roots] = i;
            _m_denominator_exponent[_m_number_of_real_roots] = exp;
            _m_number_of_real_roots++;
						timer_one_root.stop();
            return;
          }
        }
      }
			timer_one_root.stop();
      return; 
    }
		timer_one_root.stop();
  
		timer_mu_k_test.start();      
    if (mu_k_test(new_radius, midpoint, poly_circle_test, Bound(3) / Bound(2),
                  denominator_correction)) { 
      _m_precision_sufficient = false;
			timer_mu_k_test.stop();
      return;
    }
		timer_mu_k_test.stop();

    // we are left with more than one root
    // Refine the interval.
    i = Numerator(2) * i; 
    exp = exp + Denominator(1);

    Denominator new_denom_correction = denominator_correction * 
                           power<Denominator,int>(Denominator(2), _m_degree);
 
    // Consider the first half of the interval.
    Bound radius_next = radius / Bound(2);
    Bound midpoint_next = midpoint - radius_next;
    zero_one_descartes(q, i, exp, midpoint_next, radius_next,
                       new_denom_correction);
     
    // Consider the second half of the interval.
    midpoint_next = midpoint + radius_next; 
    zero_one_descartes(s, i + Numerator(1), exp, midpoint_next, radius_next,
                       new_denom_correction); 
  }

  /*!
   * computes isolating intervals by adding a computed value on both sides
   * \param p polynomial (scaled on current interval)
   * \param radius radius of current interval
   * \param midpoint midpoint of current interval
   * \param denominator_correction value with which the polynomial was multiplied previously
   * \return true if precision is sufficient
   */
	//TODO: delete midpoint from parameters?
  bool isolating_intervals(const Polynomial& p,
                           const Bound& radius, const Bound& midpoint,
                           const Denominator& denominator_correction) {

		timer_certify.start();
    
    typename FT_bound::Decompose decompose;
    
    Coefficient tau;
    Bound r;
    Numerator r_numerator; 
    Denominator r_denominator;
    Bound test_value;
    bool found_correct_bound = false;

    if (_m_degree == 1) {
      decompose(radius, r_numerator, r_denominator);
      _m_radius_numerator[_m_number_of_real_roots] = r_numerator;
      _m_radius_denominator[_m_number_of_real_roots] = r_denominator;
			timer_certify.stop();
      return true;
    }

    tau = Coefficient(5 * _m_degree + 1);
    Polynomial p_diff = CGAL::differentiate(p);
    
		//TODO: warum nicht tau*r?
    while (circle_test(Bound(3) / Bound(2), p_diff, tau)) {
      tau *= Coefficient(2);
    }

    r = Bound(2) * Bound(tau) * radius / Bound(5 * _m_degree + 1);
    test_value = r * CGAL::abs(get_coefficient(p_diff, 0)) /
        (_m_common_denom * denominator_correction * Denominator(8 * _m_degree));
    
    if (CGAL::compare(_m_precision, test_value) != CGAL::SMALLER) {
      _m_precision_sufficient = false;
      timer_certify.stop();
			return false;
    } else {
      decompose(r, r_numerator, r_denominator);
      _m_radius_numerator[_m_number_of_real_roots] = r_numerator;
      _m_radius_denominator[_m_number_of_real_roots] = r_denominator;
			timer_certify.stop();      
			return true;
    }
  }
  
  //! tests whether all real roots are covered
  void all_real_roots_covered() { 

    typename FT_bound::Decompose decompose;

    Bound left_bound_interval, right_bound_interval;
    Bound radius, midpoint;
    Numerator left_bound_num, right_bound_num;
    Denominator left_bound_den, right_bound_den;
    Bound common_denom_temp = _m_common_denom;
    for (int j = 0; j <= _m_number_of_real_roots; j++) {

      if (j == 0) {
        left_bound_interval = Bound(-1) / Bound(4);
        right_bound_interval = left_bound(j);
      } else if (j == _m_number_of_real_roots) {
        left_bound_interval = right_bound(j - 1);
        right_bound_interval = Bound(1) / Bound(4);
      } else {
        left_bound_interval = right_bound(j - 1);
        right_bound_interval = left_bound(j);
      }

      radius = (right_bound_interval - left_bound_interval) / Bound(2);
      midpoint = (right_bound_interval + left_bound_interval) / Bound(2);

      decompose(left_bound_interval, left_bound_num, left_bound_den);
      decompose(right_bound_interval, right_bound_num, right_bound_den);
       
      Numerator scale = left_bound_den * right_bound_num -
                        left_bound_num * right_bound_den;
						
	    //TODO: polynomial is scaled, save in _m_common_denom (???)
      
      Polynomial p = scale_down(_m_poly,
                                Coefficient(left_bound_den * right_bound_den));
      Polynomial q = translate(p,Coefficient(left_bound_num * right_bound_den));
      Polynomial r = scale_up(q, Coefficient(scale));
      
			//TODO: common denom *= (left_bound_den * right_bound_den)^n
//      _m_common_denom = common_denom_temp;
//      _m_common_denom = _m_common_denom /
//                     Bound(scale * left_bound_den * right_bound_den);
      _m_common_denom = common_denom_temp *
						power<Denominator,int>(left_bound_den * right_bound_den, _m_degree);

      all_real_roots_covered(r, midpoint, radius, Denominator(1));
    }        
  }
  
  /*!
   * tests whether all real roots are covered
   * recursive method to check all subintervals of given interval
   * \param p polynomial (already scaled onto current interval)
   * \param midpoint midpoint of current interval
   * \param radius radius of current interval
   * \param denominator_correction value with which the polynomial
	 *                               was multiplied previously
   */
  void all_real_roots_covered (const Polynomial& p, const Bound& midpoint,
                               const Bound& radius, 
                               const Denominator& denominator_correction) {

    typename FT_bound::Decompose decompose;
    typename PT::Evaluate evaluate;

    if (!_m_precision_sufficient) {
      return;
    }
   
    //TODO: multiply r with 1 / (b-a)
		Polynomial p_diff = CGAL::differentiate(p);
		Polynomial r = variation_transformation(p_diff);
    int descarte = sign_variations(r);

    Numerator radius_num;
    Denominator radius_den;
    decompose(radius, radius_num, radius_den);
    
    Polynomial q = scale_down(p, Coefficient(2));
    Polynomial s = translate_by_one(q);

		Polynomial poly_circle_test = s;

    if (descarte == 0 ||
	      circle_test (Bound(1), CGAL::differentiate(poly_circle_test),
										 Coefficient(1))) {
      if (CGAL::compare(CGAL::min(CGAL::abs(get_coefficient(p,0)),
           CGAL::abs(evaluate(p,1))), _m_degree * _m_precision)
                      != CGAL::LARGER) {
        _m_precision_sufficient = false;
      }
      return;
    }

    if (circle_test(Bound(3) / Bound(2), poly_circle_test, Coefficient(1))) {
      if (CGAL::compare(Bound(1) / Bound(3) *
                        CGAL::abs(get_coefficient(poly_circle_test, 0)),
                        _m_degree * _m_precision)
                        != CGAL::LARGER) {
        _m_precision_sufficient = false;
      }
      return;
    }

    if (mu_k_test (Coefficient(7 * _m_degree * _m_degree), midpoint,
                   poly_circle_test, Bound(3) / Bound(2),
                   denominator_correction)) {
      _m_precision_sufficient = false;
      return;
    }

    Bound radius_next = radius / Bound(2);
    Bound midpoint_next = midpoint - radius_next;

    Denominator new_denom_correction = denominator_correction * 
                           power<Denominator,int>(Denominator(2), _m_degree);

    // Consider the first half of the interval.
    all_real_roots_covered(q, midpoint_next, radius_next, new_denom_correction);
     
    // Consider the second half of the interval.
    midpoint_next = midpoint + radius_next; 
    all_real_roots_covered(s, midpoint_next, radius_next, new_denom_correction); 
    
  }

  template < class NT >
  /*! 
   * computes signum(value)
   * \param v value
   * \return signum(value)
   */
  NT sign_value(const NT& v) {
    CGAL::Sign s = CGAL::sign(v);
    switch (s) {
    case NEGATIVE: return NT(-1); break;
    case ZERO: return NT(0); break;
    case POSITIVE: return NT(1); break;
    }
  }
  
  template <class NTBase, class NTExp>
  /*!
   *  computes base^e
   * \param base base
   * \param e exponent
   * \return base^e
   */
  NTBase power (const NTBase base, NTExp e) {
    NTBase result = base;
    if (e == NTExp(0)) {
      return NTBase(1);
    }
    for (NTExp i = NTExp(0); i < e - NTExp(1); i++) {
      result *= base;
    }
    return result;
  }

private:
  
  //! polynomial
  Polynomial _m_poly;

  //! common denominator
  Bound _m_common_denom;

  //! number of real roots
  int _m_number_of_real_roots;   

  //! numerator
  Numerator* _m_numerator;   

  //! denominator
  Denominator* _m_denominator_exponent; 

  //! precision
  Bound _m_precision;

  //! is the precision sufficient?
  bool _m_precision_sufficient;

  //! number of round
  int _m_number_of_rounds;

  //! numerator of radius
  Numerator* _m_radius_numerator;

  //! denominator of radius
  Denominator* _m_radius_denominator;

  //! degree
  int _m_degree;

  //! various timer for benchmarking
	CGAL::Timer timer_certify;
	CGAL::Timer timer_polys;
	CGAL::Timer timer_no_root;
	CGAL::Timer timer_one_root;
	CGAL::Timer timer_mu_k_test;
  CGAL::Timer timer_mu_k_test1;
  CGAL::Timer timer_mu_k_test2;
  CGAL::Timer timer_mu_k_test3;
  CGAL::Timer timer_mu_k_test4;
  CGAL::Timer timer_mu_k_test5;
	CGAL::Timer timer_circle_test;
	CGAL::Timer timer_circle_test1;
	CGAL::Timer timer_circle_test2;
	CGAL::Timer timer_circle_test3;
	CGAL::Timer timer_circle_test4;
	CGAL::Timer timer_circle_test5;
	CGAL::Timer timer_circle_test6;
    
};

} // namespace internal

} // namespace CGAL

#endif // CGAL_STAR_NEIGHBORHOOD_DESCARTES_H



