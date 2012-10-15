// TODO: Add licence
//
// TODO: remove _m_left, _m_scale and _m_denom
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

/*!\file include/CGAL/Algebraic_kernel_d/Controlled_neighborhood_descartes.h
 * \brief defines class \c Controlled_neighborhood_descartes
 *
 * Isolates real roots of an approximated polynomial. Whenever the isolating
 * intervals get too small with respect to a given precision, it is reported.
 */


#ifndef CGAL_BITSTREAM_ISOLATORS_CONTROLLED_NEIGHBORHOOD_DESCARTES_H
#define CGAL_BITSTREAM_ISOLATORS_CONTROLLED_NEIGHBORHOOD_DESCARTES_H 1

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

#include <CGAL/Algebraic_kernel_d/univariate_polynomial_utils.h>
#include <CGAL/Algebraic_kernel_d/construct_binary.h>

#include <CGAL/Timer.h>


namespace CGAL {

namespace internal {


template <class Polynomial_, class Bound_> 
class Controlled_neighborhood_descartes {
  
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
  
  
  //!\name Constructors
  //!@{

public:

  /*!
   * Default constructor.
   */
  Controlled_neighborhood_descartes() :
    _m_poly(Polynomial(1)),
    _m_degree(1),
    _m_precision(1),
    _m_number_of_real_roots(0),
    _m_precision_sufficient(true) {

    _m_numerator = new Numerator[1]; 
    _m_denominator_exponent = new Denominator[1];
    _m_is_exact = new bool[1];
  }
  
  /*!
   * Constructor from univariate square free polynomial.
   * \param p the given polynomial
   * \param precision  
   * \pre the polynomial is square free      
   */
  Controlled_neighborhood_descartes(
      const Polynomial& p, 
      Bound precision = Bound(2)) :
    _m_poly(p),
    _m_degree(CGAL::degree(_m_poly)),
    _m_precision(precision),
    _m_number_of_real_roots(0),
    _m_precision_sufficient(true) {

    if (0 == _m_degree) { 
      if (_m_poly.is_zero()) {
        _m_number_of_real_roots = -1;
        return;
      }
    }
    
    _m_numerator = new Numerator[_m_degree]; 
    _m_denominator_exponent = new Denominator[_m_degree];
    _m_is_exact = new bool[_m_degree];

    CGAL::Timer timer_n_th_root;
    timer_n_th_root.start();
    _m_precision_n_th_root = n_th_root_of_precision();
    timer_n_th_root.stop();
    CGAL_real_roots_extended_log("\t \t Runtime n_th_root_of_precision: "
                        << timer_n_th_root.time() << std::endl)
    
    isolate();
  }
  
  /*!
   * copy constructor
   * \param cnd copy instance
   */
  Controlled_neighborhood_descartes(
      const Controlled_neighborhood_descartes& t):
    _m_poly(t._m_poly), 
    _m_degree(t._m_degree), 
    _m_number_of_real_roots(t._m_number_of_real_roots),

    _m_precision(t._m_precision),
    _m_precision_sufficient(t._m_precision_sufficient),
    _m_precision_n_th_root(t._m_precision_n_th_root) {
    
    _m_numerator = new Numerator[_m_degree]; 
    _m_denominator_exponent = new Denominator[_m_degree];
    _m_is_exact = new bool[_m_degree]; 
    
    for (int i = 0; i < _m_degree; i++) {
      _m_numerator[i] = t._m_numerator[i];
      _m_denominator_exponent[i] = t._m_denominator_exponent[i];
      _m_is_exact[i] = t._m_is_exact[i];
    }
  }
  
  //!@}

  //!\name Desctructors
  //!@{
  
  //! destructor
  ~Controlled_neighborhood_descartes() {
    delete[] _m_numerator;
    delete[] _m_denominator_exponent;
    delete[] _m_is_exact;
  }
  
  //!@}
  
  //! = operator
  void operator=(const Controlled_neighborhood_descartes& cnd) {
    delete[] _m_numerator;
    delete[] _m_denominator_exponent;
    delete[] _m_is_exact;

    _m_poly = cnd._m_poly;
    _m_degree = cnd._m_degree;
    _m_number_of_real_roots = cnd._m_number_of_real_roots;

    _m_numerator = new Numerator[_m_degree]; 
    _m_denominator_exponent = new Denominator[_m_degree];
    _m_is_exact = new bool[_m_degree];
    
    for (int i = 0; i < _m_degree; i++) {
      _m_numerator[i] = cnd._m_numerator[i];
      _m_denominator_exponent[i] = cnd._m_denominator_exponent[i];
      _m_is_exact[i] = cnd._m_is_exact[i];
    }
    
    _m_precision = cnd._m_precision;
    _m_precision_sufficient = cnd._m_precision_sufficient;
    _m_precision_n_th_root = cnd._m_precision_n_th_root;
  }   
  
public:
  
  /*!
   * returns the number of real roots
   * \return number of real roots
   */
  int number_of_real_roots() const { 
    return _m_number_of_real_roots; 
  }
  
  /*!
   * \brief returns true if the isolating interval is degenerated to a 
   *  single point.
   * 
   * If is_exact_root(i) is true, 
   * then left_bound(int i) equals  \f$root_i\f$. \n
   * If is_exact_root(i) is true, 
   * then right_bound(int i) equals  \f$root_i\f$. \n 
   * \param i number of root
   * \return true if root i is exact, false otherwise
   */
  bool is_exact_root(int i) const { 
    return _m_is_exact[i]; 
  }
  
  //! returns whether used precision sufficed while executing/constructing
  bool precision_sufficient() const { 
    return _m_precision_sufficient; 
  }
  
public:   
  
  /*!
   * computes the left bound of the isolating interval for root i
   * \param i number of root
   * \param numerator variable in which numerator of left bound is written
   * \param denominator variable in which denominator of left bound is written 
   */
  void left_bound(int i, 
                  Numerator& numerator, Denominator& denominator) const {

    CGAL_assertion(i >= 0);
    CGAL_assertion(i < _m_number_of_real_roots);

    Numerator n_th_root_num;
    Denominator n_th_root_den;

    typename FT_bound::Decompose decompose;
    decompose(_m_precision_n_th_root,n_th_root_num,n_th_root_den);

    construct_binary(_m_denominator_exponent[i], denominator);
    numerator = Numerator(2) * _m_numerator[i] + denominator;
    denominator = denominator * Denominator(4);
    numerator = numerator * n_th_root_den -
                Numerator(9) * n_th_root_num * denominator;
    denominator = denominator * n_th_root_den;
  }
  
  /*!
   * computes the right bound of the isolating interval for root i
   * \param i number of root
   * \param numerator variable in which numerator of right bound is written
   * \param denominator variable in which denominator of right bound is written 
   */
  void right_bound(int i,
                   Numerator& numerator, Denominator& denominator) const {

    CGAL_assertion(i >= 0);
    CGAL_assertion(i < _m_number_of_real_roots);

    Numerator n_th_root_num;
    Denominator n_th_root_den;

    typename FT_bound::Decompose decompose;
    decompose(_m_precision_n_th_root,n_th_root_num,n_th_root_den);
    construct_binary(_m_denominator_exponent[i], denominator);

    if (_m_is_exact[i]){
      numerator = Numerator(2) * _m_numerator[i] + denominator;
    } else {
      numerator = Numerator(2) * (_m_numerator[i]+1) + denominator;
    }
    denominator = denominator * Denominator(4);
    numerator = numerator * n_th_root_den +
                Numerator(9) * n_th_root_num * denominator;
    denominator = denominator * n_th_root_den;
  }
  
private:

  /*!
   * computes an upper bound of the n-th root of _m_precision as follows:
   * decompose precision into numerator and denominator
   * compute root1 of y1(x) = x^n - numerator
   * compute root2 of y2(x) = x^n - denominator
   * compute root1(rounded up) / root2(rounded down)
   * \return n-th root of _m_precision
   */
  Bound n_th_root_of_precision () {

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

    Numerator precision_num;
    Denominator precision_den;

    typename FT_bound::Decompose decompose;
    decompose(_m_precision,precision_num,precision_den);

    Polynomial_1 x = typename PT_1::Shift()(Polynomial_1(1),1);
    Polynomial_1 n_th_root_polynomial_num = x;
    for (int j = 0; j < _m_degree-1; j++) {
      n_th_root_polynomial_num *= x;
    }
    n_th_root_polynomial_num = n_th_root_polynomial_num - precision_num;

    Polynomial_1 n_th_root_polynomial_den = x;
    for (int j = 0; j < _m_degree-1; j++) {
      n_th_root_polynomial_den *= x;
    }
    n_th_root_polynomial_den = n_th_root_polynomial_den - precision_den;

    solve(n_th_root_polynomial_num, std::back_inserter(roots_num));
    Algebraic_real_1 n_th_root_num = roots_num.back().first;  

    solve(n_th_root_polynomial_den, std::back_inserter(roots_den));
    Algebraic_real_1 n_th_root_den = roots_den.back().first;

    int approximation_precision_exponent;
    if (_m_degree == 0) {
      approximation_precision_exponent = 0;
    }
    else {
      approximation_precision_exponent =
            CGAL::internal::ceil_log2_abs(Numerator(_m_degree)) * 7;
    }
    Approximate_relative_1 approximation;  
    std::pair<Bound,Bound> isolating_interval_num =
              approximation(n_th_root_num,approximation_precision_exponent);
    std::pair<Bound,Bound> isolating_interval_den =
              approximation(n_th_root_den,approximation_precision_exponent);
    return isolating_interval_num.second / isolating_interval_den.first;
  }

  
  /*!
   * scales and translates _m_poly onto [0,1] 
   * starts extended descartes algorithm
   */
  void isolate() {

    Polynomial p = ::CGAL::scale_down(_m_poly,Coefficient(4));
    Polynomial r = ::CGAL::translate_by_one(p);
    Polynomial q = ::CGAL::scale_up(r, Coefficient(2));

    zero_one_descartes(q,0,0);

    CGAL_real_roots_log("\t \t Runtime check length of interval: "
                        << timer_check_interval.time() << std::endl)
    CGAL_real_roots_log("\t \t Runtime transform polynomials: "
                        << timer_transform_polys.time() << std::endl)

    CGAL_real_roots_log("\t \t \t Runtime transform polynomials 1: "
                        << timer_transform_polys1.time() << std::endl)
    CGAL_real_roots_log("\t \t \t Runtime transform polynomials 2: "
                        << timer_transform_polys2.time() << std::endl)
    CGAL_real_roots_log("\t \t \t Runtime transform polynomials 3: "
                        << timer_transform_polys3.time() << std::endl)
    CGAL_real_roots_log("\t \t \t Runtime transform polynomials 4: "
                        << timer_transform_polys4.time() << std::endl)

    CGAL_real_roots_log("\t \t Runtime one root: "
                        << timer_one_root.time() << std::endl)
    CGAL_real_roots_log("\t \t Runtime refine: "
                        << timer_refine.time() << std::endl)

  }
  
  //! returns the polynomial $(1 + x)^n p(1/(1 + x))$.
  Polynomial variation_transformation(const Polynomial& p) { 
    Polynomial r = reversal(p);
    return translate_by_one(r); 
  }
  
  /*!
   * Descartes algoritm to determine isolating intervals for the roots 
   * lying in the interval (0,1).
   * The parameters $(i,exp)$ describe the interval $(i/2^exp, (i+1)/2^exp)$.
   * Here $0\leq i < 2^exp$.
   * \param p polynomial with roots in [0,1]
   * \param i numerator of interval: [i/2^exp, (i+1)/2^exp]
   * \param exp denominator of interval: [i/2^exp, (i+1)/2^exp]
   */
  void zero_one_descartes(const Polynomial& p,
                          Numerator i, Denominator exp) { 
						  
	typedef CGAL::Polynomial_traits_d<Polynomial> PT;
	typename PT::Evaluate evaluate;

    if (!_m_precision_sufficient) {
      return;
    }
    
    //TODO: length of interval computed correctly?
    timer_check_interval.start();
    Bound length_interval = 
      (Bound(5)/Bound(power<Numerator>(Numerator(2), exp))); 
    
    //TODO: correct?
	  //if length of extended interval < 19n * nth-root(precision) stop algorithm
	  //and report that precision isn't sufficient
    if (CGAL::compare(length_interval,
                      Bound(19 * _m_degree) * _m_precision_n_th_root)
                      == CGAL::SMALLER ) {
      _m_precision_sufficient = false;
      return;
    }
    timer_check_interval.stop();
    
    // Determine the number of sign variations of the transformed 
    // polynomial $(1+x)^nP(1/(1+x))$. This gives the number of 
    // roots of $P$ in $(0,1)$.
    timer_transform_polys.start();

    timer_transform_polys1.start();
    Polynomial r = variation_transformation(p);
    timer_transform_polys1.stop();

    timer_transform_polys2.start();
  	int descarte = sign_variations(r);
    timer_transform_polys2.stop();

	  //Map polynomial p onto $(-2,3)$.
	  // Determine the number of sign variations of the transformed 
    // polynomial $(1+x)^nP(1/(1+x))$.
    timer_transform_polys3.start();
    Polynomial q_extended = scale_up(p,Coefficient(-2));
    Polynomial r_extended = translate_by_one(q_extended);
    Polynomial s_extended = scale_up(r_extended,Coefficient(-5));
    Polynomial t_extended = scale_down(s_extended,Coefficient(2));
    timer_transform_polys3.stop();

    timer_transform_polys4.start();
    Polynomial u_extended = variation_transformation(t_extended);
    int descarte_extended = sign_variations(u_extended);
    timer_transform_polys4.stop();

    timer_transform_polys.stop();

    // no root
    if (descarte == 0 && descarte_extended <= 1) {
      return;
    }
    
    // exactly one root
    // Note the termination criterion $P(0)\neq 0$ and $P(1)\neq 0$.
    // This ensures that the given interval is an isolating interval.
    timer_one_root.start();
    if ( descarte == 1 
           && descarte_extended == 1
           && CGAL::compare(get_coefficient(p,0), Coefficient(0)) != CGAL::EQUAL 
           && CGAL::compare(evaluate(p,Coefficient(1)),
                            Coefficient(0)) != CGAL::EQUAL ) { 
      _m_numerator[_m_number_of_real_roots] = i;
      _m_denominator_exponent[_m_number_of_real_roots] = exp;
      _m_is_exact[_m_number_of_real_roots] = false;
      _m_number_of_real_roots++;
      timer_one_root.stop();
      return; 
    }
    timer_one_root.stop();
	
    // more than one root
    // Refine the interval.
    timer_refine.start();
    i = Numerator(2) * i; 
    exp = exp + Denominator(1);
      
    // Transform the polynomial such that the first half of the interval
    // is mapped to the unit interval.
    Polynomial q1 = scale_down(p, Coefficient(2));
	  Polynomial q2 = translate_by_one(q1);
    timer_refine.stop();

    // Consider the first half of the interval.
    zero_one_descartes(q1, i, exp);  

	  // Test if the polynomial is zero at the midpoint of the interval 
    timer_refine.start();
	  if (CGAL::compare(get_coefficient(q2,0), Coefficient(0)) == CGAL::EQUAL) { 
      _m_numerator[_m_number_of_real_roots] = i + Numerator(1);
      _m_denominator_exponent[_m_number_of_real_roots] = exp;
      _m_is_exact[_m_number_of_real_roots] = true;
      _m_number_of_real_roots++;
    }
    timer_refine.stop();
      
    // Consider the second half of the interval. 
    zero_one_descartes(q2, i + Numerator(1), exp); 
  }

  
  template <class NT>
  /*!
   *  computes base^e
   * \param base base
   * \param e exponent
   * \return base^e
   */
  NT power (const NT& base, Numerator e) {
    NT result(base);
    if (CGAL::compare(e, Numerator(0)) == CGAL::EQUAL) {
      return NT(1);
    }
    for (Numerator i = Numerator(0); i < e - Numerator(1); i++) {
      result *= base;
    }
    return result;
  }
  
private:
  
  //! Polynomial
  Polynomial _m_poly;

  //! degree of polynomial
  int _m_degree;
  
  //! number of real roots
  int _m_number_of_real_roots;
  
  //! numerator
  Numerator* _m_numerator;
  
  //! denominator exponent
  Denominator* _m_denominator_exponent; 
  
  //! exact root?
  bool* _m_is_exact;
  
  //! n-th root of precision
  Bound _m_precision_n_th_root;
  
  //! precision value
  Bound _m_precision;

  //! is precision sufficient
  bool _m_precision_sufficient;

  //! various timer for benchmarking
  CGAL::Timer timer_check_interval;
  CGAL::Timer timer_transform_polys;
  CGAL::Timer timer_transform_polys1;
  CGAL::Timer timer_transform_polys2;
  CGAL::Timer timer_transform_polys3;
  CGAL::Timer timer_transform_polys4;
  CGAL::Timer timer_one_root;
  CGAL::Timer timer_refine;

};

} // namespace internal

} // namespace CGAL

#endif // CGAL_CONTROLLED_NEIGHBORHOOD_DESCARTES_H

