// Copyright (c) 2006-2012 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//                 Eric Berberich <eric.berberich@cgal.org>

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_ANALYSIS_2_ALCIX_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_CURVE_ANALYSIS_2_ALCIX_H

#include <CGAL/config.h>
#include <CGAL/Algebraic_kernel_d/flags.h>

#include <vector>
#include <set>
#include <map>

#include <boost/bind.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <CGAL/Cache.h>
#include <CGAL/function_objects.h>
#include <CGAL/Handle_with_policy.h>
#include <CGAL/Arr_enums.h>

// TODO 2012 see question on use of an isolator below
#if 0 && CGAL_BISOLVE_ENABLE_ARCAVOID
#include <CGAL/Algebraic_kernel_d/Real_roots.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>
#include <CGAL/Algebraic_kernel_1/Arcavoid_root_isolator.h>
#endif

#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Algebraic_kernel_d/macros.h>
#include <CGAL/Algebraic_kernel_d/exceptions.h>
#include <CGAL/Algebraic_kernel_d/enums.h>
#include <CGAL/Algebraic_kernel_d/algebraic_curve_kernel_2_tools.h>
#include <CGAL/Algebraic_kernel_d/Interval_evaluate_1.h>
#include <CGAL/Algebraic_kernel_d/Status_line_CA_1.h>

#if CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE

#if CGAL_ACK_FACTORIZE_UNI_POLYNOMIALS
#include <CGAL/symbolic_exports.h>
#endif

#include <CGAL/Algebraic_kernel_d/Curve_analysis_2_geotop_line_builder.h>
#else
#include <CGAL/Algebraic_kernel_d/Curve_analysis_2_mk_line_builder.h>
#include <CGAL/Algebraic_kernel_d/Shear_controller.h>
#include <CGAL/Algebraic_kernel_d/Shear_transformation.h>
#include <CGAL/Algebraic_kernel_d/shear.h>
#endif


#ifndef Bisolve_telemetry_code
# define Bisolve_telemetry_code(x)
#endif

#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
// put includes here
#endif


#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4290)
#endif

#ifdef CGAL_WEBXACLI_MULTITHREADED
extern pthread_mutex_t CGAL_static_mtx;
#endif

namespace CGAL {

template<typename AlgebraicKernelWithAnalysis_2,
         typename Rep_>
class Curve_analysis_2;

namespace internal {

template<typename Comparable,bool has_template_typedefs>
  struct Is_derived_from_Handle_with_policy {
    typedef boost::false_type Tag;
};

template<typename Comparable>
  struct Is_derived_from_Handle_with_policy<Comparable,true> {

    typedef typename
      boost::is_base_of< CGAL::Handle_with_policy
                             < typename Comparable::T,
			       typename Comparable::Handle_policy,
                               typename Comparable::Allocator >,
                         Comparable
                       >::type Tag;
};


template<typename Comparable,typename Tag> struct Compare_for_vert_line_map_
    {
      bool operator() (const Comparable& a, const Comparable& b) const {
	return a<b;
      }
};

template<typename Comparable>
  struct Compare_for_vert_line_map_<Comparable,boost::true_type> {

    bool operator() (const Comparable& a, const Comparable& b) const {
      return CGAL::Handle_id_less_than< Comparable >()(a,b);
    }
};

template<typename Comparable> struct Compare_for_vert_line_map
  : public std::binary_function<Comparable,Comparable,bool> {

  BOOST_MPL_HAS_XXX_TRAIT_DEF(T)
  BOOST_MPL_HAS_XXX_TRAIT_DEF(Handle_policy)
  BOOST_MPL_HAS_XXX_TRAIT_DEF(Allocator)

  typedef typename CGAL::internal::Is_derived_from_Handle_with_policy
    < Comparable,
      has_T<Comparable>::value &&
      has_Handle_policy<Comparable>::value &&
      has_Allocator<Comparable>::value>::Tag Tag;

  public:

  bool operator() (const Comparable& a, const Comparable& b) const {

    return eval(a,b);
  }

  private:

  Compare_for_vert_line_map_<Comparable,Tag> eval;


};


// \brief Representation class for algebraic curves.
template< typename AlgebraicKernelWithAnalysis_2>
class Curve_analysis_2_rep {

public:
    //! this instance's template parameter
    typedef AlgebraicKernelWithAnalysis_2 Algebraic_kernel_with_analysis_2;

    //! the class itself
    typedef Curve_analysis_2_rep< Algebraic_kernel_with_analysis_2 > Self;

    //! The handle class
    typedef CGAL::Curve_analysis_2
        <Algebraic_kernel_with_analysis_2,Self> Handle;

    //protected:
public:

    CGAL_ACK_SNAP_ALGEBRAIC_CURVE_KERNEL_2_TYPEDEFS(Handle);

    //! \brief \c Vert_line specification for critical x-values
    typedef internal::Status_line_CA_1< Algebraic_kernel_with_analysis_2 > Status_line_1;

#if CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
    typedef internal::Geotop_line_builder< Algebraic_kernel_with_analysis_2 > Status_line_builder;
#else
    typedef internal::Mk_line_builder< Algebraic_kernel_with_analysis_2 > Status_line_builder;
#endif

    typedef int size_type;

    typedef std::map< Bound, Status_line_1 >
      Vert_line_at_rational_map;

    typedef
    std::map< Algebraic_real_1,
              Status_line_1,
              internal::Compare_for_vert_line_map<Algebraic_real_1> >
    Vert_line_map;

    //!\name Constructors
    //!@{

    //! Default constructor
    Curve_analysis_2_rep() :
    _m_status_line_builder(0)
    {
    }

    //! Constructor with polynomial
    Curve_analysis_2_rep(Algebraic_kernel_with_analysis_2 *kernel,
                         const Polynomial_2& poly,
                         CGAL::Degeneracy_strategy strategy) :
      _m_kernel(kernel), _m_status_line_builder(0),
      f(poly), degeneracy_strategy(strategy)
    {
    }

    //!@}

    ~Curve_analysis_2_rep() {
      delete _m_status_line_builder;
    }

private:

    typedef internal::LRU_hashed_map<
        Bound,
        std::vector<Algebraic_real_1>,
        internal::To_double_hasher > Intermediate_cache;

    Intermediate_cache intermediate_cache;

    // Internal information struct about x-coordinates
    struct Event_coordinate_1 {
        Event_coordinate_1(){} //added to solve a compilation error of gcc-3.4 (bug?)
        Algebraic_real_1 val;
        size_type mult_of_prim_res_root;
        size_type index_of_prim_res_root;
        size_type mult_of_content_root;
        size_type index_of_content_root;
        size_type mult_of_prim_lcoeff_root;
        size_type index_of_prim_lcoeff_root;
        boost::optional<Status_line_1> stack;
    };

    // Functor to get the X_coordinate of an Event_coordinate
    struct Val_functor {
        typedef Event_coordinate_1 argument_type;
        typedef Algebraic_real_1 result_type;
        result_type operator() (argument_type event) const {
            return event.val;
        }
    };


    //! The object holding the information about events, as an optional
    mutable boost::optional<std::vector<Event_coordinate_1> >
        event_coordinates;

    //! The algebraic kernel to use
    Algebraic_kernel_with_analysis_2* _m_kernel;

    //! a single status line builder for this curve
    mutable Status_line_builder*  _m_status_line_builder;

    //! The polynomial defining the curve
    boost::optional<Polynomial_2> f;

    //! How degenerate situations are handled
    CGAL::Degeneracy_strategy degeneracy_strategy;

    /*!
     * \brief The polynomial without its content (the gcd of the coeffs).
     *
     * The content is the greatest common divisor of the coefficients of \c f
     * considered as polynomial <tt>y</tt>. \c The polynomial f_primitive is
     * \c f/cont(f). The corresponding curve is equal to the curve of \c f,
     * only without vertical line components.
     */
    mutable boost::optional<Polynomial_2> f_primitive;

    //! the polynomial containing all roots of the resultant of the primitive
    //! part of f and its y-derivative
    mutable boost::optional<Polynomial_1>
        resultant_of_primitive_and_derivative_y;

    //! The content of f
    mutable boost::optional<Polynomial_1> content;

#if CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE

    //! the polynomial containing all roots of the resultant of the primitive
    //! polynomial's partial derivative wrt to x and y (without common factor)
    mutable boost::optional< Polynomial_1 >
        resultant_of_coprime_derivative_x_and_derivative_y;

    //! the factors of another resultant
    // TODO EBEB delete if new strategy is chosen
    mutable boost::optional<std::vector< std::pair< Polynomial_1, Multiplicity_type > > >
        resultant_of_coprime_derivative_x_and_derivative_y_factors;

    //! and multiplied polynomials
    mutable boost::optional<std::vector< Polynomial_1 > >
      resultant_of_coprime_derivative_x_and_derivative_y_sequence;

#else
    //! the polynomial containing all roots of the resultant of the primitive
    //! part of f and its x-derivative
    mutable boost::optional<Polynomial_1>
        resultant_of_primitive_and_derivative_x;

    //! The Sturm-Habicht polynomials of f
    mutable boost::optional<std::vector<Polynomial_2> >
        sturm_habicht_of_primitive;

    //! The non-working shear factors, as far as known
    mutable std::set<Integer> bad_shears;

    //! The already known shear factors
    mutable std::map<Integer,Handle> sheared_curves;
#endif

    //! Has the curve vertical line components
    mutable boost::optional<bool> has_vertical_component;

    //! The intermediate values
    mutable boost::optional<std::vector<boost::optional<Bound> > >
    intermediate_values;

    //! stores Y_values at rational coordinate
    mutable Vert_line_at_rational_map vert_line_at_rational_map;

    //! stores vert_lines
    mutable Vert_line_map vert_line_map;

    /**! \brief Information about whether arcs at +/- infty
     *   are asymptotic to y=beta,
     *   or go to +/- infty also in y-direction
     */
    mutable boost::optional<std::vector<CGAL::Object> >
    horizontal_asymptotes_left, horizontal_asymptotes_right;

    //! friends
    friend class ::CGAL::Curve_analysis_2<Algebraic_kernel_with_analysis_2,Self>;

}; // class Curve_analysis_2_rep
} // namespace internal


/*!
 * \brief Analysis for algebraic curves of arbitrary degree.
 *
 * This class constitutes a model for the concept
 * AlgebraicKernelWithAnalysis_d_2::CurveAnalysis_2.
 * For a square-free bivariate polynomial \c f, a topologic-geometrical
 * analysis of the algebraic curve defined by the vanishing set of \c f
 * is provided. This means, one can ask for the total number, and the position
 * of the critical x-coordinates of the curve, and for each x-coordinate,
 * geometric information about the curve can be obtained. This data
 * is capsuled into an object of type \c Curve_analysis_2::Status_line_1,
 * which is in fact a \c Status_line_CA_1 object.
 *
 * The restriction to square-free curves is a weak one, since the curves
 * can be made square-free before passed to the analysis.
 * The \c Construct_curve_2 functor of \c Algebraic_curve_kernel_2 is
 * doing so, thus it accepts arbitrary bivariate polynomials.
 *
 * The analysis is implemented in a "lazy" fashion. This means, when
 * created, the analysis delays all computations until the information
 * is queried for the first time. This means, if only parts of the curves
 * are of interest, only a partial analysis is performed.
 * We remark that nevertheless, the global \e projection \e step
 * (i.e., computing the (sub)resultants) must be done once a \c Status_line_1
 * is queried. Often, this step forms the bottleneck in the whole computation.
 *
 * For more details of the algorithm, consult the reference:
 * A.Eigenwillig, M.Kerber, N.Wolpert: Fast and Exact Geometric Analysis of
 * Real Algebraic Plane Curves. Proceedings of the International Symposium
 * on Symbolic and Algebraic Computation (ISSAC 2007), pp. 151-158
 */
template<typename AlgebraicKernelWithAnalysis_2,
  typename Rep_
   = internal::Curve_analysis_2_rep< AlgebraicKernelWithAnalysis_2>
>
class Curve_analysis_2 : public ::CGAL::Handle_with_policy< Rep_ > {

    //! \name typedefs
    //! @{

public:
    //! this instance' first template parameter
    typedef AlgebraicKernelWithAnalysis_2 Algebraic_kernel_with_analysis_2;

    //! this instance' second template parameter
    typedef Rep_ Rep;

private:

    //! The internal type for event coordinates
    typedef typename Rep::Event_coordinate_1 Event_coordinate_1;

    // Base class
    typedef ::CGAL::Handle_with_policy<Rep> Base;

    // Internal class to build status lines at
    typedef typename Rep::Status_line_builder Status_line_builder;

    // This type
    typedef CGAL::Curve_analysis_2<Algebraic_kernel_with_analysis_2,Rep> Self;

public:

    //! Indexing type
    typedef typename Rep::size_type size_type;

    CGAL_ACK_SNAP_ALGEBRAIC_CURVE_KERNEL_2_TYPEDEFS(Self);

    typedef typename Rep::Status_line_1 Status_line_1;

    //! Required by the CurveKernel_2 concept
    typedef Algebraic_real_1 Coordinate_1;

    //! Traits type for Polynomial_2
    typedef CGAL::Polynomial_traits_d<Polynomial_1> Polynomial_traits_1;

    //! Traits type for Polynomial_2
    typedef CGAL::Polynomial_traits_d<Polynomial_2> Polynomial_traits_2;

private:

    /*!
     * \brief Coercion between the coefficient type of the polynomial
     * and the bound type of the curve analysis
     *
     * Interoperability of both types is required
     */
    typedef CGAL::Coercion_traits<Bound, Coefficient> Coercion;

    /*!
     * \brief The common supertype that both the coefficient and the bound
     * type are convertible to
     */
    typedef typename Coercion::Type Coercion_type;

    //! Polynomial over the \c Coercion_type
    typedef typename CGAL::Polynomial_traits_d<Polynomial_2>
        ::template Rebind<Coercion_type,1>::Other::Type Poly_coer_1;

public:

    //! Type to represent points on curves
    typedef typename Algebraic_kernel_with_analysis_2::Algebraic_real_2
      Algebraic_real_2;

    //! Required by the CurveKernel_2 concept
    typedef Algebraic_real_2 Coordinate_2;

    //! type for horizontal asymtote values
    typedef CGAL::Object Asymptote_y;

    //! @}

private:

    //! \name Helping structs
    // @{

    struct Event_functor {
        Event_functor(const Self* curve) : curve(curve) {}
        const Self* curve;
        typedef size_type argument_type;
        typedef Status_line_1 result_type;
        result_type operator() (argument_type index) const {
            return curve->status_line_at_event(index);
        }
    };

    struct Intermediate_functor {
        Intermediate_functor(const Self* curve) : curve(curve) {}
        const Self* curve;
        typedef size_type argument_type;
        typedef Status_line_1 result_type;
        result_type operator() (argument_type index) const {
            return curve->status_line_of_interval(index);
        }
    };

#if !CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
    struct Stha_functor {
        Stha_functor(const Self* curve) : curve(curve) {}
        const Self* curve;
        typedef size_type argument_type;
        typedef Polynomial_1 result_type;
        result_type operator() (argument_type index) const {
            return curve->principal_sturm_habicht_of_primitive(index);
        }
    };
#endif

    //! @}

public:

    //! \name Iterators
    //! @{

    //! Iterator type for status lines at events
    typedef boost::transform_iterator<Event_functor,
                              boost::counting_iterator<size_type> >
    Event_line_iterator;

    //! Iterator type for status lines of intervals
    typedef boost::transform_iterator<Intermediate_functor,
                              boost::counting_iterator<size_type> >
    Intermediate_line_iterator;

#if !CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
    //! Iterator type for the principal sturm habicht coefficients of the curve
    typedef boost::transform_iterator<Stha_functor,
                              boost::counting_iterator<size_type> >
    Principal_sturm_habicht_iterator;
#endif

    //! @}

public:

    //!\name Constructors
    //!@{

    //! default Constructor
    Curve_analysis_2() {}

    /*!
     * \brief Constructs the curve analysis for the given polynomial
     *
     * Analyses the curve that is defined by the vanishing set of the
     * polynomial \c f.
     * \pre \c f is square free.
     * \param strategy The default strategy
     * (\c SHEAR_ONLY_AT_IRRATIONAL_STRATEGY)
     * is to \c shear the curve
     * if a degenerate situation is detected during the analysis,
     * except at rational x-coordinates where the curve can be analysed
     * more directly. The analysis
     * is then performed in  the sheared system, and finally translated back
     * into the original system.
     * Using \c SHEAR_STRATEGY, a shear is triggered also for degeneracies
     * at rational x-coordinate. With both strategies, it is guaranteed that
     * the analysis works successfully for any square free input curve.
     * On the other hand, the EXCEPTION_STRATEGY throws an exception of type
     * \c internal::Zero_resultant_exception<Polynomial_2>,
     * instead of performing a shear.
     *
     * \Todo Currently the defualt strategy has been changed to SHEAR_STRATEGY
     * because there exist a problem if vertical asymtotes are present at
     * the rational x-coordinate.
     */
    explicit Curve_analysis_2(Algebraic_kernel_with_analysis_2 *kernel,
                              const Polynomial_2& f,
                              CGAL::Degeneracy_strategy strategy
                                  = CGAL_ACK_DEFAULT_DEGENERACY_STRATEGY)
        throw(internal::Zero_resultant_exception<Polynomial_2>)
        : Base(Rep(kernel,f,strategy))
    {
    }

    //! \brief Copy constructor
    Curve_analysis_2(const Self& alg_curve)
        : Base(static_cast<const Base&>(alg_curve))
    {
    }


    //!@}


    //! \name Members
    //!@{

private:

    /*
     * \brief Sets all status lines at events and of intervals
     *
     * Writes the status lines of events and interval into the object.
     * The value type of both \c InputIterator1 and \c InputIterator2
     * is \c Status_line_1.
     */
    template<typename InputIterator1,typename InputIterator2>
    void set_event_lines(InputIterator1 event_begin,
                         InputIterator1 event_end,
                         InputIterator2 intermediate_begin,
                         InputIterator2 CGAL_precondition_code(intermediate_end)) const {

        if(! this->ptr()->event_coordinates) {

            std::vector<Event_coordinate_1> event_coordinate_vector;

            for(InputIterator1 it = event_begin; it != event_end; it++) {
                Event_coordinate_1 curr_event;
                curr_event.val = it->x();
                event_coordinate_vector.push_back(curr_event);
            }
            this->ptr()->event_coordinates = event_coordinate_vector;

        }

        InputIterator1 it1 = event_begin;
        for(size_type i = 0; i < number_of_status_lines_with_event() ; i++ ) {
            this->ptr()->vert_line_map[event_coordinates()[i].val] = *it1;
            event_coordinates()[i].stack = *it1;

            it1++;
        }
        CGAL_assertion(it1 == event_end);

        if(! this->ptr()->intermediate_values) {
            this->ptr()->intermediate_values
                = std::vector<boost::optional<Bound> >
                    (number_of_status_lines_with_event()+1);
        }

        InputIterator2 it2 = intermediate_begin;
        for(size_type i = 0;
            i < static_cast<int>(intermediate_values().size());
            i++,it2++) {
            CGAL_assertion(it2->x().is_rational());
            Bound q = it2->x().rational();

            intermediate_values()[i] = q;
            this->ptr()->vert_line_map[it2->x()] = *it2;
            this->ptr()->vert_line_at_rational_map[q] = *it2;

        }
        CGAL_assertion(it2 == intermediate_end);

    }

public:

    /*! \brief Returns whether the curve has a valid defining polynomial
     */
    bool has_defining_polynomial() const {
        return this->ptr()->f;
    }

public:

    /*! \brief Sets the defining polynomial.
     *
     * \pre The object has no defining polynomial yet.
     */
    void set_f(Polynomial_2 f) {
        CGAL_precondition(! has_defining_polynomial());
        if((! this->ptr()->f) || f!=this->ptr()->f.get()) {
            this->copy_on_write();
            this->ptr()->f=f;
        }
    }


public:

    /*!
     * \brief Returns whether the curve is y-regular
     *
     * A curve is called y-regular if the leading coefficient of its defining
     * polynomial wrt y is a constant, i.e., contains no x
     */
    bool is_y_regular() const {
#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_is_y_regular();
        }
#endif
        return CGAL::degree(CGAL::leading_coefficient(polynomial_2())) == 0;
    }

public:

    /*!
     * \brief returns whether the curve contains a vertical line as a component
     *
     * In algebraic terms, this methods computes whether the content
     * of its defining polynomial has a real root.
     */
    bool has_vertical_component() const {
#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_has_vertical_components();
        }
#endif
        if(is_y_regular()) {
	    this->ptr()->has_vertical_component = false;
        }
        if(! this->ptr()->has_vertical_component) {
            // This is computed as side effect
            // when the event coordinates are computed
            event_coordinates();
            CGAL_assertion(this->ptr()->has_vertical_component);
        }
        return this->ptr()->has_vertical_component.get();
    }

public:

    //! Returns the defining polynomial
    Polynomial_2 polynomial_2() const {
        CGAL_precondition(this->ptr()->f);
        return this->ptr()->f.get();
    }

public:

    /*!
     * \brief Returns the number of event lines of the curve
     *
     * Algebraically, the number of real roots of the discriminant of
     * the curve's defining equation is returned.
     */
    size_type number_of_status_lines_with_event() const {
        CGAL_precondition(this->ptr()->f);
#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_number_of_status_lines_with_event();
        }
#endif
        return static_cast<size_type>(event_coordinates().size());
    }

public:

    /*!
     * \brief Returns whether the given x-coordinate is critical for the curve
     * and which event or interval index the x-coordinate belongs to.
     *
     * \param is_event is set to \c true if the curve has an event
     * at this x-coordinate, or in other words, if the discriminant of its
     * defining polynomial vanishes at \c x
     * \param i is set to the index of the event if \c x is an event. Otherwise
     * \c i is set to the index of the interval \c x is contained in.
     */
    void x_to_index(Algebraic_real_1 x,size_type& i,bool& is_event) const {
#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_x_to_index(x,i,is_event);
        }
#endif
        CGAL_precondition(has_defining_polynomial());
        typename Rep::Val_functor xval;
        i = std::lower_bound(
                ::boost::make_transform_iterator(event_coordinates().begin(),
                                                 xval),
                ::boost::make_transform_iterator(event_coordinates().end(),
                                                 xval),
                x
        ) - ::boost::make_transform_iterator(event_coordinates().begin(),
                                             xval);
        is_event = (i < static_cast<size_type>(event_coordinates().size()) &&
                    (event_coordinates()[i].val == x) );
    }

public:

    //! Returns the status line at the <tt>i</tt>-th event of the curve.
    Status_line_1& status_line_at_event(size_type i) const {

        CGAL_precondition(has_defining_polynomial());
#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_status_line_at_event(i);
        }
#endif
        CGAL_precondition_code(
                size_type n =
                static_cast<size_type>(event_coordinates().size());
        );
        CGAL_precondition(i>=0);
        CGAL_precondition(i<n);
        if(! event_coordinates()[i].stack) {
            Status_line_1 event_line = create_status_line_at_event(i);
            this->ptr()->vert_line_map.insert(std::make_pair(event_coordinates()[i].val, event_line));
            event_coordinates()[i].stack = event_line;
        }
        CGAL_postcondition(event_coordinates()[i].stack.get().is_event());
        return event_coordinates()[i].stack.get();
    }

public:

    //! Returns a status line at the rational <tt>x</tt>-coordinate \c b
    Status_line_1& status_line_at_exact_x(Bound b) const {
#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_status_line_at_exact_x(b);
        }
#endif
        return status_line_at_exact_x(Algebraic_real_1(b));
    }

private:

    /*
     * \brief Returns a status line for an exact value \c alpha that
     * is not an event of the curve
     *
     * This function controls the internal cache that stores already created
     * status line at non-events.
     */
    Status_line_1& status_line_at_exact_non_event_x(Algebraic_real_1 alpha)
        const {

        if(alpha.is_rational()) {

            typename Rep::Vert_line_at_rational_map::iterator it =
                this->ptr()->vert_line_at_rational_map.find
                (alpha.rational());

            if (it != this->ptr()->vert_line_at_rational_map.end()) {
                CGAL_assertion(!it->second.is_event());
                return it->second;
            }
        }

        typename Rep::Vert_line_map::iterator it =
            this->ptr()->vert_line_map.find(alpha);

        if (it != this->ptr()->vert_line_map.end()) {
            CGAL_assertion(!it->second.is_event());
            return it->second;
        }


        // Not stored yet, so create it and store it
        Status_line_1 cvl
            = create_status_line_at_non_event(alpha);

        CGAL_assertion(!cvl.is_event());

        // TODO two vert line maps?
        it = this->ptr()->vert_line_map.insert(it, std::make_pair(alpha, cvl));

        if(alpha.is_rational()) {
            this->ptr()->vert_line_at_rational_map.insert(std::make_pair(alpha.rational(),cvl));
        }

        return it->second;
    }

public:

    //! Returns a vert line for the <tt>x</tt>-coordinate alpha
    Status_line_1& status_line_at_exact_x(Algebraic_real_1 alpha) const {
#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_status_line_at_exact_x(alpha);
        }
#endif
        bool is_event_value;
        size_type index;
        this->x_to_index(alpha,index,is_event_value);
        if(is_event_value) {
            return status_line_at_event(index);
        }
        else {
            return status_line_at_exact_non_event_x(alpha);
        }
    }

private:

    // Creates a status line for the curve's <tt>index</tt>th critical point
    Status_line_1 create_status_line_at_event(size_type index) const
#if !CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
        throw(CGAL::internal::Non_generic_position_exception)
#endif
    {

        Event_coordinate_1& event = event_coordinates()[index];

        Algebraic_real_1 x = event.val;

        try {

#if !CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
#if CGAL_ACK_SHEAR_ALL_NOT_Y_REGULAR_CURVES
            if(event.mult_of_prim_lcoeff_root > 0) {
                throw CGAL::internal::Non_generic_position_exception();
            }
#else
            if(event.mult_of_prim_lcoeff_root > 0) {
                if(event.mult_of_prim_lcoeff_root > 1 ||
                   event.mult_of_prim_res_root > 1) {
                    throw CGAL::internal::Non_generic_position_exception();
                }
            }

#endif
#endif

#if CGAL_ACK_DEBUG_FLAG
            double ev_approx = CGAL::to_double(x);
            CGAL_ACK_DEBUG_PRINT << (index+1) << "th line: "
                                 << std::setw(6) << std::setprecision(3)
                                 << ev_approx
                                 << ".."
                                 << std::flush;
#endif
            size_type left_arcs
                = status_line_for_x(x,CGAL::NEGATIVE).number_of_events();
            size_type right_arcs
                = status_line_for_x(x,CGAL::POSITIVE).number_of_events();

            bool root_of_resultant=(event.mult_of_prim_res_root>0);
            bool root_of_content=(event.mult_of_content_root>0);

            size_type mult_of_resultant  = event.mult_of_prim_res_root;

/*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "Event line for " << index << " "
                                 << root_of_resultant << " "
                                 << root_of_content << " "
                                 << mult_of_resultant << " "
                                 << left_arcs << " " << right_arcs
                                 << std::endl;
#endif
*/

            Status_line_1 ev_line
                = status_line_builder().create_event_line(index,
                                                         x,
                                                         left_arcs,
                                                         right_arcs,
                                                         root_of_resultant,
                                                         root_of_content,
                                                         mult_of_resultant);

            event.stack = ev_line;

#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif

            return ev_line;
        } catch(CGAL::internal::Non_generic_position_exception /* exc */) {
#if CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
          // do nothing
          CGAL_error(); // should not reach here
#else
          switch(this->ptr()->degeneracy_strategy) {
            case(CGAL::EXCEPTION_STRATEGY): {
                throw CGAL::internal::Non_generic_position_exception();
                break;
            }
	    // Feature does not working atm
            case(CGAL::SHEAR_ONLY_AT_IRRATIONAL_STRATEGY): {
	      CGAL_error_msg("Currently not supported");
	      /*
	      if(x.is_rational()) {
                    return create_non_generic_event_at_rational(x,index);
                }
                // FALL INTO NEXT CASE
	      */
            }
	    case(CGAL::SHEAR_STRATEGY): {
                return create_non_generic_event_with_shear(index);
                break;
            }
            default:{
              assert(false); // !!! Never reached
            }
            }
#endif
        }

        // !!! Never reached
        return Status_line_1(x, -1, *this, -1); // dummy

    }

private:

#if !CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
    /*
     * \brief Method to create a status line using shear and backshear
     *
     * Note that this methods creates <b>all</b> event lines of the object
     * at once, and stores them in the object.
     */
    Status_line_1 create_non_generic_event_with_shear(size_type index) const {

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Use sheared technique..." << std::endl;
#endif
        internal::Shear_controller<Integer> shear_controller;
        Integer s(0);
        while(true) {
            try {
                s = shear_controller.get_shear_factor();
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "Trying shear factor "
                                     << s << std::endl;
#endif
                // TODO: Move shear somewhere else
                Self D(kernel(),
                       CGAL::internal::shear
                           (primitive_polynomial_2(),Coefficient(s)),
                       CGAL::EXCEPTION_STRATEGY);
                Shear_transformation< Algebraic_kernel_with_analysis_2 >
                    shear_transformation(kernel());
                shear_transformation.report_sheared_disc_roots
                    (boost::make_transform_iterator(
                             event_coordinates().begin(),
                             typename Rep::Val_functor()),
                     boost::make_transform_iterator(
                             event_coordinates().end(),
                             typename Rep::Val_functor())
                    );

                // Store the sheared curve for later use
                this->ptr()->sheared_curves.insert(std::make_pair(s,D));
                shear_transformation(D,-s,(Self&)*this,false);
                set_vertical_line_components();

                break;
            }
            catch(CGAL::internal::Non_generic_position_exception /* err */) {

                shear_controller.report_failure(s);
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "Bad shear factor, retrying..."
                                     << std::endl;
#endif
            }
        }

        return status_line_at_event(index);
    }

#endif

public:

    /*!
     * \brief Returns the status line for the interval
     * preceeding the <tt>i</tt>th event
     *
     * Returns a status line for a reference x-coordinate of the <tt>i</tt>th
     * interval of the curve. If called multiple times for the same <tt>i</tt>,
     * the same status line is returned.
     */
    Status_line_1 status_line_of_interval(size_type i) const
    {
        CGAL_precondition(i >= 0);
        CGAL_precondition(i <= number_of_status_lines_with_event());

#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_status_line_of_interval(i);
        }
#endif

        Bound b = bound_value_in_interval(i);

        Status_line_1 intermediate_line
            = status_line_at_exact_non_event_x(Algebraic_real_1(b));

        CGAL_postcondition(! intermediate_line.is_event());

        return intermediate_line;
    }


public:

    /*!
     * \brief returns a status line at position \c x
     *
     * If \c x is not an event of the curve, and lies in the <tt>i</tt>th
     * interval, the result is equal to <tt>status_line_of_interval(i)</tt>.
     * Different from <tt>status_line_at_exact_x(x)</tt>
     * the status line \c s returned does not satisft <tt>s.x()==x</tt>.
     * If \c x is an event, and \c perturb is set to \c CGAL::ZERO,
     * the status line for the event is returned. Otherwise, the status line
     * for the left or right neighboring interval is returned, depending
     * on whether \c perturb is set to \c CGAL::NEGATIVE or \c CGAL::POSITIVE.
     * If \c x is not an event, \c perturb has no effect.
     */
    Status_line_1 status_line_for_x(Algebraic_real_1 x,
                                    CGAL::Sign perturb = CGAL::ZERO) const
    {
#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_status_line_for_x(x,perturb);
        }
#endif

        size_type i;
        bool is_evt;
        x_to_index(x, i, is_evt);
        if(is_evt) {
            if(perturb == CGAL::ZERO)
                return status_line_at_event(i);
            if(perturb == CGAL::POSITIVE)
                i++;
        }
        return status_line_of_interval(i);
    }


private:

    /*
     * \brief Creates an intermediate line at position \c ar.
     *
     * It is required that none of the following situations occurs at position
     * <tt>ar</tt>: singularity, vertical tangent line, vertical asymptote.\n
     * Otherwise, the method might run into an infinite loop.
     *
     * \param index if set to -1, the interval containing \c ar is computed
     * within the method, and the index of the status line is set accordingly.
     */
    Status_line_1
    create_status_line_at_non_event(Algebraic_real_1 ar, int index = -1)
        const {

      if (index==-1) {
        bool event;
        x_to_index(ar,index,event);
        CGAL_assertion(!event);
      }
      CGAL_assertion(index>=0);

      // TODO .. delay creation of refinement object
      // especially when ar is rational

      Status_line_1 status_line =
        status_line_builder().create_non_event_line(index, ar);

      // TODO debug output?

      CGAL_assertion(! status_line.is_event());

      return status_line;
    }

private:

   /*
    * \brief Returns an Status_line_builder instance
    *
    * Note: So far, a new instance is created each time the function is called
    */
    Status_line_builder& status_line_builder() const {
      if (!this->ptr()->_m_status_line_builder) {
#if CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
        this->ptr()->_m_status_line_builder = new Status_line_builder(kernel(), *this);
#else
        this->ptr()->_m_status_line_builder = new Status_line_builder(kernel(), *this, primitive_polynomial_2());
#endif
      }
      return *(this->ptr()->_m_status_line_builder);
     }

public:

    /*!
     * \brief Number of arcs over the given interval
     *
     * Shortcut for <tt>status_line_of_interval(i).number_of_events()</tt>
     */
    size_type arcs_over_interval(size_type i) const {
        CGAL_precondition(has_defining_polynomial());
#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_arcs_over_interval(i);
        }
#endif
        CGAL_assertion_code(
                size_type n
                    = static_cast<size_type>(intermediate_values().size());
        );
        CGAL_precondition(i>=0);
        CGAL_precondition(i<=n);
        return status_line_of_interval(i).number_of_events();
    }

public:

    /*!
     * \brief Rational number in the <tt>i</tt>th interval between events
     *
     * The result of this method is taken as the reference x-coordinate
     * for the status lines of intervals.
     */
    Bound bound_value_in_interval(size_type i) const {
#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_bound_value_in_interval(i);
        }
#endif
        CGAL_assertion(i>=0);
        CGAL_assertion(i < static_cast<size_type>
                           (intermediate_values().size()));
        if(! intermediate_values()[i]) {
          // Create it
            if(event_coordinates().size()==0) {
                CGAL_assertion(i==0);
                intermediate_values()[0]=Bound(0);
            } else {
                if(i==0) {
                    intermediate_values()[i]
		      = bound_left_of(kernel(),event_coordinates()[i].val);
                } else if(i == static_cast<size_type>
                              (event_coordinates().size())) {
                    intermediate_values()[i]
                        = bound_right_of
		      (kernel(),event_coordinates()[i-1].val);

                } else {
                    intermediate_values()[i]
		      = kernel()->bound_between_1_object()
		      (event_coordinates()[i-1].val,
		       event_coordinates()[i].val);
                }
            }
        }
        return intermediate_values()[i].get();
    }


public:

    /*!
     * Returns the content of the defining polynomial
     *
     * The content is the gcd of its coefficients (the polynomial is considered
     * as polynomial in \c y)
     */
    Polynomial_1 content() const {
#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_content();
        }
#endif
        if(! this->ptr()->content) {
            compute_content_and_primitive_part();
        }
        return this->ptr()->content.get();
    }

public:

    /*!
     * Returns the primitive part of the defining polynomial
     *
     * The primitive part of \c f is the \c f divided by its content.
     */
    Polynomial_2 primitive_polynomial_2() const {
#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_primitive_polynomial_2();
        }
#endif
        if(! this->ptr()->f_primitive) {
            compute_content_and_primitive_part();
        }
        return this->ptr()->f_primitive.get();
    }

    Algebraic_kernel_with_analysis_2* kernel() const {
        return this->ptr()->_m_kernel;
    }

private:


    // computes and sets the content and the primitive part for the curve
    void compute_content_and_primitive_part() const {

        CGAL_assertion(has_defining_polynomial());
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Computing the content..." << std::flush;
#endif
        this->ptr()->content
            = typename CGAL::Polynomial_traits_d< Polynomial_2 >::
                Univariate_content_up_to_constant_factor()( polynomial_2() );
        if(CGAL::degree(content())==0) {
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "no vertical lines as components"
                                 << std::endl;
#endif
            this->ptr()->f_primitive=polynomial_2();
        }
        else {
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "non-trivial content found" << std::endl;
#endif
            // Content must be square free, because the curve is square free
            CGAL_assertion( typename CGAL::Polynomial_traits_d< Polynomial_1 >
                            ::Is_square_free()(content()));
            this->ptr()->f_primitive=polynomial_2() / content();

        }

    }

#if !CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
private:

    //! Returns the Sturm-Habicht sequence of the primitive part of f
    std::vector<Polynomial_2>& sturm_habicht_of_primitive() const
    throw(internal::Zero_resultant_exception<Polynomial_2>) {
        if(! this->ptr()->sturm_habicht_of_primitive) {
            compute_sturm_habicht_of_primitive();
        }
        return this->ptr()->sturm_habicht_of_primitive.get();
    }

public:

    /*!
     * \brief Returns the <tt>i</tt>th Sturm-Habicht polynomial
     * of the primitive part of the defining polynomial
     */
    Polynomial_2 sturm_habicht_of_primitive(size_type i) const
      throw(internal::Zero_resultant_exception<Polynomial_2>) {
        CGAL_assertion(i>=0);
        CGAL_assertion(i < static_cast<size_type>
                       (sturm_habicht_of_primitive().size()));
        return sturm_habicht_of_primitive()[i];
    }

public:

    /*!
     * \brief Returns the <tt>i</tt>th principal Sturm-Habicht coefficient
     * of the primitive part of the defining polynomial
     */
    Polynomial_1 principal_sturm_habicht_of_primitive(size_type i) const
        throw(internal::Zero_resultant_exception<Polynomial_2>) {
        CGAL_assertion(i>=0);
        CGAL_assertion(i < static_cast<size_type>
                       (sturm_habicht_of_primitive().size()));

        CGAL_assertion(CGAL::degree(sturm_habicht_of_primitive()[i])<=i);
        if(CGAL::degree(sturm_habicht_of_primitive()[i]) < i) {
            return Polynomial_1(0);
        } // else:
        return sturm_habicht_of_primitive()[i][i];
    }

public:

    /*!
     * \brief Returns the <tt>i</tt>th coprincipal Sturm-Habicht coefficient
     * of the primitive part of the defining polynomial
     *
     * The coprincipal Sturm-Habicht coefficient is the coefficient
     * of <tt>y^{i-1}</tt> of the <tt>i</tt>th Sturm-Habicht polynomial
     */
    Polynomial_1 coprincipal_sturm_habicht_of_primitive(size_type i) const
        throw(internal::Zero_resultant_exception<Polynomial_2>) {
      CGAL_assertion(i>=1);
      CGAL_assertion(i < static_cast<size_type>
                       (sturm_habicht_of_primitive().size()));
        CGAL_assertion(CGAL::degree(sturm_habicht_of_primitive()[i])<=i);
        if(CGAL::degree(sturm_habicht_of_primitive()[i]) < i-1) {
            return Polynomial_1(0);
        } // else:
        return sturm_habicht_of_primitive()[i][i-1];
    }

public:

    /*!
     * \brief Returns an iterator to the principal Sturm-Habicht coefficients,
     * starting with the <tt>0</tt>th one (the resultant)
     */
    Principal_sturm_habicht_iterator principal_sturm_habicht_begin() const {
        return boost::make_transform_iterator
            (boost::counting_iterator<size_type>(0),
             Stha_functor(this));
    }

    //! Returns an iterator to the end of principal Sturm-Habicht coefficients
    Principal_sturm_habicht_iterator principal_sturm_habicht_end() const {
        return boost::make_transform_iterator
            (boost::counting_iterator<size_type>
                 (sturm_habicht_of_primitive().size()),
             Stha_functor(this));
    }

private:

    // Internal method to compute the Sturm-Habicht sequence
    void compute_sturm_habicht_of_primitive() const
        throw(internal::Zero_resultant_exception<Polynomial_2>) {

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Compute Sturm-Habicht.." << std::flush;
#endif
        std::vector<Polynomial_2> stha;

        // Fix a problem for constant primitive part.
        // In this case, the St.-Ha. sequence is never needed
        if(CGAL::degree(primitive_polynomial_2()) == 0) {
            // Set the resultant
            stha.push_back(primitive_polynomial_2());
        } else {

#if CGAL_ACK_USE_BEZOUT_MATRIX_FOR_SUBRESULTANTS
#if CGAL_AK_D_SHOW_COMPILE_OPTIONS_AS_WARNING
#warning Curve_analysis_2: Uses Bezout Matrix for subresultants
#endif
            CGAL::internal::bezout_polynomial_subresultants<Polynomial_traits_2>
                (primitive_polynomial_2(),
                 CGAL::canonicalize(CGAL::differentiate(primitive_polynomial_2())),
                 std::back_inserter(stha));
            stha.push_back(primitive_polynomial_2());
            size_type p = CGAL::degree(primitive_polynomial_2());
            CGAL_assertion(static_cast<size_type>(stha.size()) == p+1);
            for(size_type i=0;i<p; i++) {
                if((p-i)%4==0 || (p-i)%4==1) {
                    stha[i] = stha[i];
                } else {
                    stha[i] = -stha[i];
                }
            }

#else
            typename Polynomial_traits_2::Sturm_habicht_sequence()
                (primitive_polynomial_2(),std::back_inserter(stha));
#endif
        }
        // Also set the resultant, if not yet set
        if(! this->ptr()->resultant_of_primitive_and_derivative_y) {
            this->ptr()->resultant_of_primitive_and_derivative_y = stha[0][0];
            if(this->ptr()->resultant_of_primitive_and_derivative_y.
                   get().is_zero()) {
                throw internal::Zero_resultant_exception<Polynomial_2>
                    (polynomial_2());
            }
        }

        this->ptr()->sturm_habicht_of_primitive = stha;
        CGAL_assertion(CGAL::canonicalize
		       (resultant_of_primitive_and_derivative_y()) ==
                       CGAL::canonicalize
		       (principal_sturm_habicht_of_primitive(0)));
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
    }
#endif


private:

    //! Returns the resultant of the primitive part of f and its y-derivative
    Polynomial_1 resultant_of_primitive_and_derivative_y() const
        throw(internal::Zero_resultant_exception<Polynomial_2>) {
        if(! this->ptr()->resultant_of_primitive_and_derivative_y) {
            compute_resultant_of_primitive_and_derivative_y();
        }
        return this->ptr()->resultant_of_primitive_and_derivative_y.get();
    }

private:

    //! Returns the resultant of the primitive part of f with its x-derivative
    Polynomial_1 resultant_of_primitive_and_derivative_x() const
        throw(internal::Zero_resultant_exception<Polynomial_2>) {
        if(! this->ptr()->resultant_of_primitive_and_derivative_x) {
            compute_resultant_of_primitive_and_derivative_x();
        }
        return this->ptr()->resultant_of_primitive_and_derivative_x.get();
    }

#if CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
 public:
     //! Returns the resultant of the primitive part of f and its y-derivative
     Polynomial_1 resultant_of_coprime_derivative_x_and_derivative_y() const
        throw(internal::Zero_resultant_exception<Polynomial_2>) {
        if(! this->ptr()->resultant_of_coprime_derivative_x_and_derivative_y) {
            this->compute_resultant_of_coprime_derivative_x_and_derivative_y();
        }
        return this->ptr()->resultant_of_coprime_derivative_x_and_derivative_y.get();
    }

 private:

    template< class Value_ >
    struct First_equal {

      typedef Value_ Value;

      First_equal(const Value& value) :
        _m_value(value) {
      }

      template < class Second_ >
      bool operator()(const std::pair< Value, Second_ >& pair) {
        return CGAL::compare(pair.first, _m_value) == CGAL::EQUAL;
      }

    private:
      Value _m_value;

    };


 public:

    Multiplicity_type multiplicity_as_root_of_resultant_fx_fy(const Algebraic_real_1& alpha) const {

      //Bisolve_out("mult fxfy alpha: " << alpha << std::endl)
      //Bisolve_out("mult fxfy dalpha: "<< CGAL::to_double(alpha) << std::endl)
      //Bisolve_out("mult fxfy res: " << resultant_of_coprime_derivative_x_and_derivative_y() << std::endl)

      if (CGAL::is_zero(resultant_of_coprime_derivative_x_and_derivative_y())) {
        return -1; // resultant vanishes
      }

#if CGAL_ACK_FACTORIZE_UNI_POLYNOMIALS

      Bisolve_telemetry_code(t_mult_ntl_div.start();)
      NTL::ZZX p, q;

      CGAL::internal::poly2ntl(resultant_of_coprime_derivative_x_and_derivative_y(), q);
      CGAL::internal::poly2ntl(alpha.polynomial(), p);

      // remark: this test only works as NTL has been used to factor R = res(f,f_y,y) and thus
      //         alpha.polynomial() is one of the irreduciable factors
      Multiplicity_type m = 0;
      while (NTL::divide(q, q, p)) {
        m++;
      }

      Bisolve_telemetry_code(t_mult_ntl_div.stop();)
      return m;

#else
      if (!this->ptr()->resultant_of_coprime_derivative_x_and_derivative_y_factors) {
        // sqf of Q=res(f_x,f_y)
        Bisolve_telemetry_code(t_mult_sqfs.start();)
          typename Polynomial_traits_1::Square_free_factorize_up_to_constant_factor sqfactorize;
        std::vector< std::pair< Polynomial_1, Multiplicity_type > > factors;
        sqfactorize(resultant_of_coprime_derivative_x_and_derivative_y(), std::back_inserter(factors));
        // sort according to mult
        std::sort(factors.begin(), factors.begin(),
                  boost::bind(&std::pair< Polynomial_1, Multiplicity_type >::second, _1) <
                  boost::bind(&std::pair< Polynomial_1, Multiplicity_type >::second, _2));
        this->ptr()->resultant_of_coprime_derivative_x_and_derivative_y_factors = factors;
        Bisolve_telemetry_code(t_mult_sqfs.stop();)

/*    std::cerr << "alpha: " << alpha << "; res_der_xy: " << resultant_of_coprime_derivative_x_and_derivative_y() << "\n";*/

        // create sequence
        Polynomial_1 qi(1);
        size_t j = 0;
        std::vector< Polynomial_1 > qis;
        int deg = CGAL::degree(resultant_of_coprime_derivative_x_and_derivative_y());
        if(deg == 0) {
            qis.push_back(Polynomial_1(1));
        }

        for (int i = 1; i <= deg; i++) {
          if ((*this->ptr()->resultant_of_coprime_derivative_x_and_derivative_y_factors)[j].second == i) {

            if (i == 1) {
              qis.push_back((*this->ptr()->resultant_of_coprime_derivative_x_and_derivative_y_factors)[j].first);
            } else {
              qis.push_back(qis.back() *
                            (*this->ptr()->resultant_of_coprime_derivative_x_and_derivative_y_factors)[j].first);
            }
            j++;
            if (j == this->ptr()->resultant_of_coprime_derivative_x_and_derivative_y_factors->size()) {
              break;
            }
          } else {
            if (i == 1) {
              qis.push_back(Polynomial_1(1));
            } else {
              qis.push_back(qis.back());
            }
          }
        }
        this->ptr()->resultant_of_coprime_derivative_x_and_derivative_y_sequence = qis;
      }

      // map for gcds
      typedef std::map< std::pair< Polynomial_1, Polynomial_1 >, Polynomial_1 > Gcd_map;
      static Gcd_map gcd_map;

      int min, max = min = this->ptr()->resultant_of_coprime_derivative_x_and_derivative_y_sequence->size();

      bool first = true;


      while (true) {

        int mid = min + (max - min) / 2;
        Bisolve_telemetry_code(t_mult_sqfs.start();)

        std::pair< Polynomial_1, Polynomial_1 > key =
          std::make_pair(alpha.polynomial(),
                         (*this->ptr()->resultant_of_coprime_derivative_x_and_derivative_y_sequence)[mid-1]);
#ifdef CGAL_WEBXACLI_MULTITHREADED
    pthread_mutex_lock(&CGAL_static_mtx);
#endif

        typename Gcd_map::iterator git = gcd_map.find(key);
        if (git == gcd_map.end()) {
            git = gcd_map.insert(git, std::make_pair(key,
              CGAL::gcd_up_to_constant_factor(key.first, key.second)));
//           std::cout << "GCD CACHE MISS" << std::endl;
        } else {
//           std::cout << "GCD CACHE HIT" << std::endl;
        }
#ifdef CGAL_WEBXACLI_MULTITHREADED
    pthread_mutex_unlock(&CGAL_static_mtx);
#endif

        Polynomial_1 s_mid = git->second;
        Bisolve_telemetry_code(t_mult_sqfs.stop();)
        Bisolve_telemetry_code(t_mult_sign.start();)
        CGAL::Sign sl = kernel()->sign_at_1_object()(s_mid, alpha.low());
        CGAL::Sign sh = kernel()->sign_at_1_object()(s_mid, alpha.high());
        Bisolve_telemetry_code(t_mult_sign.stop();)

        CGAL::Sign sc = sl*sh;

        if (first) {
          if (sc > CGAL::ZERO) {
            // if Q_kmax has no sign change, return 0;
            // std::cout << "EARLY EXIT" << std::endl;
            return 0;
          } else {
            // start binary search
            min = 1;
            first = false;
          }
        }

        if (min == max) {
          if (sc <= CGAL::ZERO) {
            return min;
          } else {
            return min+1;
          }
        }

        // else
        if (sc > CGAL::ZERO) {
          // no sign change at mid, search in right interval
          min = mid + 1;
        } else {
          // sign change at mid, search in left interval
          max = mid;
        }

      }

      // else
      return 0;
#endif

    }
#endif

private:
    // Computes <tt>res_y(f,f_y)</tt>, where \c f is the defining polynomial
    void compute_resultant_of_primitive_and_derivative_y() const
        throw(internal::Zero_resultant_exception<Polynomial_2>) {

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Compute resultant.." << std::flush;
#endif

        CGAL_assertion(has_defining_polynomial());

#if CGAL_ACK_RESULTANT_FIRST_STRATEGY
#ifndef CGAL_ACK_RESULTANT_FIRST_STRATEGY_DEGREE_THRESHOLD
        bool speed_up = true;
#else
        bool speed_up=CGAL::degree(polynomial_2()) >=
            CGAL_ACK_RESULTANT_FIRST_STRATEGY_DEGREE_THRESHOLD;
#endif
#else
        bool speed_up=false;
#endif

        if(CGAL::degree(polynomial_2()) == 0) {
	    this->ptr()->resultant_of_primitive_and_derivative_y
                = Polynomial_1(1);
        } else {

#if !CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
            Bisolve_telemetry_code(t_res[0].start();)
            Bisolve_telemetry_code(t_res[1].start();)
            if(! speed_up) {
              CGAL_assertion(false);
              // Compute resultant using the Sturm-Habicht sequence
                this->ptr()->resultant_of_primitive_and_derivative_y
                    = principal_sturm_habicht_of_primitive(0);
            } else {
                typename Polynomial_traits_2::Differentiate diff;
                this->ptr()->resultant_of_primitive_and_derivative_y =
                    CGAL::resultant(primitive_polynomial_2(),
                    CGAL::canonicalize(diff(primitive_polynomial_2(),1)));
            }
            Bisolve_telemetry_code(t_res[1].stop();)
            Bisolve_telemetry_code(t_res[0].stop();)
#else
            this->ptr()->resultant_of_primitive_and_derivative_y =
            status_line_builder().bi_diff_vanish().res_y_f_fy();
#endif

        }

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif

        if(resultant_of_primitive_and_derivative_y().is_zero()) {
            throw internal::Zero_resultant_exception<Polynomial_2>
                (polynomial_2());
        }
    }

    // Computes <tt>res_y(f,f_x)</tt>, where \c f is the defining polynomial
    void compute_resultant_of_primitive_and_derivative_x() const
        throw(internal::Zero_resultant_exception<Polynomial_2>) {

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Compute x-resultant.." << std::flush;
#endif

        CGAL_assertion(has_defining_polynomial());

        // Transpose the polynomial
        Polynomial_2 f_yx = typename Polynomial_traits_2::Swap()
            (polynomial_2(),0,1);

        if( CGAL::degree(f_yx) == 0 ) {
            // Polynomial only consists of horizontal lines
            // primitive resultant is set to 1
            this->ptr()->resultant_of_primitive_and_derivative_x
                = Polynomial_1(1);
        } else {

            Polynomial_2 f_yx_primitive;

            Polynomial_1 content_yx
                = typename CGAL::Polynomial_traits_d< Polynomial_2 >::
                    Univariate_content_up_to_constant_factor()( f_yx );
            if(CGAL::degree(content_yx)==0) {
                f_yx_primitive=f_yx;
            }
            else {
                CGAL_assertion
                    (typename CGAL::Polynomial_traits_d< Polynomial_1 >::
                         Is_square_free()(content_yx));
                f_yx_primitive=f_yx / content_yx;

            }

            Bisolve_telemetry_code(t_res[0].stop();)
            Bisolve_telemetry_code(t_res[2].start();)

            this->ptr()->resultant_of_primitive_and_derivative_x =
            CGAL::resultant
            (typename Polynomial_traits_2::Swap() (f_yx_primitive,0,1),
             typename Polynomial_traits_2::Swap()
             (CGAL::canonicalize(CGAL::differentiate(f_yx_primitive),0,1)));
            Bisolve_telemetry_code(t_res[2].stop();)
            Bisolve_telemetry_code(t_res[0].stop();)
        }

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif

        if(resultant_of_primitive_and_derivative_x().is_zero()) {
            throw internal::Zero_resultant_exception<Polynomial_2>
                (polynomial_2());
        }
    }


#if CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
private:
    // Computes <tt>res_y(f,f_y)</tt>, where \c f is the defining polynomial
    void compute_resultant_of_coprime_derivative_x_and_derivative_y() const
        throw(internal::Zero_resultant_exception<Polynomial_2>) {

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Compute resultant(f_x, f_y) ... " << std::flush;
#endif
        CGAL_assertion(has_defining_polynomial());

        Polynomial_2 fx = CGAL::canonicalize(
            CGAL::differentiate(primitive_polynomial_2(), 0));
        Polynomial_2 fy = CGAL::canonicalize(
            CGAL::differentiate(primitive_polynomial_2(), 1));

        Polynomial_1 gres(0);

          Bisolve_telemetry_code(t_mult_xy_gcd.start();)
          // divide out common factor
          Polynomial_2 h = CGAL::gcd(fx, fy);
          if (CGAL::degree(h) > 0) {
            fx = fx/h;
            fy = fy/h;
          }
          Bisolve_telemetry_code(t_mult_xy_gcd.stop();)

        //std::cout << "fx: " << fx << std::endl;
        //std::cout << "fy: " << fy << std::endl;

        Bisolve_telemetry_code(t_mult_res.start();)
        this->ptr()->resultant_of_coprime_derivative_x_and_derivative_y = CGAL::resultant(fx, fy);
        Bisolve_telemetry_code(t_mult_res.stop();)

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
    }
#endif


private:

    // Returns the critical event coordinates
    std::vector<Event_coordinate_1>& event_coordinates() const
        throw(internal::Zero_resultant_exception<Polynomial_2>) {
        if(! this->ptr()->event_coordinates) {
            compute_event_coordinates();
        }
        return this->ptr()->event_coordinates.get();
    }

private:

    // Returns the intermediate values for intervals between events
    std::vector<boost::optional<Bound> >& intermediate_values() const
        throw(internal::Zero_resultant_exception<Polynomial_2>) {

        if(! this->ptr()->intermediate_values) {
            // This is created during event_coordiantes()
            event_coordinates();
            CGAL_assertion(this->ptr()->intermediate_values);
        }
        return this->ptr()->intermediate_values.get();
    }


private:

    /*
     * \brief Computes the event coordinates of the curve.
     *
     * This function computes the content of the defining polynomial,
     * and the roots of its discriminant. These two sets form the critical
     * x-coordinates of the curve.
     */
    void compute_event_coordinates() const
      throw(internal::Zero_resultant_exception<Polynomial_2>)  {

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "compute events..." << std::flush;
#endif

        // TODO t_events.start();

        Solve_1 solve_1;

        std::vector<std::pair<Algebraic_real_1,size_type> > content_pairs;
        std::vector<Algebraic_real_1> content_roots;
        std::vector<size_type> content_mults;
        solve_1(content(),
                std::back_inserter(content_pairs));

        for(int i=0; i < static_cast<int>(content_pairs.size()); i++ ) {
            content_roots.push_back(content_pairs[i].first);
            content_mults.push_back(content_pairs[i].second);
        }

        // Set the vertical_line_components flag as side effect
        this->ptr()->has_vertical_component = (content_roots.size() > 0);

#if CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
        std::vector<std::pair<Algebraic_real_1,size_type> > res_pairs;
        std::vector<Algebraic_real_1> res_roots;
        std::vector<size_type> res_mults;

        status_line_builder().bi_diff_vanish().copy_solutions_in_x(std::back_inserter(res_pairs));

        for(int i=0; i < static_cast<int>(res_pairs.size()); i++ ) {
            res_roots.push_back(res_pairs[i].first);
            res_mults.push_back(res_pairs[i].second);
        }
#else
        Bisolve_telemetry_code(t_solve1[0].start();)
        Bisolve_telemetry_code(t_solve1[1].start();)
        std::vector<std::pair<Algebraic_real_1,size_type> > res_pairs;
        std::vector<Algebraic_real_1> res_roots;
        std::vector<size_type> res_mults;
        Polynomial_1 R = resultant_of_primitive_and_derivative_y();
        solve_1(R,std::back_inserter(res_pairs));

        for(int i=0; i < static_cast<int>(res_pairs.size()); i++ ) {
            res_roots.push_back(res_pairs[i].first);
            res_mults.push_back(res_pairs[i].second);
        }
        Bisolve_telemetry_code(t_solve1[1].stop();)
        Bisolve_telemetry_code(t_solve1[0].stop();)
#endif

        std::vector<std::pair<Algebraic_real_1,size_type> > lcoeff_pairs;
        std::vector<Algebraic_real_1> lcoeff_roots;
        std::vector<size_type> lcoeff_mults;
        solve_1(CGAL::leading_coefficient(primitive_polynomial_2()),
                std::back_inserter(lcoeff_pairs));

        for(int i=0; i < static_cast<int>(lcoeff_pairs.size()); i++ ) {
            lcoeff_roots.push_back(lcoeff_pairs[i].first);
            lcoeff_mults.push_back(lcoeff_pairs[i].second);
        }


        //Now, merge the vertical line positions with the resultant roots
        typename
            CGAL::Real_embeddable_traits<Algebraic_real_1>::Compare compare;

        std::vector<Algebraic_real_1> event_values;
        std::vector<CGAL::internal::Three_valued> event_values_info;

        CGAL::internal::set_union_with_source
            (res_roots.begin(),
             res_roots.end(),
             content_roots.begin(),
             content_roots.end(),
             std::back_inserter(event_values),
             std::back_inserter(event_values_info),
             compare);

        // Now, build the Event_coordinate_1 entries
        // for each element of event_values
        size_type curr_res_index = 0, curr_content_index = 0,
            curr_lcoeff_index = 0;
        std::vector<Event_coordinate_1> event_coordinate_vector;

        for(size_type i = 0;
            i < static_cast<size_type>(event_values.size());
            i++ ) {

            Event_coordinate_1 curr_event;
            curr_event.val = event_values[i];
            switch(event_values_info[i]) {

            case(CGAL::internal::ROOT_OF_FIRST_SET): {
                curr_event.index_of_prim_res_root = curr_res_index;
                CGAL_expensive_assertion(res_roots[curr_res_index] ==
                                         event_values[i]);
                curr_event.mult_of_prim_res_root
                    = res_mults[curr_res_index];
                curr_res_index++;
                if(curr_lcoeff_index <
                   static_cast<size_type>(lcoeff_roots.size()) &&
                   event_values[i]==lcoeff_roots[curr_lcoeff_index]) {
                    // We have a root of the leading coefficient
                    // of the primitve polynomial
                    curr_event.index_of_prim_lcoeff_root = curr_lcoeff_index;
                    curr_event.mult_of_prim_lcoeff_root
                        = lcoeff_mults[curr_lcoeff_index];
                    curr_lcoeff_index++;
                } else {
                    curr_event.index_of_prim_lcoeff_root = -1;
                    curr_event.mult_of_prim_lcoeff_root = 0;
                }

                curr_event.index_of_content_root = -1;
                curr_event.mult_of_content_root = 0;
                break;
            }
            case(CGAL::internal::ROOT_OF_SECOND_SET): {
                curr_event.index_of_content_root = curr_content_index;
                CGAL_expensive_assertion(content_roots[curr_content_index] ==
                                         event_values[i]);
                curr_event.mult_of_content_root
                    = content_mults[curr_content_index];
                curr_content_index++;
                curr_event.index_of_prim_res_root = -1;
                curr_event.mult_of_prim_res_root = 0;
                CGAL_expensive_assertion(event_values[i]!=
                                         lcoeff_roots[curr_lcoeff_index]);
                curr_event.index_of_prim_lcoeff_root = -1;
                curr_event.mult_of_prim_lcoeff_root = 0;
                break;
            }
            case(CGAL::internal::ROOT_OF_BOTH_SETS): {
                curr_event.index_of_prim_res_root = curr_res_index;
                CGAL_expensive_assertion(res_roots[curr_res_index] ==
                                         event_values[i]);
                curr_event.mult_of_prim_res_root
                    = res_mults[curr_res_index];
                curr_res_index++;
                if(curr_lcoeff_index <
                   static_cast<size_type>(lcoeff_roots.size()) &&
                   event_values[i]==lcoeff_roots[curr_lcoeff_index]) {
                    // We have a root of the leading coefficient
                    // of the primitve polynomial
                    curr_event.index_of_prim_lcoeff_root = curr_lcoeff_index;
                    curr_event.mult_of_prim_lcoeff_root
                        = lcoeff_mults[curr_lcoeff_index];
                    curr_lcoeff_index++;
                } else {
                    curr_event.index_of_prim_lcoeff_root = -1;
                    curr_event.mult_of_prim_lcoeff_root = 0;
                }
                curr_event.index_of_content_root = curr_content_index;
                CGAL_expensive_assertion(content_roots[curr_content_index] ==
                                         event_values[i]);
                curr_event.mult_of_content_root
                    = content_mults[curr_content_index];
                curr_content_index++;
                break;
            }
            } // of switch
            /*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "Constructed event_coordinate: "
                                 << CGAL::to_double(curr_event.val) << " "
                                 << "\nmult_of_prim_res_root : "
                                 << curr_event.mult_of_prim_res_root
                                 << "\nindex_of_prim_res_root : "
                                 << curr_event.index_of_prim_res_root
                                 << "\nmult_of_content_root : "
                                 << curr_event.mult_of_content_root
                                 << "\nindex_of_content_root : "
                                 << curr_event.index_of_content_root
                                 << "\nmult_of_lcoeff_root : "
                                 << curr_event.mult_of_prim_lcoeff_root
                                 << "\nindex_of_lcoeff_root : "
                                 << curr_event.index_of_prim_lcoeff_root
                                 << std::endl;
#endif
            */
            event_coordinate_vector.push_back(curr_event);
        }


        CGAL_assertion(curr_lcoeff_index ==
                       static_cast<size_type>(lcoeff_roots.size()));
        CGAL_assertion(curr_res_index ==
                       static_cast<size_type>(res_roots.size()));
        CGAL_assertion(curr_content_index ==
                       static_cast<size_type>(content_roots.size()));

        this->ptr()->intermediate_values
            = std::vector<boost::optional<Bound> >
            (event_coordinate_vector.size()+1);
        this->ptr()->event_coordinates = event_coordinate_vector;

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif

    }

public:

#if !CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
    /*!
     * \brief Returns a \c Curve_analysis_2 object for a sheared curve.
     *
     * The shear factor is given by the integer \c s.
     * This functions only shears the primitive part of the defining equation.
     * Internal caching is used to avoid repeated shears.
     *
     * \todo The sheared curves are not inserted into the curve_cache
     * of the Algebraic_curve_kernel_2 yet.
     */
    Self& shear_primitive_part(Integer s) const
        throw(CGAL::internal::Non_generic_position_exception)
    {
        CGAL_assertion(s!=0);
#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_shear_primitive_part();
        }
#endif
        if(this->ptr()->bad_shears.find(s) !=
           this->ptr()->bad_shears.end()) {
            throw CGAL::internal::Non_generic_position_exception();
        }
        typedef typename std::map<Integer,Self>::iterator
            Map_iterator;
        Map_iterator it = this->ptr()->sheared_curves.find(s);
        if(it != this->ptr()->sheared_curves.end()) {
            return it->second;
        }
        try {
            Shear_transformation<Algebraic_kernel_with_analysis_2>
                shear_transformation(kernel());
            Self D=shear_transformation((Self&)*this, s);
            std::pair<Map_iterator,bool> insertion =
                this->ptr()->sheared_curves.insert(std::make_pair(s,D));
            CGAL_assertion(insertion.second);
            return insertion.first->second;
        }
        catch(CGAL::internal::Non_generic_position_exception /* err */) {
            this->ptr()->bad_shears.insert(s);
            throw CGAL::internal::Non_generic_position_exception();
        }
    }

public:

    //! Iterator for sheared curves
    typename std::map<Coefficient,Self>::const_iterator shear_begin() {
        return this->ptr()->sheared_curves.begin();
    }

    //! Iterator for sheared curves
    typename std::map<Coefficient,Self>::const_iterator shear_end() {
        return this->ptr()->sheared_curves.end();
    }

#endif

private:

    // Sets the flag for vertical lines in all status lines that need it
    void set_vertical_line_components() const {
        for(size_type i = 0;
            i < static_cast<size_type>(event_coordinates().size());
            i++ ) {

            if(event_coordinates()[i].mult_of_content_root > 0) {
                status_line_at_event(i)._set_v_line();
            }
        }

    }


public:

    /*!
     * \brief Increases the precision of all status lines
     *
     * For each status line at an event and each status line that represents
     * an interval, all y-coordinates are approximated such that their
     * isolating interval has absolute size smaller then \c precision.
     */
    void refine_all(Bound precision) {

#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_refine_all(precision);
        }
#endif

        for(size_type i=0;
            i<static_cast<size_type>(event_coordinates().size());
            i++) {
        /*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << i << ": " << std::flush;
#endif
        */
            Status_line_1& el = status_line_at_event(i);

            for(size_type j=0;j<el.number_of_events();j++) {
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << j << " " << std::flush;
#endif
*/
                el.refine_to(j,precision);
            }
        }
        for(size_type i=0;
            i<static_cast<size_type>(intermediate_values().size());
            i++) {
            Status_line_1 il = status_line_of_interval(i);
            for(size_type j=0;j<il.number_of_events();j++) {
                il.refine_to(j,precision);
            }
        }
    }

public:

    //! \brief Iterator for the status lines at events
    Event_line_iterator event_begin() const {
        return boost::make_transform_iterator
            (boost::counting_iterator<size_type>(0),
             Event_functor(this));
    }

    //! \brief Iterator for the status lines at events
    Event_line_iterator event_end() const {
        return boost::make_transform_iterator
            (boost::counting_iterator<size_type>
             (number_of_status_lines_with_event()),
             Event_functor(this));
    }

public:

    //! \brief Iterator for the status lines for intervals
    Intermediate_line_iterator intermediate_begin() const {
        return boost::make_transform_iterator
            (boost::counting_iterator<size_type>(0),
             Intermediate_functor(this));
    }

    //! \brief Iterator for the status lines for intervals
    Intermediate_line_iterator intermediate_end() const {
        return boost::make_transform_iterator
            (boost::counting_iterator<size_type>(intermediate_values().size()),
             Intermediate_functor(this));
    }

public:

    /*!
     * \brief Returns the limit an infinite arc converges to
     *
     * \pre <tt>loc==CGAL::ARR_LEFT_BOUNDARY ||
     *          loc==CGAL::ARR_RIGHT_BOUNDARY</tt>
     *
     * This method returns for the <tt>arcno</tt>th arc that goes to -infinity
     * or +infinity (depending on \c loc) the y-coordinate it converges to.
     * Possible values are either a \c Algebraic_real_1 object, or one of the
     * values \c CGAL::ARR_TOP_BOUNDARY, \c CGAL::ARR_BOTTOM_BOUNDARY
     * that denote that the arc is unbounded in y-direction.
     * The result is wrapped into a \c CGAL::Object object.
     */
    Asymptote_y asymptotic_value_of_arc(CGAL::Arr_parameter_space loc,
                                        size_type arcno) const {

        CGAL_precondition(loc == CGAL::ARR_LEFT_BOUNDARY ||
                          loc == CGAL::ARR_RIGHT_BOUNDARY);

#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX
        if(CGAL::degree(polynomial_2(),1)==2) {
            return this->conic_asymptotic_value_of_arc(loc,arcno);
        }
#endif

        if(loc == CGAL::ARR_LEFT_BOUNDARY) {

            if(! this->ptr()->horizontal_asymptotes_left) {
                compute_horizontal_asymptotes();
            }
            std::vector<Asymptote_y>& asym_info
                = this->ptr()->horizontal_asymptotes_left.get();
            CGAL_precondition(arcno>=0);
            CGAL_precondition(arcno<static_cast<size_type>(asym_info.size()));
            return asym_info[arcno];
        } // else loc == CGAL::ARR_RIGHT_BOUNDARY

        if(! this->ptr()->horizontal_asymptotes_right) {
            compute_horizontal_asymptotes();
        }
        std::vector<Asymptote_y>& asym_info
            = this->ptr()->horizontal_asymptotes_right.get();
        CGAL_precondition(arcno>=0);
        CGAL_precondition(arcno<static_cast<size_type>(asym_info.size()));
        return asym_info[arcno];

    }


private:

    // Internal method to compute horizontal asymptotes
    void compute_horizontal_asymptotes() const {

        // TODO: Filter out curves with no arc to +/- infty

        Solve_1 solve_1 = kernel()->solve_1_object();

        Polynomial_1 leading_coefficient_in_x
            = CGAL::leading_coefficient(typename Polynomial_traits_2::Swap()
                                        (polynomial_2(),0,1));
        std::vector<Algebraic_real_1> roots_of_lcoeff;

        solve_1(leading_coefficient_in_x,
                false,
                std::back_inserter(roots_of_lcoeff));

        std::vector<Bound> stripe_bounds;
        find_intermediate_values(kernel(),
				 roots_of_lcoeff.begin(),
                                 roots_of_lcoeff.end(),
                                 std::back_inserter(stripe_bounds));

        Bound leftmost_bound = bound_value_in_interval(0);
        Bound rightmost_bound = bound_value_in_interval(this->number_of_status_lines_with_event());

        for(size_type i=0;i<static_cast<size_type>(stripe_bounds.size());i++) {
            Bound& beta = stripe_bounds[i];
            Polynomial_1 poly_at_beta
	      = kernel()->evaluate_utcf_2_object()(this->polynomial_2(),beta);
            std::vector<Algebraic_real_1> x_coordinates_at_beta;
            solve_1(poly_at_beta, false, std::back_inserter(x_coordinates_at_beta));
            size_type number_of_roots
                = static_cast<size_type>(x_coordinates_at_beta.size());
            if(number_of_roots>0) {
                if(leftmost_bound > x_coordinates_at_beta[0].low()) {
                    leftmost_bound = x_coordinates_at_beta[0].low();
                }
                if(rightmost_bound < x_coordinates_at_beta[number_of_roots-1].high()) {
                    rightmost_bound = x_coordinates_at_beta[number_of_roots-1].high();
                }
            }
        }

        // Just to be sure...
        leftmost_bound = leftmost_bound - 1;
        rightmost_bound = rightmost_bound + 1;

        Polynomial_1 curve_at_left_end
	= kernel()->evaluate_utcf_2_object()
	  (typename Polynomial_traits_2::Swap() (this->polynomial_2(),0,1),
	   leftmost_bound);
        std::vector<Algebraic_real_1> roots_at_left_end;
        solve_1(curve_at_left_end, false, std::back_inserter(roots_at_left_end));
        size_type number_of_roots_at_left_end
            = static_cast<size_type>(roots_at_left_end.size());
        std::vector<Asymptote_y> asym_left_info;
        size_type current_stripe=0,i=0;
        while(i<number_of_roots_at_left_end) {
            if(current_stripe==static_cast<size_type>(stripe_bounds.size())) {
                asym_left_info.push_back( CGAL::make_object
                                              (CGAL::ARR_TOP_BOUNDARY) );
                i++;
                continue;
            }
            if(roots_at_left_end[i].low() > stripe_bounds[current_stripe]) {
                current_stripe++;
                continue;
            }
            if(roots_at_left_end[i].high() < stripe_bounds[current_stripe]) {
                if(current_stripe==0) {
                    asym_left_info.push_back(CGAL::make_object
                                                 (CGAL::ARR_BOTTOM_BOUNDARY));
                    i++;
                    continue;
                } else {
                    asym_left_info.push_back(CGAL::make_object
                                 (roots_of_lcoeff[current_stripe-1]));
                    i++;
                    continue;
                }
            }
            roots_at_left_end[i].refine();
        }
        this->ptr()->horizontal_asymptotes_left = asym_left_info;



        Polynomial_1 curve_at_right_end
	= kernel()->evaluate_utcf_2_object()
  	    (typename Polynomial_traits_2::Swap() (this->polynomial_2(),0,1),
             rightmost_bound);
        std::vector<Algebraic_real_1> roots_at_right_end;
        solve_1(curve_at_right_end, false, std::back_inserter(roots_at_right_end));
        size_type number_of_roots_at_right_end
            = static_cast<size_type>(roots_at_right_end.size());
        std::vector<Asymptote_y> asym_right_info;
        current_stripe=0;
        i=0;
        while(i<number_of_roots_at_right_end) {
            if(current_stripe==static_cast<size_type>(stripe_bounds.size())) {
                asym_right_info.push_back(CGAL::make_object
                                              (CGAL::ARR_TOP_BOUNDARY) );
                i++;
                continue;
            }
            if(roots_at_right_end[i].low() > stripe_bounds[current_stripe]) {
                current_stripe++;
                continue;
            }
            if(roots_at_right_end[i].high() < stripe_bounds[current_stripe]) {
                if(current_stripe==0) {
                    asym_right_info.push_back(CGAL::make_object
                                                  (CGAL::ARR_BOTTOM_BOUNDARY));
                    i++;
                    continue;
                } else {
                    asym_right_info.push_back
                        (CGAL::make_object(roots_of_lcoeff[current_stripe-1]));
                    i++;
                    continue;
                }
            }
            roots_at_right_end[i].refine();
        }
        this->ptr()->horizontal_asymptotes_right = asym_right_info;

    }

    //! @}

public:
    template<typename OutputIterator> void get_roots_at_rational
    (Bound r, OutputIterator it) const {

        typename Rep::Intermediate_cache::Find_result find_result
            = this->ptr()->intermediate_cache.find(r);

	std::vector<Algebraic_real_1> p_roots;

        if(find_result.second) {
            p_roots = find_result.first->second;
        } else {
	    Polynomial_2 swapped = typename Polynomial_traits_2::Swap()
                              	    (this->polynomial_2(), 0, 1);
	    Polynomial_1 p = kernel()->evaluate_utcf_2_object()(swapped,r);
            // TODO 2012 use an isolator here? if: rational-square-free from lifter is a candidate
#if 0 && CGAL_BISOLVE_ENABLE_ARCAVOID
            typedef CGAL::Arcavoid< CGAL::internal::Bitstream_coefficient_kernel<Coefficient>,
                                    CGAL::Arcavoid_real_root_isolator_tag > Num_isolator;

            CGAL::internal::Real_roots< Algebraic_real_1, Num_isolator > real_roots;
            real_roots(p, std::back_inserter(p_roots));
#else
            kernel()->solve_1_object()(p,false,std::back_inserter(p_roots));
#endif
            this->ptr()->intermediate_cache.insert(std::make_pair(r,p_roots));

        }
        std::copy(p_roots.begin(),p_roots.end(),it);
    }


    // \name Internal functions for Conic optimization
    //! @{

#if CGAL_ACK_USE_SPECIAL_TREATMENT_FOR_CONIX

private:

    bool conic_is_y_regular() const {
        CGAL_error_msg("Implement me");
        return false;
    }

    bool conic_has_vertical_component() const {
        CGAL_error_msg("Implement me");
        return false;
    }

    size_type conic_number_of_status_lines_with_event() const {
        CGAL_error_msg("Implement me");
        return 0;
    }

    void conic_x_to_index(Algebraic_real_1 x,size_type& i,bool& is_event) const
    {
        CGAL_error_msg("Implement me");
    }

    Status_line_1& conic_status_line_at_event(size_type i) const {
        CGAL_error_msg("Implement me");
        // Just a random status line to make compiler happy
        return this->ptr()->vert_line_at_rational_map[Bound(0)];
    }

    Status_line_1& conic_status_line_at_exact_x(Bound b) const {
        CGAL_error_msg("Implement me");
        return this->ptr()->vert_line_at_rational_map[Bound(0)];
    }

    Status_line_1& conic_status_line_at_exact_x(Algebraic_real_1 alpha) const {
        CGAL_error_msg("Implement me");
        return this->ptr()->vert_line_at_rational_map[Bound(0)];
    }

    Status_line_1 conic_status_line_of_interval(size_type i) const {
        CGAL_error_msg("Implement me");
        return this->ptr()->vert_line_at_rational_map[Bound(0)];
    }

    Status_line_1 conic_status_line_for_x
        (Algebraic_real_1 x,
         CGAL::Sign perturb = CGAL::ZERO) const {
        CGAL_error_msg("Implement me");
        return this->ptr()->vert_line_at_rational_map[Bound(0)];
    }

    size_type conic_arcs_over_interval(size_type i) const {
        CGAL_error_msg("Implement me");
        return -1;
    }

    Bound conic_bound_value_in_interval(size_type i) const {
        CGAL_error_msg("Implement me");
        return Bound(0);
    }

    Polynomial_1 conic_content() const {
        CGAL_error_msg("Implement me");
        return Polynomial_1();
    }

    Polynomial_2 conic_primitive_polynomial_2() const {
        CGAL_error_msg("Implement me");
        return Polynomial_2();
    }

    Self& conic_shear_primitive_part(Integer s) const {
        CGAL_error_msg("Implement me");
        return Self();
    }

    void conic_refine_all(Bound precision) {
        CGAL_error_msg("Implement me");
    }

    Asymptote_y conic_asymptotic_value_of_arc(CGAL::Arr_parameter_space loc,
                                              size_type arcno) const {
        CGAL_error_msg("Implement me");
        return Asymptote_y();
    }

#endif


    //! @}

    //! \name friends
    //! @{

    // friend function for id-based hashing
    friend std::size_t hash_value(const Self& x) {
        return static_cast<std::size_t>(x.id());
    }

#if !CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE
    // another friend
    friend class Shear_transformation<Algebraic_kernel_with_analysis_2>;
#endif

    // status lines are also friends
    friend class CGAL::internal::Status_line_CA_1_rep< Algebraic_kernel_with_analysis_2 > ;

    //! @}

}; // class Curve_analysis_2


//! \brief Prints the objects.
template<typename AlgebraicKernelWithAnalysis_2,
         typename Rep_>
std::ostream& operator<< (
        std::ostream& out,
        const Curve_analysis_2< AlgebraicKernelWithAnalysis_2,
        Rep_ >& curve) {

  typedef AlgebraicKernelWithAnalysis_2 Algebraic_kernel_with_analysis_2;

  typedef Rep_ Rep;

  typedef Curve_analysis_2< Algebraic_kernel_with_analysis_2, Rep > Curve;

  typedef typename Curve::size_type size_type;
  typedef typename Curve::Asymptote_y Asymptote_y;


    switch (::CGAL::get_mode(out)) {
    case ::CGAL::IO::PRETTY: {

      out << "--------------- Analysis results ---------------" << std::endl;
      out << "Number of constructed event lines: "
          << curve.number_of_status_lines_with_event()
          << std::endl;
      out << "(Horizontal) asymptotes at -infty (#arcs=" << curve.arcs_over_interval(0) << "): " << std::flush;
      for (size_type i = 0; i < curve.arcs_over_interval(0); i++) {
        const Asymptote_y& curr_asym_info_obj
          = curve.asymptotic_value_of_arc(CGAL::ARR_LEFT_BOUNDARY,i);
        typename Curve::Algebraic_real_1 curr_asym_info;
        bool is_finite = CGAL::assign(curr_asym_info,curr_asym_info_obj);
        if (!is_finite) {
          // Assignment to prevent compiler warning
          CGAL::Arr_parameter_space loc = CGAL::ARR_LEFT_BOUNDARY;
          CGAL_assertion_code(bool is_valid = )
            CGAL::assign(loc, curr_asym_info_obj);
          CGAL_assertion(is_valid);
          if (loc == CGAL::ARR_TOP_BOUNDARY) {
            out << "+infty " << std::flush;
          } else {
            CGAL_assertion(loc == CGAL::ARR_BOTTOM_BOUNDARY);
            out << "-infty " << std::flush;
          }
        } else { // is_finite
          out << CGAL::to_double(curr_asym_info)
              << " " << std::flush;
        }
      }

      out << std::endl;

      out << "Intermediate line at "
          << CGAL::to_double(curve.bound_value_in_interval(0))
          << ": " << curve.arcs_over_interval(0) << " passing arcs"
          << std::endl
          << std::endl;
      for (size_type i = 0; i < curve.number_of_status_lines_with_event();
           i++) {
        out << curve.status_line_at_event(i) << std::endl;
        out << "Intermediate line at "
            << CGAL::to_double(curve.bound_value_in_interval(i+1))
            << ": " << curve.arcs_over_interval(i+1)
            << " passing arcs" << std::endl
            << std::endl;
      }
      size_type no_events = curve.number_of_status_lines_with_event();
      out << "(Horizontal) asymptotes at +infty (#arcs=" << curve.arcs_over_interval(no_events) << "): " << std::flush;
      for (size_type i = 0; i < curve.arcs_over_interval(no_events); i++) {
        const Asymptote_y& curr_asym_info_obj
          = curve.asymptotic_value_of_arc(CGAL::ARR_RIGHT_BOUNDARY,i);
        typename Curve::Algebraic_real_1 curr_asym_info;
        bool is_finite = CGAL::assign(curr_asym_info,curr_asym_info_obj);
        if(! is_finite) {
          // Assignment to prevent compiler warning
          CGAL::Arr_parameter_space loc = CGAL::ARR_LEFT_BOUNDARY;
          CGAL_assertion_code(bool is_valid = )
            CGAL::assign(loc, curr_asym_info_obj);
          CGAL_assertion(is_valid);
          if(loc == CGAL::ARR_TOP_BOUNDARY) {
            out << "+infty " << std::flush;
          } else {
            CGAL_assertion(loc == CGAL::ARR_BOTTOM_BOUNDARY);
            out << "-infty " << std::flush;
          }
        } else { // is_finite
          out << CGAL::to_double(curr_asym_info)
              << " " << std::flush;
        }
      }

      out << std::endl;

      out << "------------------------------------------------" << std::endl;
      break;
    }
    case ::CGAL::IO::BINARY:
      std::cerr << "BINARY format not yet implemented" << std::endl;
      break;
    default:
      // ASCII
      out << curve.polynomial_2();
    }

    return out;
}

//! \brief Reads the objects from stream
template<typename AlgebraicKernelWithAnalysis_2,
         typename Rep_>
std::istream& operator>> (
    std::istream& is,
    Curve_analysis_2< AlgebraicKernelWithAnalysis_2, Rep_ >& curve) {

  CGAL_precondition(CGAL::is_ascii(is));

  typedef AlgebraicKernelWithAnalysis_2 Algebraic_kernel_with_analysis_2;

  typedef Rep_ Rep;

  typename Curve_analysis_2< Algebraic_kernel_with_analysis_2, Rep >::
    Polynomial_2 f;

  is >> f;

  // TODO is get_static_instance the right way?
  curve = Algebraic_kernel_with_analysis_2::get_static_instance().
    construct_curve_2_object()(f);

  return is;
}


} //namespace CGAL


#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif


#endif // ALGEBRAIC_CURVE_2_H
