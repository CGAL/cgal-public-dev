// Copyright (c) 2009, 2010, 2011 Max-Planck-Institut fuer Informatik (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: ??
//
//
// Author(s):  Pavel Emeliyanenko <asm@mpi-inf.mpg.de>
// ============================================================================

#ifndef CGAL_LITE_CURVE_ANALYSIS_2
#define CGAL_LITE_CURVE_ANALYSIS_2

#include <CGAL/config.h>
#include <CGAL/Arr_enums.h>

#if CGAL_LITE_CA_VERBOSE
#define LCA_out(x) std::cout << x;
#define dbl(x) CGAL::to_double(x)
#define bfi(x) CGAL::lower(CGAL::convert_to_bfi(x))
#else
#define LCA_out(x) static_cast< void >(0);
#endif

#define STILL_ALIVE std::cout << __LINE__ << "\n";


// #ifndef CGAL_BISOLVE_ENABLE_NTL_FACTORIZE
// #include <CGAL/symbolic_exports.h>
// #endif

#include <CGAL/Handle_with_policy.h>

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>


namespace CGAL {

namespace internal {


/*!\brief
 * Default representation class for Arc_2
 */
template < class LiteCurveAnalysis_2 >
class Lite_arc_2_rep {

public:
    //!\name Public types
    //!@{

    typedef LiteCurveAnalysis_2 Lite_curve_analysis_2;

    typedef typename Lite_curve_analysis_2::AK_1 AK_1;

    //! type of boundary value in x-range of an arc
    typedef typename AK_1::Bound Bound;

    //! type of boundary value in x-range of an arc
    typedef typename Lite_curve_analysis_2::Coordinate_1 Coordinate_1;

    //!@}

public:
    //!\name Constructors
    //!@{
    
    //! default constructor
    Lite_arc_2_rep() : _m_arcno(-1), _m_is_vertical(false) {
    }

    //! standard constructor
    Lite_arc_2_rep(const Coordinate_1& p, const Coordinate_1& q,
            const Lite_curve_analysis_2& c, int arcno) :
        _m_min(p), _m_max(q), _m_support(c),
        _m_arcno(arcno), _m_is_vertical(false),
        _m_loc_min(CGAL::ARR_INTERIOR), _m_loc_max(CGAL::ARR_INTERIOR) {
    }

    Lite_arc_2_rep(const Coordinate_1& point,
          CGAL::Arr_curve_end inf_end,
          const Lite_curve_analysis_2& c, int arcno) :
        _m_support(c), _m_arcno(arcno), _m_is_vertical(false) {

        if(inf_end == CGAL::ARR_MIN_END) {
            _m_loc_min = CGAL::ARR_LEFT_BOUNDARY;
            _m_loc_max = CGAL::ARR_INTERIOR;
            _m_max = point;
        } else {
            _m_loc_min = CGAL::ARR_INTERIOR;
            _m_loc_max = CGAL::ARR_RIGHT_BOUNDARY;
            _m_min = point;
        }
    }

    /*!\brief
     * Constructs a non-vertical arc with two non-interior ends at the
     * left and right boundary (branch I)
     *
     * \param c The supporting curve
     * \param arcno The arcnumber wrt to \c c in the interior of the arc
     * \return The constructed branch
     */
    Lite_arc_2_rep(const Lite_curve_analysis_2& c, int arcno) :
       _m_support(c), _m_arcno(arcno), _m_is_vertical(false),
       _m_loc_min(CGAL::ARR_LEFT_BOUNDARY),
       _m_loc_max(CGAL::ARR_RIGHT_BOUNDARY) {
    }

    //!@}

public:
    //!\name Data members
    //!@{
    
    //! minimal end-points of an arc
    Coordinate_1 _m_min;

    //! maximal end-points of an arc
    Coordinate_1 _m_max;
    
    //! supporting curve
    Lite_curve_analysis_2 _m_support;

    //! interior arcno
    int _m_arcno;

    //! indicates whether arc is vertical
    bool _m_is_vertical;

    CGAL::Arr_parameter_space _m_loc_min;
    CGAL::Arr_parameter_space _m_loc_max;

    //!@}
};

} // namespace internal 

template < class LiteCurveAnalysis_2 >
class Lite_arc_2 :
        public CGAL::Handle_with_policy<
            internal::Lite_arc_2_rep< LiteCurveAnalysis_2 > > {
  
public:
    //!\name Public types
    //!@{

    typedef LiteCurveAnalysis_2 Lite_curve_analysis_2;
    
    typedef internal::Lite_arc_2_rep< Lite_curve_analysis_2 > Rep;

    typedef typename Rep::AK_1 AK_1;

    //! type of "rational" value in x-range
    typedef typename Rep::Bound Bound;

    typedef typename Rep::Coordinate_1 Coordinate_1;

    //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;
    //!@}

    using Base::ptr;

public:
    //!\name basic constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Lite_arc_2() : 
        Base(Rep()) {   
    }

    //!@}
   
public:
    //!\name Constructors for non-vertical arcs
    //!@{
    
    /*!\brief 
     * Constructs an arc with two interior end-points (segment).
     * 
     * \param p first endpoint
     * \param q second endpoint
     * \param c The supporting curve
     * \param arcno The arcnumber wrt \c c in the interior of the arc
     */
    Lite_arc_2(const Coordinate_1& p, const Coordinate_1& q,
          const Lite_curve_analysis_2& c, int arcno) :
                Base(Rep(p, q, c, arcno)) { 
    }
      
   /*!\brief
     * Constructs an arc with one interior end-point and another end
     * at the left or right boundary of the parameter space (ray I).
     *
     * \param origin The interior end-point of the ray
     * \param inf_end Defining whether the arcs emanates from the left or right
     *        boundary
     * \param c The supporting curve
     * \param arcno The arcnumber wrt \c c in the interior of the arc
     * \param arcno_o The arcnumber wrt \c c of the arc at \c origin 
     * \return The constructed ray
     */
    Lite_arc_2(const Coordinate_1& point, CGAL::Arr_curve_end inf_end,
          const Lite_curve_analysis_2& c, int arcno) :
            Base(Rep(point, inf_end, c, arcno)) {
    }

    /*!\brief
     * Constructs a non-vertical arc with two non-interior ends at the
     * left and right boundary (branch I)
     *
     * \param c The supporting curve
     * \param arcno The arcnumber wrt to \c c in the interior of the arc
     * \return The constructed branch
     */
    Lite_arc_2(const Lite_curve_analysis_2& c, int arcno) :
        Base(Rep(c, arcno)) {
    }

    /*!\brief
     * Constructs a non-vertical arc with one interior and one asymptotic end
     */    
//     Lite_arc_2(const Coordinate_1& pt_asympt, 
//         CGAL::Arr_curve_end end_asympt, const Coordinate_1& pt_interior, 
//           const Lite_curve_analysis_2& c, int arcno) :
//             Base(Rep(pt_asympt, end_asympt, pt_interior, c, arcno)) { 
//     }

    /*!\brief
     * Constructs a non-vertical arc with two asymptotic ends
     */    
//     Lite_arc_2(const Coordinate_1& pt1_asympt, 
//         CGAL::Arr_curve_end end1_asympt, const Coordinate_1& pt2_asympt, 
//         CGAL::Arr_curve_end end2_asympt,
//             const Lite_curve_analysis_2& c, int arcno) :
//                 Base(Rep(pt1_asympt, end1_asympt, pt2_asympt, end2_asympt, 
//                     c, arcno)) { 
//     }

    /*!\brief
     * Constructs a non-vertical arc with one infinite and one asymptotic end
     */    
//     Lite_arc_2(const Coordinate_1& pt_asympt, 
//         CGAL::Arr_curve_end end_asympt, CGAL::Arr_curve_end inf_end,
//           const Lite_curve_analysis_2& c, int arcno) :
//             Base(Rep(pt_asympt, end_asympt, inf_end, c, arcno)) { 
//     }
   
    //!@}
public:
    //!\name Parameter space
    //!@{
    
    /*!\brief
     * location of arc's end
     *
     * \param ce The intended end
     * \return The location of arc's \c ce in parameterspace
     */
    CGAL::Arr_parameter_space location(CGAL::Arr_curve_end ce) const {

        if(ce == CGAL::ARR_MIN_END) {
            return ptr()->_m_loc_min;
        }
        return ptr()->_m_loc_max;
    }
    //!@}

    //!\name Access functions
    //!@{
    
    /*!\brief
     * Is a curve-end finite?
     *
     * \param ce The intended end
     * \return \c true, if finite, \c false, otherwise
     */
    bool is_finite(CGAL::Arr_curve_end ce) const {

        CGAL::Arr_parameter_space loc = (ce == CGAL::ARR_MIN_END ?
                 ptr()->_m_loc_min : ptr()->_m_loc_max);

        return (loc == CGAL::ARR_INTERIOR);
    }

    /*!\brief 
     * returns arc's interior curve end
     * 
     * \param ce The intended end
     * \return The minimal point of the arc, or the maximal point of the arc
     *
     *  \pre accessed curve end has finite coordinates
     */
//     Point_2 curve_end(CGAL::Arr_curve_end ce) const {
// 
//         const Point_2& pt = 
//             (ce == CGAL::ARR_MIN_END ? _minpoint() : _maxpoint());
//         return pt;
//     }

    /*!\brief 
     * returns x-coordinate of arc's curve end
     * 
     * \param ce The intended end
     * \return x-coordinate of arc's end at \c ce
     *
     * \pre accessed curve end has finite x-coordinate
     */
    inline
    Coordinate_1 curve_end_x(CGAL::Arr_curve_end ce) const {

        return (ce == CGAL::ARR_MIN_END ? ptr()->_m_min : ptr()->_m_max);
    }

    /*!\brief
     * supporting curve of the arc
     */
    inline
    const Lite_curve_analysis_2& curve() const {
        return ptr()->_m_support;
    }
  
    /*!\brief arc number in interior
     */
    inline
    int arcno() const { 
        return ptr()->_m_arcno;
    }

    /*!\brief
     * checks if the arc is vertical 
     */
    inline
    bool is_vertical() const {
        return ptr()->_m_is_vertical;
    }

};


namespace internal {    

//! template rep class for new bi-solve
template < class AlgebraicKernel_1 >
class Lite_curve_analysis_rep {

public:

    typedef AlgebraicKernel_1 AK_1;
    
    //! type of Coefficient
    typedef typename AK_1::Coefficient Coefficient;

    //! type for univariate polynomials
    typedef typename AK_1::Polynomial_1 Polynomial_1;

    //! type for bivariate polynomials
    typedef CGAL::Polynomial< Polynomial_1 > Polynomial_2;

    typedef CGAL::Polynomial_traits_d< Polynomial_2 > PT_2;

    //! type for algebraic reals
    typedef typename AK_1::Algebraic_real_1 Coordinate_1;

    //! type for Multiplicity
    typedef typename AK_1::Multiplicity_type Multiplicity_type;

    //! type for Bound
    typedef typename AK_1::Bound Bound;

public:
    //!\name Constructors
    //!@{
    
    //! default constructor
    Lite_curve_analysis_rep() {
    }

    Lite_curve_analysis_rep(const Polynomial_2& poly) : _m_poly(poly) {

        typename PT_2::Swap swap;
        _m_poly_swapped = swap(poly, 0, 1);
    }
    
    //!@}
public:
    //!\name Data members
    //!@{

    Polynomial_2 _m_poly;
    Polynomial_2 _m_poly_swapped;

    Polynomial_1 _m_res;

    //!@}
};

} // namespace internal

template < class AlgebraicKernel_1,
        template <class> class Rep_ = internal::Lite_curve_analysis_rep >
class Lite_curve_analysis_2 :
        public CGAL::Handle_with_policy< Rep_< AlgebraicKernel_1 > > {
  
public:
    //!\name Public types
    //!@{

    typedef AlgebraicKernel_1 AK_1;
    
    //! this instance's second template parameter
    typedef Rep_< AK_1 > Rep;

    //! this instance itself
    typedef Lite_curve_analysis_2< AK_1, Rep_ > Self;

    //! type of an x-monotone arc
    typedef Lite_arc_2< Self > Arc_2;
    
    //! type for univariate polynomials
    typedef typename Rep::Polynomial_1 Polynomial_1;

    //! type for bivariate polynomials
    typedef typename Rep::Polynomial_2 Polynomial_2;

    //! type for algebraic reals
    typedef typename Rep::Coordinate_1 Coordinate_1;

    //! type for Multiplicity
    typedef typename Rep::Multiplicity_type Multiplicity_type;

    //! type of Coefficient
    typedef typename AK_1::Coefficient Coefficient;

    //! type of Coefficient
    typedef typename AK_1::Bound Bound;

    //! arithmetic kernel
    typedef typename CGAL::Get_arithmetic_kernel< Bound >::
            Arithmetic_kernel Arithmetic_kernel;
  
    //! bigfloat interval 
    typedef typename Arithmetic_kernel::Bigfloat_interval BFI;
  
    //! our lovely bigfloats
    typedef typename CGAL::Bigfloat_interval_traits< BFI >::Bound BigFloat;  

    //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;

    using Base::ptr;
    //!@}

public:
    //!@{

    /*!\brief
     * Default constructor
     */
    Lite_curve_analysis_2() : 
        Base(Rep()) {   
    }

    Lite_curve_analysis_2(const Polynomial_2& poly) : Base(Rep(poly)) {
    }

    const Polynomial_2& polynomial_2() const {
        return ptr()->_m_poly;
    }

    void compute(std::vector< Arc_2 >& arcs) {
    
        typename CGAL::Polynomial_traits_d< Polynomial_1 >::
//             Make_square_free msf;
            Square_free_factorize_up_to_constant_factor factorize;

        ptr()->_m_res = CGAL::resultant(ptr()->_m_poly, 
                 CGAL::differentiate(ptr()->_m_poly));

        std::cout << "##### resultant computed\n ";
//              <<  ptr()->_m_res << "\n\n";

        std::vector< std::pair< Polynomial_1, Multiplicity_type > >
                res_factors;

        std::cout << "##### starting sqfree factorization \n ";
        square_free_factorize(ptr()->_m_res,
                 std::back_inserter(res_factors));

        typename AK_1::Solve_1 solve;
        typename AK_1::Bound_between_1 between_bnd;
        typename AK_1::Algebraic_real_traits::Lower_bound lower_bnd;
        typename AK_1::Algebraic_real_traits::Upper_bound upper_bnd;

        std::cout << "##### found " << res_factors.size() <<
                    " factor(s).. root solving... \n ";

        unsigned i, j, n, v1, v2;
        std::vector< Coordinate_1 > roots;

        for(i = 0; i < res_factors.size(); i++) {
            solve(res_factors[i].first, true, std::back_inserter(roots));
        }
        std::sort(roots.begin(), roots.end());

        std::cout << "##### " << roots.size() << " roots found \n ";
//             solve(ptr()->_m_res, true, std::back_inserter(roots));

        std::cout.precision(50);
        int prec = CGAL::set_precision(BFI(), 53);
        

        if(roots.size() == 0) {
            n = _n_roots_at_rational(Bound(0));
            for(i = 0; i < n; i++) {
                arcs.push_back(Arc_2(*this, i));
            }
            return;
        }

        const Coordinate_1& x0 = roots[0];
        v2 = _n_vanishing_coeffs(x0);
        Bound b = lower_bnd(x0) - Bound(1);
        unsigned deg_y = ptr()->_m_poly.degree();
        
        if(v2 != 0)
            std::cout << "### root 0 asymptote\n";

        n = _n_roots_at_rational(b);
        for(i = 0; i < n; i++) {
            arcs.push_back(Arc_2(x0, CGAL::ARR_MIN_END, *this, i));
        }

// NOTE NOTE: handle asymptotes as interior points: no distinction!!
        for(i = 0; i < roots.size()-1; i++) {
            v1 = v2;
            const Coordinate_1& x1 = roots[i], x2 = roots[i+1];
            v2 = _n_vanishing_coeffs(x2);
            if(v2 != 0)
                std::cout << "### root " << (i+1) << " asymptote\n";

//             if(x1 != x2)
//             std::cout << (x1) << " - " << (x2) << " looking for bnd\n";
            b = between_bnd(x1, x2);
            
            std::cout << "bound: " << bfi(b) << "\n";
            n = _n_roots_at_rational(b);

            for(j = 0; j < n; j++) {
                arcs.push_back(Arc_2(x1, x2, *this, j));
            }
        }
        
        const Coordinate_1& xn = roots[roots.size()-1];
        b = upper_bnd(xn) + Bound(1);
        std::cout << (i+1) << ": " << bfi(b) << "\n";
        n = _n_roots_at_rational(b);

        for(i = 0; i < n; i++) {
            arcs.push_back(Arc_2(xn, CGAL::ARR_MAX_END, *this, i));
        }
        CGAL::set_precision(BFI(), prec);

        std::cout.precision(10);
    }


protected:

    //! counts the number of vanishing leading coeffs of the curve at \c x
    unsigned _n_vanishing_coeffs(const Coordinate_1& x) {

        unsigned n = 0;

//! HACK HACK: no actual need to check for all vanishing coeffs:
//! only for vertical lines suffices
//         typename AK_1::Sign_at_1 sign_at;
//         for(int i = ptr()->_m_poly.degree(); i >= 0; i--, n++) {
//             if(sign_at(ptr()->_m_poly[i], x) != CGAL::ZERO)
//                 break;
//         }
//         std::cout << "n-vanished: " << n << "\n";
        return n;
    }

    Polynomial_1 subst_homogeneous_x(const Polynomial_2& poly, 
        Coefficient num, Coefficient denom, int in_degree) const {

        typename Polynomial_2::const_iterator ci2;
        typename Polynomial_1::Vector v(poly.degree() + 1);

        int i;
        for(ci2 = poly.begin(), i = 0; ci2 != poly.end(); ci2++, i++) {
    
            const Polynomial_1& pp = *ci2;
            typename Polynomial_1::const_iterator ci = pp.end() - 1;
        
            Coefficient z(*ci), det(denom);
            while((ci--) != pp.begin()) {
                z = z * num + (*ci) * det;
                det = det * denom;
            }
            z = z * CGAL::ipower(denom, in_degree - pp.degree());
            v[i] = z;
        }
        return Polynomial_1(v.begin(), v.end());
    }

    unsigned _n_roots_at_rational(Bound b) const {

        typedef CGAL::Fraction_traits<Bound> FT;
        typedef typename FT::Numerator_type Numerator;
        typedef typename FT::Denominator_type Denominator;
        typename FT::Decompose decompose;
    
        Coefficient num, denom;
        decompose(b, num, denom);
        Polynomial_1 fiber = subst_homogeneous_x(ptr()->_m_poly, num, denom,
                 ptr()->_m_poly_swapped.degree());

//         Polynomial_1 fiber = CGAL::evaluate_homogeneous(
//                 ptr()->_m_poly_swapped, num, denom);

//         std::cout << "bound: " << bfi(b) << "; " << fiber << "\n";

        std::vector< Coordinate_1 > roots;

        std::cout << " solving for roots ..\n";

        typename AK_1::Solve_1 solve;
        solve(fiber, true, std::back_inserter(roots));

        std::cout << roots.size() << " root(s) found ..\n";
        return roots.size();
    }

    //!@}
};

} // namespace CGAL

#endif // CGAL_LITE_CURVE_ANALYSIS_2
