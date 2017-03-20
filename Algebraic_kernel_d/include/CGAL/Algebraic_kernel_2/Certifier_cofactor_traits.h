// Copyright (c) 2010, 2011 Max-Planck-Institut fuer Informatik (Germany).
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
// $URL: svn+ssh://eric@scm.gforge.inria.fr/svn/cgal/branches/features/Symbolic-mpi/Algebraic_kernel_d/include/CGAL/Algebraic_kernel_d/Curve_analysis_2.h $
// $Id: Curve_analysis_2.h 71840 2012-08-30 11:29:52Z eric $
//
//
// Author(s): Eric Berberich <eric.berberich@cgal.org>
//            Pavel Emeliyanenko <asm@mpi-inf.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_2_CERTIFIER_COFACTOR_TRAITS_H
#define CGAL_ALGEBRAIC_KERNEL_2_CERTIFIER_COFACTOR_TRAITS_H 

/*! \file
 * Traits class implementing cofactors root certification
 */

#include <CGAL/config.h>
#include <CGAL/Algebraic_kernel_2/Bi_solve_2_flags.h>

#include <CGAL/Algebraic_kernel_2/Certifier_traits_base.h>
#include <CGAL/Algebraic_kernel_2/Range_analysis_2.h>

#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>

namespace CGAL {

namespace internal {

template < class AlgebraicKernel_d_1 >
class Certifier_cofactor_traits : public
            Certifier_traits_base < AlgebraicKernel_d_1 > {

public:

    //! type of AK_1
    typedef AlgebraicKernel_d_1 Algebraic_kernel_d_1;

    //! base class
    typedef Certifier_traits_base < Algebraic_kernel_d_1 > Base;

    //! type of univariate polynomial
    typedef typename Base::Polynomial_1 Polynomial_1;

    //! type of bivariate polynomial
    typedef typename Base::Polynomial_2 Polynomial_2;  

    //! type of algebraic real
    typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real_1;

    //! type of Multiplicity
    typedef typename Algebraic_kernel_d_1::Multiplicity_type Multiplicity_type;

    //! type of Bound
    typedef typename Base::Bound Bound;

    //! arithmetic kernel
    typedef typename Base::AK AK;

    //! bigfloat interval 
    typedef typename AK::Bigfloat_interval BFI;

    //! our lovely bigfloats
    typedef typename CGAL::Bigfloat_interval_traits< BFI >::Bound BigFloat;

    //! a root box
    typedef typename Base::Root_box Root_box;

    //! range analysis with BFI using real arithmetic
    typedef CGAL::Range_analysis_2< Polynomial_2, BFI, false > IA_real_bfi;

    struct Candidate_props {

        Candidate_props() :
            _norm_x_passed(false), _norm_y_passed(false) {
        }

        //! lower bound for resultant
        mutable boost::optional< BigFloat > lb_res_in_x;
    
        //! lower bound for resultant
        mutable boost::optional< BigFloat > lb_res_in_y;
    
        //! upper bound m_u and m_v in x
        mutable boost::optional< std::pair< BigFloat, BigFloat > > m_uv_in_x;
    
        //! upper bound m_u and m_v in y
        mutable boost::optional< std::pair< BigFloat, BigFloat > > m_uv_in_y;
    
        //! norm in x holds
        mutable bool _norm_x_passed;

        //! norm in y holds
        mutable bool _norm_y_passed;
        
    };

public:
    //!\name Constuctors
    //!@{

    //! default
    Certifier_cofactor_traits() :
        Base() {
    }

    //! standard
    Certifier_cofactor_traits
      (Algebraic_kernel_d_1 *ak,
       const Polynomial_2& f, const Polynomial_2& g,
       const Polynomial_2& ft, const Polynomial_2& gt,
       const Polynomial_1& res_x, const Polynomial_1& res_y) :
        Base(ak, f, g, ft, gt, res_x, res_y) {
    }

    //!@}
public:

    //!\name main stuff
    //!@{

    //! returns \c true if the candidate passes the norm test
    bool norm_test(const Candidate_props& ps,
            const Bound& x_0, const Bound& y_0,
            const Root_box& box_y, Multiplicity_type mul_y, long& prec) {

    typedef typename Base::Poly_bfi_2 Poly_bfi_2;
    typename CGAL::Polynomial_traits_d< Poly_bfi_2 >::Substitute
      substitute_2;

    long oldp = CGAL::set_precision(BFI(), prec);
    this->_approximate_polynomials(prec);

    BFI xf = CGAL::convert_to_bfi(x_0), yf = CGAL::convert_to_bfi(y_0);
    BFI pts[] = {xf, yf};
    BigFloat f_0, g_0;

    Bisolve_telemetry_code(t_subs.start();)
    f_0 = CGAL::upper(CGAL::abs(substitute_2(this->_m_f_bfi, pts, pts + 2)));
    g_0 = CGAL::upper(CGAL::abs(substitute_2(this->_m_g_bfi, pts, pts + 2)));
    Bisolve_telemetry_code(t_subs.stop();)
    
//     Bisolve_out("f_0 = " << " (" << (f_0) << ")" << std::endl)
//     Bisolve_out("g_0 = " << " (" << (g_0) << ")\n\n")

    // else (and final test - before refinement)
    bool& _norm_x_passed = ps._norm_x_passed;
    bool& _norm_y_passed = ps._norm_y_passed;

    CGAL_assertion(!(_norm_x_passed && _norm_y_passed));

    typename Real_embeddable_extension< BigFloat >::Ceil_log2_abs
             log2_abs;

     // compute m_u and m_v
    if (!_norm_x_passed) {

      Bisolve_telemetry_code(t_bounds.start();)
      if (!ps.m_uv_in_x) {
        ps.m_uv_in_x = _cofactors_bound(box_y, false);
      }
      // lower bound of resultant
      if (!ps.lb_res_in_x) {
        ps.lb_res_in_x =
          _lower_bound(this->_m_res_x, this->_m_mul_x, this->_m_box_x);
      }  
      Bisolve_telemetry_code(t_bounds.stop();)

      Bisolve_out("M_uv_x = " << " (" << bfi(ps.m_uv_in_x->first) << "; "
                 << bfi(ps.m_uv_in_x->second) << ")\n");

//         t_timer.start();
      // Value: this sum is compared with lower bound
      BigFloat value_in_x(ps.m_uv_in_x->first * f_0 +
                       ps.m_uv_in_x->second * g_0);

//       Bisolve_out("Value in x: (" << (value_in_x) << ")\n");
      Bisolve_out("Lower bound in x: " << *(ps.lb_res_in_x) << "\n");
        
      _norm_x_passed = (value_in_x < *(ps.lb_res_in_x));
//         t_timer.stop();

#if CGAL_BISOLVE_USE_ADJUSTABLE_PRECISION
      if(!_norm_x_passed) {
        prec = std::max(prec * 3 / 2,
                log2_abs(value_in_x) - log2_abs(*ps.lb_res_in_x) + 10);
        Bisolve_out("?????????????????? X setting precision to : " << prec <<
            "\n");
      }
#endif
    }
    
    if (!_norm_y_passed) {

      Bisolve_telemetry_code(t_bounds.start();)
      if (!ps.m_uv_in_y) {
           ps.m_uv_in_y = _cofactors_bound(box_y, true);
      }
       // lower bound of resultant
      if (!ps.lb_res_in_y) {
        ps.lb_res_in_y =
          _lower_bound(this->_m_res_y, mul_y, box_y);
      }
      Bisolve_telemetry_code(t_bounds.stop();)
      
      Bisolve_out("M_uv_y = " << " (" << bfi(ps.m_uv_in_y->first) << "; "
                 << bfi(ps.m_uv_in_y->second) << ")\n");

//       t_timer.start();  
      // Value: this sums is compared with lower bound
      BigFloat value_in_y(ps.m_uv_in_y->first * f_0 +
                       ps.m_uv_in_y->second * g_0);

//       Bisolve_out("Value in y: (" << (value_in_y) << ")\n");
      Bisolve_out("Lower bound in y: " << *(ps.lb_res_in_y) << "\n");

      _norm_y_passed = (value_in_y < *(ps.lb_res_in_y));  
//       t_timer.stop();
    
#if CGAL_BISOLVE_USE_ADJUSTABLE_PRECISION
      if(!_norm_y_passed) {
        prec = std::max(prec * 3 / 2,
                log2_abs(value_in_y) - log2_abs(*ps.lb_res_in_y) + 10);
        Bisolve_out("?????????????????? Y setting precision to : " << prec <<
            "\n");
      }
#endif
    }

#if !CGAL_BISOLVE_USE_ADJUSTABLE_PRECISION
    prec *= 2;
#endif

    CGAL::set_precision(BFI(), oldp);

    if(_norm_x_passed && _norm_y_passed) {
        Bisolve_out("Candidate certified: Extended norm test holds"
                << std::endl)
        return true;
    }
    return false;
  }  

    //!@}
protected:
    //!@{

 //! compute lower bound in x
  BigFloat _lower_bound(const Polynomial_1& res, const Multiplicity_type& mult,
        const Root_box& box) const {

    typedef typename IA_real_bfi::Poly_1 Poly_1;
    long prec = 64, oldp = CGAL::get_precision(BFI());

    // Values: choose 1/4 as dist(x(), real_xmin()) \in (r,sqrt(13)r)
    //         1/4 < 1 < 3.6 = sqrt(13) < 4
    // NEXT: r(xbar) * 1/4^mult * (3n/(3n+5))^n-mult
    // better: r(xbar) * 1/4^mult * 1/6
    // better: r(xbar) * 1/4^(mult+2)

#if 1 // TODO 2012 added flag to distinguish bound 

    CGAL::set_precision(BFI(), prec);
    BigFloat crf; // compute correction with set precision
    if(box.used_t_test()) {
#if CGAL_BISOLVE_RELAXED_RATIOS 
      crf = CGAL::ipower(BigFloat(0.5), CGAL::degree(res));
#else
      crf = CGAL::ipower(BigFloat(0.25), mult + 2);
#endif
    } else {
        crf = CGAL::ipower(BigFloat(0.25), mult) *
            CGAL::ipower(BigFloat(0.5), CGAL::degree(res)-mult);
    }

    BFI xf64, xf, bfi;
    BigFloat bfval;

    typename Root_box::Bound exarg = box.mid() - 2*box.rad();
    xf64 = CGAL::convert_to_bfi (exarg);
    bfval = crf;
    for (typename Base::Squarefree_fact::const_iterator
           facit = this->_m_res_in_x_factors_ptr->begin();
         facit != this->_m_res_in_x_factors_ptr->end();
         ++facit) {
      if (facit->second == mult) {
        typename Root_box::Bound ex = CGAL::abs (facit->first.evaluate (exarg));
        BigFloat bf
          = BigFloat (ex.numerator(), std::round_toward_neg_infinity, 53)
          / BigFloat (ex.denominator(), std::round_toward_infinity, 53);
        // TODO: make sure ipower rounds into the proper direction
        bfval *= CGAL::ipower (bf, mult);
      } else {
        bfi = CGAL::abs (facit->first.evaluate (xf64));
        while (CGAL::is_zero (CGAL::lower (bfi))) {
          prec *= 2;
          CGAL::set_precision(BFI(), prec);
          // std::cerr << "Required precision " << prec << std::endl;
          xf = CGAL::convert_to_bfi (exarg);
          bfi = CGAL::abs (facit->first.evaluate (xf));
        }
        CGAL::set_precision(BFI(), 64);
        bfval *= CGAL::ipower (CGAL::lower (bfi), facit->second);
      }
    }

    CGAL::set_precision(BFI(), oldp);
    return bfval;

    typename Root_box::Bound exact = CGAL::abs (res.evaluate(box.mid() - 2*box.rad()));
    // std::cerr << "exact " << exact << std::endl;
    BigFloat bf
      = BigFloat (exact.numerator(), std::round_toward_neg_infinity, 53)
      / BigFloat (exact.denominator(), std::round_toward_infinity, 53)
      * crf;
    CGAL::set_precision(BFI(), oldp);
    return bf;
   

    while(1) {
        CGAL::set_precision(BFI(), prec);

        // precision change, no normalization
        _m_ia_real.set_polynomial(res, true, false);

        Poly_1 poly = _m_ia_real.internal_poly_1();
        // std::cerr << "Mult " << mult << ", Rad " << box.rad() << std::endl;
        xf = CGAL::convert_to_bfi(box.left() - box.rad());
        bfi = CGAL::abs(poly.evaluate(xf)) * crf;

        bfval = CGAL::lower(bfi) / BigFloat(1); // truncate bfval's mantissa
        if(CGAL::sign(bfval) == 1) {
          // std::cerr << bfi << std::endl;
          break;
        }
        // std::cerr << bfi << std::endl;
        prec = prec * 2;
        // Bisolve_out("increasing prec to: " << prec << "\n");
    }
    CGAL::set_precision(BFI(), oldp);
    return bfval;

#else

    Bound correction(1);
    if(box.used_t_test()) {
#if CGAL_BISOLVE_RELAXED_RATIOS
      correction =
        CGAL::ipower(Bound(1)/Bound(2), CGAL::degree(res));
#else
      correction =
        CGAL::ipower(Bound(1)/Bound(4), mult + 2);
#endif
    } else {
        correction =
            CGAL::ipower(Bound(1)/Bound(4), mult) * 
            CGAL::ipower(Bound(1)/Bound(2), CGAL::degree(res)-mult);
    }

    Bound val =  correction * CGAL::abs(res.evaluate(box.left()));
    CGAL::set_precision(BFI(), prec);
    BigFloat bfval2 = CGAL::lower(CGAL::convert_to_bfi(val));

    CGAL::set_precision(BFI(), oldp);
    return bfval2;
#endif
  }

  //! computes cofactor bounds for f & g if \c transposed = false
  //! of cofactor bounds for ft & gt if \c transposed = true
  std::pair< BigFloat, BigFloat > _cofactors_bound(const Root_box& box_y,
            bool transposed) const {

        //NOTE: try higher precision ? less error 
        long oldp = CGAL::set_precision(BFI(), 100);

        const Polynomial_2& f = (transposed ? this->_m_ft : this->_m_f),
                            g = (transposed ? this->_m_gt : this->_m_g);


        BigFloat xl = CGAL::lower(CGAL::convert_to_bfi(this->_m_box_x.left())),
                xh = CGAL::upper(CGAL::convert_to_bfi(this->_m_box_x.right())),
                yl = CGAL::lower(CGAL::convert_to_bfi(box_y.left())),
                yh = CGAL::upper(CGAL::convert_to_bfi(box_y.right()));
        if(transposed) { // exchange x & y intervals
            std::swap(xl, yl);
            std::swap(xh, yh);
        }

        typedef typename IA_real_bfi::Poly_2 Poly_2;
        typedef typename IA_real_bfi::Poly_1 Poly_1;

        BigFloat sf(0), sf_full, sg(0), sg_full;
        BigFloat l, h;

        int p = CGAL::degree(f), q = CGAL::degree(g), i;
        // Bisolve_out("bounds: f: " << f << "\n g: " << g << "\n");
        // Bisolve_out("p: " << p << "\n q: " << q << "\n");

        // use precached polynomials, no normalization
        _m_ia_real.set_polynomial(f, false, false);
        Poly_2 _f = _m_ia_real.internal_poly_2();
        typename Poly_2::const_iterator fxi;
        for(fxi = _f.begin(), i = 0; fxi != _f.end(); fxi++, i++) {
            // eval_range_RT_1_aff ??
            _m_ia_real.eval_range_AF1_1(*fxi, xl, xh, l, h);
            BigFloat t = CGAL::max(CGAL::abs(l), CGAL::abs(h));
            if(fxi == _f.begin())
                sf_full = CGAL::square (t); // the sum including constant coeff
            else
                sf += CGAL::square (t);  // the sum not including constant coeff
        }
        sf_full += sf;

        _m_ia_real.set_polynomial(g, false, false);
        Poly_2 _g = _m_ia_real.internal_poly_2();
        for(fxi = _g.begin(), i = 0; fxi != _g.end(); fxi++, i++) {
            _m_ia_real.eval_range_AF1_1(*fxi, xl, xh, l, h);
            BigFloat t = CGAL::max(CGAL::abs(l), CGAL::abs(h));

            if(fxi == _g.begin())
                sg_full = CGAL::square (t); // the sum including constant coeff
            else
                sg += CGAL::square (t); // the sum not including constant coeff
        }
        sg_full += sg;

        // f(x,y)*s(x,y) + g(x,y)*t(x,y) = res(f,g)
        // deg_y s < deg_y g; deg_y t < deg_y f

        int n = std::max(p, q);
        // we need n powers of max(|(yl,yh)^i|): i=0..n-1
        std::vector< BigFloat > y_exp(n);
        y_exp[0] = BigFloat(1);

        BigFloat l_even, h_even; // even exponent bounds
        BigFloat l_odd(yl), h_odd(yh); // odd exponent bounds
        BigFloat yl_sq = yl * yl, yh_sq = yh * yh,
            yl_evel_sq, yh_even_sq;

        l_even = 0, h_even = std::max(yl_sq, yh_sq);
        if(!(yl < 0 && yh > 0)) { // if interval does not contain zero
            l_even = std::min(yl_sq, yh_sq);
        }
        // this can differ from yl_sq & yh_sq
        yl_evel_sq = l_even, yh_even_sq = h_even;
        for(i = 1; i < n; i++) {
            if(i & 1) { // odd case
                y_exp[i] = std::max(CGAL::abs(l_odd), CGAL::abs(h_odd));
                l_odd *= yl_sq, h_odd *= yh_sq;
            } else {
                y_exp[i] = std::max(CGAL::abs(l_even), CGAL::abs(h_even));
                l_even *= yl_evel_sq, h_even *= yh_even_sq;
            }
        }

        BigFloat co_s_row = sf + BigFloat(1);
        for(i = 1; i < q; i++) {
            BigFloat t = y_exp[i];
            co_s_row *= (sf_full + CGAL::square (t));
        }

        if(p > 0)
            co_s_row *= CGAL::ipower(sg_full, p - 1) * sg;
        co_s_row = CGAL::sqrt(co_s_row);

        BigFloat co_t_row = sg + BigFloat(1);
        for(i = 1; i < p; i++) {
            BigFloat t = y_exp[i];
            co_t_row *= (sg_full + CGAL::square (t));
        }

        if(q > 0)
            co_t_row *= CGAL::ipower(sf_full, q - 1) * sf;
        co_t_row = CGAL::sqrt(co_t_row);

        CGAL::set_precision(BFI(), oldp);
        return std::make_pair(co_s_row, co_t_row);

    }

    //!@}
protected:

    //!\name data
    //!@{

    mutable IA_real_bfi _m_ia_real;

    //!@}

}; // Certifier_cofactor_traits

} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_CERTIFIER_COFACTOR_TRAITS_H
// EOF
