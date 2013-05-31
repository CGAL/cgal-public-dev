// Copyright (c) 2010, 2011, 2012 Max-Planck-Institut fuer Informatik (Germany).
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
// $URL: svn+ssh://asm@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Algebraic_kernel_2/include/CGAL/Certifier_cofactor_traits.h $
// $Id: Certifier_cofactor_traits.h 57108 2010-06-25 13:54:53Z eric $
//
//
// Author(s): Eric Berberich <eric.berberich@cgal.org>
//            Pavel Emeliyanenko <asm@mpi-inf.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_2_CERTIFIER_SUBDIVISION_TRAITS_H
#define CGAL_ALGEBRAIC_KERNEL_2_CERTIFIER_SUBDIVISION_TRAITS_H

/*! \file
 * Traits class implementing subdivision root certification
 */

#include <CGAL/config.h>

#include <CGAL/Algebraic_kernel_2/Interval_nt.h>
#include <CGAL/Algebraic_kernel_2/Subdivision_knot.h>
#include <CGAL/Algebraic_kernel_2/Range_analysis_2.h>

#include <CGAL/Algebraic_kernel_2/Certifier_traits_base.h>

#include <CGAL/Algebraic_kernel_2/Bi_algebraic_real_2.h>

namespace CGAL {

namespace internal {

template < class AlgebraicKernel_d_1 >
class Certifier_subdivision_traits :
        public Certifier_traits_base < AlgebraicKernel_d_1 > {

public:

    //! template parameter
    typedef AlgebraicKernel_d_1 Algebraic_kernel_d_1;

    //! type of solution 
    typedef Bi_algebraic_real_2< Algebraic_kernel_d_1 > Algebraic_real_2;

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
    typedef typename CGAL::Get_arithmetic_kernel< Bound >::Arithmetic_kernel AK;

    //! bigfloat interval 
    typedef typename AK::Bigfloat_interval BFI;

    //! our lovely bigfloats
    typedef typename CGAL::Bigfloat_interval_traits< BFI >::Bound BigFloat;

    //! a root box
    typedef typename Base::Root_box Root_box;

public:
    //!\name Constuctors 
    //!@{

    //! default
    Certifier_subdivision_traits() :
        Base() {
    }

    //! standard
    Certifier_subdivision_traits(
                                 Algebraic_kernel_d_1 *ak,
                                 const Polynomial_2& f, const Polynomial_2& g,
                                 const Polynomial_2& ft, const Polynomial_2& gt,
                                 const Polynomial_1& res_x, const Polynomial_1& res_y) :
        Base(ak, f, g, ft, gt, res_x, res_y) {

        _m_f_ia_d.set_polynomial(f);
        _m_g_ia_d.set_polynomial(g);   
    }

    //!@}
protected:
    //!@{

    // 'Dimension_index' indicates which dimension of the boundary quad is
    // "collapsed": each quad belongs to one of 8 trees originated from the
    // respective boundary of C^2 cube. Each tree keeps an associated boundary
    // value which corresponds to collapsed dimension.
    // Boundary quads are treated uniformly as 3D cubes during subdivision.
    // The "missing" dimension is only added immediately before applying the
    // interval arithmetic.
    enum Dimension_index { 
        BOTTOM_R = 0, BOTTOM_I,
        TOP_R, TOP_I,
        LEFT_R, LEFT_I,
        RIGHT_R, RIGHT_I
    };

    //! floating-point NT to represent polynomial coefficients
    typedef CGAL::Interval_nt< true > Float_coeff;

    //! type of range analysis with double-precision
    typedef CGAL::Range_analysis_2< Polynomial_2, Float_coeff, true >
      IA_engine_double;

    typedef CGAL::Range_analysis_2< Polynomial_2, BFI, false >
        IA_real_bfi;

    //! versatile data structure
    typedef Subdivision_knot< Bound, BFI > Knot;

    struct Bound_quad { // describes a boundary quad
        Bound_quad() { }
        Bound_quad(const Knot& _1, const Knot& _2, const Knot& _3,
                   Dimension_index _idx) : 
          a(_1), b(_2), c(_3), 
          idx(_idx), norm(0) {  }
      
      Knot a, b, c; // 3 associated knots (6 numbers)

      Dimension_index idx;     // index of "collapsed" boundary
      
      BigFloat norm;           // norms previously computed
    };
    
    //! a list of quads
    typedef std::list< Bound_quad > Quad_list;
    //! list iterator
    typedef typename Quad_list::iterator QL_iterator;

    //!@}
public:
    //! \name data to be attached to the candidate
    //!@{

    struct Candidate_props {

        Candidate_props() :
            subdiv_dim(0) {
        }

        mutable Bound x0, y0;

        mutable BigFloat inorm; // norm in the interior (for convenience)

        //! a set of quads requiring further subdivision
        mutable Quad_list quads;

        mutable Knot l, r, b, t; // keep boundary apprximators

        mutable unsigned subdiv_dim; // indicates along what dimension we are
                                    // currently subdividing
    };

    //!@}
    //!\name main stuff
    //!@{

    //! returns \c true if the candidate passes the norm test
    bool norm_test(const Candidate_props& ps,
            const Bound& x_0, const Bound& y_0,
            const Root_box& box_y, Multiplicity_type mul_y, long prec) {

        prec = 100;
        long oldp = CGAL::set_precision(BFI(), prec);

        ps.x0 = x_0, ps.y0 = y_0;

        BFI x0f = CGAL::convert_to_bfi(x_0),
            y0f = CGAL::convert_to_bfi(y_0);

        typename IA_real_bfi::NT xfc(CGAL::median(x0f)),
                yfc(CGAL::median(y0f)), lc, hc;

        _m_ia_bfi.set_polynomial(this->_m_f, true, false);
        _m_ia_bfi.eval_point_AF1_2(_m_ia_bfi.internal_poly_2(),
                xfc, yfc, lc, hc);

        // we take the upper boundary because the norm test succeeds if 
        // f_boundary > f_interior

        ps.inorm = 0;
        _I_norm(lc, hc, ps.inorm);

        // use precached polynomials
        _m_ia_bfi.set_polynomial(this->_m_g, true, false);
        _m_ia_bfi.eval_point_AF1_2(_m_ia_bfi.internal_poly_2(),
            xfc, yfc, lc, hc);

        _I_norm(lc, hc, ps.inorm);

//    if (ps.inorm == 0) {
//      typename CGAL::Polynomial_traits_d< Polynomial_2 >::Substitute
//        substitute_2;
//      
//      std::list< Bound > point;
//      point.push_back(x_0);
//      point.push_back(y_0);
//      Bound f_0 = CGAL::abs(substitute_2(_m_f, point.begin(), point.end()));
//      Bound g_0 = CGAL::abs(substitute_2(_m_g, point.begin(), point.end()));
//      
//      ps.inorm = std::max(f_0, g_0);
//    }

        Bisolve_out("inorm: " << ps.inorm);

        if(_evaluate_norm(ps, box_y)) {
            std::cout << "OOOOOOOOOOOOOOOOOOOO candidate certified..\n";
            CGAL::set_precision(BFI(), oldp);
            return true;
        }

        CGAL::set_precision(BFI(), oldp);

        Bisolve_out("total " << ps.quads.size() <<
            " quads left on the list..\n");

        return false;
    }

    //! subdivides the boundary quads associated with a \c box discarding
    //! those for which the boundary test (|L; U| > \c inorm )
    //! succeeds. Returns \c true if all quads are discarded.
    bool _evaluate_norm(const Candidate_props& ps, const Root_box& box_y) {

        Bisolve_out("************* interior norms: " << ps.inorm << "\n");

        // if no quads are yet setup: subdivide original boundary & quit
        if(ps.quads.empty()) {
            _setup_boundary_quads(ps, box_y);
            // test succeeds if all initial quads proceeded norm test
            return ps.quads.empty();
        }

        // ..and here the play begins
        unsigned i = 0, nzero_quads = 0;
        for(QL_iterator qi = ps.quads.begin(); qi != ps.quads.end(); i++) {
            // test succeeds if norm on the boundary is greater
            // than the norm in the interior
            if(qi->norm <= ps.inorm) {
                unsigned n = 0;
                if(qi->norm == 0) { // subdivide only in case of zero norm
                    qi = _box_subdivide(ps, qi, n);
                    nzero_quads++;
                } else
                    qi++; // just leave this quad as it is

            } else {
                qi = ps.quads.erase(qi); // otherwise just remove the box
            }
        }

        std::cout << "quads with zero norm left: " << nzero_quads << "\n";
        ps.subdiv_dim ^= 1; // alternate subdivision direction
        // norm test succeeds if no quads left on the list
        return ps.quads.empty();
    }

    //! sets up the initial lists of boxes for subdivision
    void _setup_boundary_quads(const Candidate_props& ps,
                const Root_box& box_y) const {

        ps.quads.clear();
        QL_iterator qi = ps.quads.begin();

        typedef std::complex< Bound > C_bound;
        C_bound mx((this->_m_box_x.left + this->_m_box_x.right) / Bound(2),
               Bound(0)), my((box_y.left + box_y.right) / Bound(2), Bound(0));
        C_bound radx((this->_m_box_x.right - this->_m_box_x.left) / Bound(2),
                    this->_m_box_x.irad),
                rady((box_y.right - box_y.left) / Bound(2), box_y.irad);

    // [left x right] - [bottom x top]
    // (l_R l_I) (r_R r_I) - (b_R b_I) (t_R t_I)
    //! (l_R l_I) (r_R r_I) - (b_R b_I) (t_R b_I) BOTTOM_I
    //! keep: l, r, b_R, t_R
    //! (l_R l_I) (r_R r_I) - (b_R t_I) (t_R t_I) TOP_I
    //! keep: l, r, b_R, t_R
    //! (l_R l_I) (r_R r_I) - (b_R b_I) (b_R t_I) BOTTOM_R
    //! keep: l, r, b_I, t_I
    //! (l_R l_I) (r_R r_I) - (t_R b_I) (t_R t_I) TOP_R
    //! keep: l, r, b_I, t_I

    //! (l_R l_I) (r_R l_I) - (b_R b_I) (t_R t_I) LEFT_I
    //! keep: l_R, r_R, b, t
    //! (l_R r_I) (r_R r_I) - (b_R b_I) (t_R t_I) RIGHT_I
    //! keep: l_R, r_R, b, t
    //! (l_R l_I) (l_R r_I) - (b_R b_I) (t_R t_I) LEFT_R
    //! keep: l_I, r_I, b, t
    //! (r_R l_I) (r_R r_I) - (b_R b_I) (t_R t_I) RIGHT_R
    //! keep: l_I, r_I, b, t

        typedef Bound NT;
        NT x0 = mx.real(), y0 = my.real();
//         NT x0 = ps.x0, y0 = ps.y0;

        Knot l(mx-radx), r(mx+radx), b(my-rady), t(my+rady);

        ps.l = l, ps.r = r, ps.b = b, ps.t = t;
        Knot bt_r(b.z().real(), t.z().real()),
             bt_i(b.z().imag(), t.z().imag()),
             lr_r(l.z().real(), r.z().real()),
             lr_i(l.z().imag(), r.z().imag());

  // mx.imag() == my.imag() == 0
# if 0
        typedef Polynomial_2 Poly_2;
        typedef Polynomial_1 Poly_1;
        
        Poly_2 ff = _m_f,//_m_f_ia_d.internal_poly_2(),
                gg = _m_g;//_ia_d.internal_poly_2();

        typedef CGAL::Polynomial_traits_d< Poly_2 > PT_2;
        typename PT_2::Substitute subs;

        typedef typename Polynomial_type_generator< Bound, 2 >::Type 
            Poly_bound_2;
        typedef typename Poly_bound_2::NT Poly_bound_1;

//         typedef CGAL::Interval_nt< true> II;
        typedef Bound II;
        Poly_bound_2 pxm(Poly_bound_1(II(x0), II(1.0))),
               pym(Poly_bound_1(II(y0)),Poly_bound_1(II(1.0)));

        Poly_bound_2 terms[2] = {pxm, pym};
//         std::cout << pxm << "; " << pym << "---poly = " << poly << "\n";

        Poly_bound_2 ff_shift = subs(ff, terms, terms+2),
            gg_shift = subs(gg, terms, terms+2); 

        typedef CGAL::Fraction_traits< Poly_bound_2 > FT;

        Polynomial_2 ff2, gg2;
        typename FT::Denominator_type den;
        typename FT::Decompose()(ff_shift, ff2, den);
        typename FT::Decompose()(gg_shift, gg2, den);

//         std::cout << "====== shifted f: " << ff2 << "\n";
//         std::cout << "====== shifted g: " << gg2 << "\n";

        _m_f_ia_d.set_polynomial(ff2, false);
        _m_g_ia_d.set_polynomial(gg2, false);

        std::cout << "====== shifted f: " << _m_f_ia_d.internal_poly_2() << "\n";
        std::cout << "====== shifted g: " << _m_g_ia_d.internal_poly_2() << "\n";

        typedef typename IA_engine_double::Complex_sector Complex_sector;

        double ax, ay;
        unsigned i, j;
    
        double phi = 8*M_PI/180;

        double minr(1e10), maxr(0);

//         for(phi = 0; phi <= M_PI/2; phi += M_PI/180.0) {

        double rx = CGAL::to_double(radx.real()),
                ry = CGAL::to_double(rady.real());
        rx *= std::cos(phi), ry *= std::sin(phi);
            
        unsigned nx = 10, ny = 10;
        double dax = 2*M_PI/nx, day = 2*M_PI/ny;
        for(ax = 0.0, i = 0; i < nx; ax += dax, i++) {
            for(ay = 0.0, j = 0; j < ny; ay += day, j++) {
        

        Complex_sector sx(rx, rx, ax, ax+dax),
            sy(ry, ry, ay, ay+day), s;

         s = _m_f_ia_d.eval_boundary_sector_2(_m_f_ia_d.internal_poly_2(),
               sx, sy);
            minr = std::min(minr, s.magnitude().first);
            maxr = std::max(maxr, s.magnitude().second);
     std::cerr << "phi: " << (180/M_PI*phi) << "; ax: " << ax << "; ay: " << ay << "; " << s << "\n";
       }
        }
         }
#endif

        _add_quad(ps, qi, Bound_quad(l, r, bt_r, BOTTOM_I));
        _add_quad(ps, qi, Bound_quad(l, r, bt_r, TOP_I));
        _add_quad(ps, qi, Bound_quad(l, r, bt_i, BOTTOM_R));
        _add_quad(ps, qi, Bound_quad(l, r, bt_i, TOP_R));

        _add_quad(ps, qi, Bound_quad(b, t, lr_r, LEFT_I));
        _add_quad(ps, qi, Bound_quad(b, t, lr_r, RIGHT_I));
        _add_quad(ps, qi, Bound_quad(b, t, lr_i, LEFT_R));
        _add_quad(ps, qi, Bound_quad(b, t, lr_i, RIGHT_R));
    }
     
    //! subdivides a quad pointed to by \c qi, and replaces it by those
    //! subquads which failed the norm test. return an interator pointed
    //! to a quad immediately after \c qi (before subdivision)
    //! \c n indicates how many new quads have been added
    QL_iterator _box_subdivide(const Candidate_props& ps, QL_iterator qi,
              unsigned& n) const {

        n = 0;
        const Bound_quad& q = *qi;
        // subdivision should proceed in a uniform way: treating
        // all boxes as 3D cubes (a, b, c)
        if(ps.subdiv_dim == 0) { // do 2D subdivision along (a; b)

            Knot m((q.a.z() + q.b.z()) * Bound(0.5));
            n += _add_quad(ps, qi, Bound_quad(q.a, m, q.c, q.idx));
            n += _add_quad(ps, qi, Bound_quad(m, q.b, q.c, q.idx));
            n += _add_quad(ps, qi, Bound_quad(
                Knot(m.z().real(), q.a.z().imag()),
                Knot(q.b.z().real(), m.z().imag()), q.c, q.idx));
            n += _add_quad(ps, qi, Bound_quad(
                Knot(q.a.z().real(), m.z().imag()),
                Knot(m.z().real(), q.b.z().imag()), q.c, q.idx));

        } else { // 1D subdivision on c
            Bound m((q.c.z().real() + q.c.z().imag()) * Bound(0.5));
            // here real part is first boundary, imag - second boundary
            n += _add_quad(ps, qi, Bound_quad(q.a, q.b,
                        Knot(q.c.z().real(), m), q.idx));
            n += _add_quad(ps, qi, Bound_quad(q.a, q.b,
                        Knot(m, q.c.z().imag()), q.idx));
        }
//     Bisolve_out("%%%%%%%%% number of quads left: " << n << "\n");
        return ps.quads.erase(qi);
    }

    //! norm on interval, returns \c true if [l; h] does not contain zero
    //! (test succeeds)
    template < class NT > inline
    bool _I_norm(const NT& l, const NT& h, NT& res) const {

        int s1 = CGAL::sign(l), s2 = CGAL::sign(h);
        if(s1 * s2 <= 0) 
            return false;

        if(s1 > 0) // positive interval: |l| is minimal
            res = std::max(res, CGAL::abs(l));
        else       // negative interval: |h| is minimal
            res = std::max(res, CGAL::abs(h));
        return true;
    }

    //! checks the quad \c q and adds it to the list before \c qi if the
    //! norm test failed. returns \c true if insertion takes place
    bool _add_quad(const Candidate_props& ps, QL_iterator qi,
            const Bound_quad& q) const {
        
        typedef typename IA_engine_double::NT C_double;
        typename IA_engine_double::Quad qd;

        if(q.idx < 4) { // BOTTOM / TOP dimension

            qd.l = q.a.zd_low(), qd.r = q.b.zd_high();
            const Knot& x = (q.idx == BOTTOM_I || q.idx == BOTTOM_R ? ps.b :
                 ps.t);

            if(q.idx == BOTTOM_I || q.idx == TOP_I) {
                // q.c is (b_r; t_r): bottom rounded down, top rounded up
                qd.b = C_double(q.c.zd_low().real(), x.zd_low().imag());
                qd.t = C_double(q.c.zd_high().imag(), x.zd_high().imag());
            } else {
                // q.c is (b_i; t_i): bottom rounded down, top rounded up
                qd.b = C_double(x.zd_low().real(), q.c.zd_low().real());
                qd.t = C_double(x.zd_high().real(), q.c.zd_high().imag());
            }        
        } else { // LEFT / RIGHT dimension

            qd.b = q.a.zd_low(), qd.t = q.b.zd_high();
            const Knot& x = (q.idx == LEFT_I || q.idx == LEFT_R ? ps.l :
                 ps.r);

            if(q.idx == LEFT_I || q.idx == RIGHT_I) {
                // q.c is (l_r; r_r): left rounded down, right rounded up
                qd.l = C_double(q.c.zd_low().real(), x.zd_low().imag());
                qd.r = C_double(q.c.zd_high().imag(), x.zd_high().imag());
            } else {
                // q.c is (l_i; r_i): left rounded down, right rounded up
                qd.l = C_double(x.zd_low().real(), q.c.zd_low().real());
                qd.r = C_double(x.zd_high().real(), q.c.zd_high().imag());
            }
        }

        double q_norm(0);
        bool succeeds = false;

//         Bisolve_out("checking quad: " << qd.l << "; " << qd.b << "\n");
        Bisolve_out("checking quad (" << q.idx << ") " << (qd.r - qd.l) <<
                " x " << (qd.t - qd.b) << "\n")

        C_double l, h;
         _m_f_ia_d.eval_range_RT_2_aff(qd, l, h);
//         _m_f_ia_d.exact_range_2_complex(_m_f_ia_d.internal_poly_2(), qd, l, h);
        Bisolve_out("boundary f-int: [" << l << "; " << h << "]\n");

        succeeds |= _I_norm(l.real(), h.real(), q_norm);
        succeeds |= _I_norm(l.imag(), h.imag(), q_norm);

        _m_g_ia_d.eval_range_RT_2_aff(qd, l, h);
//         _m_g_ia_d.exact_range_2_complex(_m_g_ia_d.internal_poly_2(), qd, l, h);
        Bisolve_out("boundary g-int: [" << l << "; " << h << "]\n");

        succeeds |= _I_norm(l.real(), h.real(), q_norm);
        succeeds |= _I_norm(l.imag(), h.imag(), q_norm);

        Bisolve_out("fg-norm: " << q_norm << ", inorm: " << ps.inorm << "\n")

        if (succeeds) {
          Bisolve_out("Quad's norm non-zero\n")
          if (q_norm > ps.inorm) {// boundary norm is larger -> discard this
            Bisolve_out("Quad's norm larger than interior norm\n")
            return false;
          }
        }

        QL_iterator it = ps.quads.insert(qi, q);
        it->norm = q_norm; // save precomputed norm
        Bisolve_out("Quad will be subdivided\n")
        return true;
    }

    //!@}
protected:
    //!\name data
    //!@{

    mutable IA_engine_double _m_f_ia_d, _m_g_ia_d;

    mutable IA_real_bfi _m_ia_bfi;

    //!@}

}; // Certifier_subdivision_traits

} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_CERTIFIER_SUBDIVISION_TRAITS_H
// EOF
