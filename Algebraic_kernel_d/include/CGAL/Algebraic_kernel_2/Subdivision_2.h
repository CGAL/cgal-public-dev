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
// $URL:$
// $Id: $
// 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>

/*! \file Subdivision_2
 *  \brief 2D curve plotting using recursive subdivision
 */

#ifndef CGAL_ALGEBRAIC_KERNEL_2_SUBDIVISION_2
#define CGAL_ALGEBRAIC_KERNEL_2_SUBDIVISION_2

#include <CGAL/config.h>
#include <CGAL/Algebraic_kernel_2/Bi_solve_2_flags.h>

#include <CGAL/Algebraic_kernel_2/Range_analysis_2.h>

namespace CGAL {

template < class PolyInt_2, class FloatCoeff >
class Subdivision_2
{
public: 
    //! \name public typedefs 
    //!@{ 
    //! this instance's template argument (concise explanation ??)
    typedef PolyInt_2 Poly_int_2; // integer polynomial ??
    //! this instance's template argument ??
    typedef FloatCoeff Float_coeff; // interval of doubles or BigFloat

    //! integer number type
    typedef typename CGAL::Polynomial_traits_d< Poly_int_2 >::
            Innermost_coefficient_type Integer;

    //! arithmetic kernel
    typedef typename CGAL::Get_arithmetic_kernel< Integer >::
         Arithmetic_kernel Arithmetic_kernel;
    //! rational NT
    typedef typename Arithmetic_kernel::Rational Rational;

protected:

    //! range analysis instance (\c false = no need for complex arithmetic)
//     typedef Range_analysis_2< Poly_int_2, Float_coeff, false >
//         IA_engine_real;

    typedef Range_analysis_2< Poly_int_2, Float_coeff, true >
             IA_engine;

    //! floating-point NT for internal use
    typedef typename IA_engine::Real NT;

    //! polynomials for internal use
    typedef typename IA_engine::Poly_1 Poly_1;
    typedef typename IA_engine::Poly_2 Poly_2;

    //! a rectangular region (quad)
//     typedef typename IA_engine::Quad Quad;
    typedef internal::Quad_ < NT > Quad;

    //!@}
public:
    //! \name Constructors
    //!@{ 
    //! default constructor
    Subdivision_2() : initialized(false), polynomial_set(false) {

        z_level = 1;
    } 

    //!@}
public:
    //! \name public methods
    //!@{ 

    //! specifies drawing window and pixel resolution
    void setup(const CGAL::Bbox_2& box_,
                int res_w_, int res_h_) {

        window = box_;
        x_min = static_cast<NT>(box_.xmin());
        y_min = static_cast<NT>(box_.ymin());
        x_max = static_cast<NT>(box_.xmax());
        y_max = static_cast<NT>(box_.ymax());
        res_w = res_w_, res_h = res_h_;
    
        if(x_min >= x_max||y_min >= y_max||res_w < 4||res_h < 4||res_w > 2048||
              res_h > 2048) {
            std::cerr << "Incorrect setup parameters" << std::endl;
            initialized = false;
            return;
        }
        pixel_w = (x_max - x_min) / res_w;
        pixel_h = (y_max - y_min) / res_h;
        initialized = true;
    }
    
    //! sets the curve equation
    void set_polynomial(const Poly_int_2& poly) {
        engine.set_polynomial(poly);
        polynomial_set = true;
    }
    
    //! returns the curve equation
    Poly_int_2 polynomial_2() const {
        return engine.polynomial_2();
    }

    template < class Coord_2 >
    void draw_quadtree(std::vector< Coord_2 >& pts) {
        
        if(!initialized || !polynomial_set)
            return;

        z_level = 1;
        std::cout << "resolution: " << res_w << " x " << res_h << std::endl;
        std::cout << "box: " << window << std::endl;
        points.clear();
        quad_tree(Quad(x_min, x_max, y_min, y_max));

        for(typeof(points.begin()) it = points.begin(); 
            it != points.end(); it++) {
            pts.push_back(Coord_2(it->first, it->second));
        }
    }

    //! same as before but uses 1D subdivision along x/y-axes instead
    template < class Coord_2 >
    void draw_subdiv_1(std::vector< Coord_2 >& pts) {
        
        if(!initialized || !polynomial_set)
            return;

        z_level = 0;
        std::cout << "resolution: " << res_w << " x " << res_h << std::endl;
        std::cout << "box: " << window << std::endl;
        points.clear();

        typedef CGAL::Polynomial_traits_d< Poly_2 > PT_2;

        Poly_2 f_y = engine.internal_poly_2(),
                f_x = typename PT_2::Swap()(f_y, 0, 1);

        what_direction = 0; // first subdivide in y-dir
        pixel_sz = pixel_h; // stopping criteria: pixel size in y-dir
        for(int ix = 0; ix < res_w; ix++) {
            NT x = x_min + ix * pixel_w;
            poly_subdiv1 = subs_inner_var(f_y, x);
            //TODO we do not need to set polynomial each step - only once
            // because binomials will be reused
            engine.set_polynomial(poly_subdiv1);
            coord = x;
            subdiv_1(y_min, y_max);
/*            if(!hit && std::abs(x) < 0.24){
            std::cerr << "FATAL: no points along the line: " << x << " (" <<
                    ix << ")\n";
            }*/
        }
//         what_direction = 1; // next in x-dir
//         pixel_sz = pixel_w;
//         for(int iy = 0; iy < res_h; iy++) {
//             NT y = y_min + iy * pixel_h;
// //TODO we do not need to set polynomial each step - only once
//             // because binomials will be reused
//             poly_subdiv1 = subs_inner_var(f_x, y);
//             coord = y;
//             subdiv_1(x_min, x_max);
//         }
        for(typeof(points.begin()) it = points.begin();
            it != points.end(); it++) {
            pts.push_back(Coord_2(it->first, it->second));
        }
    }

    void debug_run(
            double& csec_l, double& csec_h, double& csec_al, double& csec_ah) {
        typedef typename IA_engine::Complex_sector Complex_sector;
        typedef typename IA_engine::Real Real;
        typedef typename IA_engine::NT Complex;

        Poly_2 f_y = engine.internal_poly_2();
        NT x = x_min + 3 * pixel_w;
        poly_subdiv1 = subs_inner_var(f_y, x);
        engine.set_polynomial(poly_subdiv1);

        Real rad(2.0), a(rad / sqrt(2.0));
        Complex_sector s(1.0, rad, M_PI/4, M_PI/2), r;
        
        r = engine.eval_range_sector_1(poly_subdiv1, s);
        csec_l = r.magnitude().first;
        csec_h = r.magnitude().second;
        csec_al = r.angle().first;
        csec_ah = r.angle().second;
                
        Complex xlc(0, 0), xhc(a, a), lc, hc, le, he;
        engine.exact_range_1_complex(poly_subdiv1, xlc, xhc, le, he);
        engine.template eval_range_AF1_RT_1_aff< 2 >(xlc, xhc, lc, hc, 2);

        std::cout << "\n resulting sector: " << r << "\n";

        std::cout << "\n exact complex: [" << le << "; " << he << "]\n";
        std::cout << "\n AF1 complex: [" << lc << "; " << hc << "]\n";

// ====== adding Sector r: []; a: [90 (1.5708) - 180.021 (3.14195)]; and Sector r: [136.172; 136.172]; a: [180 (3.14159) - 180 (3.14159)]

        Complex_sector AA(0.0324272, 0.0324388, 1.5708, 3.14195);
        Complex_sector BB(136.172, 136.172, 3.14159, 3.14159);



//         Complex_sector AA(10.0, 25.0, 0, 135*M_PI/180);
//         Complex_sector BB(1, 3.0, 150*M_PI/180, 285*M_PI/180);
// result: r: [7; 27.908]; a: [342.542; 152.457] (512.458)

//         Complex_sector AA(1.0, 2.0, M_PI/6, M_PI/3);
//         Complex_sector BB(3, 5.0, 5*M_PI/4, 11*M_PI/8);
// result:  Sector r: [1; 4.25047]; a: [199.145 (3.47573) - 288.244 (5.03081)]

/*        Complex_sector AA(3.0, 3.5, 0, M_PI*2);
            Complex_sector BB(4.5, 5.0, M_PI/2, 3*M_PI/2);*/
// result: Sector r: [1; 8.5]; a: [38.9424 (0.679674) - 321.058 (5.60351)]

        BB = AA.add(BB);
        std::cout << "\n sector: " << BB << "\n"; 

//         Complex_sector aa(10, 20, 10*M_PI/180.0, 20*M_PI/180.0),
//             bb(10, 20, 130*M_PI/180.0, 140*M_PI/180.0),
//             cc(1, 2, 250*M_PI/180.0, 260*M_PI/180.0);
// 
//         Complex_sector rr;
//         rr = aa.add(bb), rr = rr.add(cc);
//         std::cout << "\n (a + b) + c: " << rr << "\n";
// 
//         rr = aa.add(cc), rr = rr.add(bb);
//         std::cout << "\n (a + c) + b: " << rr << "\n";
// 
//         rr = bb.add(aa), rr = rr.add(cc);
//         std::cout << "\n (b + a) + c: " << rr << "\n";
// 
//         rr = bb.add(cc), rr = rr.add(aa);
//         std::cout << "\n (b + c) + a: " << rr << "\n";
// 
//         rr = cc.add(aa), rr = rr.add(bb);
//         std::cout << "\n (c + a) + b: " << rr << "\n";
// 
//         rr = cc.add(bb), rr = rr.add(aa);
//         std::cout << "\n (c + b) + a: " << rr << "\n";
    }

    NT z_level;

    //!@}
protected:
    //!@{ 

    void quad_tree(const Quad& q) {

        if(!get_range_2(q))
            return;

        if(q.r - q.l <= pixel_w && q.t - q.b <= pixel_h) {
            NT x = q.l, y = q.b;
            int pix_x = static_cast<int>(
                    CGAL::to_double((x - x_min) / pixel_w)),
                pix_y = static_cast<int>(
                    CGAL::to_double((y - y_min) / pixel_h));
            points.push_back(std::make_pair(pix_x, res_h - pix_y));
            return;            
        } 
        NT xm = (q.l + q.r) / 2, ym = (q.b + q.t) / 2;
        quad_tree(Quad(q.l, xm, q.b, ym));
        quad_tree(Quad(q.l, xm, ym, q.t));
        quad_tree(Quad(xm, q.r, ym, q.t));
        quad_tree(Quad(xm, q.r, q.b, ym));
    }

    //! evaluates polynomial over quad \c q
    bool get_range_2(const Quad& q) {

//         NT l, h;
//         engine.eval_range_AF1_2(engine.get_poly_y(), q, l, h);
//         engine.eval_range_RT_2_aff(q, l, h);
//         return !(l * h > 0);
//          engine.eval_range_AF1_2(engine.poly_y_(), q, l, h);

#if 1
        typedef typename IA_engine::Quad QuadC;
        typedef typename IA_engine::NT Complex;
        Complex lc, hc;

        QuadC qc(Complex(q.l, z_level), Complex(q.r, z_level),
                    Complex(q.b, z_level), Complex(q.t, z_level));

//         engine.eval_range_AF1_2(engine.get_poly_y(), qc, lc, hc);
        engine.eval_range_RT_2_aff(qc, lc, hc);
        if(lc.real() * hc.real() > 0)
            return false;
        return true;
#endif
    }

    void subdiv_1(const NT& l, const NT& h) {

        if(!get_range_1(l, h))
            return;

        // TODO: use pixel_h for subdiv in y-direction (y_min - y_max)
        // and pixel_w for subdiv in x-direction
        if(h - l <= pixel_sz) {

            hit = true;

            NT x = coord, y = l;
            if(what_direction == 1) //! 1 indicates that we subdivide in x-dir
                std::swap(x, y);

            int pix_x = static_cast<int>(
                    CGAL::to_double((x - x_min) / pixel_w)),
                pix_y = static_cast<int>(
                    CGAL::to_double((y - y_min) / pixel_h));
            points.push_back(std::make_pair(pix_x, res_h - pix_y));
            return;
        }
        NT m = (l + h) / 2;
        subdiv_1(l, m);
        subdiv_1(m, h);
    }

    //! substitutes the inner variable in \c poly with \c x
    Poly_1 subs_inner_var(const Poly_2& poly, const NT& x) {

    //TODO: check if it makes any difference with eval_point_AF1_1
        typename Poly_1::Vector v(poly.degree() + 1);
        typename Poly_2::const_iterator fxi;
        unsigned i;

        for(fxi = poly.begin(), i = 0; fxi != poly.end(); i++, fxi++) {
            v[i] = engine.eval_point_1(*fxi, x).real();
//             typedef typename IA_engine::NT Complex;
//             Complex xl, xh;
//             engine.eval_point_AF1_12(*fxi, Complex(x), xl, xh);
//             v[i] = Float_coeff(xl.real(), xh.real());
        }
        return Poly_1(v.begin(), v.end());
    }

    

    bool get_range_1(const NT& xl, const NT& xh) {

        typedef typename IA_engine::NT Complex;
#if 1
    z_level = 0;
        
        Complex xlc(xl, z_level-0.04f), xhc(xh, z_level+0.2f), lc, hc;
//         Complex xlc(xl, z_level), xhc(xh, z_level), lc, hc;

        engine.template eval_range_AF1_RT_1_aff< 2 >(xlc, xhc, 
                lc, hc, 2);
//         engine.eval_range_AF1_RD_1(poly_subdiv1, xlc, xhc, lc, hc);
//         engine.eval_range_AF1_1(poly_subdiv1, xlc, xhc, lc, hc);

        Complex le, he;
        engine.exact_range_1_complex(poly_subdiv1, xlc, xhc, le, he);
//         engine.eval_range_AF1_RT_1_aff(poly_subdiv1, xlc, xhc, le, he, 6);
//         engine.eval_range_AF1_1< false >(poly_subdiv1, xlc, xhc, le, he);

        if(!(le.real() >= lc.real() && he.real() <= hc.real())) {

            std::cout << "\n FATAL (REAL): [" <<
                lc << "; " << hc << "] - [" << le << "; " << he << "]\n";
        }
        if(!(le.imag() >= lc.imag() && he.imag() <= hc.imag())) {

            std::cout << "\n FATAL (IMAG): [" <<
                lc << "; " << hc << "] - [" << le << "; " << he << "]\n";
        }
#endif    

/*         Complex dl = lc - le, dh = hc - he;

        // lc < le && hc > he
        if(dl.real() <= 0 && dl.imag() <= 0 && 
                dh.real() >= 0 && dh.imag() >= 0) {

            std::cout << "\n AF1 wins: [" <<
                dl << " x " << dh << "]\n";
            std::cout << "\n info: [" <<
                lc << " x " << hc << "] and [" << le << "; " << he << "\n";

        } else if(dl.real() >= 0 && dl.imag() >= 0 && 
                dh.real() <= 0 && dh.imag() <= 0) {*/
           
//             std::cout << "\n RT2 wins: [" <<
//                 dl << " x " << dh << "] on [" << xlc << "; " << xhc << "\n";
//              std::cout << "\n info RT3: [" <<
//                 lc << " x " << hc << "] and [" << le << "; " << he << "\n";

//         } else;
//             std::cout << "\n bounds overlap";

        if(lc.real() * hc.real() >= 0)
            return false;
        return true;
    }
    //!@}
protected:
    //! \name protected properties ??
    //!@{
    CGAL::Bbox_2 window;
    IA_engine engine; //! interval analysis

    NT x_min, x_max, y_min, y_max; //! drawing window boundaries
    int res_w, res_h;              //! pixel resolution
    NT pixel_w, pixel_h;           //! pixel dimensions w.r.t. resolution

    int what_direction;    //! indicates along which direction we currently
                           //! subdivide, x or y (applies to 1D subdivision)
    NT coord, pixel_sz;    //! current coordinate (x or y) &
                           //! minimal pixel size for 1D subdivision
    Poly_1 poly_subdiv1;   //! polynomial currently being used for 1D subdiv
    bool hit;

    //! points container
    std::vector< std::pair< int, int > > points;
    bool initialized, polynomial_set;
    //!@}

}; // class Subdivision_2

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_SUBDIVISION_2
// EOF
