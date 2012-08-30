// Copyright (c) 2009, 2010, 2011, 2012 Max-Planck-Institut fuer Informatik (Germany).
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
// $URL: svn+ssh://asm@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Algebraic_kernel_2/include/CGAL/Complex_sector.h $
// $Id: Complex_sector.h 54940 2010-03-26 10:33:49Z $
//
//
// Author(s): Pavel Emeliyanenko <asm@mpi-sb.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_2_COMPLEX_SECTOR_H
#define CGAL_ALGEBRAIC_KERNEL_2_COMPLEX_SECTOR_H

/*! \file Complex_sector.h
 * arithmetic on sectors in C specified by angles and magnitudes range
 */

#include <CGAL/config.h>

namespace CGAL {

namespace internal {

template < class NT_ >
struct Complex_sector_rep {

    //! this instance's second template parameter
    typedef NT_ NT;
    //! complex NT
    typedef std::complex< NT > Complex;

    //!\name Constructors
    //!@{

    //! default constructor
    Complex_sector_rep() {
    }

    //! from real interval
    Complex_sector_rep(const NT& l_, const NT& h_) :
            l(l_), h(h_) {
        
        if(l < 0) {
            l = -l, h = -h;
            al = M_PI;
        } else 
            al = 0;
        ah = al;
    }

    // from complex interval
    Complex_sector_rep(const Complex& l_, const Complex& h_) {
        throw "NYI";
    }

    // from a point in C
    Complex_sector_rep(const Complex& z) {
        h = l = std::abs(z); // works for doubles only
        al = ah = std::atan2(z.imag(), z.real());
    }

    // from angle & magnitude
    Complex_sector_rep(const NT& rl_, const NT& rh_, const NT& al_, 
            const NT& ah_) : l(rl_), h(rh_), al(al_), ah(ah_) {
    }
            
    //!@}
    //!\name props
    //!@{

    NT l, h;  // magnitude range
    NT al, ah;  // angle range

    //!@}
}; // class Complex_sector_rep

//! represents a sector in C: <tt> [rl,rh] *  exp(2*Pi*[al,ah]) </tt>
template < class NT_ >
class Complex_sector : 
        public CGAL::Handle_with_policy< Complex_sector_rep < NT_ > > {

public:
    //!\name seemingly typedefs
    //!@{

    //! this instance's first and only template argument
    typedef NT_ NT;
    //! rep class
    typedef Complex_sector_rep< NT > Rep;
    //! base class
    typedef CGAL::Handle_with_policy< Rep > Base;

    //! complex NT
    typedef typename Rep::Complex Complex;

    //!@}
public:
    //!\name Constructors
    //!@{

    //! default constructor
    Complex_sector() : 
        Base(Rep()) {
    }

    //! from a real interval
    Complex_sector(const NT& l, const NT& h) : Base(Rep(l, h)) {
    }
  
    //! from a complex interval
    Complex_sector(const Complex& l, const Complex& h) : Base(Rep(l, h)) {
    }

    //! from a complex point
    explicit Complex_sector(const Complex& z) : Base(Rep(z)) {
    }

    //! from magnitude and angle range
    Complex_sector(const NT& rl_, const NT& rh_, const NT& al_, const NT& ah_)
             : Base(Rep(rl_, rh_, al_, ah_)) {
    }
 
    //!@}
    //!\name methods
    //!@{

    using Base::ptr;

    inline NT r2d(NT r) const {
        return r*180.0 / M_PI;
    }

    // adds two sectors 
    Complex_sector add(const Complex_sector& s) const {

        std::cout << "====== adding " << *this << "; and " << s << "\n";

        Rep const *s1 = ptr(), *s2 = s.ptr();

        if(s1->h <= s2->l)
            std::swap(s1, s2);

        if(s2->h <= NT(1e-15)) {
            Complex_sector r(s1->l, s1->h, s1->al, s1->ah);
            return r;
        }

        NT phi_l, phi_h;
        std::cout << "subtract al: [" <<
            r2d(s1->al) << "; " << r2d(s1->ah) << "] minus [" <<
            r2d(s2->al) << "; " << r2d(s2->ah) << "]\n";

        phi_l = s1->al - s2->ah;
        phi_h = s1->ah - s2->al;

        _normalize_angle(phi_l, phi_h);
        std::cout << "normalized phi: " << r2d(phi_l) << "; " <<
                r2d(phi_h) << "\n";

        bool opposite_quads, equal;
        // map angle range to opposite quadrant
        NT a1l_rev, a1h_rev;
        if(s1->ah > 2*M_PI) {
            a1l_rev = s1->al - M_PI, a1h_rev = s1->ah - M_PI;
        } else {
            a1l_rev = s1->al + M_PI, a1h_rev = s1->ah + M_PI;
        }
        std::cout << "===== a1_rev: [" << a1l_rev << "; " <<
                a1h_rev << "]\n";

        opposite_quads = _angles_overlap(a1l_rev, a1h_rev, s2->al, s2->ah,
            equal);

        NT r_min, r_max;
        _magnitude(s1, s2, phi_l, phi_h, opposite_quads|equal, r_min, r_max);
     
        NT a_min, a_max;
        // intervals overlap: check for easy case
        if(opposite_quads && s1->l < s2->h &&  s2->l < s1->h) {
            std::cout << "opposite_quads & magnitude overlap\n";
            a_min = NT(0), a_max = NT(2*M_PI);
        
        } else { // handle easy case first 

            _angle(s1, s2, a_min, a_max);
std::cout << "===== before normalize: [" << a_min << "; " << a_max << "]\n";
            _normalize_angle(a_min, a_max, true);
std::cout << "after normalize: [" << a_min << "; " << a_max << "]\n";
        }

        Complex_sector r(r_min, r_max, a_min, a_max);
        std::cout << "result: " << r << "\n";
        if(r_min == 0 && r_max == 0) {
                std::cout << "FATAL\n"; throw -1;
        }

        return r;
    }

/**
    bool _atan_min(NT y1, x1, y2, x2) {

        s1 = (y1 >= 0), s2 = (y2 >= 0);
        if(s1 != s2) {
            if(s1)
                return 1;
            return 2;
        }
        // y's have the same sign
        if(s1) { // first two quadrants
            // angles increase when x decreases => min angle with max x
            (x1 >= x2) ? x1 : x2
        } else {
            // angles increase when x increases => min angle with min x
            (x1 <= x2) ? x1 : x2
        }
    }
*/

    // muls two sectors
    Complex_sector mul(const Complex_sector& s) const {

        // magnitudes are assumed to be nonnegative => multiplication is
        // trivial
// std::cout << "mul: " << *this << "; by: " << s << "\n";

        NT rl = ptr()->l * s.ptr()->l, rh = ptr()->h * s.ptr()->h; 
        NT ral = ptr()->al + s.ptr()->al, rah = ptr()->ah + s.ptr()->ah;
        _normalize_angle(ral, rah);

        return Complex_sector(rl, rh, ral, rah);
    }

    std::pair< NT, NT > magnitude() const {
        return std::make_pair(ptr()->l, ptr()->h);
    }

    std::pair< NT, NT > angle() const {
        return std::make_pair(ptr()->al, ptr()->ah);
    }

    //! adds a sector and a real interval
    Complex_sector add_interval(const NT& l, const NT& h) const {

        bool positive = (l >= 0); // we assume that zero is not crossed
        if(l < 0 && h > 0) {
            std::cout << "FATAL: interval with zero inclusion\n";   
            throw "done";
        }

        Complex_sector s(l, h);
        std::cout << "====== adding " << *this << "; and " << s << "\n";

        Rep const *s1 = ptr(), *s2 = s.ptr();
        bool swapped = (s1->h <= s2->l);

        if(swapped)
            std::swap(s1, s2);

        if(s2->h <= NT(1e-15)) {
            Complex_sector r(s1->l, s1->h, s1->al, s1->ah);
            return r;
        }

        NT phi_l, phi_h;
        std::cout << "subtract al: [" <<
            r2d(s1->al) << "; " << r2d(s1->ah) << "] minus [" <<
            r2d(s2->al) << "; " << r2d(s2->ah) << "]\n";

        phi_l = s1->al - s2->ah;
        phi_h = s1->ah - s2->al;

        _normalize_angle(phi_l, phi_h);
        std::cout << "normalized phi: " << r2d(phi_l) << "; " <<
                r2d(phi_h) << "\n";

        bool opposite_quads, equal;
        // map angle range to opposite quadrant
        NT a1l_rev, a1h_rev;
        if(s1->ah > 2*M_PI) {
            a1l_rev = s1->al - M_PI, a1h_rev = s1->ah - M_PI;
        } else {
            a1l_rev = s1->al + M_PI, a1h_rev = s1->ah + M_PI;
        }
        std::cout << "===== a1_rev: [" << a1l_rev << "; " <<
                a1h_rev << "]\n";

        opposite_quads = _angles_overlap(a1l_rev, a1h_rev, s2->al, s2->ah,
            equal);

        NT r_min, r_max;
        _magnitude(s1, s2, phi_l, phi_h, opposite_quads|equal, r_min, r_max);
     
        NT a_min, a_max;
        // intervals overlap: check for easy case
        if(opposite_quads && s1->l < s2->h &&  s2->l < s1->h) {
            std::cout << "opposite_quads & magnitude overlap\n";
            a_min = NT(0), a_max = NT(2*M_PI);
        
        } else { // handle easy case first 

            if(!swapped) // s2 is real interval
                _angle_real_case1(s1, s2, positive, a_min, a_max);
            else    // s1 is real interval
                _angle_real_case2(s1, positive, s2, a_min, a_max);
         
            _normalize_angle(a_min, a_max, true);
        }

        Complex_sector r(r_min, r_max, a_min, a_max);
        std::cout << "result: " << r << "\n";
        return r;
    }

    //!@}
protected:
    //!\name protected methods
    //!@{

    inline void _normalize_angle(NT& a_l, NT& a_h, 
            bool clamp_range = false) const {

        while(a_l < NT(0)) {
            a_l += 2 * M_PI;
            a_h += 2 * M_PI;
        }
        while(a_l >= 2 * M_PI) {
            a_l -= 2 * M_PI;
            a_h -= 2 * M_PI;
        }

        if(a_h - a_l > 2*M_PI) {
            a_l = 0, a_h = 2*M_PI; // clamp to the final range
        if(clamp_range)
        std::cout << "FATAL: clamped angle range..\n";
        }
    }

    //! check if \c [a1l;a1h] overlaps with \c [a2l;a2h]
    //! 0 <= a1l < 2*PI; 0 <= a1h < 4*PI
    inline bool _angles_overlap(NT a1l, NT a1h, NT a2l, NT a2h,
                bool& equal) const {

        NT _l = a1l - 2*M_PI, _h = a1h - 2*M_PI;

        equal = (a1l == a2h || a2l == a1h) || (_l == a2h || a2l == _h);
        if((a1l < a2h && a2l < a1h) || (_l < a2h && a2l < _h))
            return true;
        return false;
    }

    //! check if \c [al;ah] includes \c a
    //! 0 <= a < 2*PI
    inline bool _angle_includes(Rep const *s, NT a) const {

        if(a < 0)
            a += 2*M_PI;
        else if(a >= 2*M_PI)
            a -= 2*M_PI;

        if((s->al <= a && a <= s->ah) || (s->al - 2*M_PI <= a
                    && a <= s->ah - 2*M_PI))
            return true;
        return false;
    }

    //! magnitude of the sum of two sectors
    inline void _magnitude(Rep const *s1, Rep const *s2, NT phi_l, NT phi_h,
         bool opposite_quads, NT& r_min, NT& r_max) const {

        NT l1 = s1->l, h1 = s1->h, l2 = s2->l, h2 = s2->h;
        NT l1q = l1*l1, h1q = h1*h1, l2q = l2*l2, h2q = h2*h2;
        NT sin_l = std::sin(phi_l), sin_h = std::sin(phi_h);
        NT cos_l = NT(2) * std::cos(phi_l),
           cos_h = NT(2) * std::cos(phi_h);

        r_max = NT(0);

        bool equal, overlap;
        overlap = _angles_overlap(s1->al, s1->ah, s2->al, s2->ah, equal);
        if(overlap | equal) {
            r_max = h1 + h2;
            goto Lcompute_r_min;
        } 

        {
        NT r11 = h1q + h2q, r12 = h1q + l2q, r21 = l1q + h2q;
        // TODO: can use conditions on angle without computing sin
        // TODO: use sincos
        
        NT r_max1, r_max2;
        if(sin_h >= 0 || sin_l > 0) {
            NT r113 = r11 + h1 * h2 * cos_l,
               r123 = r12 + h1 * l2 * cos_l,
               r213 = r21 + l1 * h2 * cos_l;
            r_max1 = std::max(r113, r123);
            r_max1 = std::max(r_max1, r213);
            if(sin_h >= 0) {
                r_max = r_max1;
            }
        }
        if(sin_h < 0) {
            // we have to compute all left candidates here anyway
            NT r112 = r11 + h1 * h2 * cos_h,
                r122 = r12 + h1 * l2 * cos_h,
                r212 = r21 + l1 * h2 * cos_h;
            r_max2 = std::max(r112, r122);
            r_max2 = std::max(r_max2, r212);
        //NOTE: remember: r_max1/2 are squared magnitudes !!        
            if(sin_l <= 0)
                r_max = r_max2;
            else
                r_max = std::max(r_max1, r_max2);
        }
        r_max = std::sqrt(r_max); // because this is a square
        } // end of block

Lcompute_r_min:
        NT r_min1, r_min2;

        if(opposite_quads) {
            if(l1 <= h2 && l2 <= h1) { // magnitudes overlap => 0
                r_min = 0; // r111
            } else
                r_min = std::abs(l1 - h2); // r321: cos(phi) = pi
            return;
        } 

        NT r32 = l1q + h2q, r33 = l1q + l2q;
//NOTE that here magnitudes are squared !!!
        if(sin_h <= 0 || sin_l < 0) {
            NT r323 = r32 + l1 * h2 * cos_l,
               r333 = r33 + l1 * l2 * cos_l;
            r_min1 = std::min(r323, r333);

            NT r313(0), r133(0);
            if(cos_l < 0) {
                r313 = -l1 * cos_l * NT(0.5); // cos_l == 2 * cos(phi_l)
                if(l2 <= r313 && r313 <= h2) {
                    r313 = l1 * sin_l, r313 *= r313;
                    r_min1 = std::min(r_min1, r313);
                }

                r133 = -l2 * cos_l * NT(0.5);
                if(l1 <= r133 && r133 <= h1) {
                    r133 = l2 * sin_l, r133 *= r133;
                    r_min1 = std::min(r_min1, r133);
                }
            }

    printf("r323: %f; r333: %f; r313: %f; r133: %f\n",
            r323, r333, r313, r133);

            if(sin_h <= 0) {
                r_min = r_min1;
                goto Lend;
            }
        }

        { // sin_h > 0
        NT r322 = r32 + l1 * h2 * cos_h,
            r332 = r33 + l1 * l2 * cos_h;
        r_min2 = std::min(r322, r332);
        
        NT r312(0), r132(0);
        if(cos_h < 0) {
            r312 = -l1 * cos_h * NT(0.5); // cos_h == 2 * cos(phi_h)
            if(l2 <= r312 && r312 <= h2) {
                r312 = l1 * sin_h, r312 *= r312;
                r_min2 = std::min(r_min2, r312);
            }
            r132 = -l2 * cos_h * NT(0.5);
            if(l1 <= r132 && r132 <= h1) {
                r132 = l2 * sin_h, r132 *= r132;
                r_min2 = std::min(r_min2, r132);
            }
        }

printf("r322: %f; r332: %f; r312q: %f; r132q: %f\n",
            r322, r332, r312, r132);

        // NOTE: does it really make sense to check the conditions ??
        if(sin_l >= 0)
            r_min = r_min2;
        else
            r_min = std::min(r_min1, r_min2);
        }
Lend:   r_min = std::sqrt(r_min);
    }

    //! divides [l1; h1] by [l2; h2]
    //! assumes that intervals are positive
    inline void _ia_div(Rep const *s1, Rep const *s2, NT& lx, NT& hx,
            bool& infty) const {
    
        infty = false;
        lx = s1->l / s2->h;
        if(s2->l != NT(0)) {
            hx = s1->h / s2->l;
            return;
        } 
        hx = NT(1e100);
        infty = true;
    }

#define ATAN2(z, x, y) { \
        z = std::atan2(x,y); \
        if(z < 0) \
            z += 2*M_PI; \
        }

#define ANGLE_max(z, w) \
    (z >= w ? z : z + 2*M_PI)

#define ATAN2_max(z, x, y, w) { \
        z = std::atan2(x,y); \
        if(z < 0) \
            z += 2*M_PI; \
        z = ANGLE_max(z, w); \
        }  


#if 0
    //! \c infty is true: hx is at infinity
    inline void _angle(Rep const *s1, Rep const *s2, NT& a_min, NT& a_max)
                 const {

        bool infty;
        NT lx, hx;
        _ia_div(s1, s2, lx, hx, infty);

        if(lx < NT(1e-15)) {
            std::cout << "angle indeterminate PI\n";
            a_min = NT(0), a_max = NT(M_PI*2);
            return;
        } 

        NT a1_l = s1->al, a1_h = s1->ah, a2_l = s2->al, a2_h = s2->ah;
        NT sin_1l, cos_1l, sin_1h, cos_1h, sin_2l, cos_2l, sin_2h, cos_2h;

        nt_sincos(a1_l, sin_1l, cos_1l);
        nt_sincos(a1_h, sin_1h, cos_1h);
        nt_sincos(a2_l, sin_2l, cos_2l);
        nt_sincos(a2_h, sin_2h, cos_2h);

        std::cout << "a1: [" << r2d(a1_l) << "; " << r2d(a1_h) <<
            "]; a2: [" << r2d(a2_l) << "; " << r2d(a2_h) << "]; lx = " <<
                lx << "\n";

        NT a212(-1e4), a223(-1e5), a122(-1e6), a221, a222(-1e7);

        NT s = 1e4, acos(0);
        NT sin_212, cos_212;
        nt_sincos(a1_h - a2_l, sin_212, cos_212);
        bool a212_c = (sin_212 > 0) && (lx + cos_212 < 0) &&
                (lx * cos_212  + 1 > 0);
        ATAN2(a212, lx * sin_1h + sin_2l, lx * cos_1h + cos_2l)
        if(a212_c) {
            s = a212;
        }

        if(lx >= NT(1.0)) {

            acos = std::acos(-1.0 / lx);
            NT phi2 = a1_l + acos;

            if(_angle_includes(s2, phi2)) {
//TODO: can in fact rewrite this formula without using sin/cos
                ATAN2(a223, (lx * sin_1l + std::sin(phi2)),
                    (lx * cos_1l + std::cos(phi2)))
                s = std::min(s, a223);
            }
            phi2 = a1_l - acos;
            if(_angle_includes(s2, phi2)) {
                ATAN2(a223, (lx * sin_1l + std::sin(phi2)),
                    (lx * cos_1l + std::cos(phi2)))
                s = std::min(s, a223);
            }
        }

        NT u = s, v = s;
        bool a122_c = (hx * std::cos(a1_l - a2_l) + 1 > 0);

        if(a122_c) {
            if(!infty) {
                ATAN2(a122, hx * sin_1l + sin_2l, hx * cos_1l + cos_2l)
            } else {
                ATAN2(a122, sin_1l, cos_1l) 
            }
            u = std::min(u, a122);
        }
        
        NT sin_a2h, cos_a2h, sin_a2l, cos_a2l;        
        nt_sincos(a1_l - a2_h, sin_a2h, cos_a2h);
        nt_sincos(a1_l - a2_l, sin_a2l, cos_a2l);

        bool a222_c = (lx + cos_a2l > 0) && (lx * cos_a2l + 1 > 0);
        if(a222_c) {
            ATAN2(a222, lx * sin_1l + sin_2l, lx * cos_1l + cos_2l)
            v = std::min(v, a222);
        }
            
        a_min = (sin_a2l <= 0 ? u : v);
        ATAN2(a221, lx * sin_1l + sin_2h, lx * cos_1l + cos_2h)

        if(sin_a2h <= 0) {
            bool a121_c = (hx * cos_a2h + 1 < 0);
            if(a121_c) {
                NT a121;
                if(!infty)
                    ATAN2(a121, hx * sin_1l + sin_2h, hx * cos_1l + cos_2h)
                else // NOTE: this could also be zero
                    ATAN2(a121, sin_1l, cos_1l)
                a_min = std::min(a121, a_min);
            }
        } else {

            bool a221_c = (lx + cos_a2h > 0) && (lx * cos_a2h + 1 < 0);
            if(a221_c) {
                a_min = std::min(a221, a_min);
            }
        }

        NT a213;
        s = ANGLE_max(a221, a_min);

        if(lx >= NT(1.0)) {
            NT phi2 = a1_h + acos;

            if(_angle_includes(s2, phi2)) {
                ATAN2_max(a213, (lx * sin_1h + std::sin(phi2)),
                    (lx * cos_1h + std::cos(phi2)), a_min)
                s = std::max(s, a213);
            }

            phi2 = a1_h - acos;
            if(_angle_includes(s2, phi2)) {
                ATAN2_max(a213, (lx * sin_1h + std::sin(phi2)),
                    (lx * cos_1h + std::cos(phi2)), a_min)
                s = std::max(s, a213);
            }
        }

        NT a111;
        if(!infty)
            ATAN2_max(a111, hx * sin_1h + sin_2h, hx * cos_1h + cos_2h,
                a_min)
        else 
            ATAN2_max(a111, sin_1h, cos_1h, a_min) 

        NT a211;
        ATAN2_max(a211, lx * sin_1h + sin_2h, lx * cos_1h + cos_2h, a_min)
        u = std::max(a111, s), v = std::max(a211, s);
        
        sin_a2l = std::sin(a1_h - a2_l), sin_a2h = std::sin(a1_h - a2_h);
        a_max = (sin_a2h >= 0 ? u : v); 

        if(sin_a2l >= 0) {
            NT a112;
            if(!infty)
                ATAN2_max(a112, hx * sin_1h + sin_2l, hx * cos_1h + cos_2l,
                    a_min)
            else
                ATAN2_max(a112, sin_1h, cos_1h, a_min)
            a_max = std::max(a112, a_max);
        } else {
            a_max = std::max(ANGLE_max(a212, a_min), a_max);
        }

//      NT d=0;
//     printf("a212: %f; a223: %f; a122: %f; a222: %f; a221: %f; a213: %f; a211: %f"
//     "; a111: %f; a121: %f\n",
//       a212+d, a223+d, a122+d, a222+d, a221+d, a213+d, a211+d, a111+d, a121+d);

    }
#else

// NOTE NOTE NOTE: this is an old version
 //! \c infty is true: hx is at infinity
    inline void _angle(Rep const *s1, Rep const *s2, NT& a_min, NT& a_max)
                 const {

        bool infty;
        NT lx, hx;
        _ia_div(s1, s2, lx, hx, infty);

        if(lx < NT(1e-15)) {
            std::cout << "angle indeterminate PI\n";
            a_min = NT(0), a_max = NT(M_PI*2);
            return;
        } 

        NT a1_l = s1->al, a1_h = s1->ah, a2_l = s2->al, a2_h = s2->ah;
        NT sin_1l = std::sin(a1_l), cos_1l = std::cos(a1_l);
        NT sin_2l = std::sin(a2_l), cos_2l = std::cos(a2_l);
        NT sin_1h = std::sin(a1_h), cos_1h = std::cos(a1_h);
        NT sin_2h = std::sin(a2_h), cos_2h = std::cos(a2_h);

        std::cout << "a1: [" << r2d(a1_l) << "; " << r2d(a1_h) <<
            "]; a2: [" << r2d(a2_l) << "; " << r2d(a2_h) << "]; lx = " <<
                lx << "\n";

// NOTE: remember these are tangents of actual angles !!
        NT a212(0), a223(0), a122(0), a221, a222(0);

        NT s = 10*M_PI, acos(0);
        // TODO: it may be required to test other conditions if
        // this does not suffice
        bool a212_c = (std::sin(a1_h - a2_l) > 0) &&
            (lx + std::cos(a1_h - a2_l) < 0) &&
            (lx * std::cos(a1_h - a2_l) + 1 > 0);
        ATAN2(a212, lx * sin_1h + sin_2l, lx * cos_1h + cos_2l)

        std::cout << "a212: " << a212 << "\n";
        if(a212_c) {
            s = a212;
        }
        //NT sqx = std::sqrt(lx * lx - 1);

        if(lx >= NT(1.0)) { //! skip the acos test !!!
            acos = std::acos(-1.0 / lx);

            NT phi2 = a1_l + acos;
 std::cout << "phi_2: " << phi2 << "\n";
            if(_angle_includes(s2, phi2)) {
//TODO: can in fact rewrite this formula without using sin/cos
                ATAN2(a223, (lx * sin_1l + std::sin(phi2)),
                    (lx * cos_1l + std::cos(phi2)))
                s = std::min(s, a223);
            }
            phi2 = a1_l - acos;
std::cout << "phi_2: " << phi2 << "\n";
            if(_angle_includes(s2, phi2)) {
                ATAN2(a223, (lx * sin_1l + std::sin(phi2)),
                    (lx * cos_1l + std::cos(phi2)))
                s = std::min(s, a223);
            }
        }
//         ATAN2(a223, (lx * sin_1l + std::sin(phi2)),
//                 (lx * cos_1l + std::cos(phi2)))
//         std::cout << "phi2: " << r2d(phi2) << "; a223: " <<
//             a223 << "; " << r2d(a223) << "\n";

        NT u = s, v = s;
        bool a122_c = (hx * std::cos(a1_l - a2_l) + 1 > 0);
        std::cout << "a122_c : " << a122_c << "\n";

        if(a122_c) {
            if(!infty) {
                ATAN2(a122, hx * sin_1l + sin_2l, hx * cos_1l + cos_2l)
            } else {
                ATAN2(a122, sin_1l, cos_1l) // NOTE: this could also be zero
            }
            u = std::min(u, a122);
        }

        NT sin_a2h = std::sin(a1_l - a2_h), sin_a2l = std::sin(a1_l - a2_l);
        NT cos_a2h = std::cos(a1_l - a2_h), cos_a2l = std::cos(a1_l - a2_l);

//TODO: some of these conditions might be unnecessary...
        bool a222_c = (lx + cos_a2l > 0) &&
                (lx * cos_a2l + 1 > 0);
        if(a222_c) {
            ATAN2(a222, lx * sin_1l + sin_2l, lx * cos_1l + cos_2l)
            v = std::min(v, a222);
        }
        a_min = (sin_a2l <= 0 ? u : v);

        ATAN2(a221, lx * sin_1l + sin_2h, lx * cos_1l + cos_2h)

        if(sin_a2h <= 0) {
            bool a121_c = (hx * cos_a2h + 1 < 0);
            if(a121_c) {
                NT a121;
                if(!infty)
                    ATAN2(a121, hx * sin_1l + sin_2h, hx * cos_1l + cos_2h)
                else // NOTE: this could also be zero
                    ATAN2(a121, sin_1l, cos_1l)
                a_min = std::min(a121, a_min);
            }    
        } else {

            bool a221_c = (lx + cos_a2h > 0) && (lx * cos_a2h + 1 < 0);
            if(a221_c) {
                a_min = std::min(a221, a_min);
            }
        }

        NT a213;
        s = a221;
        if(s < a_min)
            s += 2*M_PI;

        if(lx >= NT(1.0)) {
            NT phi2 = a1_h + acos;

std::cout << "phi2_1: " << r2d(phi2) << "\n";
            if(_angle_includes(s2, phi2)) {
                ATAN2_max(a213, (lx * sin_1h + std::sin(phi2)),
                    (lx * cos_1h + std::cos(phi2)), a_min)
                s = std::max(s, a213);
            }

            phi2 = a1_h - acos;
std::cout << "phi2_2: " << r2d(phi2) << "\n";
            if(_angle_includes(s2, phi2)) {
                ATAN2_max(a213, (lx * sin_1h + std::sin(phi2)),
                    (lx * cos_1h + std::cos(phi2)), a_min)
                s = std::max(s, a213);
            }
        }

        NT a111;
        if(!infty)
            ATAN2_max(a111, hx * sin_1h + sin_2h, hx * cos_1h + cos_2h,
                a_min)
        else // NOTE: this could also be zero
            ATAN2_max(a111, sin_1h, cos_1h, a_min) 

        NT a211;
        ATAN2_max(a211, lx * sin_1h + sin_2h, lx * cos_1h + cos_2h, a_min)
        u = std::max(a111, s), v = std::max(a211, s);
        
        sin_a2l = std::sin(a1_h - a2_l), sin_a2h = std::sin(a1_h - a2_h);
        a_max = (sin_a2h >= 0 ? u : v); 

        if(sin_a2l >= 0) {
            NT a112;
            if(!infty)
                ATAN2_max(a112, hx * sin_1h + sin_2l, hx * cos_1h + cos_2l,
                    a_min)
            else
                ATAN2_max(a112, sin_1h, cos_1h, a_min)

            a_max = std::max(a112, a_max);

        } else {
            if(a212 < a_min)
                a212 += 2*M_PI;
            a_max = std::max(a212, a_max);
        }


NT a112;    if(!infty)
                ATAN2(a112, hx * sin_1h + sin_2l, hx * cos_1h + cos_2l)
            else
                ATAN2(a112, sin_1h, cos_1h)

 NT d=0;
NT a121;
 if(!infty)
       ATAN2(a121, hx * sin_1l + sin_2h, hx * cos_1l + cos_2h)
    else
      ATAN2(a121, sin_1l, cos_1l) // NOTE: this could also be zero


printf("a212: %f; a223: %f; a122: %f; a222: %f; a221: %f; a213: %f; a211: %f"
    "; a111: %f; a121: %f\n",
            a212+d, a223+d, a122+d, a222+d, a221+d, a213+d, a211+d, a111+d, a121+d);

// TODO: make sure that atan(+inf) = PI/2 !!
// std::cout << "******************before atan: " << a_min << "; " << a_max << "\n\n";
//         a_min = std::atan(a_min);
//         a_max = std::atan(a_max);
    }
#endif

  //! adds \c s1 and \c s2 at angle a = (0, 180) (\c positive )
    //! assumed that \c s1->h > \c s2->l
    inline void _angle_real_case1(Rep const *s1, Rep const *s2,
            bool positive, NT& a_min, NT& a_max) const {

        bool infty;
        NT lx, hx;
        _ia_div(s1, s2, lx, hx, infty);

        if(lx < NT(1e-15)) {
            std::cout << "angle indeterminate PI\n";
            a_min = NT(0), a_max = NT(M_PI*2);
            return;
        } 

    // cases xy1 and xy2 are identical
    // a2_l == a2_h
        NT a1_l = s1->al, a1_h = s1->ah;

        NT a2 = (positive ? NT(0) : NT(M_PI));
        NT sin_1l, cos_1l, sin_1h, cos_1h;
        NT sin_2(0), cos_2 = std::cos(a2);

        nt_sincos(a1_l, sin_1l, cos_1l);
        nt_sincos(a1_h, sin_1h, cos_1h);

        std::cout << "_angle_real_case1: [" << r2d(a1_l) << "; " << r2d(a1_h)
             << "]; a2: " << r2d(a2) << "; lx = " << lx << "\n";

        NT a212(-1e4), a122(-1e6), a221, a222(-1e7);

        NT s = 1e4, acos(0);
        NT sin_212, cos_212;
        nt_sincos(a1_h - a2, sin_212, cos_212);
        bool a212_c = (sin_212 > 0) && (lx + cos_212 < 0) && 
                (lx * cos_212 + 1 > 0);
        ATAN2(a212, lx * sin_1h + sin_2, lx * cos_1h + cos_2)

        if(a212_c) {
            s = a212;
        }

        NT u = s, v = s;
        if(!infty) {
            ATAN2(a122, hx * sin_1l + sin_2, hx * cos_1l + cos_2)
        } else {
            ATAN2(a122, sin_1l, cos_1l) 
        }
        u = std::min(u, a122);
    
        NT sin_a2lh, cos_a2lh;
        nt_sincos(a1_l - a2, sin_a2lh, cos_a2lh);

        bool a222_c = (lx + cos_a2lh > 0);
        ATAN2(a222, lx * sin_1l + sin_2, lx * cos_1l + cos_2)
        if(a222_c) {
            v = std::min(v, a222);
        }
            
        a_min = (sin_a2lh <= 0 ? u : v);
        s = ANGLE_max(a222, a_min);

        NT a111;
        if(!infty)
            ATAN2_max(a111, hx * sin_1h + sin_2, hx * cos_1h + cos_2,
                a_min)
        else 
            ATAN2_max(a111, sin_1h, cos_1h, a_min) 

        u = std::max(a111, s), v = std::max(ANGLE_max(a212, a_min), s);
        
        sin_a2lh = std::sin(a1_h - a2);
        a_max = (sin_a2lh >= 0 ? u : v); 
    }

    
    //! adds \c [l1,h1] at angle a = (0, 180) (\c positive ) and \c s2 
    //! assumed that \c s1->h > \c s2->l
    inline void _angle_real_case2(Rep const *s1, bool positive,
         Rep const *s2, NT& a_min, NT& a_max) const {

        NT lx, hx;
        bool infty;
        _ia_div(s1, s2, lx, hx, infty);

        if(lx < NT(1e-15)) {
            std::cout << "angle indeterminate PI\n";
            a_min = NT(0), a_max = NT(M_PI*2);
            return;
        } 

// a1_l == a1_h
// cases x1y and x2y are identical
        NT a1 = (positive ? NT(0) : NT(M_PI));
        NT a2_l = s2->al, a2_h = s2->ah;
        NT sin_2l, cos_2l, sin_2h, cos_2h;
        NT sin_1(0), cos_1 = std::cos(a1);

        nt_sincos(a2_l, sin_2l, cos_2l);
        nt_sincos(a2_h, sin_2h, cos_2h);

        std::cout << "_angle_real_case2: a1: " << a1 <<
            "; a2: [" << r2d(a2_l) << "; " << r2d(a2_h) << "]; lx = " <<
                lx << "\n";

        NT a212(-1e4), a223(-1e5), a122(-1e6), a221, a222(-1e7);

        NT s = 1e4, acos(0);
        NT sin_212, cos_212;
        nt_sincos(a1 - a2_l, sin_212, cos_212);
        bool a212_c = (sin_212 > 0) && (lx + cos_212 < 0) && 
                (lx * cos_212 + 1 > 0);
        
        ATAN2(a212, lx * sin_1 + sin_2l, lx * cos_1 + cos_2l)
        if(a212_c) {
            s = a212;
        }

        bool a223_c = false;
        if(lx >= NT(1.0)) {

            acos = std::acos(-1.0 / lx);
            NT phi2 = a1 + acos;

            if(_angle_includes(s2, phi2)) {
                ATAN2(a223, (lx * sin_1 + std::sin(phi2)),
                    (lx * cos_1 + std::cos(phi2)))
                a223_c = true;
                s = std::min(s, a223);
            }
            phi2 = a1 - acos;
            if(_angle_includes(s2, phi2)) {
                ATAN2(a223, (lx * sin_1 + std::sin(phi2)),
                    (lx * cos_1 + std::cos(phi2)))
                a223_c = true;
                s = std::min(s, a223);
            }
        }

        NT u = s, v = s;
        bool a122_c = (hx * std::cos(a1 - a2_l) + 1 > 0);

        if(!infty) {
            ATAN2(a122, hx * sin_1 + sin_2l, hx * cos_1 + cos_2l)
        } else {
            ATAN2(a122, sin_1, cos_1) 
        }

        if(a122_c) {
            u = std::min(u, a122);
        }

        NT sin_a2h, cos_a2h, sin_a2l, cos_a2l;        
        nt_sincos(a1 - a2_h, sin_a2h, cos_a2h);
        nt_sincos(a1 - a2_l, sin_a2l, cos_a2l);

        ATAN2(a222, lx * sin_1 + sin_2l, lx * cos_1 + cos_2l)
        bool a222_c = (lx + cos_a2l > 0) && (lx * cos_a2l + 1 > 0);
        if(a222_c) {
            v = std::min(v, a222);
        }
            
        a_min = (sin_a2l <= 0 ? u : v);
        ATAN2(a221, lx * sin_1 + sin_2h, lx * cos_1 + cos_2h)

        NT a121;
        if(!infty)
            ATAN2(a121, hx * sin_1 + sin_2h, hx * cos_1 + cos_2h)
        else 
            ATAN2(a121, sin_1, cos_1)

        if(sin_a2h <= 0) {
            bool a121_c = (hx * cos_a2h + 1 < 0);
            if(a121_c) 
                a_min = std::min(a121, a_min);

        } else {
            bool a221_c = (lx + cos_a2h > 0) && (lx * cos_a2h + 1 < 0);
            if(a221_c) {
                a_min = std::min(a221, a_min);
            }
        }

        s = ANGLE_max(a221, a_min);

        if(a223_c) {
            s = std::max(s, ANGLE_max(a223, a_min));
        }

        u = std::max(ANGLE_max(a121, a_min), s);
        v = std::max(ANGLE_max(a221, a_min), s);
        a_max = (sin_a2h >= 0 ? u : v); 

        if(sin_a2l >= 0) {
            a_max = std::max(ANGLE_max(a122, a_min), a_max);
        } else {
            a_max = std::max(ANGLE_max(a222, a_min), a_max);
        }
    }

public:
    void write(std::ostream& os) const {

        NT al_grad = ptr()->al * 180 / M_PI,
            ah_grad = ptr()->ah * 180 / M_PI;
        os << "Sector r: [" << ptr()->l << ", " << ptr()->h <<
            "]; a: [" << ptr()->al << ", " << ptr()->ah << " (" << 
            al_grad << ", " << ah_grad << ")]";
    }

    //!@}
};

template < class NT >
inline void nt_sincos(const NT& a, NT& s, NT& c) {
    throw "NYI";
}

template < >
inline void nt_sincos< double >(const double& a, double& s, double& c) {
// use builtin sincos
    //::sincos(a, &s, &c);
}

//! output
template < class NT_ >
std::ostream& operator<< (std::ostream& os, const Complex_sector< NT_ >& s) {
    s.write(os);
    return os;
}

} // namespace internal

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_2_COMPLEX_SECTOR_H
// EOF
