//    (c) 2007-2009 National and Kapodistrian University of Athens
//    (c) 2009-2011 INRIA Nancy
//    (c) 2011-2012 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED

#ifndef CGAL_ELLIPSE_TRAITS_H
#define CGAL_ELLIPSE_TRAITS_H

#include<iostream>
#include<CGAL/basic.h>
#include<CGAL/GMP_arithmetic_kernel.h>
//#include<CGAL/Arithmetic_kernel.h>
#include<CGAL/Algebraic_kernel_d_1.h>
#include<CGAL/Gmpz.h>
#include<CGAL/Polynomial.h>
//#include<CGAL/Polynomial_traits_d.h>
#include<CGAL/Fraction_traits.h>
#include <CGAL/Range.h>
#include <CGAL/vorell.h>

#include<CGAL/Coercion_traits.h>

#include<CGAL/Cartesian.h>

namespace CGAL {

template<class QQ_, class IT_>
class Ellipse_traits {
public:
    typedef QQ_ QQ;
    typedef typename Fraction_traits<QQ>::Numerator_type ZZ;
    typedef Algebraic_kernel_d_1<ZZ> AK;
    typedef typename AK::Algebraic_real_1 Root;
    typedef typename AK::Polynomial_1 upolz_t;
    typedef Polynomial<upolz_t> bpolz_t;
    typedef Polynomial<bpolz_t> tpolz_t;
    typedef Polynomial<tpolz_t> qpolz_t;

    typedef Polynomial_traits_d<upolz_t> PTZ;
    typedef typename PTZ::template Rebind<QQ, 1>::Other PTQ;
    typedef typename PTQ::Polynomial_d upolq_t;

    typedef Polynomial_traits_d<bpolz_t> BPTZ;
    typedef typename BPTZ::template Rebind<upolq_t, 1>::Other BPTQ;
    typedef typename BPTQ::Polynomial_d bpolq_t;
//    typedef Fraction_traits<upolq_t> FRT;
    
    typedef Polynomial_traits_d<tpolz_t> TPTZ;
    typedef typename TPTZ::template Rebind<bpolq_t, 1>::Other TPTQ;
    typedef typename TPTQ::Polynomial_d tpolq_t;

    typedef Polynomial_traits_d<qpolz_t> QPTZ;
    typedef typename QPTZ::template Rebind<tpolq_t, 1>::Other QPTQ;
    typedef typename QPTQ::Polynomial_d qpolq_t;
    
    typedef IT_ IT;
    typedef typename Interval_traits<IT_>::Bound BT;
    typedef typename CGAL::Polynomial<IT_> upoli_t;
    typedef Polynomial<upoli_t> bpoli_t;
    typedef Polynomial<bpoli_t> tpoli_t;
    typedef Polynomial_traits_d<upoli_t> PTI;
    typedef Polynomial_traits_d<bpoli_t> BPTI;
    typedef Polynomial_traits_d<tpoli_t> TPTI;
    
//    static IT make_interval(BT lower, BT upper) {
//        return IT(std::make_pair(lower,upper));
//    }

//    static IT make_interval(BT lower, BT upper) {
//        return IT(lower, upper);
//    }


    static IT to_interval(QQ x) {
        return typename CGAL::Coercion_traits<IT, QQ>::Cast()(x);
//        return IT(x);
    }

    // ugly hack -- due to the fact that get_arithmetic_kernel thinks that only Gmpfi exists :)
    // also: precision control
//#ifdef CGAL_BOOST_INTERVAL_GMPFR_H
    //adapted from Algebraic_readl_d_1
    static IT to_interval(Root x, long prec = -1) {
        if (x.is_rational()) return to_interval(x.rational());
        if (CGAL::sign(x) == CGAL::ZERO) return IT(0);
        CGAL_postcondition(CGAL::sign(x.low()) == CGAL::sign(x.high()));

        IT bfi = CGAL::hull(to_interval(x.low()), to_interval(x.high()));
        long last_prec = CGAL::get_significant_bits(bfi);
        while( !CGAL::singleton(bfi) ) {
            if (prec > 0 && last_prec >= prec) break;
            x.refine();
            bfi = CGAL::hull(to_interval(x.low()), to_interval(x.high()));
            int new_prec = CGAL::get_significant_bits(bfi);
            if (new_prec <= last_prec && new_prec > 1) break; else last_prec = new_prec;
        }
        return bfi;
    }
//#else
//    static IT to_interval(Root x) {
//        return CGAL::convert_to_bfi(x);
//    }
//#endif

//    WARNING: this is SLOOW!
//    static IT to_interval(Root x) {
//        Root y = x;
//        typename AK::Algebraic_real_traits::Refine()(y, mpfr_get_default_prec());
//        QQ l = typename AK::Algebraic_real_traits::Lower_bound()(y);
//        QQ r = typename AK::Algebraic_real_traits::Upper_bound()(y);
//        IT li = to_interval(l);
//        IT ri = to_interval(r);
//        return IT(CGAL::lower(li), CGAL::upper(ri));
//    }

    static std::vector<Root> Real_roots(const upolz_t& poly) {
    	std::vector<std::pair<Root, int> > solm;
        std::vector<Root> sol;

        typename AK::Solve_1()(poly, std::back_inserter(solm));
        for (int i = 0; i < solm.size(); i++) sol.push_back(solm[i].first);
        return sol;
    }

    static QQ to_rational(BT x) {
        ZZ num, den;
        mp_exp_t e = mpfr_get_z_exp(num.mpz(), x.fr());
        if (e == 0) return QQ(num);
        if (e > 0) {
            mpz_mul_2exp(den.mpz(), num.mpz(), e);
            return QQ(den);
        } else {
            ZZ tmp(1);
            mpz_mul_2exp(den.mpz(), tmp.mpz(), -e);
            return QQ(num, den);
        }
    }

    typedef Cartesian<QQ> Kernel;
    typedef Cartesian<IT> Kernel_apx;
    typedef Cartesian<double> Kernel_scr;
    
    // approximate arc
    typedef VORELL::Range<BT> ARX;
    
    static ARX to_arx(VORELL::Range<IT> x, bool inwards = false) {
        IT li = x.left();
        IT ri = x.right();
        if (CGAL::overlap(li, ri)) {
            std::cerr << li << "!!" << ri << std::endl;
        }
        CGAL_assertion( ! CGAL::overlap(li,ri) );
        if (!inwards) return ARX(lower(li), upper(ri));
        else return ARX(upper(li), lower(ri));
    }

    static ARX to_arx(VORELL::Range<Root> x, bool inwards = false, long prec = -1) {
        IT li = to_interval(x.left(), prec);
        IT ri = to_interval(x.right(), prec);
        if (CGAL::overlap(li, ri)) {
            std::cerr << li << "!!" << ri << std::endl;
        }
        CGAL_assertion( ! CGAL::overlap(li,ri) );
        if (!inwards) return ARX(lower(li), upper(ri));
        else return ARX(upper(li), lower(ri));
    }
    
#if 0
    static BT midpoint(const ARX& r) {
        BT a(r.left()), b(r.right()), m;
        if (r.is_finite()) {
            BT m = a + (b - a) / BT(2);
            CGAL_assertion (m >= a && m <= b);
            return m;
        } else { // via symmetric ;-) !
            BT pa, pb, pm;
            if ((a > 0 && b > 0) || (a < 0 && b < 0)) {
                pa = b; pb = a;
            } else {
                pa = BT(-1)/a;
                pb = BT(-1)/b;
            }
            pm = pa + (pb - pa) / BT(2);
            CGAL_assertion (pm > 0 || pm < 0);
//            std::cerr << "projective : " << pa << ' ' << pm << ' ' << pb << std::endl;
            BT m = BT(-1)/pm;
            CGAL_assertion (m >= a || m <= b);
            return m;
        }
//        return (m > b)? b: m;
    }
#endif
private:    
    // map (-Inf,+Inf) to (-2,2)
    static BT arc_like(const BT t) {
        if (t <= BT(1) && t >= BT(-1)) return t;
        if (t > BT(1)) return BT(2) - BT(1) / t;
        return BT(-2) - BT(1) / t;
    }

    
    // map (-2,2) to (-Inf,+Inf)
    static BT tan_like(const BT a) {
        if (a <= BT(1) && a >= BT(-1)) return a;
        if (a > BT(1)) return BT(1) / (BT(2)-a);
        return BT(1) / (BT(-2)-a);
    }

public:
    static BT midpoint(const ARX& r) {
        BT a(r.left()), b(r.right()), m;
        if (r.is_finite() && (b-a) < BT(2)) {
            if ((a >= 0 && b >= 0) || (a <= 0 && b <= b)) 
                m = a + (b-a)/BT(2);
            else m = (a+b)/BT(2);
            CGAL_assertion (m >= a && m <= b);
            return m;
        }
        a = arc_like(a);
        b = arc_like(b);
        if ((a >= 0 && b >= 0) || (a <= 0 && b <= b)) m = a + (b-a)/BT(2);
        else m = (a+b)/BT(2);
        m = tan_like(m);
        if (r.is_finite()) {
            CGAL_assertion (m >= r.left() && m <= r.right());
            return m;
        } else {
            CGAL_assertion( m > 0 || m < 0);
            m = BT(-1)/m; 
            CGAL_assertion (m >= r.left() || m <= r.right());
            return m;
        }
    }
};

typedef CGAL::Algebraic_kernel_d_1<CGAL::Gmpz> AKz;


//template <typename T>
//struct Coercion_traits<CGAL::VORELL::Range<CGAL::Gmpfr>,
//                        CGAL::VORELL::Range<T> > {
//    typedef Tag_true  Are_explicit_interoperable;
//    typedef Tag_false Are_implicit_interoperable;

//private:

//    typedef CGAL::Gmpfi INT;
//    typedef CGAL::VORELL::Range<CGAL::Gmpfr> ARX;
//    typedef CGAL::VORELL::Range<T> Range;
//public:

//    typedef ARX Type;

//    struct Cast{
//        typedef ARX result_type;
//        ARX operator()(const ARX& x)  const { return x; }
//        ARX operator()(const Range& x, bool inwards = false) const {
//            INT li = typename CGAL::Coercion_traits<INT, T>::Cast()(x.left());
//            INT ri = typename CGAL::Coercion_traits<INT, T>::Cast()(x.right());
//            if (CGAL::overlap(li, ri)) {
//                std::cerr << li << "!!" << ri << std::endl;
//            }
//            CGAL_assertion( ! CGAL::overlap(li,ri) );
//            if (!inwards) return ARX(lower(li), upper(ri));
//            else return ARX(upper(li), lower(ri));
//        }
//    };
//};



//template <>
//struct Coercion_traits< mpfr_interval,
//                        AKz::Algebraic_real_1> {
//    typedef Tag_true  Are_explicit_interoperable;
//    typedef Tag_false Are_implicit_interoperable;
    
//private:
//    typedef CGAL::mpfr_interval INT;
//    typedef CGAL::Gmpz ZZ;
//    typedef CGAL::Gmpq QQ;
//    typedef AKz AK;
//    typedef AK::Algebraic_real_1 Root;

//public:
//    typedef INT Type;
    
//    struct Cast{
//        typedef INT result_type;
//        INT operator()(const INT& x)  const { return x; }
//        INT operator()(const Root& x) const {
//            Root y = x;
//            AK::Algebraic_real_traits::Refine()(y, mpfr_get_default_prec());
//            QQ l = AK::Algebraic_real_traits::Lower_bound()(y);
//            QQ r = AK::Algebraic_real_traits::Upper_bound()(y);
//            INT li = CGAL::Coercion_traits<INT, QQ>::Cast()(l);
//            INT ri = CGAL::Coercion_traits<INT, QQ>::Cast()(r);
//            return INT(lower(li), upper(ri));
//        }
//    };
//};


//// how can I do it with templates? (ZZ,QQ)
//template <>
//struct Coercion_traits<CGAL::Polynomial<CGAL::Gmpz>, CGAL::Polynomial<CGAL::Gmpq> >{
//    typedef Tag_true  Are_explicit_interoperable;
//    typedef Tag_false Are_implicit_interoperable;
    
//private:
//    typedef CGAL::Polynomial<CGAL::Gmpz> POLYZ;
//    typedef CGAL::Polynomial<CGAL::Gmpq> POLYQ;
////    typedef CGAL::Fraction_traits<POLYQ>::Numerator_type POLYZ;
//    typedef CGAL::Polynomial_traits_d<POLYZ> PT;

//public:
//    typedef POLYZ Type;

//    struct Cast{
//        typedef POLYZ result_type;
//        POLYZ operator()(const POLYZ& x)  const { return x;}
//        POLYZ operator()(const POLYQ& x) const {
//            POLYZ polyz;
//            Fraction_traits<POLYQ>::Denominator_type denom;
//            Fraction_traits<POLYQ>::Decompose()(x, polyz, denom);
            
//            Sign s1 = CGAL::sign(PT::Innermost_leading_coefficient()(polyz));
//            polyz = PT::Canonicalize()(polyz);
//            Sign s2 = CGAL::sign(PT::Innermost_leading_coefficient()(polyz));
//            // keep sign!
//            if (s1 != s2) return -polyz;
//            return polyz;
//        }
//    };
//};



} //namespace CGAL
#endif
