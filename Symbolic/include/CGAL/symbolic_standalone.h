// Copyright (c) 2010, 2011 Max-Planck-Institut fuer Informatik (Germany).
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
// $URL: svn+ssh://eric@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Algebraic_kernel_2/symbolic_speedups.cpp $
// $Id: symbolic_speedups.cpp 61278 2011-02-17 10:15:03Z eric $
//
//
// Author(s): Eric Berberich    <eric@mpi-inf.mpg.de>
//            Pavel Emeliyanenko <asm@mpi-inf.mpg.de>
// ============================================================================

/*! \file symbolic_standalone.h
 *  
 *  - enable/disable GPU-computations; standalone version
 * with multithreading support
 *  - conversion of polynomials to NTL format & back
 */

#include <CGAL/config.h>

#ifndef CGAL_BISOLVE_USE_BIGCD 
#define CGAL_BISOLVE_USE_BIGCD 0

#define CGAL_BIGCD_USE_SHIFT 0
#define CGAL_BIGCD_CHECK_SANITY 1

#endif

#if CGAL_BISOLVE_USE_BIGCD
#warning Using bi-gcd algorithm
#endif

#ifndef CGAL_BISOLVE_USE_GPU_RESULTANTS
#define CGAL_BISOLVE_USE_GPU_RESULTANTS 1 // default?
#define CGAL_BISOLVE_CHECK_GPU_RESULTANTS_SANITY 0 // default 0
#endif

#if CGAL_BISOLVE_USE_GPU_RESULTANTS
#warning Using gpu-resultants
#endif

#ifndef CGAL_BISOLVE_USE_GPU_GCDS
#define CGAL_BISOLVE_USE_GPU_GCDS 1  // default?
#define CGAL_BISOLVE_CHECK_GPU_GCDS_SANITY 0 // default 1
#endif

#if CGAL_BISOLVE_USE_GPU_GCDS
#warning Using gpu-gcds
#endif

#ifndef CGAL_BISOLVE_USE_NTL
#define CGAL_BISOLVE_USE_NTL  1 // default 1 ??
#endif

#define CGAL_MODULAR_FILTER_OFF

#include <CGAL/GMP_arithmetic_kernel.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>

#if (CGAL_BISOLVE_USE_GPU_RESULTANTS || CGAL_BISOLVE_USE_GPU_GCDS)
#include <CGAL/GPU_algorithms/GPU_algorithm_facade.h>
#endif

#if CGAL_BISOLVE_USE_NTL
#warning Bisolve using NTL
#include <NTL/ZZX.h>

#ifdef CGAL_WEBXACLI_MULTITHREADED
CGAL_WEBXALCI_DEFMUTEX(NTL_arithm_mtx)
#endif

#endif 

#if (CGAL_BISOLVE_USE_GPU_RESULTANTS || CGAL_BISOLVE_USE_GPU_GCDS)
#ifdef CGAL_WEBXACLI_MULTITHREADED
#warning Compiling with mutlithreading support
CGAL_WEBXALCI_DEFMUTEX(GPU_arithm_mtx)
#endif
#endif

namespace CGAL {

namespace internal {

template < class NT >
void writeout(const CGAL::Polynomial< NT >& p1,
        const CGAL::Polynomial< NT >& p2) {

    std::stringstream ss, ss2;
    CGAL::set_pretty_mode(ss);
    CGAL::set_pretty_mode(ss2);
    ss << p1 << "\n" << p2;
    std::ofstream o("ERROR_OUT"/*, std::ios::app*/);
    o << ss.str();
    o.close();
    ss2 << "F := [" << p1 << ", " << p2 << "]:";
    std::ofstream o2("ERROR_OUT_maple");
    o2 << ss2.str();
    o2.close();
}


#if CGAL_BISOLVE_USE_GPU_RESULTANTS

// repeat original version of "resultant"
template < class Coeff > 
inline Coeff resultant_cpu(
        const CGAL::Polynomial<Coeff>& F_, 
        const CGAL::Polynomial<Coeff>& G_){
    // Is it possible not to copy and paste the code of original "resultant" here?
  // make the variable to be elimnated the innermost one.
    typedef CGAL::Polynomial_traits_d<CGAL::Polynomial<Coeff> > PT;
    CGAL::Polynomial<Coeff> F = typename PT::Move()(F_, PT::d-1, 0);
    CGAL::Polynomial<Coeff> G = typename PT::Move()(G_, PT::d-1, 0);
    return internal::resultant_(F,G);
}

template < class Coeff >
inline Coeff resultant_gpu(const CGAL::Polynomial<Coeff>& F_,
                           const CGAL::Polynomial<Coeff>& G_,
                            bool zero_when_failed = false) {

    if(F_.degree() < 5 && G_.degree() < 5)
        return resultant_cpu(F_,G_);

    std::cout << "\nusing GRES\n";
//     std::cout << F_ << "\n" << G_ << "\n\n";
    bool failed;
#ifdef CGAL_WEBXACLI_MULTITHREADED
    CGAL_WEBXALCI_MUTEXLOCK(GPU_arithm_mtx)
#endif
    Coeff res = GPU_algorithm_facade::resultant(F_,G_, &failed);
#ifdef CGAL_WEBXACLI_MULTITHREADED
    CGAL_WEBXALCI_MUTEXUNLOCK(GPU_arithm_mtx)
#endif

    if(failed) {
        if(zero_when_failed)
            return Coeff(0);
        
//         writeout(F_,G_);
        std::cout << "\nGRES failed\n";
        return resultant_cpu(F_,G_);
    }
#if CGAL_BISOLVE_CHECK_GPU_RESULTANTS_SANITY
    std::cout << "\nGRES sanity check!!\n";
    Coeff check = CGAL::canonicalize(res),
          truth = CGAL::canonicalize(resultant_cpu(F_,G_));

    if(check != truth) {
        std::cout << "\nWrong GRES:\n" << check << "\n" <<
            truth << "\n";
        writeout(F_, G_);
        throw "WTF!?";
    }
#endif
    std::cout << "\nGRES done\n";
    return res;
}

// partial specialization to enable gpu'ed resultants
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL

template <> inline
CGAL::Polynomial< CGAL::Gmpz > resultant(
        const CGAL::Polynomial< CGAL::Polynomial< CGAL::Gmpz > >& F_,
        const CGAL::Polynomial< CGAL::Polynomial< CGAL::Gmpz > >& G_) {

    return resultant_gpu(F_, G_);
}
#endif

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

template <> inline
CGAL::Polynomial< CORE::BigInt > resultant(
        const CGAL::Polynomial< CGAL::Polynomial< CORE::BigInt > >& F_,
        const CGAL::Polynomial< CGAL::Polynomial< CORE::BigInt > >& G_) {

    return resultant_gpu(F_, G_);
}
#endif

#endif // CGAL_BISOLVE_USE_GPU_RESULTANTS

#if CGAL_BISOLVE_USE_GPU_GCDS

template < class NT >
Polynomial< NT > gcd_NTL(const Polynomial< NT >& p1,
        const Polynomial< NT >& p2);

template < class Poly > inline
Poly gcd_gpu(const Poly& F_, const Poly& G_) {

//     printf("\nstart GGCD ..\n");

    bool failed;
#ifdef CGAL_WEBXACLI_MULTITHREADED
    CGAL_WEBXALCI_MUTEXLOCK(GPU_arithm_mtx)
#endif
    Poly ggcd = GPU_algorithm_facade::gcd(F_, G_, &failed);
    // TODO: need to use failed flag in GPU symbolic as well !!!
#ifdef CGAL_WEBXACLI_MULTITHREADED
    CGAL_WEBXALCI_MUTEXUNLOCK(GPU_arithm_mtx)
#endif

    if(failed) {
        return gcd_NTL(F_, G_);
    }

#if CGAL_BISOLVE_CHECK_GPU_GCDS_SANITY
    //printf("\nGGCD sanity check..\n");
    Poly truth = gcd_NTL(F_, G_);

    if(truth != ggcd) {
//         Poly diff = truth + ggcd;
       writeout(F_, G_);
//         writeout(truth, ggcd);
        std::cerr << "Wrong gcd!!\n";
        throw "WTF?";
    }
//     printf("GGCD ok\n");
#endif
    return ggcd;
}

// partial specialization to enable gpu'ed gcds
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
template <> inline
CGAL::Polynomial< CGAL::Gmpz > gcd_(
        const CGAL::Polynomial< CGAL::Gmpz >& F_,
        const CGAL::Polynomial< CGAL::Gmpz >& G_) {

    return gcd_gpu(F_, G_);
}

template <> inline
CGAL::Polynomial< CGAL::Gmpz > gcd_utcf_(
        const CGAL::Polynomial< CGAL::Gmpz >& F_,
        const CGAL::Polynomial< CGAL::Gmpz >& G_) {

    return gcd_gpu(F_, G_);
}
#endif

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
template <> inline
CGAL::Polynomial< CORE::BigInt > gcd_(
        const CGAL::Polynomial< CORE::BigInt >& F_,
        const CGAL::Polynomial< CORE::BigInt >& G_) {

    return gcd_gpu(F_, G_);
}

template <> inline
CGAL::Polynomial< CORE::BigInt > gcd_utcf_(
        const CGAL::Polynomial< CORE::BigInt >& F_,
        const CGAL::Polynomial< CORE::BigInt >& G_) {

    return gcd_gpu(F_, G_);
}
#endif

#endif // CGAL_BISOLVE_USE_GPU_GCDS

#if CGAL_BISOLVE_USE_BIGCD

//! maps y to (2*x)^yexp in a bivariate polynomial \c f
//! \c degf: estimated degree of the resulting univariate polynomial
template < class Poly_1 >
Poly_1 _bi_map_degree(const CGAL::Polynomial< Poly_1 >& f, int yexp,
                        int degf) {

    typedef typename Poly_1::NT Integer;
    typename Poly_1::Vector vf(degf + 1, Integer(0));
    for(int i = 0; i <= f.degree(); i++) {
        const Poly_1& p = f[i];
        int deg = p.degree(), shift = i * yexp; // shift the polynomial

        if(deg + shift >= vf.size()) {
            vf.resize(deg + shift + 1, Integer(0));
//             printf("f resizing poly to: %d\n", vf.size());
        }
        int sh = shift * CGAL_BIGCD_USE_SHIFT;
        for(int j = 0; j <= deg; j++) {
            vf[j + shift] += (p[j] << sh);
        }
    }
    return Poly_1(vf.begin(), vf.end());
}

//! maps every non-zero monomial \c c*x^a of \c f to \c c/2^s*y^s*x^t
//! \c s = (a div yexp) and \c t = (a mod yexp)
//! finally the resulting bivariate polynomial is multiplied by \c y^trail_cf
template < class Poly_1 >
CGAL::Polynomial< Poly_1 > _bi_unmap_degree(const Poly_1& f, int yexp,
                                                int trail_cf) {

    typedef CGAL::Polynomial< Poly_1 > Poly_2;
    typedef typename Poly_1::Vector Vector;
    typedef typename Poly_1::NT Integer;

    int f_deg = f.degree(), f_ydeg = f_deg / yexp;
    //! need to multiply the result by \c gcd_trail
    typename Poly_2::Vector v(f_ydeg + 1 + trail_cf, Poly_1(Integer(0)));

    //! 1. take all coeffs with degrees smaller than yexp -> y^0
    //! 2. take all coeffs with degrees btw yexp and 2*yexp -> y^1
    //! etc..
    Vector vtmp(yexp);
    for(int i = 0; i <= f_deg; i += yexp) {
        int ydeg = (i/yexp) + trail_cf;
        int xend = std::min(i + yexp, f_deg + 1);
//          std::cout << i << " range: [" << i << "; " <<
//              xend << ") shift: " << (1 << i) << "\n";
        int sh = i * CGAL_BIGCD_USE_SHIFT;
        for(int j = i; j < xend; j++) {
            vtmp[j - i] = f[j] >> sh;
        }
        v[ydeg] = Poly_1(vtmp.begin(), vtmp.begin() + xend - i);
//         vgcd[ydeg] = Poly_1(ci, ci + xend - i);
//         ci += xend - i;
    }
    return Poly_2(v.begin(), v.end());
}

template < class Poly_2 >
Poly_2 _bi_remove_trail_cf(const Poly_2& f, int& trail_cf) {

    trail_cf = 0;
    for(; CGAL::is_zero(f[trail_cf]); trail_cf++);
    if(trail_cf == 0)
        return f;
    return Poly_2(f.begin() + trail_cf, f.end());
}

template < class Poly_2 >
Poly_2 _bi_gcd(const Poly_2& F_, const Poly_2& G_) {

    typedef typename Poly_2::NT Poly_1;
    typedef typename Poly_2::Vector Vector_2;
    typedef typename Poly_1::NT Integer;

    if(CGAL::is_zero(F_)) {
        if(CGAL::is_zero(G_))
            return Poly_2(Poly_1(1));
        else
            return G_;
    }
    if(CGAL::is_zero(G_) || F_ == G_)
        return F_;

    typedef typename CGAL::Polynomial_traits_d< Poly_2 > PT_d;
    typename PT_d::Swap swap;
    typename CGAL::internal::Degree< Poly_2 > degree;

    // 1. remove y exponent in both F and G
    // 2. compute gcd of new polynomials
    // 3. put the content back

    Poly_1 fc = F_.content(), gc = G_.content();
    Poly_2 f = F_ / fc, g = G_ / gc;
    Poly_1 gcdc = CGAL::gcd(fc, gc);

    int degfx = degree(f, 0), deggx = degree(g, 0);
    int degfy = degree(f, 1), deggy = degree(g, 1);

    printf("deg f(x/y): %d/%d; deg g(x/y): %d/%d\n",
        degfx, degfy, deggx, deggy);

    if(degfy == 0 || deggy == 0) {
//         if(gcd_trail == 0)
        return Poly_2(gcdc);
/*        Vector_2 vv(gcd_trail + 1, Poly_1(Integer(0)));
        vv[gcd_trail] = gcdc;
        return Poly_2(vv.begin(), vv.end());*/
    }

    if(degfx == 0 || deggx == 0) {
        Poly_2 res = internal::_bi_gcd(swap(f, 0, 1), swap(g, 0, 1));
        res = swap(res, 0, 1) * gcdc;
        return res;
    }

    int ftrail, gtrail, gcd_trail;
    f = _bi_remove_trail_cf(f, ftrail);
    g = _bi_remove_trail_cf(g, gtrail);
    gcd_trail = std::min(ftrail, gtrail);
    degfy -= ftrail, deggy -= gtrail;

    int minx = std::min(degfx, deggx);

    // replacing y to x^(minx + 1)
    int degf = std::max(degfx, degfy * (minx + 1));
    int degg = std::max(deggx, deggy * (minx + 1));
    int yexp = std::max(minx + 1, 2);

    printf("replacing y -> x^(minx + 1): %d %d\n", degf, degg);

//     std::cout << "f/fc: " << f << "\ng/gc: " << g << "\n";
//     std::cout << "content1: " << fc << "; content2: " <<
//             gc << "\n";

//     printf("yexp: %d\n", yexp);

    Poly_1 fx = _bi_map_degree(f, yexp, degf),
           gx = _bi_map_degree(g, yexp, degg);
    Poly_1 ggcd = CGAL::gcd(fx, gx);

//     std::cout << "fx: " << fx << "\ngx: " << gx << "\nggcdx: " <<
//          ggcd << "\n";

    Poly_2 res = _bi_unmap_degree(ggcd, yexp, gcd_trail);
    Poly_1 cnt = res.content();
    res /= cnt;
    res *= gcdc;

    return res;
}

template < class Poly_2 > inline
Poly_2 bi_gcd(const Poly_2& F_, const Poly_2& G_) {

    Poly_2 res = internal::_bi_gcd(F_, G_);
#if CGAL_BIGCD_CHECK_SANITY
    std::cout << "bigcd: computing reference solution..\n";
    Poly_2 truth = CGAL::internal::gcd_utcf_UFD(F_, G_);
    std::cout << "bigcd: done..\n";

    if(res != truth) {
        std::cout << "bigcd failed!\n";
        writeout(F_, G_);

        std::cout << "res: " << res << "\ntruth: " << truth << "\n";
        throw -1;
    }
    std::cout << "$$$$$$$$$$$$$$$$$$ bigcd done..\n";
#endif
    return res;
}

template < class Poly_2 >
Poly_2 bi_exact_div(const Poly_2& F_, const Poly_2& G_) {

    typedef typename Poly_2::NT Poly_1;
    typedef typename Poly_1::NT Integer;
    typedef typename Poly_1::Vector Vector;

    if(CGAL::is_zero(F_)) {
        return Poly_2(Poly_1(0));
    }
    if(CGAL::is_zero(G_)) {
        printf("FATAL: bi_exact_div zero divisor!\n");
        throw -1;
    }

    if(G_.degree() == 0) {
        return F_ / G_[0];
    }

    typedef typename CGAL::Polynomial_traits_d< Poly_2 > PT_d;
    typename PT_d::Swap swap;
//     typedef typename PT_d::Degree degree;
    typename CGAL::internal::Degree< Poly_2 > degree;

    Poly_2 f, g;
    int ftrail, gtrail, div_trail;
    f = _bi_remove_trail_cf(F_, ftrail);
    g = _bi_remove_trail_cf(G_, gtrail);
    div_trail = ftrail - gtrail;

    Poly_1 fc = f.content(), gc = g.content();
    f = f / fc, g = g / gc;
    Poly_1 divc = fc / gc;

//     std::cout << "F: " << F_ << "\nG: " << G_ << "\n";

    int degfx = degree(f, 0), degfy = degree(f, 1),
        deggx = degree(g, 0), deggy = degree(g, 1);

    int maxx = degfx - deggx;
        //std::max(deggx, degfx - deggx);

    if(div_trail < 0 || maxx < 0) {
        printf("FATAL: bi_exact_div polynomials do not divide!\n");
        throw -1;
    }

//     printf("bi_exact_div: deg f(x/y): %d/%d; deg g(x/y): %d/%d\n",
//         degfx, degfy, deggx, deggy);

    int degf = std::max(degfx, degfy * (maxx + 1));
    int degg = std::max(deggx, deggy * (maxx + 1));

    int yexp = std::max(maxx + 1, 2);

     Poly_1 fx = _bi_map_degree(f, yexp, degf),
           gx = _bi_map_degree(g, yexp, degg);

//     std::cout << "fx: " << fx << "\ngx: " << gx << "\n";
    Poly_1 divx = fx / gx;

    Poly_2 res = _bi_unmap_degree(divx, yexp, div_trail);
    Poly_1 cnt = res.content();
    res /= cnt;
    res *= divc;

#if CGAL_BIGCD_CHECK_SANITY
    std::cout << "bi_exact_div: checking solution..\n";

    if(F_ != res * G_) {
        std::cout << "bi_exact_div failed!\n";
        writeout(F_, G_);
        throw -1;
    }
    std::cout << "$$$$$$$$$$$$$$$$$$ bi_exact_div done..\n";
#endif

//     std::cout << "res: " << res << "\n";
//     std::cout << "div cnt: " << cnt << "\n";

    return res;
}

#endif // CGAL_BISOLVE_USE_BIGCD

#if CGAL_BISOLVE_USE_NTL

struct NTL_bigint_rep {
    long alloc;
    long size;
    mp_limb_t data;
};

template < class _ >
struct ___ {
// dummy
};

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

template < >
struct ___< CORE::BigInt > {

    typedef CORE::BigInt Integer;

    inline mpz_srcptr operator()(const Integer& x) const {
        return x.get_mp();
    }

    inline Integer operator()(mpz_srcptr z) const {
        return Integer(z);
    }
};

#endif // CGAL_HAS_CORE_ARITHMETIC_KERNEL

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL

template < >
struct ___< CGAL::Gmpz > {

    typedef CGAL::Gmpz Integer;

    inline mpz_srcptr operator()(const Integer& x) const {
        return x.mpz();
    }

    inline Integer operator()(mpz_srcptr z) const {
        return Integer(z);
    }
};

#endif // CGAL_HAS_GMP_ARITHMETIC_KERNEL

template < class NT >
void poly2ntl(const Polynomial< NT >& p, NTL::ZZX& q) {

    q.rep.SetLength(p.degree() + 1);
    ___< NT > cvt; // instantiate secret object
    
    int i;
    typename Polynomial< NT >::const_iterator pit;
    for(i = 0, pit = p.begin(); pit != p.end(); pit++, i++) {

        NTL::ZZ& zz = q.rep[i];
        mpz_srcptr tmp = cvt(*pit);
        int sz = tmp->_mp_size;
        if(sz == 0)
            continue;
        if(sz < 0)
            sz = -sz;
        
        zz.SetSize(sz);
        NTL_bigint_rep *rep = (NTL_bigint_rep *)zz.rep;
        rep->size = tmp->_mp_size;
        // copy limbs directly
        memcpy(&rep->data, tmp->_mp_d, sz*sizeof(mp_limb_t));
    }
}

template < class NT >
Polynomial< NT > ntl2poly(const NTL::ZZX& r) {

    typedef Polynomial< NT > Poly_1;

    int d = NTL::deg(r);
    typename Poly_1::Vector v(d + 1);
    
    mpz_t tmp;
    mpz_init(tmp);
    for(int i = 0; i <= d; i++) {
        
        const NTL::ZZ& zz = r.rep[i];
        if(NTL::IsZero(zz)) {
            v[i] = NT(0);
            continue;
        } 

        NTL_bigint_rep *rep = (NTL_bigint_rep *)zz.rep;
        int sz = rep->size;
        if(sz < 0)
            sz = -sz;
         
        mpz_realloc2(tmp, sz * GMP_NUMB_BITS);
        tmp->_mp_size = rep->size;
        memcpy(tmp->_mp_d, &rep->data, sz*sizeof(mp_limb_t));
         
        v[i] = NT(tmp);
    }
    mpz_clear(tmp);
    return (Poly_1(v.begin(), v.end()));
}

template < class NT >
Polynomial< NT > gcd_NTL(const Polynomial< NT >& p1,
        const Polynomial< NT >& p2) {

    typedef Polynomial< NT > Poly_1;

    if(p1 == Poly_1(0))
        return p2;
    if(p2 == Poly_1(0))
        return p1;

#ifdef CGAL_WEBXACLI_MULTITHREADED
    CGAL_WEBXALCI_MUTEXLOCK(NTL_arithm_mtx)
#endif
    NTL::ZZX r, q1, q2;
    poly2ntl(p1, q1);
    poly2ntl(p2, q2);
    
    NTL::GCD(r, q1, q2);
    
    Polynomial< NT > p = ntl2poly< NT >(r);
#ifdef CGAL_WEBXACLI_MULTITHREADED
    CGAL_WEBXALCI_MUTEXUNLOCK(NTL_arithm_mtx)
#endif

    return p;
}

#if !CGAL_BISOLVE_USE_GPU_GCDS

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

template <> inline
Polynomial< CORE::BigInt > gcd_(const Polynomial< CORE::BigInt >& p1,
            const Polynomial< CORE::BigInt >& p2) {

/*    n_ggcd_calls++;
    tm_ggcd.start();*/
    Polynomial< CORE::BigInt > ggcd = gcd_NTL(p1, p2);
//     tm_ggcd.stop();
    return ggcd;
}

template <> inline
Polynomial< CORE::BigInt > gcd_utcf_(const Polynomial< CORE::BigInt >& p1,
            const Polynomial< CORE::BigInt >& p2) {

/*    n_ggcd_calls++;
    tm_ggcd.start();*/
    Polynomial< CORE::BigInt > ggcd = gcd_NTL(p1, p2);
//     tm_ggcd.stop();
    return ggcd;
}
#endif

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL

template <> inline
Polynomial< CGAL::Gmpz > gcd_(const Polynomial< CGAL::Gmpz >& p1,
            const Polynomial< CGAL::Gmpz >& p2) {

    Polynomial< CGAL::Gmpz > ggcd = gcd_NTL(p1, p2);
    return ggcd;
}

template <> inline
Polynomial< CGAL::Gmpz > gcd_utcf_(const Polynomial< CGAL::Gmpz >& p1,
            const Polynomial< CGAL::Gmpz >& p2) {

    Polynomial< CGAL::Gmpz > ggcd = gcd_NTL(p1, p2);
    return ggcd;
}

#endif // CGAL_HAS_GMP_ARITHMETIC_KERNEL

#endif // !CGAL_BISOLVE_USE_GPU_GCDS

#endif // CGAL_BISOLVE_USE_NTL

#if CGAL_BISOLVE_USE_BIGCD
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL

template <> inline
CGAL::Polynomial< CGAL::Polynomial< CGAL::Gmpz > > gcd_(
        const CGAL::Polynomial< CGAL::Polynomial< CGAL::Gmpz > >& F_,
        const CGAL::Polynomial< CGAL::Polynomial< CGAL::Gmpz > >& G_) {

    return bi_gcd(F_, G_);
}

template <> inline
CGAL::Polynomial< CGAL::Polynomial< CGAL::Gmpz > > gcd_utcf_(
        const CGAL::Polynomial< CGAL::Polynomial< CGAL::Gmpz > >& F_,
        const CGAL::Polynomial< CGAL::Polynomial< CGAL::Gmpz > >& G_) {

    return bi_gcd(F_, G_);
}
#endif

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
template <> inline
CGAL::Polynomial< CGAL::Polynomial< CORE::BigInt > > gcd_(
        const CGAL::Polynomial< CGAL::Polynomial< CORE::BigInt > >& F_,
        const CGAL::Polynomial< CGAL::Polynomial< CORE::BigInt > >& G_) {

    return bi_gcd(F_, G_);
}

template <> inline
CGAL::Polynomial< CGAL::Polynomial< CORE::BigInt > > gcd_utcf_(
        const CGAL::Polynomial< CGAL::Polynomial< CORE::BigInt > >& F_,
        const CGAL::Polynomial< CGAL::Polynomial< CORE::BigInt > >& G_) {

    return bi_gcd(F_, G_);
}
#endif

#endif // CGAL_BISOLVE_USE_BIGCD

} // namespace internal

#if (defined CGAL_USE_GPU || defined CGAL_USE_NTL)

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
template <> inline
CGAL::Polynomial< CGAL::Gmpz > gcd(
        const CGAL::Polynomial< CGAL::Gmpz >& F_,
        const CGAL::Polynomial< CGAL::Gmpz >& G_) {
    return internal::gcd_(F_,G_);
}
#endif

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

template <> inline
CGAL::Polynomial< CORE::BigInt > gcd(
        const CGAL::Polynomial< CORE::BigInt >& F_,
        const CGAL::Polynomial< CORE::BigInt >& G_) {
    return internal::gcd_(F_,G_);
}

#endif

#endif  // (defined CGAL_USE_GPU || defined CGAL_USE_NTL)

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL

#if CGAL_BISOLVE_USE_NTL
template < > inline
void Polynomial< CGAL::Gmpz >::euclidean_division(
    const Polynomial< CGAL::Gmpz >& f, const Polynomial< CGAL::Gmpz >& g,
    Polynomial< CGAL::Gmpz >& q, Polynomial< CGAL::Gmpz >& r) {

#ifdef CGAL_WEBXACLI_MULTITHREADED
    CGAL_WEBXALCI_MUTEXLOCK(NTL_arithm_mtx)
#endif

    typedef CGAL::Gmpz NT;
    NTL::ZZX q1, r1, f1, g1;
    internal::poly2ntl(f, f1);
    internal::poly2ntl(g, g1);

    NTL::DivRem(q1, r1, f1, g1);
    
    q = internal::ntl2poly< NT >(q1);
    r = internal::ntl2poly< NT >(r1);
#ifdef CGAL_WEBXACLI_MULTITHREADED
    CGAL_WEBXALCI_MUTEXUNLOCK(NTL_arithm_mtx)
#endif
}
#endif

#if CGAL_BISOLVE_USE_BIGCD
template < > inline
void Polynomial< Polynomial< CGAL::Gmpz > >::euclidean_division(
    const Polynomial< Polynomial< CGAL::Gmpz > >& f,
    const Polynomial< Polynomial< CGAL::Gmpz > >& g,
    Polynomial< Polynomial< CGAL::Gmpz > >& q,
    Polynomial< Polynomial< CGAL::Gmpz > >& r) {

    q = internal::bi_exact_div(f, g);
    r = Polynomial< Polynomial< CGAL::Gmpz > >
        (Polynomial< CGAL::Gmpz >(CGAL::Gmpz(0)));
}
#endif

#endif // CGAL_HAS_GMP_ARITHMETIC_KERNEL

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

#if CGAL_BISOLVE_USE_NTL

template < > inline
void Polynomial< CORE::BigInt >::euclidean_division(
    const Polynomial< CORE::BigInt >& f, const Polynomial< CORE::BigInt >& g,
    Polynomial< CORE::BigInt >& q, Polynomial< CORE::BigInt >& r) {

#ifdef CGAL_WEBXACLI_MULTITHREADED
    CGAL_WEBXALCI_MUTEXLOCK(NTL_arithm_mtx)
#endif
    typedef CORE::BigInt NT;
    NTL::ZZX q1, r1, f1, g1;
    internal::poly2ntl(f, f1);
    internal::poly2ntl(g, g1);

    NTL::DivRem(q1, r1, f1, g1);
    
    q = internal::ntl2poly< NT >(q1);
    r = internal::ntl2poly< NT >(r1);
#ifdef CGAL_WEBXACLI_MULTITHREADED
    CGAL_WEBXALCI_MUTEXUNLOCK(NTL_arithm_mtx)
#endif
}
#endif // CGAL_BISOLVE_USE_NTL

#if CGAL_BISOLVE_USE_BIGCD
template < > inline
void Polynomial< Polynomial< CORE::BigInt > >::euclidean_division(
    const Polynomial< Polynomial< CORE::BigInt > >& f,
    const Polynomial< Polynomial< CORE::BigInt > >& g,
    Polynomial< Polynomial< CORE::BigInt > >& q,
    Polynomial< Polynomial< CORE::BigInt > >& r) {

    q = internal::bi_exact_div(f, g);
    r = Polynomial< Polynomial< CORE::BigInt > >
        (Polynomial< CORE::BigInt >(CORE::BigInt(0)));
}
#endif

#endif // CGAL_HAS_CORE_ARITHMETIC_KERNEL

} // namespace CGAL

