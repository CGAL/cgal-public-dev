// ============================================================================
//
// Copyright (c) 2001-2008 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// this file is not part of any library ;-)
//
// ----------------------------------------------------------------------------
//
// Library       : CUDA MP
//
// File          : 
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef GPU_ALGORITHM_TEMPLATES_H
#define GPU_ALGORITHM_TEMPLATES_H

#include <set>
#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>
#include <CGAL/convert_to_bfi.h>

#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>
/*
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d_1_generator.h>*/

// #define STILL_ALIVE std::cout << __LINE__ << "\n";

namespace CGAL {

namespace internal {

template < class _ >
struct MPZ_traits {
// dummy
};

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

template < >
struct MPZ_traits< CORE::BigInt > {

    typedef CORE::BigInt Integer;

    inline mpz_srcptr int2mpz(const Integer& x) const {
        return x.get_mp();
    }

    inline Integer mpz2int(mpz_srcptr z) const {
        return Integer(z);
    }
};

#endif // CGAL_HAS_CORE_ARITHMETIC_KERNEL

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL

template < >
struct MPZ_traits< CGAL::Gmpz > {

    typedef CGAL::Gmpz Integer;

    inline mpz_srcptr int2mpz(const Integer& x) const {
        return x.mpz();
    }

    inline Integer mpz2int(mpz_srcptr z) const {
        return Integer(z);
    }
};

#endif // CGAL_HAS_GMP_ARITHMETIC_KERNEL

//! converts a univariate polynomial to \c MPZ_vector_1
template < class NT >
void construct_mpz_vector_1(const CGAL::Polynomial< NT >& f,
         MPZ_vector_1& zv) {

    int deg = f.degree();
    zv.resize(deg + 1);
    MPZ_traits< NT > traits;

    for(unsigned i = 0; (int)i <= deg; i++) {
        zv[i] = *traits.int2mpz(f[i]);
    }
}

//! converts a bivariate polynomial to \c MPZ_vector_2 , additionally
//! polynomial's x-degree in \c deg_x
template < class NT >
void construct_mpz_vector_2(
    const CGAL::Polynomial< CGAL::Polynomial< NT > >& f,
         MPZ_vector_2& zv, unsigned& deg_x) {

    int deg = f.degree();
    zv.resize(deg + 1);
    MPZ_traits< NT > traits;

    deg_x = 0;
    for(unsigned i = 0; (int)i <= deg; i++) {

        const CGAL::Polynomial< NT >& ci = f[i];
        MPZ_vector_1& zvi = zv[i];
        zvi.resize(ci.degree() + 1);
        // to handle -1 degree correctly
        deg_x = (unsigned)std::max((int)deg_x, ci.degree());

        for(unsigned j = 0; (int)j <= ci.degree(); j++) {
            zvi[j] = *traits.int2mpz(ci[j]);
        }
    }
}

inline void dispose_mpz_vector(MPZ_vector_1& v) {
    
    for(unsigned i = 0; i < v.size(); i++) {
        mpz_clear((mpz_ptr)&v[i]);
    }
}

template < class NT >
CGAL::Polynomial< NT > construct_polynomial_from_mpz(const MPZ_vector_1& zv) {

    unsigned sz = zv.size();
    typename CGAL::Polynomial< NT >::Vector v(sz);
    MPZ_traits< NT > _;

    for(unsigned i = 0; i < sz; i++) {
        v[i] = _.mpz2int((mpz_ptr)&zv[i]);
    }
    return CGAL::Polynomial< NT >(v.begin(), v.end());
}

#define ILOG2(x) mpz_sizeinbase(_.int2mpz(x), 2)

#define ILOG2_MPZ(x) mpz_sizeinbase(x, 2)

inline unsigned mpz_vector_bitlength(const MPZ_vector_1& v) {

    unsigned bits(0), i;
    for(i = 0; i < v.size(); i++) {
      bits = std::max(bits, (unsigned)ILOG2_MPZ(&v[i]));
    }
    return bits;
}

//! computes maximal bitlength of a vector
template < class NT >
unsigned poly_bitlength(const CGAL::Polynomial< NT >& f) {

    MPZ_traits< NT > _;
    unsigned bits(0), i;
    for(i = 0; i <= f.degree(); i++) {
        bits = std::max(bits, ILOG2(f[i]));
    }
//     std::cout << "############### BITS: " << bits << "\n\n";
    return bits;
}

// log |f|_2
template < class NT >
unsigned log_poly_norm2(const CGAL::Polynomial< NT >& f) {

    typedef typename CGAL::Get_arithmetic_kernel< NT >::Arithmetic_kernel::
        Bigfloat_interval BFI;

    typedef typename CGAL::Bigfloat_interval_traits< BFI >::Bound BigFloat;
    // use low precision for estimates
    long oldp = CGAL::set_precision(BFI(), 64);
    typename Real_embeddable_extension< BigFloat >::Ceil_log2_abs log2_abs;

    MPZ_traits< NT > _;
    int log_lc = ILOG2(f.lcoeff());
    int n = f.degree(), k;
    BigFloat norm(0);
    for(k = 0; k <= n; k++) {
        BigFloat X = CGAL::upper(CGAL::convert_to_bfi(f[k]));
        norm += X*X;
    }
    unsigned lgg = log2_abs(norm);
    CGAL::set_precision(BFI(), oldp);
    return lgg / 2;
}

// computes upper bound on roots of f
template < class NT >
int log_UB1(const CGAL::Polynomial< NT >& f) {

    MPZ_traits< NT > _;
    int log_lc = ILOG2(f.lcoeff());
    int n = f.degree(), k, log_bnd(-(1<<30));
    for(k = 1; k <= n; k++) {

        NT X = CGAL::abs(f[n-k]);
        int logc = ILOG2(X);
        logc = (logc - log_lc) / k;
//             std::cout << " logc: " << logc << "; " << log_lc << "\n";
        log_bnd = std::max(log_bnd, logc);
    }
    log_bnd += 1; // mul by 2
    return log_bnd;
}

// computes upper bound on roots of f (Cauchy bound)
template < class NT >
int log_UB2(const CGAL::Polynomial< NT >& f) {

    MPZ_traits< NT > _;
    int log_lc = ILOG2(f.lcoeff());
    int n = f.degree(), k, log_bnd(-(1<<30));
    for(k = 1; k < n; k++) {

        NT X = CGAL::abs(f[n-k]);
        int logc = ILOG2(X);
        logc = (logc - log_lc) / k;
        log_bnd = std::max(log_bnd, logc);
    }
    // add the last term f[0] / (f[n] * 2)
    int logc = ILOG2(f[0]);
    logc = (logc - log_lc - 1) / n;
    log_bnd = std::max(log_bnd, logc) + 1; // and mul by 2
    return log_bnd;
}

// degree-aware bound for poly factors
template < class NT >
unsigned poly_deg_bound(const CGAL::Polynomial< NT >& f, unsigned gcd_deg,
            unsigned lcgcd_bits) {

//     typedef typename CGAL::Get_arithmetic_kernel< NT >::Arithmetic_kernel::
//        Rational Rational;
//     typedef typename CGAL::Algebraic_kernel_d_1_generator
//         < NT, Rational>::Algebraic_kernel_with_qir_and_descartes_1 AK1;
// 
//     typedef typename AK1::Algebraic_real_1 Algebraic_real_1;
//     std::vector< std::pair< Algebraic_real_1, int > > r_f;
//     typename AK1::Solve_1 solve;

//     solve(f, std::back_inserter(r_f));

    std::cout << "\npoly bitlength: " << poly_bitlength(f) << "\n";
    int logUB1 = log_UB1(f), logUB2 = log_UB2(f),
            min_logUB = std::min(logUB1, logUB2);
    std::cout << "\nlogUB1: " << logUB1 << "; logUB2: " << logUB2 << "\n";
//! from "INEQUALITIES ON POLYNOMIAL HEIGHTS"
//! can use Proposition 2.3 which computes bounds using upper bound on
//! polynomial roots

//     std::cout << "\nroots: ";
//     double fmax(0);
//     for(unsigned i = 0; i < r_f.size(); i++) {
//         double root = CGAL::to_double(r_f[i].first);
//         fmax = std::max(fmax, std::abs(root));
// //         std::cout << root << "  ";
//     }
//     std::cout << "\n###### max_root: " << fmax << "\n";
    if(min_logUB <= 0)
        min_logUB = 1;
    unsigned binom_bnd = gcd_deg * min_logUB;
//                     log2(fmax);
    unsigned lognorm2 = log_poly_norm2(f);
    std::cout << "binom_bnd bits: " << (binom_bnd + lcgcd_bits) << "\n";
    std::cout << "log_poly_norm: " << lognorm2 <<
        "; mignotte_bnd bits: " << (lognorm2 + gcd_deg) << "\n";

    return 0;
}

template < class NT >
void gcd_alt_bounds(const CGAL::Polynomial< NT >& f,
        const CGAL::Polynomial< NT >& g, unsigned gcd_deg, unsigned& bits) {

    MPZ_traits< NT > _;

    CGAL::Polynomial< NT > gcd_truth = CGAL::gcd(f,g);
    gcd_deg = gcd_truth.degree();
    std::cout <<  "gcd_fg_deg: " << gcd_deg <<
            ":\n GCD bitlength: " <<  poly_bitlength(gcd_truth) << "\n";

    NT lcgcd = CGAL::gcd(f.lcoeff(), g.lcoeff());
    unsigned lcgcd_bits = ILOG2(lcgcd);
    std::cout << "lcgcd bits: " << lcgcd_bits << ":\n";

    std::cout << "\n====================== f ========================\n";
    poly_deg_bound(f, gcd_deg, lcgcd_bits);

    std::cout << "\n====================== g ========================\n";
    poly_deg_bound(g, gcd_deg, lcgcd_bits);

    CGAL::Polynomial< NT > frev = f, grev = g;
    frev.reversal(); grev.reversal();
    NT lcgcd_rev = CGAL::gcd(frev.lcoeff(), grev.lcoeff());
    unsigned lcgcd_rev_bits = ILOG2(lcgcd_rev);
    std::cout << "lcgcd_rev bits: " << lcgcd_bits << ":\n";
    
    std::cout << "\n====================== f rev ========================\n";
    poly_deg_bound(frev, gcd_deg, lcgcd_rev_bits);

    std::cout << "\n====================== g rev ========================\n";
    poly_deg_bound(grev, gcd_deg, lcgcd_rev_bits);
}

template < class NT >
void compute_gcd_bitlength(
        const CGAL::Polynomial< NT >& f, const CGAL::Polynomial< NT >& g,
        unsigned& bits) {

    MPZ_traits< NT > _;

    unsigned n = f.degree(), m = g.degree(), i, bits_f(0), bits_g(0);
    for(i = 0; i <= n; i++) {
        unsigned b = ILOG2(f[i]);
        bits_f = std::max(bits_f, b);
    }
    for(i = 0; i <= m; i++) {
        unsigned b = ILOG2(g[i]);
        bits_g = std::max(bits_g, b);
    }
//     min of two bitlength might not be enough
    bits = std::max(bits_f, bits_g);
    bits = std::max(bits, 80u); // at least 80 bits

//     printf("bits_f: %d; bits_g: %d ##gcdbits: %d\n", bits_f, bits_g, bits);
}

#undef ILOG2

//! computes Hamadard's bounds on resultant bitlength, lowest and highest
//! degree
template < class NT >
void compute_resultant_bounds(
        const CGAL::Polynomial< CGAL::Polynomial< NT > >& f,
        const CGAL::Polynomial< CGAL::Polynomial< NT > >& g,
        unsigned& low_deg, unsigned& high_deg, unsigned& bits) {

    typedef typename CGAL::Get_arithmetic_kernel< NT >::Arithmetic_kernel::
        Bigfloat_interval BFI;

    typedef typename CGAL::Bigfloat_interval_traits< BFI >::Bound BigFloat;
    typename Real_embeddable_extension< BigFloat >::Ceil_log2_abs log2_abs;

    unsigned n = f.degree(), m = g.degree(), i, j, k;
    // collects sums of absolute values of polynomial coeffs
    std::vector< BigFloat > sum_f(n + 1, BigFloat(0)),
            sum_g(m + 1, BigFloat(0));

    //! assume n >= m ?
    BigFloat s(0), t;
    long rows_bits(0), cols_bits(0);
    unsigned cols_min_deg(0);
    unsigned rows_deg(0), cols_deg(0), degx_f(0), degx_g(0);

// use low precision for estimates
//     t = CGAL::upper(CGAL::convert_to_bfi(NT(1)));
    long oldp = CGAL::set_precision(BFI(), 53);
    CGAL::set_precision(BFI(), 64);

    unsigned t_degf(0), t_degg(0); // total degrees of f/g
    for(i = 0; i <= n; i++) {
        const CGAL::Polynomial< NT >& poly = f[i];
        degx_f = std::max(degx_f, (unsigned)poly.degree());
        t_degf = std::max(t_degf, i + poly.degree());

        NT cf(0);
        for(k = poly.degree(); (int)k >= 0 ; k--) {
            cf += CGAL::abs(poly[k]);
        }
        t = CGAL::upper(CGAL::convert_to_bfi(cf)), s += t * t;
        sum_f[i] = t;
    }

    rows_bits = log2_abs(s) * m; // the # of bits in fp-number
    rows_deg = degx_f * m;
    s = 0;
    for(i = 0; i <= m; i++) {
        const CGAL::Polynomial< NT >& poly = g[i];
        degx_g = std::max(degx_g, (unsigned)poly.degree());

        t_degg = std::max(t_degg, i + poly.degree());
        NT cf(0);
        for(k = 0; (int)k <= poly.degree(); k++) { 
            cf += CGAL::abs(poly[k]);
        }
        t = CGAL::upper(CGAL::convert_to_bfi(cf)), s += t * t;
        sum_g[i] = t;
    }
    if(!CGAL::is_zero(s))
        rows_bits += log2_abs(s) * n;
    rows_deg += degx_g * n;

    unsigned minf, ming, deg(0), deg2(0);
    for(i = 0, minf = n, ming = m; i < n + m; i++, minf--, ming--) {
        s = 0, j = 0, k = minf, deg = 0, deg2 = -1u;
        if((int)k < 0) {
            j -= k, k = 0;
        }
        for(; j < m && k <= n; j++, k++) {
            const CGAL::Polynomial< NT >& poly = f[k];
            deg = std::max(deg, (unsigned)poly.degree());

            unsigned l = 0;
            while((int)l <= poly.degree()) {
                if(poly[l] != NT(0))
                    break;
                l++;
            }
            deg2 = std::min(deg2, l); // deg2 is a minimal resultant degree
            t = sum_f[k], s += t * t;
        }

        j = 0, k = ming;
        if((int)k < 0) {
            j -= k, k = 0;
        }
        for(; j < n && k <= m; j++, k++) {
            const CGAL::Polynomial< NT >& poly = g[k];
            deg = std::max(deg, (unsigned)poly.degree());

            unsigned l = 0;
            while((int)l <= poly.degree()) {
                if(poly[l] != NT(0))
                    break;
                l++;
            }
            deg2 = std::min(deg2, l); // deg2 is a minimal resultant degree
            // NOTE: what if 's' does not fit a double value ??
            t = sum_g[k], s += t * t;
        }
        if(!CGAL::is_zero(s))
            cols_bits += log2_abs(s);  
        cols_deg += deg;
        cols_min_deg += deg2;
    }
    rows_deg = std::min(rows_deg, cols_deg);
    rows_bits = std::min(rows_bits, cols_bits) / 2;

    unsigned bez_deg = t_degf * t_degg,
            new_deg = (n + m) * std::max(degx_f, degx_g);
    printf("rows/cols deg: %d bez_deg: %d new_deg: %d\n", rows_deg, bez_deg,
                new_deg);
    rows_deg = std::min(rows_deg, bez_deg);
    rows_deg = std::min(rows_deg, new_deg);

    CGAL::set_precision(BFI(), oldp);
    printf("resultant bounds: degree = %d; height = %d\n\n",
            rows_deg, rows_bits);

    bits = (unsigned)rows_bits, low_deg = cols_min_deg, high_deg = rows_deg;

    if(high_deg >= 1024*4 - 30)
        high_deg = 1024*4 - 30;

} // compute_resultant_bounds

extern gmp_randstate_t rands;

template < class NT >
struct Rand_coeff {
};

template < >
struct Rand_coeff< CORE::BigInt > {

    void operator()(CORE::BigInt& x, unsigned nbits) {
        CORE::BigInt y, maxx(1);
        mpz_urandomb(y.get_mp(), rands, nbits);
        y -= (maxx << (nbits-1));
        x = y;
    }
};


//! sparse bivariate polynomial: has at most n non-zero entries
//! \c nonzero_coeffs - generate with non-zero leading/trailing coefficients
//! \c RandCoeff: an object generating random coefficients of type \c NT
//! with magnitude / bitlength given by \c coeff_mod
template < class NT, class RandCoeff >
CGAL::Polynomial < CGAL::Polynomial < NT > >
    generate_sparse_random_poly2(unsigned deg_y, unsigned deg_x,
         unsigned n, unsigned coeff_mod, bool nonzero_coeffs = true) {

    typedef CGAL::Polynomial< NT > Poly_1;
    typedef CGAL::Polynomial< Poly_1 > Poly_2;

    typedef CGAL::Polynomial_traits_d< Poly_2 > PT;
    typename PT::Construct_polynomial cons; 

    typedef std::pair< CGAL::Exponent_vector, NT > Monomial;
    typedef std::vector< Monomial > MVector;

    MVector mv;
    RandCoeff rand_coeff;

    std::set< unsigned > unique;
    CGAL::Exponent_vector ev(0, 0);
    if(nonzero_coeffs) {
        NT cf;
        rand_coeff(cf, coeff_mod);
        ev[0] = deg_x; ev[1] = deg_y;
        mv.push_back(Monomial(ev, (cf == NT(0) ? NT(1) : cf)));
        unique.insert(deg_x * (deg_y + 1) + deg_y);

        rand_coeff(cf, coeff_mod);
        ev[0] = 0; ev[1] = 0;
        mv.push_back(Monomial(ev, (cf == NT(0) ? NT(1) : cf)));
        n--;
        unique.insert(0);
    }

    for(int i = 0; i < (int)n; i++) {

        unsigned ix = lrand48() % (deg_x + 1),
                 iy = lrand48() % (deg_y + 1);

        unsigned key = ix * (deg_y + 1) + iy;
        if(!unique.insert(key).second) {
            i--; continue;
        }

        ev[0] = ix; ev[1] = iy;

        NT cf;
        rand_coeff(cf, coeff_mod);
//    std::cout << "ix: " << ix << "; iy: " << iy << "; cf: " << cf << "\n";
        mv.push_back(Monomial(ev, cf));
    }
    return cons(mv.begin(), mv.end());
}

template < class NT, class RandCoeff >
CGAL::Polynomial< NT > generate_sparse_random_poly1(unsigned deg_x,
         unsigned n, unsigned coeff_mod) {

    typedef CGAL::Polynomial < NT > Poly_1;
    typename Poly_1::Vector pv(deg_x + 1, NT(0));
    RandCoeff rand_coeff;

    NT cf;
    rand_coeff(cf, coeff_mod);
    pv[deg_x] = (cf == NT(0) ? NT(1) : cf);

    rand_coeff(cf, coeff_mod);
    pv[0] = (cf == NT(0) ? NT(1) : cf);

    for(unsigned i = 0; i < n - 1; i++) {

        unsigned ix = lrand48() % (deg_x + 1);
        if(ix == deg_x) {
            i--; continue;
        }
        rand_coeff(cf, coeff_mod);
        pv[ix] = cf;
    }

    return Poly_1(pv.begin(), pv.end());
}

} // namespace internal

} // namespace CGAL

#endif // GPU_ALGORITHM_TEMPLATES_H
