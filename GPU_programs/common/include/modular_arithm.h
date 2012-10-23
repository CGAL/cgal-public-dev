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

#ifndef MODULAR_ARITHM_H
#define MODULAR_ARITHM_H

// #include <CGAL/config.h>
#include <gmp.h>
#include <iostream>
#include <vector>
#include <stdlib.h>

#include <string.h> // memcpy is defined here ??
#include <include/MOD_entry.h>
#include <include/zmod.h>
#include <math.h>

namespace CGAL {

namespace internal {

// twice the size of array for primality test (even numbers omitted)
const ulong small_prime_limit = 1024*1024;

typedef std::vector< bool > Bit_vector;
typedef std::pair< umod_t, umod_t > Int_pair;
typedef std::vector< Int_pair > Pair_vector;

static Bit_vector oddprime_array;

#if 0

void make_oddprime_bitarray(ulong n, Bit_vector& ba);
// Return next prime >= n.
// Return zero if table of primes is not big enough.
ulong next_small_prime(ulong n);

bool is_prime(umod_t x);

void get_next_mod(unsigned& startm, unsigned nskips);

#endif

//! generates moduli table of size \c nmods and saves it to the file \c fname
//! reports error in case moduli exhausted
void generate_mod_table(unsigned nmods, const char *fname);

//! returns x ^ -1 mod m
inline unsigned mod_inverse(unsigned x, unsigned m) {

//NOTE: m must be set using zmod::set_mod() prior to this operation !!!  
    return pow_mod(x, m - 2, m);
}

//! computes x ^ -1 mod m using Montgomery inverse algorithm
//! \c mu = -m ^ -1 mod beta (beta = 2^24 or 2^32 depending on mod bitlength)
unsigned mod_inverse_montgomery(unsigned x, unsigned m, unsigned mu);

unsigned load_moduli_set(std::vector< MOD_entry >& mod_table,
         std::vector< double >& mod_bits, const char *fname);

//! returns # of moduli required to recover \c nbits
//! returns 0 if moduli set exhausted
unsigned bits_to_moduli(unsigned nbits, const std::vector< double >& mod_bits);

//! takes a sequence of moduli \c mods and computes decreasing modular
//! inverses, i.e., for n_mods = 4: m4^-1 mod m3, (m4*m3)^-1 mod m2,
//! (m4*m3*m2)^-1 mod m1
//! \c r must have enough space for \c n_mods-1 residues
//! \c stride - memory stride to access moduli
//! precondition: m[0] > m[1] > m[2] > .. > m[n-1]
void mod_inverse_seq(unsigned *r, const unsigned *mods, unsigned n,
        unsigned stride = 1);

//! adds residue \c x modulo \c m to final result \c r
//! \c nr & \c nprod - length of r and prod resp. (in words)
//! \c inv_prod - \c prod^-1 mod m
//! attention: \c r & \c prod must have enough for residue and modulo resp.
void incremental_CRT(limb *r, unsigned& nr, limb *prod, unsigned& nprod,
        limb x, limb m, limb inv_prod, unsigned n_max);
#if 0
void CRA_bases(unsigned *r, unsigned r_stride,  const unsigned *mods,
        unsigned n_mods, unsigned m_stride = 1);

//! computes weights \c w for core function \c CM 
void core_weights(unsigned *w, unsigned w_stride, unsigned CM,
    const unsigned *bases, const unsigned *mods, unsigned n_mods,
    unsigned m_stride = 1);

//! computes evaluates core function \c CM for residue set \c r ,
//! CRA \c bases and set of core \c weights 
void core_function(const unsigned *r, unsigned r_stride, unsigned CM,
    const unsigned *weights, const unsigned *bases,
    const unsigned *mods, unsigned n_mods, unsigned m_stride = 1);

//!\c bases - CRA bases: (M / m[i])^(-1) mod m[i]
void RNS_magnitude(const unsigned *r, unsigned r_stride,
    const unsigned *bases, const unsigned *mods, unsigned n_mods,
        unsigned m_stride = 1);
#endif

void convert_to_RNS(unsigned *r, unsigned r_stride, const mpz_t mp_in,
        const unsigned *mods, unsigned n_mods, unsigned m_stride = 1);

//! reduces input \c in by \c n moduli \c mods , writes results to \c r
//! \c r must have enough space for \c n_mods residues
void convert_to_RNS_pure(unsigned *r, unsigned r_stride, const limb *in,
        unsigned n, const unsigned *mods, unsigned n_mods,
                unsigned m_stride = 1);
// given n residues modulo mods[i], computes associated mixed radix digits
// c[i] = (m[0]*m[1]*..*m[i-1])^-1 mod m[i]
// assume: mods[0] < mods[1] < .. < mods[n-1]
//! returns in \c y an MRS digit representation of a large integer given by its residues \c x ; \c x and \c y can overlap
void compute_MR_digits(unsigned *y, const unsigned *x, const unsigned *mods,
        const unsigned *c, unsigned n);

// the same as above but everything is stored in reversed order
//! x_stride: in-out data stride
//! m_stride: moduli stride
void compute_MR_digits_rev(unsigned *y, const unsigned *x, unsigned x_stride,
        const unsigned *mods, unsigned m_stride, const unsigned *c,
            unsigned n);


// extern gmp_randstate_t rands;

template < class NT >
void rand_coeff(NT&, unsigned) {
}

template < class NT >
void print_vector(const std::vector< NT >& v, bool no_sign = false) {
    printf("NYI\n");
}

template < >
inline void rand_coeff(zmod& x,  unsigned coeff_mod) {
    x = zmod(lrand48() % coeff_mod) - zmod(coeff_mod/2);
}

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

extern gmp_randstate_t rands;

template < >
inline void rand_coeff(CORE::BigInt& x, unsigned nbits) {

    CORE::BigInt y, maxx(1);
    mpz_urandomb(y.get_mp(), rands, nbits);
    y -= (maxx << (nbits-1));
    x = y;
}

template < >
void print_vector< CORE::BigInt >(const std::vector< CORE::BigInt >& v,
        bool hex);
#endif 

template < >
void print_vector< zmod >(const std::vector< zmod >& v, bool no_sign);

template < class NT >
std::ostream& operator <<(std::ostream& os, const std::vector< NT >& p) {

    for(unsigned i = p.size() - 1; i != -1u; i--) {
        NT t = p[i];
        if(t == NT(0)) 
            continue;
        
        if(i != p.size() - 1) {
            if(t > 0)
                os << " + ";
            else 
                os << " ";
        }
        os << t;

        if(i >= 1) {
            os << "*x";
            if(i > 1)
                os << "^" << i;
        }
    }
    return os;
}

std::ostream& operator <<(std::ostream& os, const std::vector< zmod >& p);

} // namespace internal

} // namespace CGAL

#endif // MODULAR_ARITHM_H
