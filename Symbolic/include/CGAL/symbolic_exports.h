
#ifndef CGAL_SYMBOLIC_EXPORTS_H
#define CGAL_SYMBOLIC_EXPORTS_H

#if (defined CGAL_USE_GPU) || (defined CGAL_USE_NTL)
#warning Switching off CGAL`s modular filter !!
#define CGAL_MODULAR_FILTER_OFF     // do not use modular filter
#endif

#include <CGAL/config.h>

#include <cstring>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>

#if (defined CGAL_USE_NTL)
#include <NTL/ZZX.h>
#include <NTL/ZZXFactoring.h>
#endif 

#if (defined CGAL_USE_GPU)
#include <CGAL/GPU_algorithms/GPU_algorithm_facade.h>
#endif

namespace CGAL {

namespace internal {

#if (defined CGAL_USE_GPU)

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL

template <> 
CGAL::Polynomial< CGAL::Gmpz > resultant(
        const CGAL::Polynomial< CGAL::Polynomial< CGAL::Gmpz > >& F_,
        const CGAL::Polynomial< CGAL::Polynomial< CGAL::Gmpz > >& G_);
#endif

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

template <> 
CGAL::Polynomial< CORE::BigInt > resultant(
        const CGAL::Polynomial< CGAL::Polynomial< CORE::BigInt > >& F_,
        const CGAL::Polynomial< CGAL::Polynomial< CORE::BigInt > >& G_);
#endif

#endif // CGAL_USE_GPU

#if (defined CGAL_USE_GPU) || (defined CGAL_USE_NTL)

// partial specialization to enable gpu'ed gcds
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL

template <> 
CGAL::Polynomial< CGAL::Gmpz > gcd_(
        const CGAL::Polynomial< CGAL::Gmpz >& F_,
        const CGAL::Polynomial< CGAL::Gmpz >& G_);

template <> 
CGAL::Polynomial< CGAL::Gmpz > gcd_utcf_(
        const CGAL::Polynomial< CGAL::Gmpz >& F_,
        const CGAL::Polynomial< CGAL::Gmpz >& G_);
#endif

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

template <>
CGAL::Polynomial< CORE::BigInt > gcd_(
        const CGAL::Polynomial< CORE::BigInt >& F_,
        const CGAL::Polynomial< CORE::BigInt >& G_);

template <> 
CGAL::Polynomial< CORE::BigInt > gcd_utcf_(
        const CGAL::Polynomial< CORE::BigInt >& F_,
        const CGAL::Polynomial< CORE::BigInt >& G_);

#endif

#endif // CGAL_USE_GPU or CGAL_USE_NTL

} // namespace internal

#if (defined CGAL_USE_NTL)

namespace internal {

struct NTL_bigint_rep {
  long alloc;
  long size;
  mp_limb_t data;
};

template < class Integer_ >
struct Integer_vs_mpz {
  // dummy
  typedef Integer_ Integer;
};

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

template < >
struct Integer_vs_mpz< CORE::BigInt > {

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
struct Integer_vs_mpz< CGAL::Gmpz > {

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
    Integer_vs_mpz< NT > cvt; // instantiate secret object
    
    int i;
    typename Polynomial< NT >::const_iterator pit;
    for(i = 0, pit = p.begin(); pit != p.end(); pit++, i++) {

        // TODO avoid double code (see other conversions)

        NTL::ZZ& zz = q.rep[i];
        mpz_srcptr tmp = cvt(*pit);
        int sz = tmp->_mp_size;
        if(sz == 0)
            continue;
        if(sz < 0)
            sz = -sz;
        
        zz.SetSize(sz);
        NTL_bigint_rep *rep = reinterpret_cast< NTL_bigint_rep* >(zz.rep);
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
      
        // TODO avoid double code (see other conversions)

        const NTL::ZZ& zz = r.rep[i];
        if(NTL::IsZero(zz)) {
            v[i] = NT(0);
            continue;
        } 

        NTL_bigint_rep *rep = reinterpret_cast< NTL_bigint_rep* >(zz.rep);
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

} // namespace internal

//!\brief factorizes an integral polynomial (CGAL::Gmpz, CORE::BigInt) with NTL's factor method
template < class Coeff, class OutputIterator >
inline Coeff factorize_NTL(const Polynomial< Coeff >& p, 
                           OutputIterator factors) {
  
  //std::cout << "FACTORIZE NTL ..." << std::flush;

  NTL::ZZX np;
  internal::poly2ntl(p, np);
  
  NTL::vec_pair_ZZX_long nfactors;
  NTL::ZZ nc;

  NTL::factor(nc, nfactors, np);

  for (int i = 0; i < nfactors.length(); i++) {
    *factors++ = std::make_pair(internal::ntl2poly< Coeff >(nfactors[i].a), nfactors[i].b);
  }

  //std::cout << "done." << std::endl;

  // TODO avoid double code (see other conversions)

  if (NTL::IsZero(nc)) {
    return Coeff(0);
  }
  
  mpz_t tmp;
  mpz_init(tmp);
  
  internal::NTL_bigint_rep *rep = reinterpret_cast< internal::NTL_bigint_rep* >(nc.rep);
  int sz = rep->size;
  if(sz < 0)
    sz = -sz;
  
  mpz_realloc2(tmp, sz * GMP_NUMB_BITS);
  tmp->_mp_size = rep->size;
  memcpy(tmp->_mp_d, &rep->data, sz*sizeof(mp_limb_t));
  
  return Coeff(tmp);

}


#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
  
template <>
void Polynomial< CGAL::Gmpz >::euclidean_division(
    const Polynomial< CGAL::Gmpz >& f, const Polynomial< CGAL::Gmpz >& g,
    Polynomial< CGAL::Gmpz >& q, Polynomial< CGAL::Gmpz >& r);

#endif // CGAL_HAS_GMP_ARITHMETIC_KERNEL

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

template <>
void Polynomial< CORE::BigInt >::euclidean_division(
    const Polynomial< CORE::BigInt >& f, const Polynomial< CORE::BigInt >& g,
    Polynomial< CORE::BigInt >& q, Polynomial< CORE::BigInt >& r);

#endif

#endif // CGAL_USE_NTL

// get rid of modular stuff
// #if (defined CGAL_USE_GPU || defined CGAL_USE_NTL)
// 
// #ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
// template <> inline
// CGAL::Polynomial< CGAL::Gmpz > gcd(
//         const CGAL::Polynomial< CGAL::Gmpz >& F_,
//         const CGAL::Polynomial< CGAL::Gmpz >& G_) {
//     return gcd_(F_,G_);
// }
// #endif
// 
// #ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
// 
// template <> inline
// CGAL::Polynomial< CORE::BigInt > gcd(
//         const CGAL::Polynomial< CORE::BigInt >& F_,
//         const CGAL::Polynomial< CORE::BigInt >& G_) {
//     return gcd_(F_,G_);
// }
// 
// #endif
// 
// #endif  // CGAL_USE_??

} // namespace CGAL

#endif // CGAL_SYMBOLIC_EXPORTS_H
