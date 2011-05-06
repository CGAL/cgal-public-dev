

#include <CGAL/config.h>

#define CGAL_BISOLVE_CHECK_GPU_RESULTANTS_SANITY 0 // default 0
#define CGAL_BISOLVE_CHECK_GPU_GCDS_SANITY 0 // default 1

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>

#include <CGAL/symbolic_exports.h>

#if (defined CGAL_USE_GPU)
#include <CGAL/GPU_algorithms/GPU_algorithm_facade.h>
#endif

#if (defined CGAL_USE_NTL)
#include <NTL/ZZX.h>
#endif 

namespace CGAL {

namespace internal {

#if (defined CGAL_USE_GPU)

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
                           const CGAL::Polynomial<Coeff>& G_) {

    if(F_.degree() < 5 && G_.degree() < 5)
        return resultant_cpu(F_,G_);

    bool failed;
    Coeff res = GPU_algorithm_facade::resultant(F_,G_, &failed);

    if(failed) {
/*        writeout(F_,G_);
        throw "WTF!?";*/
        std::cout << "\nGRES failed!!\n";
        return resultant_cpu(F_,G_);
    }

#if CGAL_BISOLVE_CHECK_GPU_RESULTANTS_SANITY

    std::cout << "\nGRES sanity check!!\n";
    Coeff check = CGAL::canonicalize(res),
          truth = CGAL::canonicalize(resultant_cpu(F_,G_));

    if(check != truth) {
        std::cout << "\nWrong GRES:\n" << check << "\n" <<
            truth << "\n";
//         writeout(F_, G_);
        throw "WTF!?";
    }
#endif
    return res;
}

// partial specialization to enable gpu'ed resultants
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL

#warning !!!!!!!!!!!!!!!! gmp is present

template <>
CGAL::Polynomial< CGAL::Gmpz > resultant(
        const CGAL::Polynomial< CGAL::Polynomial< CGAL::Gmpz > >& F_,
        const CGAL::Polynomial< CGAL::Polynomial< CGAL::Gmpz > >& G_) {

    return resultant_gpu(F_, G_);
}
#endif

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

#warning !!!!!!!!!!!!!!!! core is present

template <>
CGAL::Polynomial< CORE::BigInt > resultant(
        const CGAL::Polynomial< CGAL::Polynomial< CORE::BigInt > >& F_,
        const CGAL::Polynomial< CGAL::Polynomial< CORE::BigInt > >& G_) {

    return resultant_gpu(F_, G_);
}
#endif

template < class Poly >
Poly gcd_gpu(const Poly& F_, const Poly& G_) {

    std::cout << "--------- Using our lovely libs!! -------------\n";

//     n_ggcd_calls++;
//     tm_ggcd.start();
    Poly ggcd = GPU_algorithm_facade::gcd(F_,G_);
//     tm_ggcd.stop();

#if CGAL_BISOLVE_CHECK_GPU_GCDS_SANITY
    printf("\nGGCD sanity check..\n");
    Poly truth = /*gcd_NTL(F_, G_);*/
                    modular_gcd_utcf_dfai(F_, G_);
    if(truth != ggcd) {
        writeout(F_, G_);
        std::cerr << "Wrong gcd!!\n";
        throw "WTF?";
    }
    printf("GGCD ok\n");
#endif
    return ggcd;
}

// partial specialization to enable gpu'ed gcds
#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
template <> 
CGAL::Polynomial< CGAL::Gmpz > gcd_(
        const CGAL::Polynomial< CGAL::Gmpz >& F_,
        const CGAL::Polynomial< CGAL::Gmpz >& G_) {

    return gcd_gpu(F_, G_);
}

template <> 
CGAL::Polynomial< CGAL::Gmpz > gcd_utcf_(
        const CGAL::Polynomial< CGAL::Gmpz >& F_,
        const CGAL::Polynomial< CGAL::Gmpz >& G_) {

    return gcd_gpu(F_, G_);
}
#endif

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

template <>
CGAL::Polynomial< CORE::BigInt > gcd_(
        const CGAL::Polynomial< CORE::BigInt >& F_,
        const CGAL::Polynomial< CORE::BigInt >& G_) {

    return gcd_gpu(F_, G_);
}

template <> 
CGAL::Polynomial< CORE::BigInt > gcd_utcf_(
        const CGAL::Polynomial< CORE::BigInt >& F_,
        const CGAL::Polynomial< CORE::BigInt >& G_) {

    return gcd_gpu(F_, G_);
}
#endif

#endif // CGAL_USE_GPU

#if (defined CGAL_USE_NTL)

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

template < class NT > inline
Polynomial< NT > gcd_NTL(const Polynomial< NT >& p1,
                         const Polynomial< NT >& p2) {

    typedef Polynomial< NT > Poly_1;

    if(p1 == Poly_1(0))
        return p2;
    if(p2 == Poly_1(0))
        return p1;

    NTL::ZZX r, q1, q2;
    poly2ntl(p1, q1);
    poly2ntl(p2, q2);

    NTL::GCD(r, q1, q2);
    return ntl2poly< NT >(r);
}

template < class NT > inline
Polynomial< NT > div_NTL(
        const Polynomial< NT >& f, const Polynomial< NT >& g,
        Polynomial< NT >& q, Polynomial< NT >& r) {

    NTL::ZZX q1, r1, f1, g1;
    internal::poly2ntl(f, f1);
    internal::poly2ntl(g, g1);

    NTL::DivRem(q1, r1, f1, g1);

    q = internal::ntl2poly< NT >(q1);
    r = internal::ntl2poly< NT >(r1);
}

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL

#if !(defined CGAL_USE_GPU) // give preference to gpu impl
template <> 
Polynomial< CGAL::Gmpz > gcd_(const Polynomial< CGAL::Gmpz >& p1,
                              const Polynomial< CGAL::Gmpz >& p2) {

    return gcd_NTL(p1, p2);
}
#endif

#endif // CGAL_HAS_GMP_ARITHMETIC_KERNEL

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

#if !(defined CGAL_USE_GPU) // give preference to gpu impl
template <> 
Polynomial< CORE::BigInt > gcd_(const Polynomial< CORE::BigInt >& p1,
                                const Polynomial< CORE::BigInt >& p2) {

    return gcd_NTL(p1, p2);
}
#endif

#endif // CGAL_HAS_CORE_ARITHMETIC_KERNEL

#endif // CGAL_USE_NTL

} // namespace internal

#if (defined CGAL_USE_NTL)

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL

template <>
void Polynomial< CGAL::Gmpz >::euclidean_division(
    const Polynomial< CGAL::Gmpz >& f, const Polynomial< CGAL::Gmpz >& g,
    Polynomial< CGAL::Gmpz >& q, Polynomial< CGAL::Gmpz >& r) {

  internal::div_NTL(f, g, q, r);
}

#endif // CGAL_HAS_GMP_ARITHMETIC_KERNEL

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

template <>
void Polynomial< CORE::BigInt >::euclidean_division(
    const Polynomial< CORE::BigInt >& f, const Polynomial< CORE::BigInt >& g,
    Polynomial< CORE::BigInt >& q, Polynomial< CORE::BigInt >& r) {

  internal::div_NTL(f, g, q, r);
}

#endif // CGAL_HAS_CORE_ARITHMETIC_KERNEL

#endif // defined CGAL_USE_NTL

} // namespace CGAL

