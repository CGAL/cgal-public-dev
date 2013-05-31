

#include <CGAL/config.h>

#define CGAL_BISOLVE_CHECK_GPU_RESULTANTS_SANITY 0 // default 0
#define CGAL_BISOLVE_CHECK_GPU_GCDS_SANITY 0 // default 1

#include <CGAL/symbolic_exports.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>

#if (defined CGAL_USE_GPU)
#include <CGAL/GPU_algorithms/GPU_algorithm_facade.h>
#endif

#if (defined CGAL_USE_NTL)
#include <NTL/ZZX.h>
#include <NTL/ZZXFactoring.h>
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
                           const CGAL::Polynomial<Coeff>& G_,
                            bool zero_when_failed = false) {

    if(F_.degree() < 5 && G_.degree() < 5)
        return resultant_cpu(F_,G_);

    bool failed;
    Coeff res = GPU_algorithm_facade::resultant(F_,G_, &failed);

    if(failed) {
        if(zero_when_failed)
            return Coeff(0);

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
	writeout(F_, G_);
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

#if (defined CGAL_USE_NTL)
template < class NT > inline
Polynomial< NT > gcd_NTL(const Polynomial< NT >& p1,
                         const Polynomial< NT >& p2);
#endif

template < class Poly >
Poly gcd_gpu(const Poly& F_, const Poly& G_) {
    bool failed;
    Poly ggcd = GPU_algorithm_facade::gcd(F_, G_, &failed);

    if(failed) {
        std::cerr << "GGCD failed!\n";
        return gcd_utcf_UFD(F_, G_);
    }

#if CGAL_BISOLVE_CHECK_GPU_GCDS_SANITY
    std::cerr << "\nGGCD sanity check..\n";
    Poly truth = gcd_NTL(F_, G_);
//                     modular_gcd_utcf_dfai(F_, G_);
    if(truth != ggcd) {
        std::cerr << "truth: " << truth << "\n\nggcd: " << ggcd << "\n";
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
void div_NTL(
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

template <> 
Polynomial< CGAL::Gmpz > gcd_utcf_(const Polynomial< CGAL::Gmpz >& p1,
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

template <> 
Polynomial< CORE::BigInt > gcd_utcf_(const Polynomial< CORE::BigInt >& p1,
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

