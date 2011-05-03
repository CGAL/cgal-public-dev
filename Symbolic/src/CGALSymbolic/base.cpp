

#include <CGAL/config.h>

#define CGAL_BISOLVE_CHECK_GPU_RESULTANTS_SANITY 0 // default 0
#define CGAL_BISOLVE_CHECK_GPU_GCDS_SANITY 0 // default 1

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>

#include <CGAL/symbolic_exports.h>
#include <CGAL/GPU_algorithms/GPU_algorithm_facade.h>

static int _VERBOSE_Symbolic_ = 1;

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

#endif // CGAL_USE_GPU

} // namespace internal

} // namespace CGAL

