
#include <CGAL/config.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Polynomial.h>

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

#endif // CGAL_USE_GPU

} // namespace internal

// get rid of modular stuff
#if (defined CGAL_USE_GPU || defined CGAL_USE_NTL)

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
template <> inline
CGAL::Polynomial< CGAL::Gmpz > gcd(
        const CGAL::Polynomial< CGAL::Gmpz >& F_,
        const CGAL::Polynomial< CGAL::Gmpz >& G_) {
    return gcd_(F_,G_);
}
#endif

#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL

template <> inline
CGAL::Polynomial< CORE::BigInt > gcd(
        const CGAL::Polynomial< CORE::BigInt >& F_,
        const CGAL::Polynomial< CORE::BigInt >& G_) {
    return gcd_(F_,G_);
}

#endif

#endif  // CGAL_USE_??

} // namespace CGAL
