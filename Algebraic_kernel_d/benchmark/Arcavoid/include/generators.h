/*
** generators.h
** 
** Made by Alexander Kobel
** Login   <perpeduumimmobile@lidinoid>
** 
** Started on  Tue Apr 20 16:02:30 2010 Alexander Kobel
** Last update Tue Apr 20 16:02:30 2010 Alexander Kobel
*/

#ifndef   	GENERATORS_H_
# define   	GENERATORS_H_

#include <cmath>
#include <vector>

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/GMP_arithmetic_kernel.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Cartesian_complex.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace CGAL {

template< typename NT >
class Random_polynomial_generator {
};

#ifdef CGAL_USE_GMP

template<>
class Random_polynomial_generator< Gmpz > {
public:
  typedef Gmpz                           NT;
  typedef Polynomial< NT >               Poly;
  typedef Polynomial_traits_d< Poly >    PT;
  typedef std::vector< NT >              Vector;
  typedef Polynomial< Poly >             Poly_2;
  typedef Polynomial_traits_d< Poly_2 >  PT_2;
  typedef std::vector< Poly >            Poly_vector;
  typedef std::vector< Vector >          Vector_2;

  typedef Polynomial< Poly_2 >           Poly_3;
  typedef Polynomial_traits_d< Poly_3 >  PT_3;

  typedef Cartesian_complex< NT >        CNT;
  typedef Polynomial< CNT >              CPoly;
  typedef Polynomial_traits_d< CPoly >   CPT;
  typedef std::vector< CNT >             CVector;
  typedef Polynomial< CPoly >            CPoly_2;
  typedef Polynomial_traits_d< CPoly_2 > CPT_2;
  typedef std::vector< CPoly >           CPoly_vector;
  typedef std::vector< CVector >         CVector_2;

  typedef std::vector< size_t > Size_t_vector;
  typedef unsigned long         mp_bitcnt_t;

  gmp_randstate_t rstate;

  typedef boost::mt19937 BaseRNG;
  typedef boost::normal_distribution< double >
  DoubleNormalDistribution;
	typedef boost::variate_generator< BaseRNG, DoubleNormalDistribution >
  DoubleNormalGenerator;

  BaseRNG base_rng;
  DoubleNormalGenerator normal_gen;

  Random_polynomial_generator (boost::mt19937::result_type seed = 42)
    : base_rng (seed),
      normal_gen(base_rng, DoubleNormalDistribution()) {
    gmp_randinit_mt (rstate);
    gmp_randseed_ui (rstate, seed);
  }

  const int sign() {
    return (gmp_urandomb_ui (rstate, 1) == 1 ? 1 : -1);
  }

  const size_t log2 (size_t n) {
    size_t r = 0;
    while (n /= 2) ++r;
    return r;
  }

  const Size_t_vector choose (size_t n, size_t k) {
    Size_t_vector a (n);

    for (size_t i = 0; i < n; ++i)
      a[i] = i;

    for (size_t i = n-1; i > 0; --i) {
      size_t j = gmp_urandomm_ui (rstate, n);
      std::swap (a[i], a[j]);
    }

    return Size_t_vector (a.begin(), a.begin()+k);
  }
  
  const Poly uniform (size_t n, mp_bitcnt_t L) {
    Vector fc (n+1);
    for (size_t i = 0; i <= n; ++i) fc[i] = 0;

    if (L < 2) L = 2;

    for (size_t i = 0; i <= n; ++i) {
      mpz_urandomb (fc[i].mpz(), rstate, L);
      if (sign() != 1) fc[i] = - fc[i];
    }

    while (fc[n] == 0)
      mpz_urandomb (fc[n].mpz(), rstate, L);
    if (sign() != 1) fc[n] = - fc[n];

    return CGAL::canonicalize (Poly (fc.begin(), fc.end()));
  }

  const Poly_2 uniform_bivariate (size_t n, mp_bitcnt_t L) {
    Poly_vector fc (n+1);

    for (size_t i = 0; i <= n; ++i)
      fc[i] = uniform (n-i, L);
    
    return CGAL::canonicalize (Poly_2 (fc.begin(), fc.end()));
  }

	const Poly fahramand (size_t n, mp_bitcnt_t L) {
    Vector fc (n+1, 0);
    double scale = std::pow (2., L);
    
    for (size_t i = 0; i <= n; ++i) fc[i] = 0;

    double bin = n;
    for (size_t i = 0; i <= n; ++i) {
			mpz_set_d (fc[i].mpz(), scale * normal_gen() * std::sqrt (bin / (i+1)));
			bin *= n-i;
			bin /= i+1;
    }

		return CGAL::canonicalize (Poly (fc.begin(), fc.end()));
	}

  const Poly uniform_sparse (size_t n, mp_bitcnt_t L) {
    Vector fc (n+1, 0);
    for (size_t i = 0; i <= n; ++i) fc[i] = 0;

    size_t k = log2(n) * 2;
    if (k > n) k = n;
    Size_t_vector a = choose (n+1, k);
    a.push_back (0);

    if (L < 2) L = 2;
    
    for (size_t i = 0; i <= k; ++i) {
      mpz_urandomb (fc[a[i]].mpz(), rstate, L);
      if (sign() != 1) fc[a[i]] = - fc[a[i]];
    }
    
    while (fc[n] == 0)
      mpz_urandomb (fc[n].mpz(), rstate, L);
    if (sign() != 1) fc[n] = - fc[n];
    
    return CGAL::canonicalize (Poly (fc.begin(), fc.end()));
  }

  const CPoly complex_uniform_sparse (size_t n, mp_bitcnt_t L) {
    const Poly f = uniform_sparse (n, L);
    return CPoly (f.begin(), f.end());
  }

  const Poly_2 uniform_sparse_bivariate (size_t n, mp_bitcnt_t L) {
    Poly_vector fc (n+1);

    for (size_t i = 0; i <= n; ++i)
      fc[i] = uniform_sparse (n-i, L);
    
    return CGAL::canonicalize (Poly_2 (fc.begin(), fc.end()));
  }
  const CPoly_2 complex_uniform_sparse_bivariate (size_t n, mp_bitcnt_t L) {
    CPoly_vector fc (n+1);

    for (size_t i = 0; i <= n; ++i)
      fc[i] = complex_uniform_sparse (n-i, L);
    
    return CPoly_2 (fc.begin(), fc.end());
  }

  const Poly monic_uniform (size_t n, mp_bitcnt_t L) {
    Vector fc (n+1);
    for (size_t i = 0; i <= n; ++i) fc[i] = 0;

    for (size_t i = 0; i < n; ++i) {
      mpz_urandomb (fc[i].mpz(), rstate, L);
      if (sign() != 1) fc[i] = - fc[i];
    }
    fc[n] = NT(1);

    return CGAL::canonicalize (Poly (fc.begin(), fc.end()));
  }

  const Poly monic_uniform_sparse (size_t n, mp_bitcnt_t L) {
    Vector fc (n+1, 0);
    for (size_t i = 0; i <= n; ++i) fc[i] = 0;

    size_t k = log2(n) * 2;
    if (k > n) k = n;
    Size_t_vector a = choose (n+1, k);
    a.push_back (0);

    if (L < 2) L = 2;
    
    for (size_t i = 0; i <= k; ++i) {
      mpz_urandomb (fc[a[i]].mpz(), rstate, L);
      if (sign() != 1) fc[a[i]] = - fc[a[i]];
    }
    fc[n] == 1;
    
    return CGAL::canonicalize (Poly (fc.begin(), fc.end()));
  }

  const Poly mignotte (size_t n, mp_bitcnt_t L) {
    Vector fc (n+1);
    for (size_t i = 0; i <= n; ++i) fc[i] = NT(0);

    mpz_urandomb (fc[1].mpz(), rstate, L / 2);
    fc[n] = NT(1);
    fc[2] = -CGAL::square (fc[1]);
    fc[1] *= NT(2);
    fc[0] = -NT(1);

    return CGAL::canonicalize (Poly (fc.begin(), fc.end()));
  }

  const Poly k_mignotte (size_t n, mp_bitcnt_t L, int k) {
    Vector fc (n+1);
    for (size_t i = 0; i <= n; ++i) fc[i] = NT(0);
    fc[n] = NT(1);

    Poly f (fc.begin(), fc.end());
    
    Poly m (1);
    for (int i = 0; i < k; ++i) {
      Gmpz a = 0;
      mpz_urandomb (a.mpz(), rstate, L / (2*k));
      m *= CGAL::ipower (CGAL::ipower (Poly (NT(0), NT(1)), k) - a, 2);
    }

    return CGAL::canonicalize (f + m);
  }

  const Poly bad_mignotte (size_t n, mp_bitcnt_t L) {
    const size_t nf = n/2;
    const size_t ng = n - nf;

    const mp_bitcnt_t Lf = L / 4;
    const mp_bitcnt_t Lg = L - Lf;

    Vector fc (nf+1, NT(0));
    for (size_t i = 0; i <= nf; ++i) fc[i] = 0;
    Vector gc (ng+1, NT(0));
    for (size_t i = 0; i <= ng; ++i) gc[i] = 0;

    mpz_urandomb (fc[1].mpz(), rstate, Lf);
    fc[nf] = NT(1);
    fc[2] = CGAL::square (fc[1]) * NT(2);
    fc[1] *= NT(4);
    fc[0] = NT(2);

    mpz_urandomb (gc[1].mpz(), rstate, Lg);
    gc[ng] = NT(1);
    gc[2] = CGAL::square (gc[1]) * NT(3);
    gc[1] *= NT(6);
    gc[0] = NT(3);

    return CGAL::canonicalize (Poly (fc.begin(), fc.end())
                               * Poly (gc.begin(), gc.end()));
  }

  const Poly resultant_of_two (size_t n, mp_bitcnt_t L) {
    std::vector< size_t > factors;
    size_t rem = n;
    for (size_t i = 2; i <= rem; ++i) {
      while (rem % i == 0) {
        rem /= i;
        factors.push_back (i);
      }
    }

    size_t nf = 1;
    size_t ng = n;

    for (size_t i = 0; i < factors.size() && nf < ng; ++i) {
      nf *= factors[i];
      ng /= factors[i];
    }

    mp_bitcnt_t Lf = L / 2 / nf;
    mp_bitcnt_t Lg = L / 2 / ng;

    const Poly_2 f = uniform_bivariate (nf, Lf);
    const Poly_2 g = uniform_bivariate (ng, Lg);

    return CGAL::resultant (f, g);
  }

  const Poly resultant_of_two_sparse (size_t n, mp_bitcnt_t L) {
    std::vector< size_t > factors;
    size_t rem = n;
    for (size_t i = 2; i <= rem; ++i) {
      while (rem % i == 0) {
        rem /= i;
        factors.push_back (i);
      }
    }

    size_t nf = 1;
    size_t ng = n;

    for (size_t i = 0; i < factors.size() && nf < ng; ++i) {
      nf *= factors[i];
      ng /= factors[i];
    }

    mp_bitcnt_t Lf = L / 2 / nf;
    mp_bitcnt_t Lg = L / 2 / ng;

    const Poly_2 f = uniform_sparse_bivariate (nf, Lf);
    const Poly_2 g = uniform_sparse_bivariate (ng, Lg);

    return CGAL::resultant (f, g);
  }

  const Poly_2 mandelbrot (size_t n) const {
    CPT_2::Shift cshift;
    const CPoly_2 x = cshift (CPoly_2 (1), 1, 0);
    const CPoly_2 y = cshift (CPoly_2 (1), 1, 1);

    std::vector< CPoly_2 > f (n+1);
    std::vector< CPoly_2 > g (n+1);

    f[0] = x + y * CNT (0,1);
    g[0] = x - y * CNT (0,1);

    for (size_t i = 0; i < n; ++i) {
      f[i+1] = f[i] * f[i] + f[0];
      g[i+1] = g[i] * g[i] + g[0];
    }

    const CPoly_2 m = f[n] * g[n];
    std::list< std::pair< Exponent_vector, NT > > mc;

    mc.push_back (std::make_pair (Exponent_vector (0, 0), -NT(n+1)));
    for (size_t i = 0; i <= m.degree(); ++i)
      for (size_t j = 0; j <= m[i].degree(); ++j)
        mc.push_back (std::make_pair (Exponent_vector ((int)j, (int)i), m[i][j].real()));

    return PT_2::Construct_polynomial() (mc.begin(), mc.end());
  }

  const Poly mandelbrot_resultant (size_t n) const {
    // const Poly_2 f = mandelbrot (n);
    const Poly r = CGAL::resultant (mandelbrot (n), mandelbrot (n-1));
    return CGAL::make_square_free (r);
    // return CGAL::resultant (f, CGAL::differentiate (f));
  }

  const Poly_2 lemniscate (size_t n, mp_bitcnt_t L) {
    const CPoly_2 f = complex_uniform_sparse_bivariate (n/2, L/2);

    const CPoly_2 ff = CGAL::scale (f, CNT::I()) * CGAL::scale (f, -CNT::I());
    std::list< std::pair< Exponent_vector, NT > > ffc;

    ffc.push_back (std::make_pair (Exponent_vector (0, 0), -NT(n+1)));
    for (size_t i = 0; i <= ff.degree(); ++i)
      for (size_t j = 0; j <= ff[i].degree(); ++j)
        ffc.push_back (std::make_pair (Exponent_vector ((int)j, (int)i), ff[i][j].real()));

    return PT_2::Construct_polynomial() (ffc.begin(), ffc.end());
  }

  const Poly lemniscate_resultant (size_t n, mp_bitcnt_t L) {
    std::vector< size_t > factors;
    size_t rem = n;
    for (size_t i = 2; i <= rem; ++i) {
      while (rem % i == 0) {
        rem /= i;
        factors.push_back (i);
      }
    }

    size_t nf = 1;
    size_t ng = n;

    for (size_t i = 0; i < factors.size() && nf < ng; ++i) {
      nf *= factors[i];
      ng /= factors[i];
    }

    mp_bitcnt_t Lf = L / 2 / nf;
    mp_bitcnt_t Lg = L / 2 / ng;

    const Poly_2 f = lemniscate (nf, Lf);
    const Poly_2 g = lemniscate (ng, Lg);

    return CGAL::make_square_free (CGAL::resultant (f, g));
  }

  const Poly exp_taylor (size_t n) const {
    Vector fc (n+1);
    for (size_t i = 0; i <= n; ++i) fc[i] = 0;

    NT i_fac = 1;
    NT n_fac = 1;
    mpz_fac_ui (n_fac.mpz(), n);

    for (size_t i = 0; i <= n; ++i) {
      mpz_fac_ui (i_fac.mpz(), i);
      fc[i] = n_fac / i_fac;
    }

    return CGAL::canonicalize (Poly (fc.begin(), fc.end()));
  }

  const Poly bezier_grid (size_t n, size_t L) {
    std::vector< NT > cp (n+1);
    NT range = L;
    // std::cout << "CP: [";
    for (size_t i = 0; i <= n; ++i) {
      cp[i] = 0;
      mpz_urandomm (cp[i].mpz(), rstate, range.mpz());
      if (sign() != 1) cp[i] = - cp[i];
      // std::cout << " " << cp[i];
    }
    // std::cout << " ]" << std::endl;
    
    Poly ONE = Poly(NT(1));
    Poly t = PT::Shift() (ONE,1,0);
    Poly bi = PT::Shift() (ONE,n,0);
    NT bin = NT(1);

    Poly f = bi * cp[n];
    for (size_t i = 0; i < n; ++i) {
      mpz_bin_uiui (bin.mpz(), n, i);
      bi = PT::Pseudo_division_quotient() (bi, t);
      bi *= (ONE - t);
      f += bi * bin * cp[i];
    }

    return f;
  }

  const Poly_2 bezier_2 (size_t n, size_t L) {
    // CGAL::set_pretty_mode (std::cout);
    Poly x = bezier_grid (n, L);
    // std::cout << x << std::endl;
    Poly y = bezier_grid (n, L);
    // std::cout << y << std::endl;
    Poly_3 xt (x.begin(), x.end());
    xt -= PT_3::Shift() (Poly_3(1),1,0);
    Poly_3 yt (y.begin(), y.end());
    yt -= PT_3::Shift() (Poly_3(1),1,1);
    // std::cout << xt << std::endl;
    // std::cout << yt << std::endl;
    // std::cout << CGAL::canonicalize (CGAL::resultant (xt, yt)) << std::endl;
    return CGAL::canonicalize (CGAL::resultant (xt, yt));
  }

  const Poly bezier_resultant (size_t n, size_t L) {
    const Poly_2 f = bezier_2 (n, L);
    const Poly_2 g = bezier_2 (n, L);
    return CGAL::canonicalize (CGAL::resultant (f, g));
  }

  const Poly unit (size_t n) {
    std::vector< NT > fc (n+1, NT(1));
    return Poly (fc.begin(), fc.end());
  }

  const Poly wilkinson (size_t n) {
    Poly f (NT(1));
    const Poly x = PT::Shift() (Poly (1),1,0);
    for (size_t i = 1; i <= n; ++i)
      f *= x - Poly (NT(i));
    return f;
  }
};

} // namespace CGAL

#endif // CGAL_HAS_GMP

#endif 	    /* !GENERATORS_H_ */
