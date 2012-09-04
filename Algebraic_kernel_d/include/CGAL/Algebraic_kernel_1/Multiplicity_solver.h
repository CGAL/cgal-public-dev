#ifndef   	MULTIPLICITY_SOLVER_H_
# define   	MULTIPLICITY_SOLVER_H_

#include <iostream>
#include <fstream>
#include <numeric>
#include <algorithm>

#include <boost/random.hpp>

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_1/Cartesian_complex.h>

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

#include <boost/lexical_cast.hpp>

#include <CGAL/Algebraic_kernel_1/Bitsize.h>
#include <CGAL/Algebraic_kernel_1/Bigfloat_traits.h>

namespace CGAL {

template< class AK_ >
class Multiplicity_solver {
public:
  typedef AK_ AK;
  typedef typename AK::Integer    IT;
  typedef typename AK::Bigfloat   RT;
  typedef Cartesian_complex< RT > CT;
  typedef typename Polynomial_type_generator< IT, 1 >::Type IPoly;
  typedef typename Polynomial_type_generator< RT, 1 >::Type RPoly;
  typedef typename Polynomial_type_generator< CT, 1 >::Type CPoly;
  // TODO: Do we really need complex-coefficient polynomials?
  //       Make sure simple complex (op) real operations are used if not
  typedef Polynomial_traits_d< IPoly > IPT;
  typedef Polynomial_traits_d< RPoly > RPT;
  typedef Polynomial_traits_d< CPoly > CPT;

private:
  typedef Bigfloat_traits< RT > BFT;
  typename BFT::Integral_pow2 ipow2;

public:
  typedef std::vector< CT > CVector;
  typedef typename std::vector< CT >::iterator CVectorIterator;
  typedef typename std::vector< CT >::const_iterator CVectorConstIterator;

  struct Cluster {
    CT m;
    RT r;
    int k;
    std::vector <int> is;
  };
  std::vector< Cluster > clusters;
  std::vector< int >     in_cluster;

  std::ofstream asy;

private:
  /* MEMBER VARIABLES */
  IPoly f_int, df_int;
  std::vector< IPoly > ddfs_int;
  CPoly f, df;
  std::vector< CPoly > ddfs;
  int n, L;
  RT nn;
  CVector roots, oldroots;
  std::vector< bool > verification;

  /* MEMBER FUNCTIONS */
  void init () {
    df_int = differentiate (f_int);

    ddfs_int.push_back (f_int);
    ddfs_int.push_back (df_int);
    while (ddfs_int.back().degree() > 0)
      ddfs_int.push_back (differentiate (ddfs_int.back()));

    f = CPoly (f_int.begin(), f_int.end());
    df = CPoly (df_int.begin(), df_int.end());
    for (int i = 0; i < ddfs_int.size(); ++i)
      ddfs.push_back (CPoly (ddfs_int[i].begin(), ddfs_int[i].end()));

    n = f_int.degree();
    L = CGAL::bitsize (f_int);

    nn = square (RT (n));
    in_cluster = std::vector< int > (n, -1);

    asy.open ("clusters.asy", std::ios::out | std::ios::trunc);
    std::ifstream header ("/home/perpeduumimmobile/master/asy/clusters.asy");
    std::string line;
    while (std::getline (header, line))
      asy << line << std::endl;
    header.close();
  }

public:
  ~Multiplicity_solver () {
    asy.close();
  }

private:
  void init_roots () {
    roots = CVector (n, CT::ZERO());
    oldroots = CVector (n, CT::ZERO());
    verification = std::vector< bool > (n, false);

    const CT z = CT (3,4) / RT(5);
    roots[0] = z;
    for (int i = 1; i < n; ++i)
      roots[i] = roots[i-1] * z;
    for (int i = 1; i < n; ++i)
      oldroots[i] = roots[i];
  }

  const bool verify_root (size_t i) {
    /****************************************************************
     * Checks if simplified T_1-Test fails on a disc of radius
     *   r = delta (z_i) / 2 / n.
     * If so, we know that the disc of radius n*r = delta (z_i) / 2
     * contains a root of f.
     ****************************************************************/

    RT RR = (roots[(i+1)%n] - roots[i]).squared_norm_2();
    for (size_t j = (i+2)%n; j != i; j = (j+1)%n)
      RR = min (RR, (roots[j] - roots[i]).squared_norm_2());

    // TODO: make sure that discs around z_i and z_j do not overlap
    RR /= RT (4); // TODO: Shift here
    const RT rr = RR / nn; // TODO: Round down here

    const CT f_at_m = evaluate (f, roots[i]);
    const CT df_at_m = evaluate (df, roots[i]);

    verification[i] = (f_at_m.squared_norm_2() <= df_at_m.squared_norm_2() * rr);
    return verification[i];
  }

  const bool verify_all_roots () {
    bool result = true;
    for (size_t i = 0; i < n; ++i)
      result = result && verify_root (i);
    return result;
  }

  void print_verification () const {
    for (size_t i = 0; i < n; ++i)
      std::cerr << (verification[i] ? '+' : '-');
    std::cerr << std::endl;
  }

  const RT durand_kerner_gauss () {
    RT max_diff = 0;

    for (int i = 0; i < n; ++i) {
      const CT f_at_z_i = f.evaluate (roots[i]);
      CT den = f.lcoeff();
      for (int j = 0; j < n; ++j) {
        if (i == j) continue;
        den *= roots[i] - roots[j];
      }
      const CT difference = f_at_z_i / den;
      roots[i] = roots[i] - difference;

      if (is_zero (max_diff))
        max_diff = difference.norm_inf();
      else
        max_diff = max (max_diff, difference.norm_inf());
    }

    return max_diff;
  }

  const RT durand_kerner_jacobi () {
    // NYETWORKING ?!
    RT max_diff = 0;

    for (int i = 0; i < n; ++i)
      oldroots[i] = roots[i];

    for (int i = 0; i < n; ++i) {
      const CT f_at_z_i = f.evaluate (oldroots[i]);
      CT den = f.lcoeff();
      for (int j = 0; j < n; ++j) {
        if (i == j) continue;
        den *= oldroots[i] - oldroots[j];
      }
      const CT difference = f_at_z_i / den;
      roots[i] = oldroots[i] - difference;

      if (is_zero (max_diff))
        max_diff = difference.norm_inf();
      else
        max_diff = max (max_diff, difference.norm_inf());
    }

    return max_diff;
  }

  const RT durand_kerner_mult_old () {
    RT max_diff = 0;

    for (int i = 0; i < n; ++i) {
      CT f_at_z_i;
      CT den = f.lcoeff();
      if (in_cluster[i] == -1) {
        f_at_z_i = f.evaluate (roots[i]);
        for (int j = 0; j < n; ++j) {
          if (i == j) continue;
          den *= roots[i] - roots[j];
        }
      } else {
        Cluster c = clusters[in_cluster[i]];
        den *= c.k;
        f_at_z_i = f.evaluate (roots[i]);
        //f_at_z_i = ddfs[c.k-1].evaluate (roots[i]);
        for (int j = 0; j < n; ++j) {
          //          if (i == j || in_cluster[i] == in_cluster[j]) continue;
          if (i == j) continue;
          den *= roots[i] - roots[j];
        }
      }
      const CT difference = f_at_z_i / den;
      roots[i] = roots[i] - difference;

      if (is_zero (max_diff))
        max_diff = difference.norm_inf();
      else
        max_diff = max (max_diff, difference.norm_inf());
    }

    return max_diff;
  }

  const RT durand_kerner_mult_2 () {
    for (int i = 0; i < n; ++i) {
      if (in_cluster[i] == -1) {
        const CT f_at_z_i = f.evaluate (roots[i]);
        CT den = f.lcoeff();
        for (int j = 0; j < n; ++j) {
          if (i == j) continue;
          den *= roots[i] - roots[j];
        }

        const CT difference = f_at_z_i / den;
        roots[i] = roots[i] - difference;
      } else {
        const Cluster c = clusters[in_cluster[i]];
        const int k = c.k;

        const CT f_at_z_i = ddfs[k-1].evaluate (roots[i]);
        CT den = f.lcoeff();
        for (int j = 0; j < n; ++j) {
          if (in_cluster[i] == in_cluster[j]) continue;
          den *= roots[i] - roots[j];
        }

        const CT difference = f_at_z_i / den;
        roots[i] = roots[i] - difference;
      }
    }

    for (int i = 0; i < n; ++i) {
      if (in_cluster[i] == -1) continue;
      Cluster c = clusters[in_cluster[i]];
      std::cerr << c.k << std::endl;
      roots[i] = c.m + (roots[i] - c.m) / c.k / n;
    }

    RT max_diff = (roots[0] - oldroots[0]).norm_inf();
    for (int i = 1; i < n; ++i)
      max_diff = max (max_diff, (roots[i] - oldroots[i]).norm_inf());

    for (int i = 0; i < n; ++i)
      oldroots[i] = roots[i];

    return max_diff;
  }

  const RT durand_kerner_mult () {
    for (int i = 0; i < n; ++i) {
      const CT f_at_z_i = f.evaluate (roots[i]);
      CT den = f.lcoeff();
      for (int j = 0; j < n; ++j) {
        if (i == j) continue;
        den *= roots[i] - roots[j];
      }

      const CT difference = f_at_z_i / den;
      roots[i] = roots[i] - difference;
    }

    for (int i = 0; i < n; ++i) {
      if (in_cluster[i] == -1) continue;
      Cluster c = clusters[in_cluster[i]];
      roots[i] = c.m + (roots[i] - c.m) / c.k / n;
    }

    RT max_diff = (roots[0] - oldroots[0]).norm_inf();
    for (int i = 1; i < n; ++i)
      max_diff = max (max_diff, (roots[i] - oldroots[i]).norm_inf());

    for (int i = 0; i < n; ++i)
      oldroots[i] = roots[i];

    return max_diff;
  }

  struct Less_distance {
    int base;
    CVector &roots;
    Less_distance (int base, CVector &roots) : base(base), roots(roots) {}
    bool operator() (int i, int j) {
      return (roots[base] - roots[i]).squared_norm_2() < (roots[base] - roots[j]).squared_norm_2();
    }
  };

  void find_clusters (const RT &r) {
    for (int i = 0; i < n; ++i)
      in_cluster[i] = -1;
    clusters.clear();

    for (int i = 0; i < n; ++i) {
      std::vector< int > sorted;
      for (int j = 0; j < n; ++j) {
        if (in_cluster[j] < 0)
          sorted.push_back (j);
      }
      std::sort (sorted.begin(), sorted.end(), Less_distance (i, roots));

      int j = 1;
      while (j < sorted.size()
             && (roots[i] - roots[sorted[j]]).squared_norm_2() < r)
        ++j;

      if (j > 1) {
        Cluster c;
        c.m = 0;
        for (int k = 0; k < j; ++k) {
          c.is.push_back (sorted[k]);
          c.m += roots[sorted[k]];
          in_cluster[sorted[k]] = clusters.size();
        }
        c.m /= j;
        c.k = j;
        c.r = r;
        clusters.push_back (c);
      }
    }
  }

  void print_clusters () {
    for (int i = 0; i < n; ++i)
      std::cerr << (char)('A' + in_cluster[i]);
    std::cerr << std::endl;
    for (int k = 2; k <= n; ++k) {
      int nr = 0;
      for (int i = 0; i < clusters.size(); ++i)
        if (clusters[i].k == k)
          ++nr;
      if (nr != 0)
        std::cerr << "k=" << k << ": " << nr << "\t";
    }
    std::cerr << std::endl;
  }

  void asy_clusters () {
    for (int i = 0; i < n; ++i)
      asy << "zs[" << i << "] = " << roots[i] << ";" << std::endl;
    for (int i = 0; i < clusters.size(); ++i)
      asy << "ms[" << i << "] = " << clusters[i].m << ";" << std::endl
          << "rs[" << i << "] = " << clusters[i].r << ";" << std::endl;
    asy << "draw_all (0.2);" << std::endl
        << "zs.delete(); ms.delete(); rs.delete();" << std::endl;
  }

public:
  Multiplicity_solver (const IPoly &f_int) : f_int (f_int) { init(); }
  Multiplicity_solver (std::istream &in) : f_int() {
    std::string line;
    while (line.empty() || line[0] == '#')
      std::getline (in, line);
    f_int = boost::lexical_cast< IPoly > (line);

    init();
  }

  void solve () {
    init_roots();

    int count = 0;

    long prec = 53;
    std::cerr << "Try precision " << prec << std::endl;
    leda::bigfloat::set_precision (prec);
    RT eps = leda_bigfloat (1, - static_cast< long >(prec / 2));

    while (true) {
      ++count;
      const RT diff = durand_kerner_mult();
      //      std::cerr << diff << std::endl;

      if (verify_all_roots()) break;

      if (count % 5 == 0)
        {
          find_clusters(ipow2 (min (-5L, - static_cast< long >(count / 2 - 5))));
          //find_clusters (eps * nn);
          //print_clusters();
          asy_clusters();
        }

      if (diff > eps && count < prec) continue;

      print_verification();
      prec = prec * 2 + 1;
      std::cerr << "# iterations: " << count << std::endl;
      std::cerr << "Try precision " << prec << std::endl;
      leda::bigfloat::set_precision (prec);
      eps = leda_bigfloat (1, - static_cast< long >(prec / 2));
    }

    print_verification();
    std::cerr << "Total # iterations: " << count << std::endl;
  }

  const IPoly & polynomial () const { return f_int; }
  const int degree() const { return n; }
  const int bitsize() const { return L; }
  const int number_of_real_roots () const { return 0; }
  const int number_of_complex_roots () const { return n; }

  const CT & operator[] (int n) const { return roots[n]; }
  const CVectorConstIterator complex_roots_begin () const { return roots.begin(); }
  const CVectorConstIterator complex_roots_end ()   const { return roots.end(); }
};

} // namespace CGAL

#endif 	    /* !MULTIPLICITY_SOLVER_H_ */
