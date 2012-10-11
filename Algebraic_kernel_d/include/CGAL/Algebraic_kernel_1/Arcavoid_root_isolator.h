#ifndef CGAL_ARCAVOID_ROOT_ISOLATOR_H
#define CGAL_ARCAVOID_ROOT_ISOLATOR_H

#include <CGAL/Algebraic_kernel_1/Arcavoid.h>
#include <boost/detail/algorithm.hpp>

#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel.h>
// For Bitstream-like interface. We offer
//   Square_free_descartes_tag
//   M_k_descartes_tag

namespace CGAL {

class Arcavoid_real_root_isolator_tag {};
class Arcavoid_complex_root_isolator_tag {};


namespace internal {
namespace STL_extension {
// borrowed from Boost 1.51.0: <boost/detail/is_sorted.hpp>
// this is just because older Boost releases don't provide is_sorted.
// TODO : eventually remove this, and replace by C++11 STL extension's is_sorted
//        or use Boost directly
template<class Iterator, class Comp>
inline Iterator is_sorted_until (Iterator first, Iterator last, Comp c) {
  if (first == last)
    return last;

  Iterator it = first; ++it;

  for (; it != last; first = it, ++it)
    if (c(*it, *first))
      return it;

  return it;
}

template<class Iterator>
inline Iterator is_sorted_until (Iterator first, Iterator last) {
  typedef typename boost::detail::iterator_traits<Iterator>::value_type
    value_type;

  typedef std::less<value_type> c;

  return ::boost::detail::is_sorted_until(first, last, c());
}

template<class Iterator, class Comp>
inline bool is_sorted (Iterator first, Iterator last, Comp c) {
  return ::boost::detail::is_sorted_until(first, last, c) == last;
}

template<class Iterator>
inline bool is_sorted (Iterator first, Iterator last) {
  return ::boost::detail::is_sorted_until(first, last) == last;
}
} // namespace STL_extension
} // namespace internal


//typedef internal::Square_free_descartes_tag Square_free_arcavoid_tag;
struct Square_free_arcavoid_tag {};
struct M_arcavoid_tag {};
struct K_arcavoid_tag {};
//typedef internal::M_k_descartes_tag M_k_arcavoid_tag;
struct M_k_arcavoid_tag {};

#define CGAL_ARCAVOID_ROOT_ISOLATOR_COMMON_TYPEDEFS                     \
  public:                                                               \
  typedef BitstreamCoefficientKernel Bitstream_coefficient_kernel;      \
  typedef Arcavoid< Bitstream_coefficient_kernel,                       \
                    Arcavoid_complex_root_isolator_tag >                \
  Complex_root_isolator;                                                \
  typedef Arcavoid< Bitstream_coefficient_kernel,                       \
                    Arcavoid_real_root_isolator_tag >                   \
  Real_root_isolator;                                                   \
protected:                                                              \
 typedef internal::Arcavoid_list< Bitstream_coefficient_kernel > List;  \
public:                                                                 \
 typedef typename List::Coefficient                  Coefficient;       \
 typedef typename List::Arithmetic_kernel            Arithmetic_kernel; \
protected:                                                              \
 typedef typename List::Bigfloat                     Bigfloat;          \
 typedef typename List::Bigfloat_interval            Bigfloat_interval; \
 typedef typename List::Approximation                Approximation;     \
 typedef typename List::Approximation_iterator       Approximation_iterator; \
 typedef typename List::Approximation_const_iterator Approximation_const_iterator; \
 typedef typename List::Cluster                      Cluster;           \
 typedef typename List::Cluster_iterator             Cluster_iterator;  \
 typedef typename List::Cluster_const_iterator       Cluster_const_iterator; \
 typedef typename List::Cluster_range                Cluster_range;     \
public:                                                                 \
 typedef Bigfloat Bound;                                                \
 typedef Bound Boundary;                                                \
protected:                                                              \
 typedef typename CGAL::Polynomial_type_generator< Coefficient, 1 >::Type Input_polynomial; \
 typedef Bigfloat              BF;                                      \
 typedef Bigfloat_interval     BFI;                                     \
 typedef Bigfloat                      RR;                              \
 typedef CGAL::Cartesian_complex< RR > CC;                              \
 typedef typename CGAL::Polynomial_type_generator< RR, 1 >::Type RRPoly; \
 typedef typename CGAL::Polynomial_type_generator< CC, 1 >::Type CCPoly; \
 typedef boost::shared_ptr< List >     List_handle;                     \
 typedef std::vector< Algebraic_complex_1 >                    Algebraic_complex_1_vector; \
 typedef typename Algebraic_complex_1_vector::iterator         Algebraic_complex_1_iterator; \
 typedef typename Algebraic_complex_1_vector::const_iterator   Algebraic_complex_1_const_iterator
// end of CGAL_ARCAVOID_ROOT_ISOLATOR_COMMON_TYPEDEFS

template< class BitstreamCoefficientKernel >
class Arcavoid< BitstreamCoefficientKernel, Arcavoid_complex_root_isolator_tag > {
public:
  class Algebraic_complex_1;
  CGAL_ARCAVOID_ROOT_ISOLATOR_COMMON_TYPEDEFS;

  friend class Arcavoid< Bitstream_coefficient_kernel, Arcavoid_real_root_isolator_tag >;

public:
  class Algebraic_complex_1 {
    friend class Arcavoid< Bitstream_coefficient_kernel, Arcavoid_real_root_isolator_tag >;
    friend class Arcavoid< Bitstream_coefficient_kernel, Arcavoid_complex_root_isolator_tag >;
  public:
    mutable List_handle list_ptr;
    mutable Cluster_iterator cit;

  public:
    const size_t multiplicity () const {
      return cit->multiplicity();
    }
    const CC & center () const {
      return cit->center();
    }
    const RR & radius () const {
      return cit->radius();
    }
    void refine () const {
      const Cluster_range range =  list_ptr->real_subdivide_cluster (cit);
      cit = range.first;
    }

    friend
    std::ostream & operator<< (std::ostream & out, const Algebraic_complex_1 &c) {
      // CGAL_precondition (std::distance (c.cit->begin(), c.cit->end()) == 1
      //                    || (CGAL::is_zero (c.center()) && CGAL::is_zero (c.radius())));
      return out << *(c.cit);
    }

    struct Compare_real
      : public std::binary_function< bool, Algebraic_complex_1, Algebraic_complex_1 > {
      const bool operator() (const Algebraic_complex_1 &lhs,
                             const Algebraic_complex_1 &rhs) const {
        return typename Cluster::Compare_real() (*(lhs.cit), *(rhs.cit));
      }
    };
  };

private:
  List_handle list_ptr;
  Algebraic_complex_1_vector roots;
  int N, mult_zero, n;
  // N: degree
  // n: nr roots
  // mult_zero: multiplicity of zero

public:
  // squarefree constructors
  Arcavoid (const Input_polynomial &input)
    : list_ptr (new List (Bitstream_coefficient_kernel(), input)),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (mult_zero <= 1 ? N : N - mult_zero + 1) {
    solve_to_n_clusters();
  }

  Arcavoid (const Bitstream_coefficient_kernel &bck, const Input_polynomial &input)
    : list_ptr (new List (bck, input)),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (mult_zero <= 1 ? N : N - mult_zero + 1) {
    solve_to_n_clusters();
  }

  template< class InputIterator >
  Arcavoid (const Bitstream_coefficient_kernel &bck,
            const InputIterator first, const InputIterator beyond)
    : list_ptr (new List (bck, first, beyond)),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (mult_zero <= 1 ? N : N - mult_zero + 1) {
    solve_to_n_clusters();
  }

  Arcavoid (Square_free_arcavoid_tag,
            const Input_polynomial &input,
            const Bitstream_coefficient_kernel &bck = Bitstream_coefficient_kernel())
    : list_ptr (new List (bck, input)),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (mult_zero <= 1 ? N : N - mult_zero + 1) {
    solve_to_n_clusters();
  }

  template< class InputIterator >
  Arcavoid (Square_free_arcavoid_tag,
            const InputIterator first, const InputIterator beyond,
            const Bitstream_coefficient_kernel &bck = Bitstream_coefficient_kernel())
    : list_ptr (new List (bck, first, beyond)),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (mult_zero <= 1 ? N : N - mult_zero + 1) {
    solve_to_n_clusters();
  }

  // k-constructors (known number of complex roots)
  Arcavoid (K_arcavoid_tag,
            const Input_polynomial &input,
            int deg_gcd,
            const Bitstream_coefficient_kernel &bck = Bitstream_coefficient_kernel())
    : list_ptr (new List (bck, input)),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (N - deg_gcd) {
    solve_to_n_clusters();
  }

  template< class InputIterator >
  Arcavoid (K_arcavoid_tag,
            const InputIterator first, const InputIterator beyond,
            int deg_gcd,
            const Bitstream_coefficient_kernel &bck = Bitstream_coefficient_kernel())
    : list_ptr (new List (bck, first, beyond)),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (N - deg_gcd) {
    solve_to_n_clusters();
  }

private:
  void solve_to_n_clusters() {
    // BFS subdivision
    int nr_clusters = std::distance (list_ptr->begin(), list_ptr->end());

    while (nr_clusters < n) {
      Cluster_iterator it = list_ptr->begin();
      while (it != list_ptr->end() && nr_clusters < n) {
        if (it->multiplicity() > 1 && ! CGAL::is_zero (it->radius())) {
          Cluster_range range = list_ptr->subdivide_cluster (it);
          it = range.second;

          nr_clusters += std::distance (range.first, range.second) - 1;
          CGAL_assertion (nr_clusters == std::distance (list_ptr->begin(), list_ptr->end()));
        } else
          ++it;
      }
    }

    for (Cluster_iterator it = list_ptr->begin(); it != list_ptr->end(); ++it) {
      Algebraic_complex_1 c;
      c.list_ptr = list_ptr;
      c.cit = it;
      roots.push_back (c);
    }
  }

public:
  const Algebraic_complex_1_iterator complex_roots_begin () {
    return roots.begin();
  }
  const Algebraic_complex_1_iterator complex_roots_end () {
    return roots.end();
  }
  const Algebraic_complex_1_const_iterator complex_roots_begin () const {
    return roots.begin();
  }
  const Algebraic_complex_1_const_iterator complex_roots_end () const {
    return roots.end();
  }
  const size_t number_of_complex_roots () const {
    return roots.size();
  }
  const Algebraic_complex_1 & operator[] (int i) const {
    return roots[i];
  }
};

template< class BitstreamCoefficientKernel >
class Arcavoid< BitstreamCoefficientKernel, Arcavoid_real_root_isolator_tag > {
public:
  typedef typename Arcavoid< BitstreamCoefficientKernel, Arcavoid_complex_root_isolator_tag >
  ::Algebraic_complex_1 Algebraic_complex_1;
  CGAL_ARCAVOID_ROOT_ISOLATOR_COMMON_TYPEDEFS;

  // TODO: polynomial() function? What is the return type? Univariate? Bivariate + alpha?

private:
  List_handle list_ptr;
  Algebraic_complex_1_vector real_roots;
  enum { SQUARE_FREE, M_ARCAVOID, K_ARCAVOID, M_K_ARCAVOID } variant;
  int N, mult_zero, n, n_real, nr_multiple_real_root_clusters;

public:
  // Squarefree constructors
  Arcavoid (const Input_polynomial &input)
    : list_ptr (new List (Bitstream_coefficient_kernel(), input)),
      variant (SQUARE_FREE),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (mult_zero <= 1 ? N : N - mult_zero + 1) {
    solve_squarefree();
  }

  Arcavoid (const Bitstream_coefficient_kernel &bck, const Input_polynomial &input)
    : list_ptr (new List (bck, input)),
      variant (SQUARE_FREE),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (mult_zero <= 1 ? N : N - mult_zero + 1) {
    solve_squarefree();
  }

  template< class InputIterator >
  Arcavoid (const Bitstream_coefficient_kernel &bck,
            const InputIterator first, const InputIterator beyond)
    : list_ptr (new List (bck, first, beyond)),
      variant (SQUARE_FREE),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (mult_zero <= 1 ? N : N - mult_zero + 1) {
    solve_squarefree();
  }

  Arcavoid (Square_free_arcavoid_tag,
            const Input_polynomial &input,
            const Bitstream_coefficient_kernel &bck = Bitstream_coefficient_kernel())
    : list_ptr (new List (bck, input)),
      variant (SQUARE_FREE),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (mult_zero <= 1 ? N : N - mult_zero + 1) {
    solve_squarefree();
  }

  template< class InputIterator >
  Arcavoid (Square_free_arcavoid_tag,
            const InputIterator first, const InputIterator beyond,
            const Bitstream_coefficient_kernel &bck = Bitstream_coefficient_kernel())
    : list_ptr (new List (bck, first, beyond)),
      variant (SQUARE_FREE),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (mult_zero <= 1 ? N : N - mult_zero + 1) {
    solve_squarefree();
  }

  // k-constructors (known number of complex roots)
  Arcavoid (K_arcavoid_tag,
            const Input_polynomial &input,
            int deg_gcd,
            const Bitstream_coefficient_kernel &bck = Bitstream_coefficient_kernel())
    : list_ptr (new List (bck, input)),
      variant (K_ARCAVOID),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (N - deg_gcd) {
    solve_to_n_complex_clusters();
  }

  template< class InputIterator >
  Arcavoid (K_arcavoid_tag,
            const InputIterator first, const InputIterator beyond,
            int deg_gcd,
            const Bitstream_coefficient_kernel &bck = Bitstream_coefficient_kernel())
    : list_ptr (new List (bck, first, beyond)),
      variant (K_ARCAVOID),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (N - deg_gcd) {
    solve_to_n_complex_clusters();
  }

  // m-constructors (known number of real roots)
  Arcavoid (M_arcavoid_tag,
            const Input_polynomial &input,
            int nr_real_roots,
            const Bitstream_coefficient_kernel &bck = Bitstream_coefficient_kernel())
    : list_ptr (new List (bck, input)),
      variant (M_ARCAVOID),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (mult_zero <= 1 ? N : N - mult_zero + 1),
      n_real (nr_real_roots) {
    solve_to_m_real_clusters();
  }

  template< class InputIterator >
  Arcavoid (M_arcavoid_tag,
            const InputIterator first, const InputIterator beyond,
            int nr_real_roots,
            const Bitstream_coefficient_kernel &bck = Bitstream_coefficient_kernel())
    : list_ptr (new List (bck, first, beyond)),
      variant (M_ARCAVOID),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (mult_zero <= 1 ? N : N - mult_zero + 1),
      n_real (nr_real_roots) {
    solve_to_m_real_clusters();
  }

  // m-k-constructors (known number of real and complex roots)
  Arcavoid (M_k_arcavoid_tag,
            const Input_polynomial &input,
            int nr_real_roots, int deg_gcd,
            const Bitstream_coefficient_kernel &bck = Bitstream_coefficient_kernel())
    : list_ptr (new List (bck, input)),
      variant (M_K_ARCAVOID),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (N - deg_gcd),
      n_real (nr_real_roots) {
    solve_to_m_real_clusters();
  }

  template< class InputIterator >
  Arcavoid (M_k_arcavoid_tag,
            const InputIterator first, const InputIterator beyond,
            int nr_real_roots, int deg_gcd,
            const Bitstream_coefficient_kernel &bck = Bitstream_coefficient_kernel())
    : list_ptr (new List (bck, first, beyond)),
      variant (M_K_ARCAVOID),
      N (list_ptr->degree()),
      mult_zero (list_ptr->multiplicity_of_zero()),
      n (N - deg_gcd),
      n_real (nr_real_roots) {
    solve_to_m_real_clusters();
  }

private:
  void solve_squarefree () {
    if (list_ptr->begin() == list_ptr->end())
      return;

    list_ptr->subdivide_cluster (list_ptr->begin());
    // Precondition for termination:
    // f has only simple roots on the real axis (but, possibly, zero)

    // BFS subdivision
    bool all_real_have_mult_1 = false;
    while (! all_real_have_mult_1) {
      all_real_have_mult_1 = true;

      Cluster_iterator it = list_ptr->begin();
      do {
        if (it->multiplicity() > 1
            && ! CGAL::is_zero (it->radius())
            && it->touch_real()) {
          all_real_have_mult_1 = false;
          it = list_ptr->subdivide_cluster (it).second;
        } else
          ++it;
      } while (it != list_ptr->end());
    }

    const Cluster_iterator beyond_reals = list_ptr->partition_and_sort_reals();

    for (Cluster_iterator it = list_ptr->begin(); it != beyond_reals; ++it) {
      CGAL_assertion (it->touch_real());

      Algebraic_complex_1 c;
      c.list_ptr = list_ptr;
      c.cit = it;
      real_roots.push_back (c);
    }

    CGAL_postcondition (CGAL::internal::STL_extension::is_sorted (real_roots.begin(), real_roots.end(),
                                                                  typename Algebraic_complex_1::Compare_real()));
  }

  void solve_to_n_complex_clusters() {
    if (list_ptr->begin() == list_ptr->end())
      return;

    list_ptr->subdivide_cluster (list_ptr->begin());
    // BFS subdivision
    int nr_clusters = std::distance (list_ptr->begin(), list_ptr->end());

    while (nr_clusters < n) {
      bool all_real_have_mult_1 = true;

      Cluster_iterator it = list_ptr->begin();
      while (it != list_ptr->end() && nr_clusters < n) {
        if (it->multiplicity() > 1 && ! CGAL::is_zero (it->radius())) {
          if (it->touch_real())
            all_real_have_mult_1 = false;

          Cluster_range range = list_ptr->subdivide_cluster (it);
          it = range.second;

          nr_clusters += std::distance (range.first, range.second) - 1;
          CGAL_assertion (nr_clusters == std::distance (list_ptr->begin(), list_ptr->end()));
        } else
          ++it;
      }

      if (all_real_have_mult_1)
        break;
    }

    const Cluster_iterator beyond_reals = list_ptr->partition_and_sort_reals();

    for (Cluster_iterator it = list_ptr->begin(); it != beyond_reals; ++it) {
      CGAL_assertion (it->touch_real());

      Algebraic_complex_1 c;
      c.list_ptr = list_ptr;
      c.cit = it;
      real_roots.push_back (c);
    }

    CGAL_postcondition (CGAL::internal::STL_extension::is_sorted (real_roots.begin(), real_roots.end(),
                                                                  typename Algebraic_complex_1::Compare_real()));
  }

  void solve_to_m_real_clusters() {
    if (list_ptr->begin() == list_ptr->end())
      return;

    list_ptr->subdivide_cluster (list_ptr->begin());
    // BFS subdivision
    Cluster_iterator beyond_reals = list_ptr->partition_and_sort_reals();
    int nr_real_clusters = std::distance (list_ptr->begin(), beyond_reals);
    int nr_clusters = std::distance (list_ptr->begin(), list_ptr->end());

    //DBG_ARCA (std::cerr << "M_k_arcavoid. GOAL: #real (m) = " << n_real << " #total (N-k) = " << n << std::endl);

    while (nr_real_clusters < n_real
           && nr_clusters < n) {

      //DBG_ARCA (std::cerr << "now: #real = " << nr_real_clusters << " #total = " << nr_clusters << std::endl);

      bool all_real_have_mult_1 = true;

      Cluster_iterator it = list_ptr->begin();
      while (it != beyond_reals
             && nr_real_clusters < n_real
             && nr_clusters < n) {
        CGAL_assertion (it->touch_real());

        if (it->multiplicity() > 1
            && ! CGAL::is_zero (it->radius())) {
          all_real_have_mult_1 = false;

          Cluster_range range = list_ptr->subdivide_cluster (it);
          it = range.second;

          nr_clusters += std::distance (range.first, range.second) - 1;

          --nr_real_clusters;
          for (Cluster_iterator jt = range.first; jt != range.second; ++jt)
            if (jt->touch_real())
              ++nr_real_clusters;

          CGAL_assertion (nr_clusters == std::distance (list_ptr->begin(), list_ptr->end()));
        } else
          ++it;
      }

      beyond_reals = list_ptr->partition_and_sort_reals();
      CGAL_assertion (nr_real_clusters == std::distance (list_ptr->begin(), beyond_reals));

      if (all_real_have_mult_1)
        break;
    }

    nr_multiple_real_root_clusters = 0;
    for (Cluster_iterator it = list_ptr->begin(); it != beyond_reals; ++it) {
      CGAL_assertion (it->touch_real());

      if (it->multiplicity() > 1) ++nr_multiple_real_root_clusters;

      Algebraic_complex_1 c;
      c.list_ptr = list_ptr;
      c.cit = it;
      real_roots.push_back (c);
    }

    CGAL_postcondition (CGAL::internal::STL_extension::is_sorted (real_roots.begin(), real_roots.end(),
                                                                  typename Algebraic_complex_1::Compare_real()));
  }

public:
  // Bitstream_descartes facade
  void refine_interval(int i) const { real_roots[i].refine(); }

  const bool is_certainly_simple_root (int i) const {
    switch (variant) {
    case SQUARE_FREE:
      CGAL_assertion (real_roots[i].multiplicity() == 1);
      return true;
    case M_ARCAVOID:
    case K_ARCAVOID:
    case M_K_ARCAVOID:
      return real_roots[i].multiplicity() == 1;
    default:
      CGAL_assertion (variant == SQUARE_FREE
                      || variant == M_ARCAVOID
                      || variant == K_ARCAVOID
                      || variant == M_K_ARCAVOID);
    }
    return false;
  }
  const bool is_certainly_multiple_root (int i) const {
    switch (variant) {
    case SQUARE_FREE:
      CGAL_assertion (real_roots[i].multiplicity() == 1);
      return false;
    case M_ARCAVOID:
    case K_ARCAVOID:
    case M_K_ARCAVOID:
      return (nr_multiple_real_root_clusters == 1
              && real_roots[i].multiplicity() > 1);
    default:
      CGAL_assertion (variant == SQUARE_FREE
                      || variant == M_ARCAVOID
                      || variant == K_ARCAVOID
                      || variant == M_K_ARCAVOID);
    }
    return false;
  }
  const int get_upper_bound_for_multiplicity (int i) const {
    return real_roots[i].multiplicity();
  }

  const Input_polynomial & polynomial () const {
    return list_ptr->polynomial();
  }
  const Algebraic_complex_1_iterator real_roots_begin () {
    return real_roots.begin();
  }
  const Algebraic_complex_1_iterator real_roots_end () {
    return real_roots.end();
  }
  const Algebraic_complex_1_const_iterator real_roots_begin () const {
    return real_roots.begin();
  }
  const Algebraic_complex_1_const_iterator real_roots_end () const {
    return real_roots.end();
  }
  const size_t number_of_real_roots () const {
    return real_roots.size();
  }
  const RR left_bound (int i) const {
    return real_roots[i].cit->center().real() - real_roots[i].cit->radius();
  }
  const RR right_bound (int i) const {
    return real_roots[i].cit->center().real() + real_roots[i].cit->radius();
  }
  const bool is_exact_root (int i) const {
    return CGAL::is_zero (real_roots[i].cit->radius());
  }
  const Algebraic_complex_1 & operator[] (int i) const {
    return real_roots[i];
  }
};

} // namespace CGAL

#endif // CGAL_ARCAVOID_ROOT_ISOLATOR_H
