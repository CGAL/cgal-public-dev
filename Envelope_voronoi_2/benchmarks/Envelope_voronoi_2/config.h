#ifndef CONFIG_H
#define CONFIG_H


#include <CGAL/Gmpq.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>

#if defined(CGAL_TRAITS_COUNTING)
#include <CGAL/Arr_counting_traits_2.h>
#endif
#if defined(CGAL_TRAITS_TRACING)
#include <CGAL/Arr_tracing_traits_2.h>
#endif

#include <CGAL/Envelope_voronoi_2/Voronoi_diagram_2.h>
#include <CGAL/envelope_3.h>
#include <CGAL/Envelope_voronoi_traits_2/Point_diagram_traits_2.h>
#include <CGAL/Envelope_voronoi_traits_2/Apollonius_diagram_traits_2.h>
#include <CGAL/Envelope_surface_id_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

// CGAL apolloinus
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>

#ifdef USE_EXACUS

// EXACUS apolloinus
#include <CGAL/EXACUS_apollonius_traits_2.h>

// EXACUS apollonius using Qdx
#include <QdX/P_quadric_3.h>
#include <QdX/P_quadric_3_envelope_traits.h>

#endif

// Voronoi on sphere
#include <CGAL/Envelope_voronoi_traits_2/Spherical_voronoi_diagram_traits_2.h>
#include <CGAL/algorithm.h>

// envelope voronoi
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

#if defined(FILTERED_ZONE)
using Envelope_voronoi_traits_2 =
  CGAL::Linear_filtered_voronoi_traits<CGAL::Point_diagram_traits_2<Kernel>>;
#else
using Envelope_voronoi_traits_2 = CGAL::Point_diagram_traits_2<Kernel>;
#endif
using Point_2 = Kernel::Point_2;
using Plane_3 = Kernel::Plane_3;
using Envelope_diagram_2 =
  CGAL::Envelope_voronoi_2::Voronoi_diagram_2<Envelope_voronoi_traits_2>;

// triangulation voronoi
using Triangulation = CGAL::Delaunay_triangulation_2<Kernel>;
using Adaptation_traits =
  CGAL::Delaunay_triangulation_adaptation_traits_2<  Triangulation>;
using Adaptation_policy =
  CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2
  <Triangulation>;
using Triangulation_voronoi =
  CGAL::Voronoi_diagram_2<Triangulation, Adaptation_traits, Adaptation_policy>;

// envelope apollonius
using Nt_traits = CGAL::CORE_algebraic_number_traits;
using Rational = Nt_traits::Rational;
using Algebraic = Nt_traits::Algebraic;
using Rat_kernel = CGAL::Cartesian<Rational>;
using Alg_kernel = CGAL::Cartesian<Algebraic>;
#if defined(FILTERED_ZONE)
using Env_apollonius_traits =
  CGAL::Filtered_voronoi_conic_traits
  <CGAL::Apollonius_diagram_traits_2<Rat_kernel, Alg_kernel, Nt_traits>>;
#else
using Env_apollonius_traits =
  CGAL::Apollonius_diagram_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;
#endif

using Apollonius_env_diagram =
  CGAL::Envelope_voronoi_2::Voronoi_diagram_2<Env_apollonius_traits>;

// CGAL Apollonius
// Apollonius doesn't work with lazy kernel.
// using Apollonius_kernel =
//  CGAL::Filtered_kernel<CGAL::Simple_cartesian<CGAL::Lazy_exact_nt<CGAL::Gmpq >>>;
using Apollonius_kernel = CGAL::Simple_cartesian<CORE::Expr>;
using Apollonius_traits =
  CGAL::Apollonius_graph_traits_2<Apollonius_kernel, CGAL::Field_with_sqrt_tag>;
//  Apollonius_kernel>                                    Apollonius_traits;
using Apollonius_graph = CGAL::Apollonius_graph_2<Apollonius_traits>;

#ifdef USE_EXACUS

// EXACUS apollonius
using AT = NiX::Arithmetic_traits;
using IDK = CnX::Integer_descartes_kernel< AT >;

#if defined(FILTERED_ZONE)
using EXACUS_apollonius_traits =
  CGAL::Basic_filtered_voronoi_traits
  <CGAL::Envelope_surface_id_traits_2<CGAL::EXACUS_apollonius_traits_2<IDK>>>;
#else
using EXACUS_apollonius_traits = CGAL::EXACUS_apollonius_traits_2<IDK>;
#endif

using EXACUS_apollonius_diagram =
  CGAL::Envelope_diagram_2<EXACUS_apollonius_traits>;

// EXACUS using Quadrics
#if defined(FILTERED_ZONE)
using EXACUS_Qdx_traits =
  CGAL::Basic_filtered_voronoi_traits
  <CGAL::Envelope_surface_id_traits_2<QdX::P_quadric_3_envelope_traits<AT>>>;
#else
using EXACUS_Qdx_traits = QdX::P_quadric_3_envelope_traits< AT >;
#endif
using EXACUS_Qdx_apollonius_diagram =
  CGAL::Envelope_diagram_2<EXACUS_Qdx_traits>;

using P_quadric_3 = QdX::P_quadric_3<AT>;
#endif // #ifdef USE_EXACUS

// Voronoi on sphere
#if defined(FILTERED_ZONE)
using Sphere_traits =
  CGAL::Filtered_voronoi_sphere_traits
  <CGAL::Spherical_voronoi_diagram_traits_2<Kernel>>;
#else
using Sphere_traits = CGAL::Spherical_voronoi_diagram_traits_2<Kernel>;
#endif
using Envelope_sphere_2 =
  CGAL::Envelope_voronoi_2::Spherical_voronoi_diagram_2<Sphere_traits>;

template <typename Site>
class Basic_voronoi {
public:
  using Sites = std::vector<Site>;

  Basic_voronoi() : m_filename(0), m_verbose_level(0), m_postscript(false) {}

  virtual ~Basic_voronoi() {}

  int init() {
    // temporary until I'll now how to initialize a benchmark.
    std::ifstream in(m_filename);
    if (in.good() == false) return -1;
    int n_sites;
    in >> n_sites;
    CGAL::copy_n(std::istream_iterator<Site>(in), n_sites,
                 std::back_inserter(_sites));
/*
      int n_sites;
      in >> n_sites;
      for(int i = 0; i < n_sites; ++i)
      {
        Kernel::FT x, y;
        in >> x >> y;

        _sites.push_back(Point_2(x, y));
      }
*/
    if (m_verbose_level > 0)
      std::cout << _sites.size() << " sites" << std::endl;
    return 0;
  }

  void clean() { _sites.clear(); }
  void sync() {}

  void set_file_name(const char * filename, int file_index=0)
  { m_filename = filename; }

  void set_verbose_level(const unsigned int level) { m_verbose_level = level; }
  void set_postscript(const bool postscript) { m_postscript = postscript; }

protected:
  const char * m_filename;
  Sites _sites;

  unsigned int m_verbose_level;
  bool m_postscript;
};

#endif
