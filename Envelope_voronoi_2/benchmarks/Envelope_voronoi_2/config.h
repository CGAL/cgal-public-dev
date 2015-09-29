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
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

#if defined(FILTERED_ZONE)
typedef CGAL::Linear_filtered_voronoi_traits<
    CGAL::Point_diagram_traits_2<Kernel> >                  Envelope_voronoi_traits_2;
#else
typedef CGAL::Point_diagram_traits_2<Kernel>               
  Envelope_voronoi_traits_2;
#endif
typedef Kernel::Point_2                                   Point_2;
typedef Kernel::Plane_3                                   Plane_3;
typedef CGAL::Envelope_voronoi_2::Voronoi_diagram_2<
Envelope_voronoi_traits_2>                                Envelope_diagram_2;


// triangulation voronoi
typedef CGAL::Delaunay_triangulation_2<Kernel>            Triangulation;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<
  Triangulation>                                          Adaptation_traits;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<
  Triangulation>                                          Adaptation_policy;
typedef CGAL::Voronoi_diagram_2<Triangulation, 
                                Adaptation_traits, 
                                Adaptation_policy>        Triangulation_voronoi;

// envelope apollonius
typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
#if defined(FILTERED_ZONE)
typedef CGAL::Filtered_voronoi_conic_traits<
  CGAL::Apollonius_diagram_traits_2<Rat_kernel, 
                                  Alg_kernel,
                                  Nt_traits> >        Env_apollonius_traits;
#else
typedef CGAL::Apollonius_diagram_traits_2<Rat_kernel, 
				      Alg_kernel,
				      Nt_traits>        Env_apollonius_traits;
#endif

typedef CGAL::Envelope_voronoi_2::Voronoi_diagram_2<Env_apollonius_traits>
                                                        Apollonius_env_diagram;



// CGAL Apollonius
// Apollonius doesn't work with lazy kernel.
// typedef CGAL::Filtered_kernel<CGAL::Simple_cartesian<CGAL::Lazy_exact_nt<CGAL::Gmpq > > > Apollonius_kernel;
typedef CGAL::Simple_cartesian<CORE::Expr> Apollonius_kernel;
typedef CGAL::Apollonius_graph_traits_2<
  Apollonius_kernel, CGAL::Field_with_sqrt_tag>           Apollonius_traits;
//  Apollonius_kernel>                                    Apollonius_traits;
typedef CGAL::Apollonius_graph_2<Apollonius_traits>     Apollonius_graph;

#ifdef USE_EXACUS

// EXACUS apollonius
typedef NiX::Arithmetic_traits                            AT;
typedef CnX::Integer_descartes_kernel< AT >               IDK;

#if defined(FILTERED_ZONE)
typedef CGAL::Basic_filtered_voronoi_traits< CGAL::Envelope_surface_id_traits_2<
CGAL::EXACUS_apollonius_traits_2<IDK> > >
EXACUS_apollonius_traits;
#else
typedef CGAL::EXACUS_apollonius_traits_2<IDK>             
EXACUS_apollonius_traits;
#endif

typedef CGAL::Envelope_diagram_2<EXACUS_apollonius_traits> 
EXACUS_apollonius_diagram;

// EXACUS using Quadrics
#if defined(FILTERED_ZONE)
typedef CGAL::Basic_filtered_voronoi_traits< CGAL::Envelope_surface_id_traits_2<
QdX::P_quadric_3_envelope_traits< AT > > >     EXACUS_Qdx_traits;
#else
typedef QdX::P_quadric_3_envelope_traits< AT > EXACUS_Qdx_traits;
#endif
typedef CGAL::Envelope_diagram_2< EXACUS_Qdx_traits > 
EXACUS_Qdx_apollonius_diagram;

typedef QdX::P_quadric_3< AT >             P_quadric_3;
#endif // #ifdef USE_EXACUS

// Voronoi on sphere
#if defined(FILTERED_ZONE)
typedef CGAL::Filtered_voronoi_sphere_traits<
    CGAL::Spherical_voronoi_diagram_traits_2<Kernel> >                   Sphere_traits;
#else
typedef CGAL::Spherical_voronoi_diagram_traits_2<Kernel>                  Sphere_traits;
#endif
//typedef CGAL::Envelope_diagram_2<
typedef CGAL::Envelope_voronoi_2::Spherical_voronoi_diagram_2<
              Sphere_traits>                              Envelope_sphere_2;



template <class Site>
class Basic_voronoi
{
public:
  typedef std::vector<Site>        Sites;

  Basic_voronoi() :
    m_filename(0), m_verbose_level(0), m_postscript(false)
    { }
  
  virtual ~Basic_voronoi() {}
  
  int init()
    {
      
      // temporary until I'll now how to initialize a benchmark.
      std::ifstream in(m_filename);
      
      if (in.good() == false)
        return -1;

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
    
  void clean()
    {
      _sites.clear();
    }
  void sync(){}

  void set_file_name(const char * filename, int file_index=0) 
    { 
      m_filename = filename; 
    }
  
  void set_verbose_level(const unsigned int level) { m_verbose_level = level; }
  void set_postscript(const bool postscript) { m_postscript = postscript; }
  
protected:
  const char * m_filename;
  Sites   _sites;

  unsigned int m_verbose_level;
  bool m_postscript;
};


#endif
