#ifndef BENCH_CLASSES_H
#define BENCH_CLASSES_H

#include "config.h"
#include "Option_parser.hpp"

class Bench_envelope_voronoi : public Basic_voronoi<Point_2> {
public:
  void op();
};

// cgal triangulation
class Bench_triangulation_voronoi : public Basic_voronoi<Point_2> {
public:
  void op();
};

// apollonius using CGAL conic class
class Bench_envelope_apollonius :
  public Basic_voronoi<Env_apollonius_traits::Site_2> {
public:
  void op();
};

// CGAL apollonius benchmark
class Bench_apollonius : public Basic_voronoi<Apollonius_traits::Site_2> {
public:
  void op();
};

#ifdef USE_EXACUS

// EXACUS apollonius
class Bench_EXACUS_apollonius :
  public Basic_voronoi<EXACUS_apollonius_traits::Surface_3> {
public:
  void op();
};

// EXACUS apollonius using Qdx cones
class Bench_EXACUS_Qdx_apollonius :
  public Basic_voronoi<EXACUS_apollonius_traits::Surface_3> {
  using Base = Basic_voronoi<EXACUS_apollonius_traits::Surface_3>;
  using Base_sites = Base::Sites;
  using Rational = AT::Rational;
  using Quadrics = std::vector< P_quadric_3 >;

public:
  Bench_EXACUS_Qdx_apollonius()
  {
    EXACUS_Qdx_traits::toggle_to_apollonius_diagram_mode();
  }

  int init();

  void op();

protected:
  void convert_sites_to_quadrics();

  Quadrics _quadric_sites;
};

#endif // #ifdef USE_EXACUS

// bench voronoi diagram on sphere
class Bench_sphere : public Basic_voronoi<Sphere_traits::Site_2> {
public:
  void op();
};


#endif
