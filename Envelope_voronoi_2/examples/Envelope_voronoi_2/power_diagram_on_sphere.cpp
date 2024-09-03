//! \file examples/Envelope_3/ex_envelope_planes.cpp
// Constructing the lower and the upper envelope of a set of planes.

#include <iostream>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/envelope_voronoi_2.h>
#include <CGAL/Envelope_voronoi_traits_2/Spherical_power_diagram_traits_2.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Plane_3 = Kernel::Plane_3;
using Traits_2 = CGAL::Spherical_power_diagram_traits_2<Kernel>;
using Voronoi_diagram_2 =
  CGAL::Envelope_voronoi_2::Spherical_voronoi_diagram_2<Traits_2>;

int main(int argc, char* argv[]) {
  std::vector<Traits_2::Site_2> sites;
  sites.push_back(Point_3());
  sites.push_back(Plane_3());
  Voronoi_diagram_2 VD;
  voronoi_2(sites.begin(), sites.end(), VD);
  return 0;
}
