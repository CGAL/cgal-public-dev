/*! \file   voronoi_of_points_on_sphere.cpp
 *\brief  An example for a Voronoi diagram of points embedded on the sphere.
 */

#include <fstream>
#include <iostream>
#include <vector>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/envelope_voronoi_2.h>
#include <CGAL/envelope_voronoi_2.h>
#include <CGAL/Envelope_voronoi_traits_2/Spherical_voronoi_diagram_traits_2.h>

#include <read_objects.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Voronoi_traits = CGAL::Spherical_voronoi_diagram_traits_2<Kernel>;
using Voronoi_diagram =
  CGAL::Envelope_voronoi_2::Spherical_voronoi_diagram_2<Voronoi_traits>;
using Site_2 = Voronoi_traits::Site_2;

int main(int argc, char* argv[]) {
  const char* filename = (argc > 1) ? argv[1] : "data/voronoi_of_points_on_sphere/voronoi.in";
  std::vector<Kernel::Point_3> points;
  read_objects<Kernel::Point_3>(filename, std::back_inserter(points));
  std::vector<Site_2> sites;
  sites.reserve(points.size());
  Voronoi_traits traits;
  auto ctr_pnt = traits.construct_point_2_object();
  for (const auto& p : points) {
    Kernel::Vector_3 v(p, CGAL::ORIGIN);
    sites.push_back(ctr_pnt(v.direction()));
  }
  std::cout << "Number of sites: " << sites.size() << std::endl;
  for (const auto& s : sites) std::cout << s << std::endl;
  Voronoi_diagram diagram;
  CGAL::voronoi_2(sites.begin(), sites.end(), diagram);
  std::cout << "Voronoi diagram:" << std::endl
            << "V = " << diagram.number_of_vertices() << ", E = "
            << diagram.number_of_edges() << ", F = "
            << diagram.number_of_faces() <<
    std::endl;
  return 0;
}
