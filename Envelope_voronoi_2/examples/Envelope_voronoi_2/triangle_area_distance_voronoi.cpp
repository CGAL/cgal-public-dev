/*! \file triangle_area_distance_voronoi.cpp
 * \brief Example of the 2-point triangle area distance function.
 */

#include <iostream>
#include <list>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/envelope_voronoi_2.h>
#include <CGAL/Envelope_voronoi_2/Voronoi_diagram_2.h>
#include <CGAL/Envelope_voronoi_traits_2/Triangle_area_distance_traits_2.h>

#include "cartesian_product.h"
#include "read_objects.h"

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using VD_traits = CGAL::Triangle_area_distance_traits_2<Kernel>;
using Triangle_VD = CGAL::Envelope_voronoi_2::Voronoi_diagram_2<VD_traits>;
using Point_2 = Kernel::Point_2;
using Site_2 = VD_traits::Site_2;

int main(int argc, char* argv[]) {
  const char* filename = (argc > 1) ? argv[1] : "tri_dist_sites.txt";
  std::list<Point_2> sites;
  read_objects<Point_2>(filename, std::back_inserter(sites));
  std::list<Site_2> site_pairs;
  cartesian_square_without_order(sites.begin(), sites.end(),
                                 std::back_inserter(site_pairs));
  std::cout << "Number of sites: " << site_pairs.size() << std::endl;
  Triangle_VD diagram;
  CGAL::voronoi_2(site_pairs.begin(), site_pairs.end(), diagram);
  std::cout << "Triangle area VD:" << std::endl
            << "V = " << diagram.number_of_vertices()
            << ", E = " << diagram.number_of_edges()
            << ", F = " << diagram.number_of_faces() << std::endl;
  return 0;
}
