/*! \file   mobius_diagram.cpp
 * \brief  An example file for the Mobius diagram of weighted points.
 */

#include <iostream>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/envelope_voronoi_2.h>
#include <CGAL/Envelope_voronoi_2/Voronoi_diagram_2.h>
#include <CGAL/Envelope_voronoi_traits_2/Moebius_diagram_traits_2.h>

#include "read_objects.h"

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using VD_traits = CGAL::Moebius_diagram_traits_2<Kernel>;
using Mobius_diagram = CGAL::Envelope_voronoi_2::Voronoi_diagram_2<VD_traits>;
using Point_2 = Kernel::Point_2;
using Site_2 = VD_traits::Site_2;

int main(int argc, char* argv[]) {
  const char* filename = (argc > 1) ? argv[1] : "mobius_sites.txt";
  std::list<Site_2> sites;
  read_objects<Site_2>(filename, std::back_inserter(sites));
  std::cout << "Number of sites: " << sites.size() << std::endl;
  Mobius_diagram diagram;
  CGAL::voronoi_2(sites.begin(), sites.end(), diagram);
  std::cout << "Mobius diagram:" << std::endl
            << "V = " << diagram.number_of_vertices()
            << ", E = " << diagram.number_of_edges()
            << ", F = " << diagram.number_of_faces() << std::endl;
  return 0;
}
