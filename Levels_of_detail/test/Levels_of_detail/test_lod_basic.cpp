// CGAL includes.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

// Local includes.
#include "include/wrappers/Mylod.h"

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using LOD = CGAL::Levels_of_detail::Mylod<Kernel>;

int main(int argc, char **argv) {

  const std::string logs_path = "/Users/danisimo/Documents/lod/logs/";
  LOD lod(argc, argv, logs_path);

  lod.run();
  return EXIT_SUCCESS;
}
