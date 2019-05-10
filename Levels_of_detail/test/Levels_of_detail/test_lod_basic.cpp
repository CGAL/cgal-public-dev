// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Local includes.
#include "Wrapper.h"

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Wrapper = CGAL::Levels_of_detail::Wrapper<Kernel>;

int main(int argc, char **argv) {

  const std::string path_to_save = "/Users/monet/Documents/lod/logs/";
  Wrapper wrapper(argc, argv, path_to_save);

  wrapper.execute();
  return EXIT_SUCCESS;
}
