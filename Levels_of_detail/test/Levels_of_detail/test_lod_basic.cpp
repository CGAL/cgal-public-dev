// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Local includes.
#include "Wrapper.h"

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using Wrapper = CGAL::Levels_of_detail::Wrapper<Kernel>;

/*
Usage:

cmake
 make
./test_lod_basic -data path_to_ply e.g. data.ply -load_params path_to_params e.g. data_params.lod -make_buildings

All the params are in the class Wrapper.h.

Extra options:
-make_buildings - makes LOD01 buildings
-make_trees     - makes LOD012 trees
-make_ground    - makes LOD012 ground
-make_lod01     - makes LOD012 lods
-make_lod2      - makes LOD2 buildings

Do not forget to change the path below to the folder for saving results.
Inside this folder, you need to create 4 folders: buildings/ground/lods/trees.

All intermediate steps are saved in these folders. The final result is saved
in the directory, where these folders are placed.
*/

int main(int argc, char **argv) {

  const std::string path_to_save = "/Users/monet/Documents/gf/lod/logs/";
  Wrapper wrapper(argc, argv, path_to_save);

  wrapper.execute();
  return EXIT_SUCCESS;
}
