// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Local includes.
#include "include/wrappers/Mylod.h"

namespace LOD = CGAL::Level_of_detail;

using Kernel      = CGAL::Exact_predicates_inexact_constructions_kernel;
using LOD_wrapper = LOD::Mylod<Kernel>;

int main(int argc, char** argv) {
   	
    /*
    std::string logs_path;
    if (std::getenv("LOD_LOGS_PATH")) logs_path = static_cast<std::string>(std::getenv("LOD_LOGS_PATH"));
    else logs_path = "/Users/danisimo/Documents/pipeline-clean/logs/"; */

    const std::string logs_path = "/Users/danisimo/Documents/pipeline-clean/logs/";
    LOD_wrapper lod_wrapper(argc, argv, logs_path);
    
    lod_wrapper.run_lod_pipeline();
    return EXIT_SUCCESS;
}
