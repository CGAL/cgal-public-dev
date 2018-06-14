// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Local includes.
#include "../include/wrappers/Myquality.h"

namespace LOD = CGAL::Level_of_detail;

/*
using Kernel      = CGAL::Exact_predicates_inexact_constructions_kernel;
using LOD_quality = LOD::Myquality<Kernel>; */

int main(int /* argc */, char** /* argv */) {
   	
    /*
    std::string logs_path;
    if (std::getenv("LOD_LOGS_PATH")) logs_path = static_cast<std::string>(std::getenv("LOD_LOGS_PATH"));
    else logs_path = "/Users/danisimo/Documents/pipeline-clean/logs/"; */

    /*
    const std::string logs_path = "/Users/danisimo/Documents/pipeline-clean/logs/";
    LOD_quality lod_quality(argc, argv, logs_path); lod_quality.run(); */

    return EXIT_SUCCESS;
}