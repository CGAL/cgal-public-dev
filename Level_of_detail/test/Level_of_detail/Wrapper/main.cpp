// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>

// CGAL new includes.
#include <CGAL/Level_of_detail/Level_of_detail_traits.h>

// Local includes.
#include "../include/Mywrapper.h"

// using Kernel = CGAL::Simple_cartesian<double>;
// using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using Kernel      = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT 		  = Kernel::FT;
using Point       = Kernel::Point_3;
using Container   = CGAL::Point_set_3<Point>;
using Lod_traits  = CGAL::Level_of_detail::Level_of_detail_traits<Kernel, Container>;
using Lod_wrapper = CGAL::Level_of_detail::Mywrapper<Lod_traits>;

int main(int argc, char** argv) {
   	Lod_wrapper lod_wrapper(argc, argv);
    
    lod_wrapper.run_lod_pipeline();
    return EXIT_SUCCESS;
}