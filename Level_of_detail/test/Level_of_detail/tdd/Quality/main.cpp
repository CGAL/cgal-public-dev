// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>

// CGAL new includes.
#include <CGAL/Level_of_detail_traits.h>
#include <CGAL/Tools/Level_of_detail_quality.h>
#include <CGAL/Tools/Level_of_detail_quality_estimator.h>

// using Kernel     = CGAL::Simple_cartesian<double>;
// using Kernel     = CGAL::Exact_predicates_exact_constructions_kernel;

using Kernel                = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT 		            = Kernel::FT;
using Point                 = Kernel::Point_3;
using Container             = CGAL::Point_set_3<Point>;
using Lod_traits            = CGAL::LOD::Level_of_detail_traits<Kernel, Container>;
using Lod_quality           = CGAL::LOD::Level_of_detail_quality<Lod_traits>;
using Lod_quality_estimator = CGAL::LOD::Level_of_detail_quality_estimator<Lod_quality>;

int main(int argc, char** argv) {
   	Lod_quality_estimator lod_quality_estimator(argc, argv);
    
    lod_quality_estimator.run_quality_test();
    return EXIT_SUCCESS;
}