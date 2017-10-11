// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <string>

// CGAL includes.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

// CGAL new includes.
#include <CGAL/Level_of_detail_traits.h>
#include <CGAL/Base/Level_of_detail_base.h>

using namespace testing;

class LOD_BaseTest: public Test {

public:
	using FT = double;

	// using Kernel = CGAL::Simple_cartesian<FT>;
	
	using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
	using Point     = Kernel::Point_3;
	using Container = CGAL::Point_set_3<Point>;
	using LodTraits = CGAL::LOD::Level_of_detail_traits<Kernel, Container>;
	using LodBase   = CGAL::LOD::Level_of_detail_base<LodTraits>;

	LodBase lodBase;
};

TEST_F(LOD_BaseTest, Compiles) {

	// Empty test.
}

/*
TEST_F(LOD_BaseTest, RunsCreateLodsFirstVersionFunction) {
	lodBase.create_lods_a();
} */

TEST_F(LOD_BaseTest, RunsCreateLodsSecondVersionFunction) {
	lodBase.create_lods_b();
}