// Google test includes.
#include "gmock/gmock.h"

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

// CGAL new includes.
#include <CGAL/Level_of_detail_traits_stub.h>
#include <CGAL/Base/Level_of_detail_base.h>

using namespace testing;

class LOD_BaseTest: public Test {

public:
	using FT = double;

	using Kernel    = CGAL::Simple_cartesian<FT>;
	using Point     = Kernel::Point_3;
	using Container = CGAL::Point_set_3<Point>;
	using LodTraits = CGAL::LOD::Level_of_detail_traits_stub<Kernel, Container>;
	using LodBase   = CGAL::LOD::Level_of_detail_base<LodTraits>;

	LodBase lodbase;
};

TEST_F(LOD_BaseTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_BaseTest, RunsCreateLodZeroFunction) {

	std::vector<Point> result;
	lodbase.create_lod_0("internallyUnusedFilePath", std::back_inserter(result));
}