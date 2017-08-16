// Google test includes.
#include "gmock/gmock.h"

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

// CGAL new includes.
#include <CGAL/Loader/Level_of_detail_loader.h>

using namespace testing;

class LOD_LoaderTest: public Test {

public:
	using FT = double;
	
	using Traits    = CGAL::Simple_cartesian<FT>;
	using Point     = Traits::Point_3;
	using Container = CGAL::Point_set_3<Point>;
	using LodLoader = CGAL::LOD::Level_of_detail_loader<Traits, Container>;

	LodLoader lodloader;
};

TEST_F(LOD_LoaderTest, Compiles) {
   
	// Empty test.
}