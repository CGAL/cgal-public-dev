// Google test includes.
#include "gmock/gmock.h"

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

// New CGAL includes.
#include <CGAL/Selector/Level_of_detail_selector.h>

using namespace testing;

class LOD_SelectorTest: public Test {

public:
	using FT = double;
	
	using Traits    = CGAL::Simple_cartesian<FT>;
	using Point     = Traits::Point_3;
	using Container = CGAL::Point_set_3<Point>;
	using LodSelector = CGAL::LOD::Level_of_detail_selector<Traits>;

	LodSelector lodSelector;
};

TEST_F(LOD_SelectorTest, Compiles) {
   
	// Empty test.
}