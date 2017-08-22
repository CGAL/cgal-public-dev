// Google test includes.
#include "gmock/gmock.h"

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

// New CGAL includes.
#include <CGAL/Regularizer/Level_of_detail_regularizer.h>

using namespace testing;

class LOD_RegularizerTest: public Test {
	
public:
	using FT = double;
	
	using Traits         = CGAL::Simple_cartesian<FT>;
	using Point          = Traits::Point_3;
	using Container      = CGAL::Point_set_3<Point>;
	using LodRegularizer = CGAL::LOD::Level_of_detail_regularizer<Traits, Container>;

	LodRegularizer lodRegularizer;
};

TEST_F(LOD_RegularizerTest, Compiles) {
   
	// Empty test.
}
