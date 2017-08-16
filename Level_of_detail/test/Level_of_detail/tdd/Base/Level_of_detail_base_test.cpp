// Google test includes.
#include "gmock/gmock.h"

// CGAL includes.
#include <CGAL/Base/Level_of_detail_base.h>

using namespace testing;

class LOD_BaseTest: public Test {

public:
	using LodBase = CGAL::Level_of_detail_base;
	LodBase lodbase;

};

TEST_F(LOD_BaseTest, Compiles) {
   
	// Empty test.
}