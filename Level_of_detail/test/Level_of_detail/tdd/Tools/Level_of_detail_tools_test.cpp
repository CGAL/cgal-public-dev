// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <string>

// CGAL includes.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>

// CGAL new includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Tools/Level_of_detail_tools.h>

using namespace testing;

class LOD_ToolsTest: public Test {

public:
	using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
	using Point     = Kernel::Point_3;
	using Container = CGAL::Point_set_3<Point>;
	
	using Log = CGAL::LOD::Mylog;
};

TEST_F(LOD_ToolsTest, Compiles) {

	// Empty test.
}