// Google test includes.
#include "gmock/gmock.h"

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

// CGAL new includes.
#include <CGAL/Loader/Level_of_detail_loader_stub.h>

using namespace testing;

class LOD_LoaderTestStub: public Test {

public:
	using FT = double;
	
	using Traits    = CGAL::Simple_cartesian<FT>;
	using Point     = Traits::Point_3;      
	using Container = CGAL::Point_set_3<Point>;
	using LodLoader = CGAL::LOD::Level_of_detail_loader_stub<Traits, Container>;

	LodLoader lodLoader;
};

TEST_F(LOD_LoaderTestStub, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_LoaderTestStub, ReturnsNonEmptyContainer) {

	Container container;
	const auto filePath = "internallyUnusedFilePath";

	lodLoader.load(filePath, container); 

	ASSERT_THAT(static_cast<int>(container.number_of_points()), Ne(0));
}