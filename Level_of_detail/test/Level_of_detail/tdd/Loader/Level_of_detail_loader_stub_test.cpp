// Google test includes.
#include "gmock/gmock.h"

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

// New CGAL includes.
#include <CGAL/Loader/Level_of_detail_loader_stub.h>

using namespace testing;

class LOD_LoaderTestStub: public Test {

public:
	using FT = double;
	
	using Traits    = CGAL::Simple_cartesian<FT>;
	using Point     = Traits::Point_3;
	using Vector    = Traits::Vector_3;
	using Container = CGAL::Point_set_3<Point>;
	using LodLoader = CGAL::LOD::Level_of_detail_loader_stub<Traits, Container>;

	LodLoader lodLoader;
	Container input;
};

TEST_F(LOD_LoaderTestStub, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_LoaderTestStub, ReturnsNonEmptyContainer) {

	const auto filePath = "internallyUnusedFilePath";

	lodLoader.get_data(filePath, input); 

	ASSERT_THAT(static_cast<int>(input.number_of_points()), Ne(0));
}

TEST_F(LOD_LoaderTestStub, ReturnsDefaultTestPoint) {

	const auto filePath = "internallyUnusedFilePath";	

	lodLoader.get_data(filePath, input);

	const auto defaultTestPoint = Point(0, 0, 0);

	ASSERT_THAT(input.point(0), Eq(defaultTestPoint));
}

TEST_F(LOD_LoaderTestStub, InsertsOnePoint) {

	const auto newPoint = Point(1, 1, 1);
	lodLoader.insert_point(newPoint, input);

	ASSERT_THAT(input.point(0), Eq(newPoint));
}

TEST_F(LOD_LoaderTestStub, ReturnsBasicMock) {

	const auto filePath = "internallyUnusedFilePath";

	const auto examplePoint  = Point(0, 0, 0);
	const auto exampleNormal = Vector(1, 1, 1);

	lodLoader.set_mock_data_type(CGAL::LOD::Mock_data_type::BASIC);
	lodLoader.get_data(filePath, input);

	ASSERT_THAT(input.point(0) , Eq(examplePoint));
	ASSERT_THAT(input.normal(0), Eq(exampleNormal));
}