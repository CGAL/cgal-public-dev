// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <string>

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/array.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Loader/Level_of_detail_loader_stub.h>

using namespace testing;

class LOD_LoaderTestStub: public Test {

public:
	using FT = double;
	
	using Traits    = CGAL::Simple_cartesian<FT>;
	using Point     = Traits::Point_3;
	using Normal    = Traits::Vector_3;
	using Container = CGAL::Point_set_3<Point>;
	using LodLoader = CGAL::LOD::Level_of_detail_loader_stub<Traits, Container>;

	using Log  = CGAL::LOD::Mylog;

	LodLoader lodLoader;
	Container input;

	const std::string internallyUnusedFilePath = "internallyUnusedFilePath";

	void getDefaultData(Container &input) {
		const auto filePath = internallyUnusedFilePath;	

		lodLoader.set_mock_data_type(CGAL::LOD::Mock_data_type::DEFAULT);
		lodLoader.get_data(filePath, input);
	}

	void getBasicData(Container &input) {
		const auto filePath = internallyUnusedFilePath;

		lodLoader.set_mock_data_type(CGAL::LOD::Mock_data_type::BASIC);
		lodLoader.get_data(filePath, input);
	}
};

TEST_F(LOD_LoaderTestStub, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_LoaderTestStub, ReturnsNonEmptyContainer) {

	getDefaultData(input);

	ASSERT_THAT(static_cast<int>(input.number_of_points()), Ne(0));
}

TEST_F(LOD_LoaderTestStub, ReturnsDefaultTestPoint) {

	getDefaultData(input);
	const auto defaultTestPoint = Point(0, 0, 0);

	ASSERT_THAT(input.point(0), Eq(defaultTestPoint));
}

TEST_F(LOD_LoaderTestStub, InsertsOnePoint) {

	const auto newPoint = Point(1, 1, 1);
	lodLoader.insert_point(newPoint, input);

	ASSERT_THAT(input.point(0), Eq(newPoint));
}

TEST_F(LOD_LoaderTestStub, ReturnsBasicMock) {

	getBasicData(input);	

	const auto examplePoint  = Point(0.7, 1.0, 0.8);
	const auto exampleNormal = Normal(-0.16, 0.0, 0.0);

	ASSERT_THAT(input.point(4) , Eq(examplePoint));
	ASSERT_THAT(input.normal(6), Eq(exampleNormal));
}

TEST_F(LOD_LoaderTestStub, SavesLogWithBasicMockData) {

	getBasicData(input);

	Log log; const auto saveWithExtraProperties = true;

	log.save_ply<Traits, Container>(input, "basic_mock", saveWithExtraProperties);
}
