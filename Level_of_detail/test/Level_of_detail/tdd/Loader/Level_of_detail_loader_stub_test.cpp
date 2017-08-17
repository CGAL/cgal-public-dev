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
	using Vector    = Traits::Vector_3;
	using Container = CGAL::Point_set_3<Point>;
	using LodLoader = CGAL::LOD::Level_of_detail_loader_stub<Traits, Container>;

	using Color = CGAL::cpp11::array<unsigned char, 3>;
	using Label = int;
	using Plane = Traits::Plane_3;

	using Color_map = Container:: template Property_map<Color>;
	using Label_map = Container:: template Property_map<Label>;
	using Plane_map = Container:: template Property_map<Plane>;

	using Log = CGAL::LOD::Mylog;

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
	const auto exampleNormal = Vector(-0.16, 0.0, 0.0);

	ASSERT_THAT(input.point(4) , Eq(examplePoint));
	ASSERT_THAT(input.normal(6), Eq(exampleNormal));
}

TEST_F(LOD_LoaderTestStub, SavesLogWithBasicMockData) {

	getBasicData(input);

	Log log;

		log.out << "ply \n"       << 
		"format ascii 1.0 \n"     << 
		"element vertex " << input.number_of_points() << " \n" << 
		"property double x \n"    << 
		"property double y \n"    << 
		"property double z \n"    <<
		"property double nx \n"   <<
		"property double ny \n"   <<
		"property double nz \n"   << 
		"property uchar red \n"   << 
		"property uchar green \n" <<
		"property uchar blue \n"  <<
		"property label \n"       <<
		"property plane \n"       <<
		"end_header \n";

		Color_map color;
		boost::tie(color, boost::tuples::ignore) = input.property_map<Color>("color");

		Label_map label;
		boost::tie(label, boost::tuples::ignore) = input.property_map<Label>("label");

		Plane_map plane;
		boost::tie(plane, boost::tuples::ignore) = input.property_map<Plane>("plane");

	for (Container::const_iterator it = input.begin(); it != input.end(); ++it) {

		log.out 
		<< input.point(*it)    << " " 
		<< input.normal(*it)   << " "
		<< (int) color[*it][0] << " " << (int) color[*it][1] << " " << (int) color[*it][2] << " "
		<< label[*it]          << " "
		<< plane[*it]          << "\n";
	}

	log.save("basic_mock");
}



