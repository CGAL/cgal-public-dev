// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <map>
#include <vector>
#include <string>

// CGAL includes.
#include <CGAL/Simple_cartesian.h>

// CGAL new includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Clutter/Level_of_detail_clutter_processor.h>

using namespace testing;

class LOD_ClutterProcessorTest: public Test {

public:
	using FT = double;

	using Kernel  = CGAL::Simple_cartesian<FT>;
	using Point_2 = Kernel::Point_2;

	using Boundary_data    = std::map<int, std::vector<int> >;
	using Projected_points = std::map<int, Point_2>;

	using LodClutterProcessor = CGAL::LOD::Level_of_detail_clutter_processor<Kernel, Boundary_data, Projected_points>;
	using Log = CGAL::LOD::Mylog;

	LodClutterProcessor lodClutterProcessor;

	Boundary_data boundary_clutter; Projected_points boundary_clutter_projected;

	LOD_ClutterProcessorTest() {
		create_data();
	}

	void create_data() {

		boundary_clutter.clear();
		boundary_clutter_projected.clear();

		set_basic_input();
	}

	void set_basic_input(){

		std::vector<int> idxs(24);

		boundary_clutter_projected[0]  = Point_2( 0.05,  0.04); idxs[0]  =  0;
		boundary_clutter_projected[1]  = Point_2( 0.11, -0.04); idxs[1]  =  1;
		boundary_clutter_projected[2]  = Point_2( 0.25,  0.05); idxs[2]  =  2;
		boundary_clutter_projected[3]  = Point_2( 0.28, -0.06); idxs[3]  =  3;
		boundary_clutter_projected[4]  = Point_2( 0.39,  0.05); idxs[4]  =  4;
		boundary_clutter_projected[5]  = Point_2( 0.49, -0.06); idxs[5]  =  5;
		boundary_clutter_projected[6]  = Point_2(-0.06, -0.06); idxs[6]  =  6;
		boundary_clutter_projected[7]  = Point_2(-0.08,  0.04); idxs[7]  =  7;
		boundary_clutter_projected[8]  = Point_2(-0.14, -0.03); idxs[8]  =  8;
		boundary_clutter_projected[9]  = Point_2( 0.03,  0.16); idxs[9]  =  9;
		boundary_clutter_projected[10] = Point_2(-0.05,  0.28); idxs[10] = 10;
		boundary_clutter_projected[11] = Point_2( 0.02,  0.44); idxs[11] = 11;
		boundary_clutter_projected[12] = Point_2( 0.00,  0.53); idxs[12] = 12;
		boundary_clutter_projected[13] = Point_2( 0.35,  0.60); idxs[13] = 13;
		boundary_clutter_projected[14] = Point_2( 0.42,  0.56); idxs[14] = 14;
		boundary_clutter_projected[15] = Point_2( 0.47,  0.66); idxs[15] = 15;
		boundary_clutter_projected[16] = Point_2( 0.40,  0.51); idxs[16] = 16;
		boundary_clutter_projected[17] = Point_2( 0.48,  0.48); idxs[17] = 17;
		boundary_clutter_projected[18] = Point_2( 0.51,  0.58); idxs[18] = 18;
		boundary_clutter_projected[19] = Point_2( 0.56,  0.53); idxs[19] = 19;
		boundary_clutter_projected[20] = Point_2( 0.43,  0.44); idxs[20] = 20;
		boundary_clutter_projected[21] = Point_2( 0.51,  0.39); idxs[21] = 21;
		boundary_clutter_projected[22] = Point_2( 0.59,  0.44); idxs[22] = 22;
		boundary_clutter_projected[23] = Point_2( 0.61,  0.36); idxs[23] = 23;

		boundary_clutter[0] = idxs;

		Log log; 
		log.export_projected_points_as_xyz("tmp/clutter_input", boundary_clutter_projected, "unused path");
	}
};

TEST_F(LOD_ClutterProcessorTest, Compiles) {

	// Empty test.
}

/*
TEST_F(LOD_ClutterProcessorTest, RunsThinning) {

	lodClutterProcessor.set_number_of_neighbours(3);
	const auto number_of_processed_points = lodClutterProcessor.apply_thinning(boundary_clutter, boundary_clutter_projected);

	ASSERT_THAT(number_of_processed_points, Eq(24));
}

TEST_F(LOD_ClutterProcessorTest, RunsGridSimplify) {

	lodClutterProcessor.set_grid_cell_length(0.05);
	const auto number_of_removed_points = lodClutterProcessor.apply_grid_simplify(boundary_clutter, boundary_clutter_projected);

	ASSERT_THAT(number_of_removed_points, Eq(10));
} */

TEST_F(LOD_ClutterProcessorTest, Runs) {

	lodClutterProcessor.set_number_of_neighbours(3);
	lodClutterProcessor.set_grid_cell_length(0.05);

	const auto number_of_removed_points = lodClutterProcessor.process(boundary_clutter, boundary_clutter_projected);

	ASSERT_THAT(number_of_removed_points, Eq(2));
}