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

	Boundary_data bd_thinning, bd_grid_simplify;
	Projected_points pp_thinning, pp_grid_simplify;
	

	LOD_ClutterProcessorTest() {
		create_data();
	}

	void create_data() {

		bd_thinning.clear();
		bd_grid_simplify.clear();

		pp_thinning.clear();
		pp_grid_simplify.clear();

		lodClutterProcessor.set_fitter_type(LodClutterProcessor::Fitter_type::LINE);
		lodClutterProcessor.set_new_point_type(LodClutterProcessor::New_point_type::CENTROID);

		set_thinning_input();
		set_grid_simplify_input();
	}

	void set_thinning_input() {

		std::vector<int> idxs(24);

		pp_thinning[0]  = Point_2( 0.05,  0.04); idxs[0]  =  0;
		pp_thinning[1]  = Point_2( 0.11, -0.04); idxs[1]  =  1;
		pp_thinning[2]  = Point_2( 0.25,  0.05); idxs[2]  =  2;
		pp_thinning[3]  = Point_2( 0.28, -0.06); idxs[3]  =  3;
		pp_thinning[4]  = Point_2( 0.39,  0.05); idxs[4]  =  4;
		pp_thinning[5]  = Point_2( 0.49, -0.06); idxs[5]  =  5;
		pp_thinning[6]  = Point_2(-0.06, -0.06); idxs[6]  =  6;
		pp_thinning[7]  = Point_2(-0.08,  0.04); idxs[7]  =  7;
		pp_thinning[8]  = Point_2(-0.14, -0.03); idxs[8]  =  8;
		pp_thinning[9]  = Point_2( 0.03,  0.16); idxs[9]  =  9;
		pp_thinning[10] = Point_2(-0.05,  0.28); idxs[10] = 10;
		pp_thinning[11] = Point_2( 0.02,  0.44); idxs[11] = 11;
		pp_thinning[12] = Point_2( 0.00,  0.53); idxs[12] = 12;
		pp_thinning[13] = Point_2( 0.35,  0.60); idxs[13] = 13;
		pp_thinning[14] = Point_2( 0.42,  0.56); idxs[14] = 14;
		pp_thinning[15] = Point_2( 0.47,  0.66); idxs[15] = 15;
		pp_thinning[16] = Point_2( 0.40,  0.51); idxs[16] = 16;
		pp_thinning[17] = Point_2( 0.48,  0.48); idxs[17] = 17;
		pp_thinning[18] = Point_2( 0.51,  0.58); idxs[18] = 18;
		pp_thinning[19] = Point_2( 0.56,  0.53); idxs[19] = 19;
		pp_thinning[20] = Point_2( 0.43,  0.44); idxs[20] = 20;
		pp_thinning[21] = Point_2( 0.51,  0.39); idxs[21] = 21;
		pp_thinning[22] = Point_2( 0.59,  0.44); idxs[22] = 22;
		pp_thinning[23] = Point_2( 0.61,  0.36); idxs[23] = 23;

		bd_thinning[0] = idxs;

		Log log; 
		log.export_projected_points_as_xyz("tmp/thinning_input", pp_thinning, "unused path");
	}

	void set_grid_simplify_input() {

		std::vector<int> idxs(13);

		// First cell.
		pp_grid_simplify[0] = Point_2(0.20, 0.40); idxs[0] = 0;
		pp_grid_simplify[1] = Point_2(0.20, 0.09); idxs[1] = 1;
		pp_grid_simplify[2] = Point_2(0.69, 0.90); idxs[2] = 2;
		pp_grid_simplify[3] = Point_2(0.91, 0.91); idxs[3] = 3;
		pp_grid_simplify[4] = Point_2(0.69, 0.59); idxs[4] = 4;
		
		// Second cell.
		pp_grid_simplify[5] = Point_2(1.10, 1.81); idxs[5] = 5;
		pp_grid_simplify[6] = Point_2(1.61, 1.31); idxs[6] = 6;
		pp_grid_simplify[7] = Point_2(1.80, 1.10); idxs[7] = 7;

		// Negative cell.
		pp_grid_simplify[8]  = Point_2(-0.77, -0.26); idxs[8]  =  8;		
		pp_grid_simplify[9]  = Point_2(-0.61, -0.82); idxs[9]  =  9;		
		pp_grid_simplify[10] = Point_2(-0.11, -0.59); idxs[10] = 10;
		pp_grid_simplify[11] = Point_2(-0.19, -0.30); idxs[11] = 11;
		pp_grid_simplify[12] = Point_2(-0.33, -0.08); idxs[12] = 12;		

		bd_grid_simplify[0] = idxs;

		Log log; 
		log.export_projected_points_as_xyz("tmp/grid_simplify_input", pp_grid_simplify, "unused path");
	}
};

TEST_F(LOD_ClutterProcessorTest, Compiles) {

	// Empty test.
}

TEST_F(LOD_ClutterProcessorTest, RunsThinning) {

	lodClutterProcessor.set_number_of_neighbours(3);
	const auto number_of_processed_points = lodClutterProcessor.apply_thinning(bd_thinning, pp_thinning);

	ASSERT_THAT(number_of_processed_points, Eq(24));
}

TEST_F(LOD_ClutterProcessorTest, RunsGridSimplifyWithCellLengthOne) {

	lodClutterProcessor.set_grid_cell_length(1.0);
	const auto number_of_removed_points = lodClutterProcessor.apply_grid_simplify(bd_grid_simplify, pp_grid_simplify);

	ASSERT_THAT(number_of_removed_points, Eq(10));

	ASSERT_THAT(pp_grid_simplify[0], Eq(Point_2(-0.5, -0.5)));
	ASSERT_THAT(pp_grid_simplify[1], Eq(Point_2( 0.5,  0.5)));
	ASSERT_THAT(pp_grid_simplify[2], Eq(Point_2( 1.5,  1.5)));
}

TEST_F(LOD_ClutterProcessorTest, RunsGridSimplifyWithCellLengthHalf) {

	lodClutterProcessor.set_grid_cell_length(0.5);
	const auto number_of_removed_points = lodClutterProcessor.apply_grid_simplify(bd_grid_simplify, pp_grid_simplify);

	ASSERT_THAT(number_of_removed_points, Eq(5));
}