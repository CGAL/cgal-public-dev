// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <map>
#include <vector>
#include <string>

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

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
	using Point_3 = Kernel::Point_3;
	using Normal  = Kernel::Vector_3; 

	using Boundary_data    = std::map<int, std::vector<int> >;
	using Projected_points = std::map<int, Point_2>;
	using Container 	   = CGAL::Point_set_3<Point_3>;

	using LodClutterProcessor = CGAL::LOD::Level_of_detail_clutter_processor<Kernel, Boundary_data, Projected_points, Container>;
	using LodClutterThinning  = CGAL::LOD::Level_of_detail_thinning<Kernel, Boundary_data, Projected_points, Container>;
	
	using Log = CGAL::LOD::Mylog;
	
	LodClutterProcessor lodClutterProcessor;
	LodClutterThinning  lodClutterThinning;


	// Simple data.
	Container 		 in_stub;
	Boundary_data    bd_thinning_naive, bd_grid_simplify;
	Projected_points pp_thinning_naive, pp_grid_simplify;
	

	// Complex data.
	Boundary_data 	 bd_complex_test, bd_complex;
	Projected_points pp_complex_test, pp_complex;
	Container 		 in_complex_test, in_complex;


	LOD_ClutterProcessorTest() {
		create_data();
	}


	void create_data() {

		bd_thinning_naive.clear();
		bd_grid_simplify.clear();

		pp_thinning_naive.clear();
		pp_grid_simplify.clear();

		set_thinning_naive_input();
		set_grid_simplify_input();

		/*
		clear_complex_data();
		set_complex_input(); */
	}


	// Simple.
	void set_thinning_naive_input() {

		std::vector<int> idxs(24);

		pp_thinning_naive[0]  = Point_2( 0.05,  0.04); idxs[0]  =  0;
		pp_thinning_naive[1]  = Point_2( 0.11, -0.04); idxs[1]  =  1;
		pp_thinning_naive[2]  = Point_2( 0.25,  0.05); idxs[2]  =  2;
		pp_thinning_naive[3]  = Point_2( 0.28, -0.06); idxs[3]  =  3;
		pp_thinning_naive[4]  = Point_2( 0.39,  0.05); idxs[4]  =  4;
		pp_thinning_naive[5]  = Point_2( 0.49, -0.06); idxs[5]  =  5;
		pp_thinning_naive[6]  = Point_2(-0.06, -0.06); idxs[6]  =  6;
		pp_thinning_naive[7]  = Point_2(-0.08,  0.04); idxs[7]  =  7;
		pp_thinning_naive[8]  = Point_2(-0.14, -0.03); idxs[8]  =  8;
		pp_thinning_naive[9]  = Point_2( 0.03,  0.16); idxs[9]  =  9;
		pp_thinning_naive[10] = Point_2(-0.05,  0.28); idxs[10] = 10;
		pp_thinning_naive[11] = Point_2( 0.02,  0.44); idxs[11] = 11;
		pp_thinning_naive[12] = Point_2( 0.00,  0.53); idxs[12] = 12;
		pp_thinning_naive[13] = Point_2( 0.35,  0.60); idxs[13] = 13;
		pp_thinning_naive[14] = Point_2( 0.42,  0.56); idxs[14] = 14;
		pp_thinning_naive[15] = Point_2( 0.47,  0.66); idxs[15] = 15;
		pp_thinning_naive[16] = Point_2( 0.40,  0.51); idxs[16] = 16;
		pp_thinning_naive[17] = Point_2( 0.48,  0.48); idxs[17] = 17;
		pp_thinning_naive[18] = Point_2( 0.51,  0.58); idxs[18] = 18;
		pp_thinning_naive[19] = Point_2( 0.56,  0.53); idxs[19] = 19;
		pp_thinning_naive[20] = Point_2( 0.43,  0.44); idxs[20] = 20;
		pp_thinning_naive[21] = Point_2( 0.51,  0.39); idxs[21] = 21;
		pp_thinning_naive[22] = Point_2( 0.59,  0.44); idxs[22] = 22;
		pp_thinning_naive[23] = Point_2( 0.61,  0.36); idxs[23] = 23;

		bd_thinning_naive[0] = idxs;

		Log log; 
		log.export_projected_points_as_xyz("tmp/thinning_naive_input", pp_thinning_naive, "unused path");
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


	// Complex.
	void clear_complex_data() {

		bd_complex_test.clear();
		pp_complex_test.clear();
		in_complex_test.clear();

		bd_complex.clear();
		pp_complex.clear();
		in_complex.clear();
	}

	void set_complex_input() {

		/* create_complex_test(); */
		/* create_complex_data(); */
	}

	void create_complex_test() {

		std::vector<int> idxs(3);

		pp_complex_test[0] = Point_2(0.2, 0.0); idxs[0] = 0;
		pp_complex_test[1] = Point_2(0.4, 0.0); idxs[1] = 1;
		pp_complex_test[2] = Point_2(0.8, 0.0); idxs[2] = 2;

		in_complex_test.add_normal_map();
		const Normal normal = Normal(0.0, 0.22, 0.0); 

		in_complex_test.insert(Point_3(0.2, 0.0, 0.3), normal);
		in_complex_test.insert(Point_3(0.2, 0.0, 0.3), normal);
		in_complex_test.insert(Point_3(0.2, 0.0, 0.3), normal);

		bd_complex_test[0] = idxs;

		Log log; 
		log.export_projected_points_as_xyz("tmp/complex_test", pp_complex_test, "unused path");
	}

	void create_complex_data() {

		lodClutterThinning.set_fuzzy_radius(0.6);
		lodClutterThinning.set_thinning_type(CGAL::LOD::Thinning_type::COMPLEX);
		lodClutterThinning.set_fitter_type(CGAL::LOD::Thinning_fitter_type::LINE);
		lodClutterThinning.set_neighbour_search_type(CGAL::LOD::Neighbour_search_type::CIRCLE);

		in_complex.add_normal_map();

		std::vector<int> idxs;
		bd_complex[0] = idxs;

		// Uncomment the one you need!
		// add_complex_data_1();
		// add_complex_data_2();
		// add_complex_data_3();
		// add_complex_data_4();
		// add_complex_data_5();
		// add_complex_data_6();

		// Log log; 
		// log.export_projected_points_as_xyz("tmp/complex", pp_complex, "unused path");
	}

	// Dense with dim = 0. Tested!
	void add_complex_data_1() {

		lodClutterThinning.set_dim_0_max(0.05);

		pp_complex[0] = Point_2(0.50, 0.50); bd_complex[0].push_back(0);

		pp_complex[1] = Point_2(0.51, 0.49); bd_complex[0].push_back(1);
		pp_complex[2] = Point_2(0.51, 0.51); bd_complex[0].push_back(2);
		pp_complex[3] = Point_2(0.49, 0.51); bd_complex[0].push_back(3);
		pp_complex[4] = Point_2(0.49, 0.49); bd_complex[0].push_back(4);

		const Normal normal = Normal(0.0, 1.0, 0.0);

		in_complex.insert(Point_3(0.50, 0.50, 0.0), normal);

		in_complex.insert(Point_3(0.51, 0.49, 0.0), normal);
		in_complex.insert(Point_3(0.51, 0.51, 0.0), normal);
		in_complex.insert(Point_3(0.49, 0.51, 0.0), normal);
		in_complex.insert(Point_3(0.49, 0.49, 0.0), normal);
	}

	// Sparse with dim = 2. Tested!
	void add_complex_data_2() {

		lodClutterThinning.set_dim_2_sparsity_num_cells(6.0);
		lodClutterThinning.set_dim_2_sparsity_max(0.25);

		pp_complex[0]  = Point_2(1.50, 0.50); bd_complex[0].push_back(0);

		pp_complex[1]  = Point_2(1.15, 0.85); bd_complex[0].push_back(1);
		pp_complex[2]  = Point_2(1.15, 0.55); bd_complex[0].push_back(2);
		pp_complex[3]  = Point_2(1.15, 0.05); bd_complex[0].push_back(3);
		pp_complex[4]  = Point_2(1.25, 0.26); bd_complex[0].push_back(4);
		pp_complex[5]  = Point_2(1.45, 0.75); bd_complex[0].push_back(5);
		pp_complex[6]  = Point_2(1.54, 0.95); bd_complex[0].push_back(6);
		pp_complex[7]  = Point_2(1.55, 0.15); bd_complex[0].push_back(7);
		pp_complex[8]  = Point_2(1.75, 0.35); bd_complex[0].push_back(8);
		pp_complex[9]  = Point_2(1.84, 0.85); bd_complex[0].push_back(9);
		pp_complex[10] = Point_2(1.85, 0.64); bd_complex[0].push_back(10);
		pp_complex[11] = Point_2(1.83, 0.06); bd_complex[0].push_back(11);

		const Normal normal = Normal(0.0, 1.0, 0.0);

		in_complex.insert(Point_3(1.50, 0.50, 0.0), normal);

		in_complex.insert(Point_3(1.15, 0.85, 0.0), normal);
		in_complex.insert(Point_3(1.15, 0.55, 0.0), normal);
		in_complex.insert(Point_3(1.15, 0.05, 0.0), normal);
		in_complex.insert(Point_3(1.25, 0.26, 0.0), normal);
		in_complex.insert(Point_3(1.45, 0.75, 0.0), normal);
		in_complex.insert(Point_3(1.54, 0.95, 0.0), normal);
		in_complex.insert(Point_3(1.55, 0.15, 0.0), normal);
		in_complex.insert(Point_3(1.75, 0.35, 0.0), normal);
		in_complex.insert(Point_3(1.84, 0.85, 0.0), normal);
		in_complex.insert(Point_3(1.85, 0.64, 0.0), normal);
		in_complex.insert(Point_3(1.83, 0.06, 0.0), normal);
	}

	// Multiple clusters with dim = 2. Tested!
	void add_complex_data_3() {

		pp_complex[0]  = Point_2(2.50, 0.50); bd_complex[0].push_back(0);

		pp_complex[1]  = Point_2(2.15, 0.75); bd_complex[0].push_back(1);
		pp_complex[2]  = Point_2(2.15, 0.65); bd_complex[0].push_back(2);
		pp_complex[3]  = Point_2(2.25, 0.65); bd_complex[0].push_back(3);
		pp_complex[4]  = Point_2(2.25, 0.75); bd_complex[0].push_back(4);

		pp_complex[5]  = Point_2(2.45, 0.25); bd_complex[0].push_back(5);
		pp_complex[6]  = Point_2(2.45, 0.15); bd_complex[0].push_back(6);
		pp_complex[7]  = Point_2(2.55, 0.15); bd_complex[0].push_back(7);
		pp_complex[8]  = Point_2(2.55, 0.25); bd_complex[0].push_back(8);

		pp_complex[9]  = Point_2(2.65, 0.85); bd_complex[0].push_back(9);
		pp_complex[10] = Point_2(2.65, 0.75); bd_complex[0].push_back(10);
		pp_complex[11] = Point_2(2.75, 0.75); bd_complex[0].push_back(11);
		pp_complex[12] = Point_2(2.75, 0.85); bd_complex[0].push_back(12);

		const Normal normal = Normal(0.0, 1.0, 0.0);

		in_complex.insert(Point_3(2.50, 0.50, 0.0), normal);

		in_complex.insert(Point_3(2.15, 0.75, 0.0), normal);
		in_complex.insert(Point_3(2.15, 0.65, 0.0), normal);
		in_complex.insert(Point_3(2.25, 0.65, 0.0), normal);
		in_complex.insert(Point_3(2.25, 0.75, 0.0), normal);

		in_complex.insert(Point_3(2.45, 0.25, 0.0), normal);
		in_complex.insert(Point_3(2.45, 0.15, 0.0), normal);
		in_complex.insert(Point_3(2.55, 0.15, 0.0), normal);
		in_complex.insert(Point_3(2.55, 0.25, 0.0), normal);

		in_complex.insert(Point_3(2.65, 0.85, 0.0), normal);
		in_complex.insert(Point_3(2.65, 0.75, 0.0), normal);
		in_complex.insert(Point_3(2.75, 0.75, 0.0), normal);
		in_complex.insert(Point_3(2.75, 0.85, 0.0), normal);
	}

	// Line with dim = 1. Untested!
	void add_complex_data_4() {

		pp_complex[0]  = Point_2(0.50, 1.50); bd_complex[0].push_back(0);

		pp_complex[1]  = Point_2(0.05, 1.51); bd_complex[0].push_back(1);
		pp_complex[2]  = Point_2(0.15, 1.49); bd_complex[0].push_back(2);
		pp_complex[3]  = Point_2(0.25, 1.51); bd_complex[0].push_back(3);
		pp_complex[4]  = Point_2(0.35, 1.49); bd_complex[0].push_back(4);
		pp_complex[5]  = Point_2(0.45, 1.51); bd_complex[0].push_back(5);
		pp_complex[6]  = Point_2(0.55, 1.49); bd_complex[0].push_back(6);
		pp_complex[7]  = Point_2(0.65, 1.51); bd_complex[0].push_back(7);
		pp_complex[8]  = Point_2(0.75, 1.49); bd_complex[0].push_back(8);
		pp_complex[9]  = Point_2(0.85, 1.51); bd_complex[0].push_back(9);
		pp_complex[10] = Point_2(0.95, 1.49); bd_complex[0].push_back(10);		

		const Normal normal = Normal(0.0, 1.0, 0.0);

		in_complex.insert(Point_3(0.50, 1.50, 0.0), normal);

		in_complex.insert(Point_3(0.05, 1.51, 0.0), normal);
		in_complex.insert(Point_3(0.15, 1.49, 0.0), normal);
		in_complex.insert(Point_3(0.25, 1.51, 0.0), normal);
		in_complex.insert(Point_3(0.35, 1.49, 0.0), normal);
		in_complex.insert(Point_3(0.45, 1.51, 0.0), normal);
		in_complex.insert(Point_3(0.55, 1.49, 0.0), normal);
		in_complex.insert(Point_3(0.65, 1.51, 0.0), normal);
		in_complex.insert(Point_3(0.75, 1.49, 0.0), normal);
		in_complex.insert(Point_3(0.85, 1.51, 0.0), normal);
		in_complex.insert(Point_3(0.95, 1.49, 0.0), normal);
	}

	// Corner with dim = 1. Untested!
	void add_complex_data_5() {

		// Horizontal line.
		pp_complex[0]  = Point_2(1.50, 1.50); bd_complex[0].push_back(0);

		pp_complex[1]  = Point_2(1.05, 1.51); bd_complex[0].push_back(1);
		pp_complex[2]  = Point_2(1.15, 1.49); bd_complex[0].push_back(2);
		pp_complex[3]  = Point_2(1.25, 1.51); bd_complex[0].push_back(3);
		pp_complex[4]  = Point_2(1.35, 1.49); bd_complex[0].push_back(4);
		pp_complex[5]  = Point_2(1.45, 1.51); bd_complex[0].push_back(5);
		pp_complex[6]  = Point_2(1.55, 1.49); bd_complex[0].push_back(6);
		pp_complex[7]  = Point_2(1.65, 1.51); bd_complex[0].push_back(7);
		pp_complex[8]  = Point_2(1.75, 1.49); bd_complex[0].push_back(8);
		pp_complex[9]  = Point_2(1.85, 1.51); bd_complex[0].push_back(9);
		pp_complex[10] = Point_2(1.95, 1.49); bd_complex[0].push_back(10);

		Normal normal = Normal(0.0, 1.0, 0.0);

		in_complex.insert(Point_3(1.50, 1.50, 0.0), normal);

		in_complex.insert(Point_3(1.05, 1.51, 0.0), normal);
		in_complex.insert(Point_3(1.15, 1.49, 0.0), normal);
		in_complex.insert(Point_3(1.25, 1.51, 0.0), normal);
		in_complex.insert(Point_3(1.35, 1.49, 0.0), normal);
		in_complex.insert(Point_3(1.45, 1.51, 0.0), normal);
		in_complex.insert(Point_3(1.55, 1.49, 0.0), normal);
		in_complex.insert(Point_3(1.65, 1.51, 0.0), normal);
		in_complex.insert(Point_3(1.75, 1.49, 0.0), normal);
		in_complex.insert(Point_3(1.85, 1.51, 0.0), normal);
		in_complex.insert(Point_3(1.95, 1.49, 0.0), normal);


		// Vertical line.
		pp_complex[11] = Point_2(1.25, 1.95); bd_complex[0].push_back(11);
		pp_complex[12] = Point_2(1.25, 1.85); bd_complex[0].push_back(12);
		pp_complex[13] = Point_2(1.25, 1.75); bd_complex[0].push_back(13);
		pp_complex[14] = Point_2(1.25, 1.65); bd_complex[0].push_back(14);
		pp_complex[15] = Point_2(1.25, 1.55); bd_complex[0].push_back(15);
		pp_complex[16] = Point_2(1.25, 1.45); bd_complex[0].push_back(16);
		pp_complex[17] = Point_2(1.25, 1.35); bd_complex[0].push_back(17);
		pp_complex[18] = Point_2(1.25, 1.25); bd_complex[0].push_back(18);
		pp_complex[19] = Point_2(1.25, 1.15); bd_complex[0].push_back(19);
		pp_complex[20] = Point_2(1.25, 1.05); bd_complex[0].push_back(20);

		normal = Normal(1.0, 0.0, 0.0);

		in_complex.insert(Point_3(1.25, 1.95, 0.0), normal);
		in_complex.insert(Point_3(1.25, 1.85, 0.0), normal);
		in_complex.insert(Point_3(1.25, 1.75, 0.0), normal);
		in_complex.insert(Point_3(1.25, 1.65, 0.0), normal);
		in_complex.insert(Point_3(1.25, 1.55, 0.0), normal);
		in_complex.insert(Point_3(1.25, 1.45, 0.0), normal);
		in_complex.insert(Point_3(1.25, 1.35, 0.0), normal);
		in_complex.insert(Point_3(1.25, 1.25, 0.0), normal);
		in_complex.insert(Point_3(1.25, 1.15, 0.0), normal);
		in_complex.insert(Point_3(1.25, 1.05, 0.0), normal);
	}

	// Arch with dim = 1. Untested!
	void add_complex_data_6() {

		assert(!"The arch data set is not yet prepared!");
	}
};


// Tests.

// Simple tests.
TEST_F(LOD_ClutterProcessorTest, Compiles) {

	// Empty test.
}

TEST_F(LOD_ClutterProcessorTest, RunsNaiveThinning) {

	lodClutterProcessor.set_number_of_neighbours(3);
	lodClutterProcessor.set_thinning_type(CGAL::LOD::Thinning_type::NAIVE);
	lodClutterProcessor.set_fitter_type(CGAL::LOD::Thinning_fitter_type::LINE);
	lodClutterProcessor.set_neighbour_search_type(CGAL::LOD::Neighbour_search_type::KNN);

	const auto number_of_processed_points = lodClutterProcessor.apply_thinning(bd_thinning_naive, pp_thinning_naive, in_stub);
	ASSERT_THAT(number_of_processed_points, Eq(24));
}

TEST_F(LOD_ClutterProcessorTest, RunsGridSimplifyWithCellLengthOne) {

	lodClutterProcessor.set_grid_cell_length(1.0);
	lodClutterProcessor.set_new_point_type(CGAL::LOD::Grid_new_point_type::CENTROID);

	const auto number_of_removed_points = lodClutterProcessor.apply_grid_simplify(bd_grid_simplify, pp_grid_simplify);
	ASSERT_THAT(number_of_removed_points, Eq(10));

	ASSERT_THAT(pp_grid_simplify[0], Eq(Point_2(-0.5, -0.5)));
	ASSERT_THAT(pp_grid_simplify[1], Eq(Point_2( 0.5,  0.5)));
	ASSERT_THAT(pp_grid_simplify[2], Eq(Point_2( 1.5,  1.5)));
}

TEST_F(LOD_ClutterProcessorTest, RunsGridSimplifyWithCellLengthHalf) {

	lodClutterProcessor.set_grid_cell_length(0.5);
	lodClutterProcessor.set_new_point_type(CGAL::LOD::Grid_new_point_type::BARYCENTRE);

	const auto number_of_removed_points = lodClutterProcessor.apply_grid_simplify(bd_grid_simplify, pp_grid_simplify);
	ASSERT_THAT(number_of_removed_points, Eq(5));
}


// Complex tests.
/*
TEST_F(LOD_ClutterProcessorTest, RunsComplexThinningTest) {

	lodClutterProcessor.set_fuzzy_radius(0.25);
	lodClutterProcessor.set_thinning_type(CGAL::LOD::Thinning_type::COMPLEX);
	lodClutterProcessor.set_fitter_type(CGAL::LOD::Thinning_fitter_type::LINE);
	lodClutterProcessor.set_neighbour_search_type(CGAL::LOD::Neighbour_search_type::CIRCLE);

	const auto number_of_processed_points = lodClutterProcessor.apply_thinning(bd_complex_test, pp_complex_test, in_complex_test);
	ASSERT_THAT(number_of_processed_points, Eq(-1));
}

TEST_F(LOD_ClutterProcessorTest, RunsComplexThinning) {

	const auto number_of_processed_points = lodClutterThinning.process(bd_complex, pp_complex, in_complex);
	ASSERT_THAT(number_of_processed_points, Eq(-1));
} */