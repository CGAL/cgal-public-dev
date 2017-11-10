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
#include <CGAL/Region_growing_2/Level_of_detail_region_growing_2.h>

using namespace testing;

class LOD_RegionGrowingTest: public Test {

public:
	using FT = double;

	using Kernel  = CGAL::Simple_cartesian<FT>;
	using Point_2 = Kernel::Point_2;
	using Point_3 = Kernel::Point_3;
	using Normal  = Kernel::Vector_3; 

	using Boundary_data    = std::map<int, std::vector<int> >;
	using Projected_points = std::map<int, Point_2>;
	using Container 	   = CGAL::Point_set_3<Point_3>;

	using LodRegionGrowing = CGAL::LOD::Level_of_detail_region_growing_2<Kernel, Boundary_data, Projected_points, Container>;
	using Log 			   = CGAL::LOD::Mylog;
	

	LodRegionGrowing lodRegionGrowing;

	Boundary_data 	 boundary_clutter, building_boundaries;
	Projected_points boundary_clutter_projected, building_boundaries_projected;
	Container 		 input;


	LOD_RegionGrowingTest() {
		create_data();
	}


	void create_data() {

		boundary_clutter.clear(); building_boundaries.clear();
		boundary_clutter_projected.clear(); building_boundaries_projected.clear();
		
		input.clear();
		input.add_normal_map();

		set_new_data();
	}

	void set_new_data() {

		lodRegionGrowing.set_epsilon(0.02);
		lodRegionGrowing.set_normal_threshold(0.9);
		lodRegionGrowing.set_minimum_shape_points(10);
		lodRegionGrowing.save_info(true);

		// Horizontal line.
		boundary_clutter_projected[0]  = Point_2(1.50, 1.50); boundary_clutter[0].push_back(0); // query point in the centre

		boundary_clutter_projected[1]  = Point_2(1.05, 1.51); boundary_clutter[0].push_back(1);
		boundary_clutter_projected[2]  = Point_2(1.15, 1.49); boundary_clutter[0].push_back(2);
		boundary_clutter_projected[3]  = Point_2(1.25, 1.51); boundary_clutter[0].push_back(3);
		boundary_clutter_projected[4]  = Point_2(1.35, 1.49); boundary_clutter[0].push_back(4);
		boundary_clutter_projected[5]  = Point_2(1.45, 1.51); boundary_clutter[0].push_back(5);
		boundary_clutter_projected[6]  = Point_2(1.55, 1.49); boundary_clutter[0].push_back(6);
		boundary_clutter_projected[7]  = Point_2(1.65, 1.51); boundary_clutter[0].push_back(7);
		boundary_clutter_projected[8]  = Point_2(1.75, 1.49); boundary_clutter[0].push_back(8);
		boundary_clutter_projected[9]  = Point_2(1.85, 1.51); boundary_clutter[0].push_back(9);
		boundary_clutter_projected[10] = Point_2(1.95, 1.49); boundary_clutter[0].push_back(10);

		Normal normal = Normal(0.0, 2.0, 0.0);

		input.insert(Point_3(1.50, 1.50, 0.0), normal);

		input.insert(Point_3(1.05, 1.51, 0.0), normal);
		input.insert(Point_3(1.15, 1.49, 0.0), normal);
		input.insert(Point_3(1.25, 1.51, 0.0), normal);
		input.insert(Point_3(1.35, 1.49, 0.0), normal);
		input.insert(Point_3(1.45, 1.51, 0.0), normal);
		input.insert(Point_3(1.55, 1.49, 0.0), normal);
		input.insert(Point_3(1.65, 1.51, 0.0), normal);
		input.insert(Point_3(1.75, 1.49, 0.0), normal);
		input.insert(Point_3(1.85, 1.51, 0.0), normal);
		input.insert(Point_3(1.95, 1.49, 0.0), normal);


		// Vertical line.
		boundary_clutter_projected[11] = Point_2(1.25, 1.95); boundary_clutter[0].push_back(11);
		boundary_clutter_projected[12] = Point_2(1.25, 1.85); boundary_clutter[0].push_back(12);
		boundary_clutter_projected[13] = Point_2(1.25, 1.75); boundary_clutter[0].push_back(13);
		boundary_clutter_projected[14] = Point_2(1.25, 1.65); boundary_clutter[0].push_back(14);
		boundary_clutter_projected[15] = Point_2(1.25, 1.55); boundary_clutter[0].push_back(15);
		boundary_clutter_projected[16] = Point_2(1.25, 1.45); boundary_clutter[0].push_back(16);
		boundary_clutter_projected[17] = Point_2(1.25, 1.35); boundary_clutter[0].push_back(17);
		boundary_clutter_projected[18] = Point_2(1.25, 1.25); boundary_clutter[0].push_back(18);
		boundary_clutter_projected[19] = Point_2(1.25, 1.15); boundary_clutter[0].push_back(19);
		boundary_clutter_projected[20] = Point_2(1.25, 1.05); boundary_clutter[0].push_back(20);

		normal = Normal(1.0, 0.0, 0.0);

		input.insert(Point_3(1.25, 1.95, 0.0), normal);
		input.insert(Point_3(1.25, 1.85, 0.0), normal);
		input.insert(Point_3(1.25, 1.75, 0.0), normal);
		input.insert(Point_3(1.25, 1.65, 0.0), normal);
		input.insert(Point_3(1.25, 1.55, 0.0), normal);
		input.insert(Point_3(1.25, 1.45, 0.0), normal);
		input.insert(Point_3(1.25, 1.35, 0.0), normal);
		input.insert(Point_3(1.25, 1.25, 0.0), normal);
		input.insert(Point_3(1.25, 1.15, 0.0), normal);
		input.insert(Point_3(1.25, 1.05, 0.0), normal);
	}
};


TEST_F(LOD_RegionGrowingTest, Compiles) {

	// Empty test.
}

TEST_F(LOD_RegionGrowingTest, DetectesOneLine) {

	lodRegionGrowing.set_cluster_epsilon(0.15);

	const auto number_of_detected_lines = lodRegionGrowing.detect(
		boundary_clutter, boundary_clutter_projected,
		building_boundaries, building_boundaries_projected,
		input);

	ASSERT_THAT(number_of_detected_lines, Eq(1));
}

TEST_F(LOD_RegionGrowingTest, DetectesTwoLines) {

	lodRegionGrowing.set_cluster_epsilon(0.2);

	const auto number_of_detected_lines = lodRegionGrowing.detect(
		boundary_clutter, boundary_clutter_projected,
		building_boundaries, building_boundaries_projected,
		input);

	ASSERT_THAT(number_of_detected_lines, Eq(2));
}