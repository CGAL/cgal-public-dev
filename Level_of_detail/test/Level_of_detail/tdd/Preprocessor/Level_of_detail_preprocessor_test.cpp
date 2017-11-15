// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

// New CGAL includes.
#include <CGAL/Preprocessor/Level_of_detail_preprocessor.h>

using namespace testing;

class LOD_PreprocessorTest: public Test {

public:
	using FT = double;
	
	using Traits          = CGAL::Simple_cartesian<FT>;
	using Point_2         = Traits::Point_2;
	using Point_3         = Traits::Point_3;
	using Container       = CGAL::Point_set_3<Point_3>;
	using LodPreprocessor = CGAL::LOD::Level_of_detail_preprocessor<Traits, Container>;
	using Planes          = std::map<int, std::vector<int> >;
	using Boundary_data   = Planes;

	using Projected_points = std::map<int, Point_2>;

	// Here Index is the index of a plane from the input!
	using Index = int; 
	using Index_map = Container:: template Property_map<Index>; 
	using Iter = Container::iterator;

	using Point_index = Container::Index;
	using Indices = std::vector<Point_index>;


	LodPreprocessor lodPreprocessor;

	void set_indices_property(Container &input, Index_map &indices) const {
		
		auto success = false;
		boost::tie(indices, success)  = input. template add_property_map<Index>("index", -1);
		assert(success);
	}


	void get_simple_input(Container &input) const {

		Index_map indices;
		set_indices_property(input, indices);

		Iter 
		it = input.insert(Point_3(0, 0, 0)); indices[*it] =  0;
		it = input.insert(Point_3(0, 0, 1)); indices[*it] =  1;
		it = input.insert(Point_3(1, 0, 0)); indices[*it] =  0;
		it = input.insert(Point_3(1, 0, 1)); indices[*it] =  1;
		it = input.insert(Point_3(0, 1, 0)); indices[*it] =  0;
		it = input.insert(Point_3(2, 0, 0)); indices[*it] =  1;
		it = input.insert(Point_3(3, 0, 0)); indices[*it] = -1;
	}

	void get_simple_indices(Indices &mapping) const {

		mapping.clear();
		mapping.resize(4);

		mapping[0] = 1;
		mapping[1] = 3;
		mapping[2] = 5;
		mapping[3] = 6;
	}


	void get_alpha_shapes_input(Container &input) const {

		Index_map indices;
		set_indices_property(input, indices);

		// Boundary.
		Iter 
		it = input.insert(Point_3(0.0, 0.0, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(0.5, 0.0, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(1.0, 0.0, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(1.0, 0.5, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(1.0, 1.0, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(0.5, 1.0, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(0.0, 0.5, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(0.0, 1.0, 1.0)); indices[*it] = -1;

		// Interior.
		it = input.insert(Point_3(0.25, 0.25, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(0.50, 0.25, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(0.75, 0.25, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(0.75, 0.50, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(0.75, 0.75, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(0.50, 0.75, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(0.25, 0.75, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(0.25, 0.50, 1.0)); indices[*it] = -1;
		it = input.insert(Point_3(0.50, 0.50, 1.0)); indices[*it] = -1;
	}

	void get_alpha_shapes_indices(Indices &mapping) const {

		mapping.clear();
		mapping.resize(16);

		mapping[0] = 0; mapping[1] = 1; mapping[2] = 2; mapping[3] = 3;
		mapping[4] = 4; mapping[5] = 5; mapping[6] = 6; mapping[7] = 7;

		mapping[8]  =  8; mapping[9]  =  9; mapping[10] = 10; 
		mapping[11] = 11; mapping[12] = 12; mapping[13] = 13; 
		mapping[14] = 14; mapping[15] = 15; /* mapping[16] = 16; */
	}


	void get_projected_points(Projected_points &projected_boundaries, Boundary_data &boundaries) const {

		// First boundary.
		projected_boundaries[0] = Point_2(-0.4, 0.0); // outlier
		projected_boundaries[1] = Point_2(0.00, 0.0);
		projected_boundaries[2] = Point_2(0.10, 0.0);
		projected_boundaries[3] = Point_2(0.20, 0.0);

		boundaries[0].push_back(0);
		boundaries[0].push_back(1);
		boundaries[0].push_back(2);
		boundaries[0].push_back(3);

		// Second boundary.
		projected_boundaries[4] = Point_2(0.75, 0.0);
		projected_boundaries[5] = Point_2(1.00, 0.0);
		projected_boundaries[6] = Point_2(1.50, 0.0); // outlier

		boundaries[1].push_back(4);
		boundaries[1].push_back(5);
		boundaries[1].push_back(6);
	}
};


TEST_F(LOD_PreprocessorTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_PreprocessorTest, ReturnsTwoPlanes) {

	Container input;
	get_simple_input(input);

	Planes planes;
	const auto number_of_planes = lodPreprocessor.get_planes(input, planes);

	ASSERT_THAT(number_of_planes, Eq(2));
}

TEST_F(LOD_PreprocessorTest, FindsTwoPlanesInArbitraryOrder) {

	Container input;
	get_simple_input(input);

	Planes planes;
	lodPreprocessor.get_planes(input, planes);

	ASSERT_THAT(planes[0][0], Eq(0));
	ASSERT_THAT(planes[0][1], Eq(2));
	ASSERT_THAT(planes[0][2], Eq(4));

	ASSERT_THAT(planes[1][0], Eq(1));
	ASSERT_THAT(planes[1][1], Eq(3));
	ASSERT_THAT(planes[1][2], Eq(5));
}

TEST_F(LOD_PreprocessorTest, RejectsNegativeIndices) {

	Container input;
	Index_map indices;

	set_indices_property(input, indices);

	Iter
	it = input.insert(Point_3(0, 0, 0)); indices[*it] = -1;
	it = input.insert(Point_3(1, 1, 1)); indices[*it] = -1;

	Planes planes;
	const auto number_of_planes = lodPreprocessor.get_planes(input, planes);

	ASSERT_THAT(number_of_planes, Eq(0));
}

TEST_F(LOD_PreprocessorTest, ReturnsOnePlaneUsingGivenIndices) {

	Container input;
	get_simple_input(input);

	Indices mapping, stub;
	get_simple_indices(mapping);

	Boundary_data boundaries, boundary_clutter;

	const bool with_shape_detection = true;
	const auto number_of_boundaries = lodPreprocessor.get_boundary_points(input, mapping, stub, with_shape_detection, boundaries, boundary_clutter);

	ASSERT_THAT(number_of_boundaries, Eq(1));
	ASSERT_THAT(static_cast<int>((*boundaries.begin()).second.size()), Eq(3));
	ASSERT_THAT(static_cast<int>(boundary_clutter.at(0).size()), Eq(1));
}

TEST_F(LOD_PreprocessorTest, RemovesTwoOutliers) {

	Projected_points projected_boundaries; Boundary_data boundaries;
	get_projected_points(projected_boundaries, boundaries);

	const auto number_of_outliers = lodPreprocessor.clean_projected_points(projected_boundaries, boundaries);

	ASSERT_THAT(number_of_outliers, Eq(2));
	
	ASSERT_THAT(static_cast<int>(boundaries.at(0).size()), Eq(3));
	ASSERT_THAT(static_cast<int>(boundaries.at(1).size()), Eq(2));
}


TEST_F(LOD_PreprocessorTest, AppliesAlphaShapesWithTwoBoundaries) {

	Container input;
	get_alpha_shapes_input(input);

	Indices stub, mapping;
	get_alpha_shapes_indices(mapping);

	Boundary_data boundaries, boundary_clutter;

	const FT alpha = 0.05;
	lodPreprocessor.set_alpha(alpha);

	const bool with_shape_detection = false;
	const auto number_of_boundaries = lodPreprocessor.get_boundary_points(input, stub, mapping, with_shape_detection, boundaries, boundary_clutter);

	ASSERT_THAT(number_of_boundaries, Eq(-1));
	ASSERT_THAT(static_cast<int>(boundary_clutter.size()), Eq(1));
	ASSERT_THAT(static_cast<int>(boundary_clutter.at(0).size()), Eq(16));
}

TEST_F(LOD_PreprocessorTest, AppliesAlphaShapesWithOneBoundary) {

	Container input;
	get_alpha_shapes_input(input);

	Indices stub, mapping;
	get_alpha_shapes_indices(mapping);

	Boundary_data boundaries, boundary_clutter;

	const FT alpha = 1.5;
	lodPreprocessor.set_alpha(alpha);

	const bool with_shape_detection = false;
	const auto number_of_boundaries = lodPreprocessor.get_boundary_points(input, stub, mapping, with_shape_detection, boundaries, boundary_clutter);

	ASSERT_THAT(number_of_boundaries, Eq(-1));
	ASSERT_THAT(static_cast<int>(boundary_clutter.size()), Eq(1));
	ASSERT_THAT(static_cast<int>(boundary_clutter.at(0).size()), Eq(8));
}