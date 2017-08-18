// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <map>

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
	using Point           = Traits::Point_3;
	using Container       = CGAL::Point_set_3<Point>;
	using LodPreprocessor = CGAL::LOD::Level_of_detail_preprocessor<Traits>;
	using Planes          = std::map<int, std::vector<int> >;

	// Here Index is the index of a plane from the input!
	using Index = int; 
	using Index_map = Container:: template Property_map<Index>; 
	using Iter = Container::iterator;

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
		it = input.insert(Point(0, 0, 0)); indices[*it] = 0;
		it = input.insert(Point(0, 0, 1)); indices[*it] = 1;
		it = input.insert(Point(1, 0, 0)); indices[*it] = 0;
		it = input.insert(Point(1, 0, 1)); indices[*it] = 1;
		it = input.insert(Point(0, 1, 0)); indices[*it] = 0;
		it = input.insert(Point(2, 0, 0)); indices[*it] = 1;
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
	it = input.insert(Point(0, 0, 0)); indices[*it] = -1;
	it = input.insert(Point(1, 1, 1)); indices[*it] = -1;

	Planes planes;
	const auto number_of_planes = lodPreprocessor.get_planes(input, planes);

	ASSERT_THAT(number_of_planes, Eq(0));
}
