// Google test includes.
#include "gmock/gmock.h"

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

// New CGAL includes.
#include <CGAL/Selector/Level_of_detail_selector.h>
#include <CGAL/Selector/Level_of_detail_selection_strategy.h>

using namespace testing;

class LOD_SelectorTest: public Test {

public:
	using FT = double;
	
	using Traits      = CGAL::Simple_cartesian<FT>;
	using Point       = Traits::Point_3;
	using Container   = CGAL::Point_set_3<Point>;
	using Iterator    = Container::iterator;
	using Index       = Container::Index;

	using ClutterStrategy 		   = CGAL::LOD::Level_of_detail_clutter<Traits, Container>;
	using GroundStrategy  		   = CGAL::LOD::Level_of_detail_ground<Traits, Container>; 
	using BuildingBoundaryStrategy = CGAL::LOD::Level_of_detail_building_boundary<Traits, Container>;
	using BuildingInteriorStrategy = CGAL::LOD::Level_of_detail_building_interior<Traits, Container>;

	using ClutterSelector 		   = CGAL::LOD::Level_of_detail_selector<Traits, ClutterStrategy>;
	using GroundSelector  		   = CGAL::LOD::Level_of_detail_selector<Traits, GroundStrategy>; 
	using BuildingBoundarySelector = CGAL::LOD::Level_of_detail_selector<Traits, BuildingBoundaryStrategy>;
	using BuildingInteriorSelector = CGAL::LOD::Level_of_detail_selector<Traits, BuildingInteriorStrategy>;

	using Label     = int; 
	using Label_map = Container:: template Property_map<Label>; 

	ClutterSelector 		 clutterSelector;
	GroundSelector  		 groundSelector;
	BuildingBoundarySelector buildingBoundarySelector;
	BuildingInteriorSelector buildingInteriorSelector;

	void set_labels_property(Container &input, Label_map &labels) {
		
		auto success = false;
		boost::tie(labels, success) = input. template add_property_map<Label>("label", -1);
		assert(success);
	}

	void set_basic_input(Container &input) {

		Label_map labels;
		set_labels_property(input, labels);

		Iterator 
		it = input.insert(Point(0, 0, 0)); labels[*it] = 0; 
		it = input.insert(Point(1, 1, 1)); labels[*it] = 1;
		it = input.insert(Point(2, 2, 2)); labels[*it] = 2;
		it = input.insert(Point(3, 3, 3)); labels[*it] = 3;
	}
};

TEST_F(LOD_SelectorTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_SelectorTest, ReturnsClutterWithOneElement) {

	std::vector<Index> output;
	
	Container input;
	set_basic_input(input);

	clutterSelector.select_elements(input, std::back_inserter(output));

	ASSERT_THAT(static_cast<int>(output.size()), Eq(1));
	ASSERT_THAT(static_cast<int>(output[0]), Eq(3));
}

TEST_F(LOD_SelectorTest, ReturnsGroundWithOneElement) {

	std::vector<Index> output;
	
	Container input;
	set_basic_input(input);

	groundSelector.select_elements(input, std::back_inserter(output));

	ASSERT_THAT(static_cast<int>(output.size()), Eq(1));
	ASSERT_THAT(static_cast<int>(output[0]), Eq(0));
}

TEST_F(LOD_SelectorTest, ReturnsBuildingBoundaryWithOneElement) {

	std::vector<Index> output;
	
	Container input;
	set_basic_input(input);

	buildingBoundarySelector.select_elements(input, std::back_inserter(output));

	ASSERT_THAT(static_cast<int>(output.size()), Eq(1));
	ASSERT_THAT(static_cast<int>(output[0]), Eq(1));
}

TEST_F(LOD_SelectorTest, ReturnsBuildingInteriorWithOneElement) {

	std::vector<Index> output;
	
	Container input;
	set_basic_input(input);

	buildingInteriorSelector.select_elements(input, std::back_inserter(output));

	ASSERT_THAT(static_cast<int>(output.size()), Eq(1));
	ASSERT_THAT(static_cast<int>(output[0]), Eq(2));
}
