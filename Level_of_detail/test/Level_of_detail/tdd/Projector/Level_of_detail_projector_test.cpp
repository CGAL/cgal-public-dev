// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

// New CGAL includes.
#include <CGAL/Projector/Level_of_detail_projector.h>

using namespace testing;

class LOD_ProjectorTest: public Test {
	
public:
	using FT = double;
	
	using Traits = CGAL::Simple_cartesian<FT>;
	
	using Point_3 = Traits::Point_3;
	using Point_2 = Traits::Point_2;
	using Normal  = Traits::Vector_3;
	using Plane   = Traits::Plane_3;

	using Planes         = std::map<int, std::vector<int> >;
	using Container      = CGAL::Point_set_3<Point_3>;
	using Output		 = std::map<int, Point_2>;
	using LodProjector   = CGAL::LOD::Level_of_detail_simple_projector<Traits, Container, Planes, Output>;

	using Iterator  = Container::iterator;

	LodProjector lodProjector;

	void create_test_input(Container &input, Planes &indices, Plane &ground) const {

		Iterator 
		it = input.insert(Point_3(0.2, 0.3, 0.22));
		it = input.insert(Point_3(0.4, 0.6, 0.22));
		it = input.insert(Point_3(0.8, 0.1, 0.22));

		it = input.insert(Point_3(1.00, 0.2, 0.9));
		it = input.insert(Point_3(1.05, 0.5, 0.3));
		it = input.insert(Point_3(1.00, 0.8, 0.9));

		indices[0].push_back(0); indices[0].push_back(1); indices[0].push_back(2);
		indices[1].push_back(3); indices[1].push_back(4); indices[1].push_back(5);

		ground = Plane(0.0, 0.0, 2.56, 0.0);
	}
};

TEST_F(LOD_ProjectorTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_ProjectorTest, ReturnsSixProjectedPoints) {

	Container input; Planes indices; Plane ground; Output output;

	create_test_input(input, indices, ground);
	const auto number_of_projected_points = lodProjector.project(input, indices, ground, output);

	ASSERT_THAT(number_of_projected_points, Eq(6));
}

TEST_F(LOD_ProjectorTest, MakesZCoordinateZero) {

	Container input; Planes indices; Plane ground; Output output;

	create_test_input(input, indices, ground);
	lodProjector.project(input, indices, ground, output);

	ASSERT_THAT(output[indices[0][0]], Eq(Point_2(0.2, 0.3)));
	ASSERT_THAT(output[indices[1][0]], Eq(Point_2(1.0, 0.2)));
}
