// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <map>
#include <vector>
#include <memory>

// CGAL includes.
#include <CGAL/Simple_cartesian.h>

// New CGAL includes.
#include <CGAL/Structuring_2/Level_of_detail_structuring_2.h>

using namespace testing;

class LOD_StructuringTest: public Test {
	
public:
	using FT = double;

	using Traits = CGAL::Simple_cartesian<FT>;
	
	using Point   = Traits::Point_2;
	using Line    = Traits::Line_2;
	using Segment = Traits::Segment_2;

	using Points     = std::map<int, Point>;
	using Components = std::map<int, std::vector<int> >;
	using Lines      = std::vector<Line>;
	using Segments   = std::vector<Segment>;

	using LodStructuring = CGAL::LOD::Level_of_detail_structuring_2<Traits>;

	std::unique_ptr<LodStructuring> lodStructuring;
	Points points; Components components; Lines lines;

	void create_test_input() {

		points.clear(); components.clear(); lines.clear();

		// First segment.
		points[0] = Point(0.0, 0.0);    components[0].push_back(0);
		points[1] = Point(0.2, 0.05);	components[0].push_back(1);
		points[2] = Point(0.4, -0.025); components[0].push_back(2);
		points[3] = Point(0.6, 0.025);  components[0].push_back(3);
		points[4] = Point(0.8, -0.5);   components[0].push_back(4);
		points[5] = Point(1.0, 0.0);    components[0].push_back(5);

		lines.push_back(Line(points[0], points[5]));

		// Second segment.


		lodStructuring = std::make_unique<LodStructuring>(points, components, lines);
	}
};

TEST_F(LOD_StructuringTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_StructuringTest, SavesLogFiles) {

	Segments segments;
	create_test_input();

	lodStructuring->structure_point_set(segments);
}
