// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <map>
#include <vector>
#include <memory>

// CGAL includes.
#include <CGAL/Simple_cartesian.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>
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

	using Structured_points  = std::vector< std::vector<Point> >; 			  
	using Structured_labels  = std::vector< std::vector<CGAL::LOD::Structured_label> >;  
	using Structured_anchors = std::vector< std::vector<std::vector<int> > >;

	std::unique_ptr<LodStructuring> lodStructuring;
	Points points; Components components; Lines lines;

	void create_test_input() {

		points.clear(); components.clear(); lines.clear();

		// First segment.
		points[0] = Point(0.0, 0.0);     components[0].push_back(0);
		points[1] = Point(0.2, 0.05);	 components[0].push_back(1);
		points[2] = Point(0.4, -0.025);  components[0].push_back(2);
		points[3] = Point(0.6, 0.025);   components[0].push_back(3);
		points[4] = Point(0.8, -0.05);   components[0].push_back(4);
		points[5] = Point(1.0, 0.0);     components[0].push_back(5);

		lines.push_back(Line(points[0], points[5]));

		// Second segment.
		points[6]  = Point(-0.01, 0.01); components[1].push_back(6);
		points[7]  = Point(-0.1, 0.1);   components[1].push_back(7);
		points[8]  = Point(-0.3, 0.3);   components[1].push_back(8);
		points[9]  = Point(-0.5, 0.5);   components[1].push_back(9);
		points[10] = Point(-0.7, 0.7);   components[1].push_back(10);
		points[11] = Point(-0.8, 0.8);   components[1].push_back(11);

		lines.push_back(Line(points[6], points[11]));

		// Third segment.
		points[12] = Point(-0.8, 0.4);   components[2].push_back(12);
		points[13] = Point(-0.8, 0.6);   components[2].push_back(13); 
		points[14] = Point(-0.8, 0.79);  components[2].push_back(14);

		lines.push_back(Line(points[12], points[14]));

		// Fourth segment.
		points[15] = Point(-0.2, 0.8);   components[3].push_back(15);
		points[16] = Point(-0.3, 0.8);   components[3].push_back(16);
		points[17] = Point(-0.4, 0.8);   components[3].push_back(17);
		points[18] = Point(-0.79, 0.8);  components[3].push_back(18);

		lines.push_back(Line(points[15], points[18]));

		// Fifth segment.
		points[19] = Point(1.0, 0.4);    components[4].push_back(19);
		points[20] = Point(1.0, 0.42);   components[4].push_back(20);
		points[21] = Point(1.0, 0.63);   components[4].push_back(21);
		points[22] = Point(1.0, 0.77);   components[4].push_back(22);
		points[23] = Point(1.0, 0.8);    components[4].push_back(23);

		lines.push_back(Line(points[19], points[23]));

		lodStructuring = std::make_unique<LodStructuring>(points, components, lines);
		lodStructuring->save_log(true);
	}

	LOD_StructuringTest() {
		
		create_test_input();
		lodStructuring->structure_point_set();
	}
};

TEST_F(LOD_StructuringTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_StructuringTest, VerifiesStructuredPoints) {

	const Structured_points &str_points = lodStructuring->get_structured_points();

	const double eps = 0.000001;
	ASSERT_THAT(str_points[0][0], Eq(Point(0.0, 0.0)));

	ASSERT_THAT(str_points[2][str_points[2].size() - 1].x() + 0.8, Lt(eps));
	ASSERT_THAT(str_points[2][str_points[2].size() - 1].y() - 0.8, Lt(eps));
}

TEST_F(LOD_StructuringTest, VerifiesStructuredLabels) {

	const Structured_labels &str_labels = lodStructuring->get_structured_labels();

	ASSERT_THAT(str_labels[0][0], Eq(CGAL::LOD::Structured_label::CORNER));
	ASSERT_THAT(str_labels[0][str_labels[0].size() - 1], Eq(CGAL::LOD::Structured_label::LINEAR));
}

TEST_F(LOD_StructuringTest, VerifiesStructuredAnchors) {

	const Structured_anchors &str_anchors = lodStructuring->get_structured_anchors();

	ASSERT_THAT(str_anchors[3][2].size(), Eq(static_cast<size_t>(1)));
	ASSERT_THAT(str_anchors[3][2][0], Eq(3));

	ASSERT_THAT(str_anchors[2][str_anchors[2].size() - 1].size(), Eq(static_cast<size_t>(2)));

	ASSERT_THAT(str_anchors[2][str_anchors[2].size() - 1][0], Eq(2));
	ASSERT_THAT(str_anchors[2][str_anchors[2].size() - 1][1], Eq(1));
}