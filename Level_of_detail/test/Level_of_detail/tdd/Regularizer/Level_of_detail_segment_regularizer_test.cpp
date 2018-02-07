// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <map>
#include <list>
#include <vector>

// CGAL includes.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/property_map.h>

// New CGAL includes.
#include <CGAL/Regularizer/Level_of_detail_line_regularizer_jean_philippe.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_2.h>

using namespace testing;

class LOD_SegmentRegularizerTest: public Test {
	
public:
	
	// Basic.
	using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
    using FT      = typename Kernel::FT;
    using Point   = typename Kernel::Point_2;
	using Segment = typename Kernel::Segment_2;
	using Line    = typename Kernel::Line_2;

    using Segments   = std::vector<Segment>;
	using SegmentMap = CGAL::Identity_property_map<Segment>;
	using Lines		 = std::vector<Line>;

	// Segment regularizer.
	using LodSegmentRegularizer = CGAL::LOD::Level_of_detail_segment_regularizer_2<Kernel>;
	LodSegmentRegularizer lodSegmentRegularizer;


	// Line regularizer Jean Philippe.
	using Boundary_data_stub    = std::vector<FT>;
	using Projected_points_stub = std::vector<FT>;

	using LodLineRegularizerJeanPhilippe = CGAL::LOD::Level_of_detail_line_regularizer_jean_philippe<Kernel, Boundary_data_stub, Projected_points_stub>;
	LodLineRegularizerJeanPhilippe lodLineRegularizerJeanPhilippe;


	// Functions.
	void create_simple_test_segments(Segments &segments) const {
		segments.clear();

		segments.push_back(Segment(Point(0.00, 0.00), Point(1.0, 0.0)));
		segments.push_back(Segment(Point(0.05, 0.05), Point(1.0, 0.1)));

		segments.push_back(Segment(Point(0.0, 0.0), Point(0.0, 1.00)));
		segments.push_back(Segment(Point(0.0, 1.0), Point(1.0, 0.95)));

		segments.push_back(Segment(Point(0.69, 0.51), Point(0.7, 0.68)));
		segments.push_back(Segment(Point(0.75, 0.70), Point(1.0, 0.70)));

		segments.push_back(Segment(Point(0.38, 0.65), Point(0.4, 0.8)));
		segments.push_back(Segment(Point(0.50, 0.65), Point(0.5, 0.8)));
	}
};

TEST_F(LOD_SegmentRegularizerTest, ExecutesSegmentRegularizer) {
	
	Segments segments;
	create_simple_test_segments(segments);

	lodSegmentRegularizer.regularize(segments, SegmentMap());
}

TEST_F(LOD_SegmentRegularizerTest, ExecutesLineRegularizerJeanPhilippe) {
	
	Segments segments;
	create_simple_test_segments(segments);

	Boundary_data_stub    boundary_data_stub;
	Projected_points_stub projected_points_stub;

	Lines lines;
	lodLineRegularizerJeanPhilippe.process(boundary_data_stub, projected_points_stub, segments, lines);
}