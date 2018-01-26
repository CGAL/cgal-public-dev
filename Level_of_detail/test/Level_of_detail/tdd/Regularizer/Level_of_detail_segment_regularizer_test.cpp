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
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_2.h>

using namespace testing;

class LOD_SegmentRegularizerTest: public Test {
	
public:
	using KernelTraits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using FT           = typename KernelTraits::FT;
    using Point        = typename KernelTraits::Point_2;
	using Segment      = typename KernelTraits::Segment_2;

    using Segments = std::list<Segment>;
	using LodSegmentRegularizer = CGAL::LOD::Level_of_detail_segment_regularizer_2<KernelTraits>;
	using SegmentMap = CGAL::Identity_property_map<Segment>;

	LodSegmentRegularizer lodSegmentRegularizer;

	void create_simple_test_segments(Segments &segments) const {
		segments.clear();

		segments.push_back(Segment(Point(0.00, 0.00), Point(1.0, 0.0)));
		segments.push_back(Segment(Point(0.05, 0.05), Point(1.0, 0.1)));

		segments.push_back(Segment(Point(0.0, 0.0), Point(0.0, 1.0)));
		segments.push_back(Segment(Point(0.0, 1.0), Point(1.0, 1.0)));

		segments.push_back(Segment(Point(0.69, 0.51), Point(0.7, 0.68)));
		segments.push_back(Segment(Point(0.75, 0.70), Point(1.0, 0.70)));

		segments.push_back(Segment(Point(0.38, 0.65), Point(0.4, 0.8)));
		segments.push_back(Segment(Point(0.50, 0.65), Point(0.5, 0.8)));
	}
};

TEST_F(LOD_SegmentRegularizerTest, Compiles) {
	
	Segments segments;
	create_simple_test_segments(segments);

	lodSegmentRegularizer.regularize(segments, SegmentMap());
}