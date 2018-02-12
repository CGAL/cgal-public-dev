// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <map>
#include <list>
#include <vector>

// CGAL includes.
#include <CGAL/number_utils.h>
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

	// Line regularizer Jean Philippe - original method.
	using Boundary_data_stub    = std::vector<FT>;
	using Projected_points_stub = std::vector<FT>;

	using LodLineRegularizerJeanPhilippe = CGAL::LOD::Level_of_detail_line_regularizer_jean_philippe<Kernel, Boundary_data_stub, Projected_points_stub>;
	LodLineRegularizerJeanPhilippe lodLineRegularizerJeanPhilippe;

	// Segment regularizer - updated method.
	using LodSegmentRegularizer = CGAL::LOD::Level_of_detail_segment_regularizer_2<Kernel>;
	LodSegmentRegularizer lodSegmentRegularizer;

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

	void load_segments_from_file(Segments &segments) {
		segments.clear();

		const std::string path = "/Users/danisimo/Documents/pipeline/logs/regularizer-data/segments.data";
		std::ifstream loader(path.c_str(), std::ios_base::in);

        if (!loader) {
            std::cerr << std::endl << std::endl << "ERROR: Error loading file with segments!" << std::endl << std::endl;
            exit(EXIT_FAILURE);
        }

		size_t num_segments;
		loader >> num_segments;

		segments.resize(num_segments);

		FT x1, y1, x2, y2;
		for (size_t i = 0; i < num_segments; ++i) {

            loader >> x1 >> y1 >> x2 >> y2;
			segments[i] = Segment(Point(x1, y1), Point(x2, y2));
		}
	}

	void create_segments(Segments &original, Segments &updated) {

		create_simple_test_segments(original);
		create_simple_test_segments(updated);
	}

	void load_segments(Segments &original, Segments &updated) {

		load_segments_from_file(original);
		load_segments_from_file(updated);
	}

	void apply_original_method(Segments &segments) {
		
		Boundary_data_stub    boundary_data_stub;
		Projected_points_stub projected_points_stub;
		Lines 				  lines_stub;
	
		lodLineRegularizerJeanPhilippe.process(boundary_data_stub, projected_points_stub, segments, lines_stub);
		segments = lodLineRegularizerJeanPhilippe.get_regularized_segments();
	}

	void apply_updated_method(Segments &segments) {
		lodSegmentRegularizer.regularize(segments, SegmentMap());
	}
};

TEST_F(LOD_SegmentRegularizerTest, IsEqualToOriginalMethodSmall) {
	
	Segments original, updated;
	create_segments(original, updated);
	
	apply_original_method(original);
	apply_updated_method(updated);

	ASSERT_THAT(original.size(), Eq(updated.size()));

	const FT eps = FT(1) / FT(1000000);
	for (size_t i = 0; i < original.size(); ++i) {
		ASSERT_LT(CGAL::abs(updated[i].source().x() - original[i].source().x()), eps);
		ASSERT_LT(CGAL::abs(updated[i].source().y() - original[i].source().y()), eps);

		ASSERT_LT(CGAL::abs(updated[i].target().x() - original[i].target().x()), eps);
		ASSERT_LT(CGAL::abs(updated[i].target().y() - original[i].target().y()), eps);
	}
}

TEST_F(LOD_SegmentRegularizerTest, IsEqualToOriginalMethodBig) {
	
	Segments original, updated;
	load_segments(original, updated);
	
	apply_original_method(original);
	apply_updated_method(updated);

	ASSERT_THAT(original.size(), Eq(updated.size()));

	const FT eps = FT(1) / FT(1000000);
	for (size_t i = 0; i < original.size(); ++i) {
		ASSERT_LT(CGAL::abs(updated[i].source().x() - original[i].source().x()), eps);
		ASSERT_LT(CGAL::abs(updated[i].source().y() - original[i].source().y()), eps);

		ASSERT_LT(CGAL::abs(updated[i].target().x() - original[i].target().x()), eps);
		ASSERT_LT(CGAL::abs(updated[i].target().y() - original[i].target().y()), eps);
	}
}