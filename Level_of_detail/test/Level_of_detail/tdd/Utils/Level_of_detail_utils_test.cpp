#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\" 
#else 
#define PSR "/" 
#endif 

// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <string>
#include <vector>

// CGAL includes.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/intersections.h>

// CGAL new includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Utils/Level_of_detail_utils.h>

using namespace testing;

class LOD_UtilsTest: public Test {

public:
	using FT = double;

	using Kernel      = CGAL::Exact_predicates_inexact_constructions_kernel;
	using Point_2     = Kernel::Point_2;
	using Point_3     = Kernel::Point_3;
	using Line_2      = Kernel::Line_2;
	using Segment_2   = Kernel::Segment_2;
	using Intersect_2 = Kernel::Intersect_2; 

	using Container = CGAL::Point_set_3<Point_3>;

	using Str_label = CGAL::LOD::Structured_label;

	using My_vertex_info = CGAL::LOD::My_vertex_info<Str_label>; 
	using My_face_info   = CGAL::LOD::My_face_info<FT>;

	using VB 		   = CGAL::Triangulation_vertex_base_with_info_2<My_vertex_info, Kernel>;
	using FB_with_info = CGAL::Triangulation_face_base_with_info_2<My_face_info, Kernel>;
	using FB 		   = CGAL::Constrained_triangulation_face_base_2<Kernel, FB_with_info>;

	using EPT = CGAL::Exact_predicates_tag;
	using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
	using CDT = CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, EPT>;

	using LodUtils = CGAL::LOD::Level_of_detail_utils<Kernel, Container, CDT>;
	using Log = CGAL::LOD::Mylog;

	using Mapping 		   = std::map<int, std::vector<int> >;
	using Projected_points = std::map<int, Point_2>;
	using Lines 		   = std::vector<Line_2>;
	using Segments 		   = std::vector<Segment_2>;
	using Boxes 		   = std::vector< std::vector<Point_2> >;

	LodUtils lodUtils;

	void load_data(Projected_points &points, Mapping &mapping, Boxes &boxes) {

		points.clear();
		mapping.clear();
		boxes.clear();

		std::string file_path = "/Users/danisimo/Documents/pipeline/bugs/data1.pwn";
		std::ifstream loader(file_path.c_str(), std::ios_base::in);

        if (!loader) {
            std::cerr << std::endl << std::endl << "ERROR: Error loading file with data!" << std::endl << std::endl;
            exit(EXIT_FAILURE);
        }

        const size_t num_points = 713;
		FT x, y, z;

		const FT big_value = FT(1000000000000);
		FT minx =  big_value, miny =  big_value;
		FT maxx = -big_value, maxy = -big_value;

        for (size_t i = 0; i < num_points; ++i) {
        	loader >> x >> y >> z;

			points[i] = Point_2(x, y);
			mapping[0].push_back(i);

			minx = CGAL::min(minx, x);
			miny = CGAL::min(miny, y);

			maxx = CGAL::max(maxx, x);
			maxy = CGAL::max(maxy, y);
		}
		loader.close();

		boxes.resize(1);
		boxes[0].push_back(Point_2(minx, miny));
		boxes[0].push_back(Point_2(maxx, miny));
		boxes[0].push_back(Point_2(maxx, maxy));
		boxes[0].push_back(Point_2(minx, maxy));
	}

	void intersect_lines_with_boxes(const Lines &lines, const Boxes &boxes, 
	Segments &segments) {
		
		assert(lines.size() == boxes.size());
		
		segments.clear();
		segments.resize(lines.size());

		for (size_t i = 0; i < lines.size(); ++i) {
			const Line_2 &line = lines[i];

			const std::vector<Point_2> &box = boxes[i];
			std::vector<Point_2> intersections;

			for (size_t j = 0; j < box.size(); ++j) {
				const size_t jp = (j + 1) % box.size();

				const Point_2 &p1 = box[j];
				const Point_2 &p2 = box[jp];

				const Segment_2 edge(p1, p2);
				typename CGAL::cpp11::result_of<Intersect_2(Line_2, Segment_2)>::type result = CGAL::intersection(line, edge);

				if (result)
					intersections.push_back(boost::get<Point_2>(*result));
			}
			segments[i] = Segment_2(intersections[0], intersections[1]);
		}
	}

	void create_segments_from_lines(const Projected_points &points, const Mapping &mapping, const Lines &lines, 
	Segments &segments) {
		
		if (lines.size() == 0) return;

		segments.clear();
		lodUtils.create_segments_from_lines(points, mapping, lines, segments);
		
		// Segments new_segments;
		// lodUtils.create_segments_from_lines(points, mapping, segments, new_segments);
		// segments = new_segments;
	}
};

TEST_F(LOD_UtilsTest, Compiles) {

	Projected_points points;
	Mapping mapping;
	Boxes boxes;

	load_data(points, mapping, boxes);

	Lines lines;
	lodUtils.fit_lines_to_projected_points(points, mapping, lines);

	// intersect_lines_with_boxes(lines, boxes, segments);
	
	Segments segments;
	create_segments_from_lines(points, mapping, lines, segments);

	std::string stub = "";
	Log segments_exporter;
	segments_exporter.export_segments_as_obj("tmp" + std::string(PSR) + "segments", segments, stub);
}