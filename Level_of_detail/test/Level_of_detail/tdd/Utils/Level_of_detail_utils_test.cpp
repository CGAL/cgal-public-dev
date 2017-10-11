// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <string>

// CGAL includes.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// CGAL new includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Utils/Level_of_detail_utils.h>

using namespace testing;

class LOD_UtilsTest: public Test {

public:
	using FT = double;

	using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
	using Point     = Kernel::Point_3;
	using Container = CGAL::Point_set_3<Point>;

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

	LodUtils lodUtils;
};

TEST_F(LOD_UtilsTest, Compiles) {

	// Empty test.
}