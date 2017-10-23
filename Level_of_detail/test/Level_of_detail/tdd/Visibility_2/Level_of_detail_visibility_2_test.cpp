// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <map>
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Visibility_2/Level_of_detail_visibility_2.h>

using namespace testing;

class LOD_VisibilityTest: public Test {

public:
	using FT = double;
	
	using Traits  = CGAL::Simple_cartesian<FT>;
	using Point_2 = Traits::Point_2;

	using Label     = int; 
	using Container = std::vector< std::pair<Point_2, Label> >;

	using My_face_info = CGAL::LOD::My_face_info<FT>;

	using VB 		   = CGAL::Triangulation_vertex_base_2<Traits>;
	using FB_with_info = CGAL::Triangulation_face_base_with_info_2<My_face_info, Traits>;
	using FB 		   = CGAL::Constrained_triangulation_face_base_2<Traits, FB_with_info>;

	using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
	using CDT = CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS>;

	using LodVisibilityWithClassification = CGAL::LOD::Level_of_detail_visibility_from_classification_2<Traits, Container, CDT>;
	using LodVisibilityWithRayShooting    = CGAL::LOD::Level_of_detail_visibility_ray_shooting_2<Traits, Container, CDT>;
	using LodVisibilityWithBlend          = CGAL::LOD::Level_of_detail_visibility_blend_2<Traits, Container, CDT>;

	using Vertex_handle = CDT::Vertex_handle;
	using Log = CGAL::LOD::Mylog;

	LodVisibilityWithClassification lodVisibilityCL;
	LodVisibilityWithRayShooting    lodVisibilityRS;
	LodVisibilityWithBlend          lodVisibilityBL;

	CDT cdt; Container input;

	LOD_VisibilityTest() {
		create_data();
	}

	void create_data() {

		cdt.clear();
		input.clear();		

		set_basic_input(cdt, input);
	}

	void set_basic_input(CDT &cdt, Container &input) {

		const Label ground     = 0;
		const Label facade     = 1;
		const Label roof       = 2;
		const Label vegetation = 3; 
 
		// IN
		input.push_back(std::make_pair(Point_2(0.25, 0.25), roof));

		// IN
		input.push_back(std::make_pair(Point_2(0.70, 0.70), roof));
		input.push_back(std::make_pair(Point_2(0.70, 0.90), ground));
		input.push_back(std::make_pair(Point_2(0.50, 0.80), roof));

		// OUT
		input.push_back(std::make_pair(Point_2(1.30, 0.30), ground));
		input.push_back(std::make_pair(Point_2(1.30, 0.10), facade));

		// OUT
		input.push_back(std::make_pair(Point_2(1.60, 0.60), vegetation));

		// UNKNOWN
		input.push_back(std::make_pair(Point_2(1.00, 0.80), facade));

		Vertex_handle va = cdt.insert(Point_2(0, 0));
		Vertex_handle vb = cdt.insert(Point_2(1, 0));
		Vertex_handle vc = cdt.insert(Point_2(0, 1));
		Vertex_handle vd = cdt.insert(Point_2(1, 1));
		
		cdt.insert(Point_2(2, 0));
		cdt.insert(Point_2(2, 1));
		cdt.insert(Point_2(3, 0));

		// Building.
		cdt.insert_constraint(va, vb);
		cdt.insert_constraint(vb, vd);
		cdt.insert_constraint(vd, vc);
		cdt.insert_constraint(vc, va);
	}
};

TEST_F(LOD_VisibilityTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_VisibilityTest, WithClassification) {

	lodVisibilityCL.compute(input, cdt);
	Log log; log.save_visibility_eps(cdt, "tmp/visibility_classification");

	auto face = cdt.finite_faces_begin();

	ASSERT_GT((++face)->info().in, 0.5); ++face;
	ASSERT_LT((face)->info().in, 0.5);
}

TEST_F(LOD_VisibilityTest, WithRayShooting) {

	lodVisibilityRS.compute(input, cdt);
	Log log; log.save_visibility_eps(cdt, "tmp/visibility_ray_shooting");

	auto face = cdt.finite_faces_begin();

	ASSERT_GT((++face)->info().in, 0.5); ++face;
	ASSERT_LT((face)->info().in, 0.5);
}

TEST_F(LOD_VisibilityTest, WithBlend) {

	lodVisibilityBL.compute(input, cdt);
	Log log; log.save_visibility_eps(cdt, "tmp/visibility_blend");

	auto face = cdt.finite_faces_begin();

	ASSERT_GT((++face)->info().in, 0.5); ++face;
	ASSERT_LT((face)->info().in, 0.5);
}