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

	using LodVisibility = CGAL::LOD::Level_of_detail_visibility_from_classification_2<Traits, Container, CDT>;

	using Vertex_handle = CDT::Vertex_handle;
	using Log = CGAL::LOD::Mylog;

	LodVisibility lodVisibility;
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
		input.push_back(std::make_pair(Point_2(1.30, 0.10), roof));

		// OUT
		input.push_back(std::make_pair(Point_2(1.60, 0.60), vegetation));

		// UNKNOWN
		input.push_back(std::make_pair(Point_2(1.00, 0.80), facade));

		Vertex_handle va = cdt.insert(Point_2(0, 0));
		Vertex_handle vb = cdt.insert(Point_2(1, 0));
		Vertex_handle vc = cdt.insert(Point_2(0, 1));
		Vertex_handle vd = cdt.insert(Point_2(1, 1));
		Vertex_handle ve = cdt.insert(Point_2(2, 0));
		Vertex_handle vf = cdt.insert(Point_2(2, 1));
		Vertex_handle vg = cdt.insert(Point_2(3, 0));

		cdt.insert_constraint(va, vb);
		cdt.insert_constraint(vb, vc);
		cdt.insert_constraint(vc, va);

		cdt.insert_constraint(vb, vd);
		cdt.insert_constraint(vd, vc);
		cdt.insert_constraint(vc, vb);

		cdt.insert_constraint(vb, ve);
		cdt.insert_constraint(ve, vd);
		cdt.insert_constraint(vd, vb);

		cdt.insert_constraint(ve, vf);
		cdt.insert_constraint(vf, vd);
		cdt.insert_constraint(vd, ve);		

		cdt.insert_constraint(ve, vg);
		cdt.insert_constraint(vg, vf);
		cdt.insert_constraint(vf, ve);		
	}
};

TEST_F(LOD_VisibilityTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_VisibilityTest, SavesVisibility) {

	lodVisibility.compute(input, cdt);

	Log log;
	log.save_visibility_eps(cdt);
}

TEST_F(LOD_VisibilityTest, VerifiesLabels) {

	lodVisibility.compute(input, cdt);

	auto face = cdt.finite_faces_begin();

	++face; ++face;
	ASSERT_LT((++face)->info().in, 0.5);
	ASSERT_THAT((++face)->info().in, Eq(0.5));
}
