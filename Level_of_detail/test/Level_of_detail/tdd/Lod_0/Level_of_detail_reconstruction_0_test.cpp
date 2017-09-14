// Google test includes.
#include "gmock/gmock.h"

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// New CGAL includes.
#include <CGAL/Lod_0/Level_of_detail_reconstruction_0.h>
#include <CGAL/Visibility_2/Level_of_detail_visibility_2.h>
#include <CGAL/Structuring_2/Level_of_detail_structuring_2.h>
#include <CGAL/Mylog/Mylog.h>

using namespace testing;

class LOD_ReconstructionTest: public Test {

public:
	using FT = double;
	
	using Traits    = CGAL::Simple_cartesian<FT>;
	using Point_2   = Traits::Point_2;
	using Point_3   = Traits::Point_3;
	using Container = CGAL::Point_set_3<Point_3>;
	using Iterator  = Container::iterator;

	using Structured_label  = typename CGAL::LOD::Level_of_detail_structuring_2<Traits>::Structured_label;
	using Structured_labels = std::vector<std::vector<Structured_label> >;
	using Structured_points = std::vector<std::vector<Point_2> >;

	using VB = CGAL::Triangulation_vertex_base_with_info_2<Structured_label, Traits>;
	using FB = CGAL::Constrained_triangulation_face_base_2<Traits>;

	using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
	using CDT = CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS>;

	using LodVisibility = CGAL::LOD::Level_of_detail_visibility_from_classification_2<Traits, Container, CDT>;
	using Visibility_result = std::map<int, LodVisibility::Visibility_label>;
	using VL = LodVisibility::Visibility_label;

	using Lod_0 = CGAL::LOD::Level_of_detail_reconstruction_0<Traits, CDT, Visibility_result, Structured_labels>;
	using Vertex_handle = CDT::Vertex_handle;

	using Log = CGAL::LOD::Mylog;

	using Segment = Traits::Segment_2;
	using Lod_0_result = std::vector<Segment>;

	Lod_0 lod_0;
	CDT cdt; Visibility_result visibility; 
	Structured_points str_points; Structured_labels str_labels;

	LOD_ReconstructionTest() {
		create_data();
	}

	void create_data() {

		cdt.clear();
		visibility.clear();

		str_points.clear();
		str_labels.clear();

		set_basic_input();
	}

	void set_basic_input() {

		// First building.
		Vertex_handle va = cdt.insert(Point_2(0.0, 0.0)); va->info() = Structured_label::CORNER;
		Vertex_handle vb = cdt.insert(Point_2(0.5, 0.0)); vb->info() = Structured_label::LINEAR;
		Vertex_handle vc = cdt.insert(Point_2(1.0, 0.0)); vc->info() = Structured_label::CORNER;
		Vertex_handle vd = cdt.insert(Point_2(1.0, 0.5)); vd->info() = Structured_label::LINEAR;
		Vertex_handle ve = cdt.insert(Point_2(1.0, 1.0)); ve->info() = Structured_label::CORNER;
		Vertex_handle vf = cdt.insert(Point_2(0.5, 1.0)); vf->info() = Structured_label::LINEAR;
		Vertex_handle vg = cdt.insert(Point_2(0.0, 1.0)); vg->info() = Structured_label::CORNER;
		Vertex_handle vh = cdt.insert(Point_2(0.0, 0.5)); vh->info() = Structured_label::LINEAR;

		cdt.insert_constraint(va, vb);
		cdt.insert_constraint(vb, vc);
		cdt.insert_constraint(vc, vd);
		cdt.insert_constraint(vd, ve);
		cdt.insert_constraint(ve, vf);
		cdt.insert_constraint(vf, vg);
		cdt.insert_constraint(vg, vh);
		cdt.insert_constraint(vh, va);

		// Second building.
		va = cdt.insert(Point_2(2.0, 0.0)); va->info() = Structured_label::CORNER;
		vb = cdt.insert(Point_2(2.5, 0.0)); vb->info() = Structured_label::LINEAR;
		vc = cdt.insert(Point_2(3.0, 0.0)); vc->info() = Structured_label::CORNER;
		vd = cdt.insert(Point_2(3.0, 0.5)); vd->info() = Structured_label::LINEAR;
		ve = cdt.insert(Point_2(3.0, 1.0)); ve->info() = Structured_label::CORNER;
		vf = cdt.insert(Point_2(2.5, 1.0)); vf->info() = Structured_label::LINEAR;
		vg = cdt.insert(Point_2(2.0, 1.0)); vg->info() = Structured_label::CORNER;
		vh = cdt.insert(Point_2(2.0, 0.5)); vh->info() = Structured_label::LINEAR;

		cdt.insert_constraint(va, vb);
		cdt.insert_constraint(vb, vc);
		cdt.insert_constraint(vc, vd);
		cdt.insert_constraint(vd, ve);
		cdt.insert_constraint(ve, vf);
		cdt.insert_constraint(vf, vg);
		cdt.insert_constraint(vg, vh);
		cdt.insert_constraint(vh, va);

		// Visibility predictions.
		visibility[0]  = VL::IN;  // 2 3 4
		visibility[1]  = VL::IN;  // 6 2 4
		visibility[2]  = VL::IN;  // 5 6 4
		visibility[3]  = VL::IN;  // 6 7 8
		visibility[4]  = VL::OUT; // 9 16 4
		visibility[5]  = VL::IN;  // 8 1 2
		visibility[6]  = VL::IN;  // 6 8 2
		visibility[7]  = VL::OUT; // 3 9 4
		visibility[8]  = VL::IN;  // 10 12 14
		visibility[9]  = VL::IN;  // 11 12 10
		visibility[10] = VL::IN;  // 12 13 14
		visibility[11] = VL::OUT; // 16 5 4
		visibility[12] = VL::IN;  // 16 14 15
		visibility[13] = VL::OUT; // 16 15 5
		visibility[14] = VL::IN;  // 10 16 9
		visibility[15] = VL::IN;  // 10 14 16

		// Structured points.
		str_points.resize(8);

		str_points[0].push_back(Point_2(0.0, 0.0));
		str_points[0].push_back(Point_2(1.0, 0.0));

		str_points[1].push_back(Point_2(1.0, 0.0));
		str_points[1].push_back(Point_2(1.0, 1.0));
		
		str_points[2].push_back(Point_2(1.0, 1.0));
		str_points[2].push_back(Point_2(0.0, 1.0));

		str_points[3].push_back(Point_2(0.0, 1.0));
		str_points[3].push_back(Point_2(0.0, 0.0));

		str_points[4].push_back(Point_2(2.0, 0.0));
		str_points[4].push_back(Point_2(3.0, 0.0));

		str_points[5].push_back(Point_2(3.0, 0.0));
		str_points[5].push_back(Point_2(3.0, 1.0));

		str_points[6].push_back(Point_2(3.0, 1.0));
		str_points[6].push_back(Point_2(2.0, 1.0));

		str_points[7].push_back(Point_2(2.0, 1.0));
		str_points[7].push_back(Point_2(2.0, 0.0));

		// Log log;
		// Container tmp;
		// log.save_visibility_eps(cdt, visibility, tmp, str_points);
	}
};

TEST_F(LOD_ReconstructionTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_ReconstructionTest, RunsProgram) {

	lod_0.max_flow(cdt, visibility);
}