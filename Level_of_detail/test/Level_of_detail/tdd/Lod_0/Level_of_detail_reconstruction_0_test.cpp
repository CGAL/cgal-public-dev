// Google test includes.
#include "gmock/gmock.h"

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Lod_0/Level_of_detail_reconstruction_0.h>

using namespace testing;

class LOD_ReconstructionTest: public Test {

public:
	using FT = double;
	
	using Traits    = CGAL::Simple_cartesian<FT>;
	using Point_2   = Traits::Point_2;

	using Label     = int; 
	using Container = std::vector< std::pair<Point_2, Label> >;

	using Str_label = CGAL::LOD::Structured_label;

	using My_vertex_info = CGAL::LOD::My_vertex_info<Str_label>; 
	using My_face_info   = CGAL::LOD::My_face_info<FT>;

	using VB 		   = CGAL::Triangulation_vertex_base_with_info_2<My_vertex_info, Traits>;
	using FB_with_info = CGAL::Triangulation_face_base_with_info_2<My_face_info, Traits>;
	using FB 		   = CGAL::Constrained_triangulation_face_base_2<Traits, FB_with_info>;

	using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
	using CDT = CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS>;

	using Lod_0 = CGAL::LOD::Level_of_detail_reconstruction_0<Traits, CDT>;

	using Vertex_handle = CDT::Vertex_handle;
	using Face_iterator = CDT::Finite_faces_iterator;

	using Log = CGAL::LOD::Mylog;

	Lod_0 lod_0;
	CDT cdt;

	LOD_ReconstructionTest() {
		create_data();
	}

	void create_data() {

		cdt.clear();
		set_basic_input();
	}

	void set_basic_input() {

		// First building.
		Vertex_handle va = cdt.insert(Point_2(0.0, 0.0)); va->info().label = Str_label::CORNER;
		Vertex_handle vb = cdt.insert(Point_2(0.5, 0.0)); vb->info().label = Str_label::LINEAR;
		Vertex_handle vc = cdt.insert(Point_2(1.0, 0.0)); vc->info().label = Str_label::CORNER;
		Vertex_handle vd = cdt.insert(Point_2(1.0, 0.5)); vd->info().label = Str_label::LINEAR;
		Vertex_handle ve = cdt.insert(Point_2(1.0, 1.0)); ve->info().label = Str_label::CORNER;
		Vertex_handle vf = cdt.insert(Point_2(0.5, 1.0)); vf->info().label = Str_label::LINEAR;
		Vertex_handle vg = cdt.insert(Point_2(0.0, 1.0)); vg->info().label = Str_label::CORNER;
		Vertex_handle vh = cdt.insert(Point_2(0.0, 0.5)); vh->info().label = Str_label::LINEAR;

		cdt.insert_constraint(va, vb);
		cdt.insert_constraint(vb, vc);
		cdt.insert_constraint(vc, vd);
		cdt.insert_constraint(vd, ve);
		cdt.insert_constraint(ve, vf);
		cdt.insert_constraint(vf, vg);
		cdt.insert_constraint(vg, vh);
		cdt.insert_constraint(vh, va);

		// Second building.
		va = cdt.insert(Point_2(2.0, 0.0)); va->info().label = Str_label::CORNER;
		vb = cdt.insert(Point_2(2.5, 0.0)); vb->info().label = Str_label::LINEAR;
		vc = cdt.insert(Point_2(3.0, 0.0)); vc->info().label = Str_label::CORNER;
		vd = cdt.insert(Point_2(3.0, 0.5)); vd->info().label = Str_label::LINEAR;
		ve = cdt.insert(Point_2(3.0, 1.0)); ve->info().label = Str_label::CORNER;
		vf = cdt.insert(Point_2(2.5, 1.0)); vf->info().label = Str_label::LINEAR;
		vg = cdt.insert(Point_2(2.0, 1.0)); vg->info().label = Str_label::CORNER;
		vh = cdt.insert(Point_2(2.0, 0.5)); vh->info().label = Str_label::LINEAR;

		cdt.insert_constraint(va, vb);
		cdt.insert_constraint(vb, vc);
		cdt.insert_constraint(vc, vd);
		cdt.insert_constraint(vd, ve);
		cdt.insert_constraint(ve, vf);
		cdt.insert_constraint(vf, vg);
		cdt.insert_constraint(vg, vh);
		cdt.insert_constraint(vh, va);

		Face_iterator fh = cdt.finite_faces_begin();
		fh->info().in = 1.0; ++fh; // 2 3 4
		fh->info().in = 1.0; ++fh; // 6 2 4
		fh->info().in = 1.0; ++fh; // 5 6 4
		fh->info().in = 1.0; ++fh; // 6 7 8
		fh->info().in = 0.0; ++fh; // 9 16 4
		fh->info().in = 1.0; ++fh; // 8 1 2
		fh->info().in = 1.0; ++fh; // 6 8 2
		fh->info().in = 0.0; ++fh; // 3 9 4
		fh->info().in = 1.0; ++fh; // 10 12 14
		fh->info().in = 1.0; ++fh; // 11 12 10
		fh->info().in = 1.0; ++fh; // 12 13 14
		fh->info().in = 0.0; ++fh; // 16 5 4
		fh->info().in = 1.0; ++fh; // 16 14 15
		fh->info().in = 0.0; ++fh; // 16 15 5
		fh->info().in = 1.0; ++fh; // 10 16 9
		fh->info().in = 1.0; ++fh; // 10 14 16
	}
};

TEST_F(LOD_ReconstructionTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_ReconstructionTest, RunsProgram) {

	lod_0.max_flow(cdt);
}