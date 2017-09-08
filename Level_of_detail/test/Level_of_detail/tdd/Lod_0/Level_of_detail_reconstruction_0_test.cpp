// Google test includes.
#include "gmock/gmock.h"

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>

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

	using VB = CGAL::Triangulation_vertex_base_2<Traits>;
	using FB = CGAL::Constrained_triangulation_face_base_2<Traits>;

	using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
	using CDT = CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS>;

	using LodVisibility = CGAL::LOD::Level_of_detail_visibility_from_classification_2<Traits, Container, CDT>;
	using Visibility_result = std::map<int, LodVisibility::Visibility_label>;

	using Structured_label  = typename CGAL::LOD::Level_of_detail_structuring_2<Traits>::Structured_label;
	using Structured_labels = std::vector<std::vector<Structured_label> >;

	using Lod_0 = CGAL::LOD::Level_of_detail_reconstruction_0<Traits, CDT, Visibility_result, Structured_labels>;
	using Vertex_handle = CDT::Vertex_handle;

	using Log = CGAL::LOD::Mylog;

	using Segment = Traits::Segment_2;
	using Lod_0_result = std::vector<Segment>;

	Lod_0 lod_0;
	CDT cdt; Visibility_result visibility; Structured_labels str_labels;

	LOD_ReconstructionTest() {
		create_data();
	}

	void create_data() {

		cdt.clear();
		visibility.clear();
		str_labels.clear();

		set_basic_input();
	}

	void set_basic_input() {

		Vertex_handle va = cdt.insert(Point_2(0, 0));
		Vertex_handle vb = cdt.insert(Point_2(1, 0));
		Vertex_handle vc = cdt.insert(Point_2(0, 1));
		Vertex_handle vd = cdt.insert(Point_2(1, 1));
		Vertex_handle ve = cdt.insert(Point_2(2, 0));
		Vertex_handle vf = cdt.insert(Point_2(2, 1));
		Vertex_handle vg = cdt.insert(Point_2(3, 0));
		Vertex_handle vk = cdt.insert(Point_2(3, 1));

		// Bounding box.
		cdt.insert_constraint(va, vg);
		cdt.insert_constraint(vg, vk);
		cdt.insert_constraint(vk, vc);
		cdt.insert_constraint(vc, va);

		// Triangle faces of my mesh.
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

		// Visibility predictions.
		visibility[0] = LodVisibility::Visibility_label::IN;
		visibility[1] = LodVisibility::Visibility_label::IN;

		visibility[2] = LodVisibility::Visibility_label::OUT;
		visibility[3] = LodVisibility::Visibility_label::OUT;

		visibility[4] = LodVisibility::Visibility_label::UNKNOWN;

		// Vertex labels.
		// to be added
	}
};

TEST_F(LOD_ReconstructionTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_ReconstructionTest, RunsProgram) {

	Lod_0_result lod_0_result;
	lod_0.reconstruct(cdt, visibility, str_labels, lod_0_result);
}