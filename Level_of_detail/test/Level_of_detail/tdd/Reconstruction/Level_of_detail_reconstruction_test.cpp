// Google test includes.
#include "gmock/gmock.h"

// STL includes.
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Polyhedron_3.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Reconstruction/Level_of_detail_reconstruction.h>

using namespace testing;

class LOD_ReconstructionTest: public Test {

public:
	using FT = double;
	
	using Traits  = CGAL::Simple_cartesian<FT>;
	using Point_2 = Traits::Point_2;

	using Str_label = CGAL::LOD::Structured_label;

	using My_vertex_info = CGAL::LOD::My_vertex_info<Str_label>; 
	using My_face_info   = CGAL::LOD::My_face_info<FT>;

	using VB 		   = CGAL::Triangulation_vertex_base_with_info_2<My_vertex_info, Traits>;
	using FB_with_info = CGAL::Triangulation_face_base_with_info_2<My_face_info, Traits>;
	using FB 		   = CGAL::Constrained_triangulation_face_base_2<Traits, FB_with_info>;

	using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
	using CDT = CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS>;

	using Vertex_handle = CDT::Vertex_handle;
	using Face_handle   = CDT::Face_handle;
	using Face_iterator = CDT::Finite_faces_iterator;

	using Building  = CGAL::LOD::Building<FT, Vertex_handle, Face_handle>;
	using Buildings = std::map<int, Building>;

	using Mesh = CGAL::Polyhedron_3<Traits>;

	using Lods = CGAL::LOD::Level_of_detail_reconstruction<Traits, CDT, Buildings, Mesh>;
	using Mesh_facet_colors = Lods::Mesh_facet_colors;

	using Log = CGAL::LOD::Mylog;

	using Ground_point = Lods::Point;
	using Ground = Lods::Ground;

	CDT cdt; Buildings buildings; Ground ground; Mesh mesh; Mesh_facet_colors mesh_facet_colors;
	Lods lods;

	LOD_ReconstructionTest() { 
		create_data();
	}

	void create_data() {

		cdt.clear();
		buildings.clear();
		ground.clear();

		set_basic_input();
	}

	void set_basic_input() {

		// First building.
		const Vertex_handle va1 = cdt.insert(Point_2(0.0, 0.0)); va1->info().label = Str_label::CORNER;
		const Vertex_handle vb1 = cdt.insert(Point_2(1.0, 0.0)); vb1->info().label = Str_label::CORNER;
		const Vertex_handle vc1 = cdt.insert(Point_2(1.0, 1.0)); vc1->info().label = Str_label::CORNER;
		const Vertex_handle vd1 = cdt.insert(Point_2(0.0, 1.0)); vd1->info().label = Str_label::CORNER;

		cdt.insert_constraint(va1, vb1); 
		cdt.insert_constraint(vb1, vc1);
		cdt.insert_constraint(vc1, vd1);
		cdt.insert_constraint(vd1, va1);


		// Second building.
		const Vertex_handle vb2 = cdt.insert(Point_2(1.0, 0.0)); vb2->info().label = Str_label::CORNER;
		const Vertex_handle ve2 = cdt.insert(Point_2(1.5, 0.0)); ve2->info().label = Str_label::CORNER;
		const Vertex_handle vf2 = cdt.insert(Point_2(1.5, 0.5)); vf2->info().label = Str_label::CORNER;
		const Vertex_handle vg2 = cdt.insert(Point_2(1.0, 0.5)); vg2->info().label = Str_label::CORNER;

		cdt.insert_constraint(vb2, ve2);
		cdt.insert_constraint(ve2, vf2);
		cdt.insert_constraint(vf2, vg2);
		cdt.insert_constraint(vg2, vb2);


		// Third building.
		const Vertex_handle vf3 = cdt.insert(Point_2(1.5, 0.5)); vf3->info().label = Str_label::CORNER;
		const Vertex_handle vh3 = cdt.insert(Point_2(2.0, 0.5)); vh3->info().label = Str_label::CORNER;
		const Vertex_handle vi3 = cdt.insert(Point_2(2.0, 1.0)); vi3->info().label = Str_label::CORNER;
		const Vertex_handle vj3 = cdt.insert(Point_2(1.5, 1.0)); vj3->info().label = Str_label::CORNER;

		cdt.insert_constraint(vf3, vh3);
		cdt.insert_constraint(vh3, vi3);
		cdt.insert_constraint(vi3, vj3);
		cdt.insert_constraint(vj3, vf3);


		// Fourth building - exterior.
		const Vertex_handle vk4 = cdt.insert(Point_2(2.5, 0.0)); vk4->info().label = Str_label::CORNER;
		const Vertex_handle vl4 = cdt.insert(Point_2(3.5, 0.0)); vl4->info().label = Str_label::CORNER;
		const Vertex_handle vm4 = cdt.insert(Point_2(3.5, 1.0)); vm4->info().label = Str_label::CORNER;
		const Vertex_handle vn4 = cdt.insert(Point_2(2.5, 1.0)); vn4->info().label = Str_label::CORNER;

		cdt.insert_constraint(vk4, vl4);
		cdt.insert_constraint(vl4, vm4);
		cdt.insert_constraint(vm4, vn4);
		cdt.insert_constraint(vn4, vk4);


		// Fourth building - interior.
		const Vertex_handle vo4 = cdt.insert(Point_2(2.8, 0.2)); vo4->info().label = Str_label::CORNER;
		const Vertex_handle vp4 = cdt.insert(Point_2(3.2, 0.2)); vp4->info().label = Str_label::CORNER;
		const Vertex_handle vq4 = cdt.insert(Point_2(3.2, 0.8)); vq4->info().label = Str_label::CORNER;
		const Vertex_handle vr4 = cdt.insert(Point_2(2.8, 0.8)); vr4->info().label = Str_label::CORNER;

		cdt.insert_constraint(vo4, vp4);
		cdt.insert_constraint(vp4, vq4);
		cdt.insert_constraint(vq4, vr4);
		cdt.insert_constraint(vr4, vo4);


		// Fifth building - with the hole inside.
		const Vertex_handle vs5 = cdt.insert(Point_2(4.0, 0.0)); vs5->info().label = Str_label::CORNER;
		const Vertex_handle vt5 = cdt.insert(Point_2(5.0, 0.0)); vt5->info().label = Str_label::CORNER;
		const Vertex_handle vu5 = cdt.insert(Point_2(5.0, 1.0)); vu5->info().label = Str_label::CORNER;
		const Vertex_handle vv5 = cdt.insert(Point_2(4.0, 1.0)); vv5->info().label = Str_label::CORNER;

		const Vertex_handle vw5 = cdt.insert(Point_2(4.3, 0.3)); vw5->info().label = Str_label::CORNER;
		const Vertex_handle vz5 = cdt.insert(Point_2(4.7, 0.3)); vz5->info().label = Str_label::CORNER;
		const Vertex_handle va5 = cdt.insert(Point_2(4.7, 0.7)); va5->info().label = Str_label::CORNER;
		const Vertex_handle vb5 = cdt.insert(Point_2(4.3, 0.7)); vb5->info().label = Str_label::CORNER;

		cdt.insert_constraint(vs5, vt5);
		cdt.insert_constraint(vt5, vu5);
		cdt.insert_constraint(vu5, vv5);
		cdt.insert_constraint(vv5, vs5);

		cdt.insert_constraint(vw5, vz5);
		cdt.insert_constraint(vz5, va5);
		cdt.insert_constraint(va5, vb5);
		cdt.insert_constraint(vb5, vw5);


		// Face and building information.
		const CGAL::Color    r(169, 0  ,   0);
		const CGAL::Color    g(0  , 169,   0);
		const CGAL::Color    b(0  , 0  , 169);
		const CGAL::Color    f(169, 0  , 169);
		const CGAL::Color    o(255, 169,   0);
		const CGAL::Color    p(255, 20 , 169);
		const CGAL::Color grey(169, 169, 169);	

		Face_iterator fh = cdt.finite_faces_begin();

		fh->info().bu =  0; fh->info().bu_color = g; buildings[0].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  0; fh->info().bu_color = g; buildings[0].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  1; fh->info().bu_color = r; buildings[1].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu = -1; fh->info().bu_color = grey; ++fh;														   // outside
		fh->info().bu =  0; fh->info().bu_color = g; buildings[0].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  1; fh->info().bu_color = r; buildings[1].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  2; fh->info().bu_color = b; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu = -1; fh->info().bu_color = grey; ++fh;														   // outside
		fh->info().bu = -1; fh->info().bu_color = grey; ++fh;														   // outside
		fh->info().bu = -1; fh->info().bu_color = grey; ++fh;														   // outside
		fh->info().bu =  2; fh->info().bu_color = b; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu = -1; fh->info().bu_color = grey; ++fh;														   // outside
		fh->info().bu =  4; fh->info().bu_color = o; buildings[4].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu = -1; fh->info().bu_color = grey; ++fh;														   // outside
		fh->info().bu =  3; fh->info().bu_color = f; buildings[3].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  3; fh->info().bu_color = f; buildings[3].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  3; fh->info().bu_color = f; buildings[3].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  3; fh->info().bu_color = f; buildings[3].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  3; fh->info().bu_color = f; buildings[3].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  3; fh->info().bu_color = f; buildings[3].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  4; fh->info().bu_color = o; buildings[4].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  3; fh->info().bu_color = f; buildings[3].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  3; fh->info().bu_color = f; buildings[3].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu = -1; fh->info().bu_color = grey; ++fh;														   // outside
		fh->info().bu = -1; fh->info().bu_color = grey; ++fh;														   // outside
		fh->info().bu = -1; fh->info().bu_color = grey; ++fh;														   // outside
		fh->info().bu =  5; fh->info().bu_color = p; buildings[5].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  5; fh->info().bu_color = p; buildings[5].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  5; fh->info().bu_color = p; buildings[5].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  5; fh->info().bu_color = p; buildings[5].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  5; fh->info().bu_color = p; buildings[5].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  5; fh->info().bu_color = p; buildings[5].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu = -1; fh->info().bu_color = grey; ++fh;														   // outside
		fh->info().bu =  5; fh->info().bu_color = p; buildings[5].faces.push_back(static_cast<Face_handle>(fh)); ++fh;
		fh->info().bu =  5; fh->info().bu_color = p; buildings[5].faces.push_back(static_cast<Face_handle>(fh)); ++fh;


		// Boundaries and other information.
		buildings[0].color = g; buildings[0].height = 1.0; buildings[0].is_oriented = true;
		buildings[1].color = r; buildings[1].height = 0.5; buildings[1].is_oriented = true;
		buildings[2].color = b; buildings[2].height = 0.8; buildings[2].is_oriented = true;
		buildings[3].color = f; buildings[3].height = 1.0; buildings[3].is_oriented = true;
		buildings[4].color = o; buildings[4].height = 1.5; buildings[4].is_oriented = true;
		buildings[5].color = p; buildings[5].height = 0.7; buildings[5].is_oriented = true;

		// First building.
		buildings[0].boundaries.resize(1);
		buildings[0].boundaries[0].resize(4);
		buildings[0].boundaries[0][0] = va1; buildings[0].boundaries[0][1] = vb1; buildings[0].boundaries[0][2] = vc1; buildings[0].boundaries[0][3] = vd1;

		// Second building.
		buildings[1].boundaries.resize(1);
		buildings[1].boundaries[0].resize(4);
		buildings[1].boundaries[0][0] = vb2; buildings[1].boundaries[0][1] = ve2; buildings[1].boundaries[0][2] = vf2; buildings[1].boundaries[0][3] = vg2;

		// Third building.
		buildings[2].boundaries.resize(1);
		buildings[2].boundaries[0].resize(4);
		buildings[2].boundaries[0][0] = vf3; buildings[2].boundaries[0][1] = vh3; buildings[2].boundaries[0][2] = vi3; buildings[2].boundaries[0][3] = vj3;

		// Fourth building - exterior.
		buildings[3].boundaries.resize(1);
		buildings[3].boundaries[0].resize(4);
		buildings[3].boundaries[0][0] = vk4; buildings[3].boundaries[0][1] = vl4; buildings[3].boundaries[0][2] = vm4; buildings[3].boundaries[0][3] = vn4;

		// - interior.
		buildings[4].boundaries.resize(1);
		buildings[4].boundaries[0].resize(4);
		buildings[4].boundaries[0][0] = vo4; buildings[4].boundaries[0][1] = vp4; buildings[4].boundaries[0][2] = vq4; buildings[4].boundaries[0][3] = vr4;

		// Fifth building.
		buildings[5].boundaries.resize(2);

		// - exterior boundary.
		buildings[5].boundaries[0].resize(4);
		buildings[5].boundaries[0][0] = vs5; buildings[5].boundaries[0][1] = vt5; buildings[5].boundaries[0][2] = vu5; buildings[5].boundaries[0][3] = vv5;		

		// - interior boundary.
		buildings[5].boundaries[1].resize(4);
		buildings[5].boundaries[1][0] = vw5; buildings[5].boundaries[1][1] = vz5; buildings[5].boundaries[1][2] = va5; buildings[5].boundaries[1][3] = vb5;


		// Set ground.
		ground.resize(4);

		ground[0] = Ground_point(0.0, 0.0, 0.0);
		ground[1] = Ground_point(5.0, 0.0, 0.0);
		ground[2] = Ground_point(5.0, 1.0, 0.0);
		ground[3] = Ground_point(0.0, 1.0, 0.0);

		lods.use_boundaries(true);

		// Log log; 
		// log.save_cdt_ply(cdt, "tmp/cdt_lod1", "bu");
		// log.save_buildings_info(cdt, buildings, "tmp/buildings_lod1");
	}
};

TEST_F(LOD_ReconstructionTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_ReconstructionTest, RunsReconstructionLOD0) {

	mesh.clear();
	mesh_facet_colors.clear();

	lods.reconstruct_lod0(cdt, buildings, ground, mesh, mesh_facet_colors);

	Log log;
	log.save_mesh_as_ply(mesh, mesh_facet_colors, "LOD0");
}

TEST_F(LOD_ReconstructionTest, RunsReconstructionLOD1) {

	mesh.clear();
	mesh_facet_colors.clear();

	lods.reconstruct_lod1(cdt, buildings, ground, mesh, mesh_facet_colors);

	Log log;
	log.save_mesh_as_ply(mesh, mesh_facet_colors, "LOD1");
}