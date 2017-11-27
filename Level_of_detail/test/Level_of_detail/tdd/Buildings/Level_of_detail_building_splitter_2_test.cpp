#if defined(WIN32) || defined(_WIN32) 
#define PS "\\" 
#else 
#define PS "/" 
#endif 

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

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Buildings/Level_of_detail_building_splitter_2.h>

using namespace testing;

class LOD_BuildingSplitterTest: public Test {

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

	using LodBuildingSplitter = CGAL::LOD::Level_of_detail_building_splitter_2<Traits, CDT>;

	using Vertex_handle = CDT::Vertex_handle;
	using Face_handle   = CDT::Face_handle;
	using Face_iterator = CDT::Finite_faces_iterator;

	using Log = CGAL::LOD::Mylog;

	using Building  = CGAL::LOD::Building<FT, Vertex_handle, Face_handle>;
	using Buildings = std::map<int, Building>;

	CDT cdt;
	LodBuildingSplitter lodBuildingSplitter;

	LOD_BuildingSplitterTest() { 
		create_data();
	}

	void create_data() {

		cdt.clear();
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
		const Vertex_handle vf2 = cdt.insert(Point_2(2.0, 0.0)); vf2->info().label = Str_label::CORNER;
		const Vertex_handle vg2 = cdt.insert(Point_2(2.0, 0.6)); vg2->info().label = Str_label::CORNER;
		const Vertex_handle ve2 = cdt.insert(Point_2(1.0, 0.6)); ve2->info().label = Str_label::LINEAR;

		cdt.insert_constraint(vb2, vf2);
		cdt.insert_constraint(vf2, vg2);
		cdt.insert_constraint(vg2, ve2);
		cdt.insert_constraint(ve2, vb2);

		// Third building - exterior.
		const Vertex_handle vh3 = cdt.insert(Point_2(2.6, 0.0)); vh3->info().label = Str_label::CORNER;
		const Vertex_handle vk3 = cdt.insert(Point_2(4.0, 0.0)); vk3->info().label = Str_label::CORNER;
		const Vertex_handle vj3 = cdt.insert(Point_2(4.0, 1.0)); vj3->info().label = Str_label::CORNER;
		const Vertex_handle vi3 = cdt.insert(Point_2(2.6, 1.0)); vi3->info().label = Str_label::CORNER;

		cdt.insert_constraint(vh3, vk3);
		cdt.insert_constraint(vk3, vj3);
		cdt.insert_constraint(vj3, vi3);
		cdt.insert_constraint(vi3, vh3);

		// Third building - interior.
		const Vertex_handle vm3 = cdt.insert(Point_2(3.0, 0.2)); vm3->info().label = Str_label::CORNER;
		const Vertex_handle vn3 = cdt.insert(Point_2(3.8, 0.2)); vn3->info().label = Str_label::CORNER;
		const Vertex_handle vo3 = cdt.insert(Point_2(3.8, 0.8)); vo3->info().label = Str_label::CORNER;
		const Vertex_handle vl3 = cdt.insert(Point_2(3.0, 0.8)); vl3->info().label = Str_label::CORNER;

		cdt.insert_constraint(vm3, vn3);
		cdt.insert_constraint(vn3, vo3);
		cdt.insert_constraint(vo3, vl3);
		cdt.insert_constraint(vl3, vm3);

		Face_iterator fh = cdt.finite_faces_begin();
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 0.0; fh->info().in_color = CGAL::Color(255, 51, 51); ++fh; // 	OUTSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 0.0; fh->info().in_color = CGAL::Color(255, 51, 51); ++fh; // 	OUTSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 0.0; fh->info().in_color = CGAL::Color(255, 51, 51); ++fh; // 	OUTSIDE
		fh->info().in = 0.0; fh->info().in_color = CGAL::Color(255, 51, 51); ++fh; // 	OUTSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE
		fh->info().in = 1.0; fh->info().in_color = CGAL::Color(51, 255, 51); ++fh; // INSIDE

		// for (typename CDT::Finite_edges_iterator eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit)
		//	if (cdt.is_constrained(*eit)) std::cout << "constrained edge" << std::endl;

		// Log log; log.save_cdt_ply(cdt, "tmp" + std::string(PS) + "cdt_splitter");
	}
};

TEST_F(LOD_BuildingSplitterTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_BuildingSplitterTest, ReturnsFourBuildings) {
   
   	Buildings buildings;
	const int number_of_buildings = lodBuildingSplitter.split(cdt, buildings);

	ASSERT_THAT(number_of_buildings, Eq(4));
	ASSERT_THAT(static_cast<int>(buildings.size()), Eq(number_of_buildings));
}