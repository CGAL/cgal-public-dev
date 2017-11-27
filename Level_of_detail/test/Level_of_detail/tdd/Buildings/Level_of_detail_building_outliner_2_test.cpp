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
#include <CGAL/Buildings/Level_of_detail_building_outliner_2.h>

using namespace testing;

class LOD_BuildingOutlinerTest: public Test {

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

	using LodBuildingOutliner = CGAL::LOD::Level_of_detail_building_outliner_2<Traits, CDT>;

	using Vertex_handle = CDT::Vertex_handle;
	using Face_handle   = CDT::Face_handle;
	using Face_iterator = CDT::Finite_faces_iterator;

	using Log = CGAL::LOD::Mylog;

	using Building  = CGAL::LOD::Building<FT, Vertex_handle, Face_handle>;
	using Buildings = std::map<int, Building>;

	using Boundaries = std::vector< std::vector<Vertex_handle> >;

	CDT cdt; Boundaries bound; Buildings buildings;
	LodBuildingOutliner lodBuildingOutliner;

	LOD_BuildingOutlinerTest() { 
		create_data();
	}

	void create_data() {

		cdt.clear();
		bound.clear();
		buildings.clear();

		lodBuildingOutliner.save_info(true);
		
		lodBuildingOutliner.set_max_outer_iterations(1000000);
		lodBuildingOutliner.set_max_inner_iterations(1000);

		set_basic_input();
	}

	void set_basic_input() {

		bound.resize(1);
		bound[0].resize(4);


		// First building.
		const Vertex_handle va1 = cdt.insert(Point_2(0.0, 0.0)); va1->info().label = Str_label::CORNER; bound[0][0] = va1;
		const Vertex_handle vb1 = cdt.insert(Point_2(1.0, 0.0)); vb1->info().label = Str_label::CORNER; bound[0][1] = vb1;
		const Vertex_handle vc1 = cdt.insert(Point_2(1.0, 1.0)); vc1->info().label = Str_label::CORNER; bound[0][2] = vc1;
		const Vertex_handle vd1 = cdt.insert(Point_2(0.0, 1.0)); vd1->info().label = Str_label::CORNER; bound[0][3] = vd1;

		cdt.insert_constraint(va1, vb1); 
		cdt.insert_constraint(vb1, vc1);
		cdt.insert_constraint(vc1, vd1);
		cdt.insert_constraint(vd1, va1);


		// Second building - exterior.
		const Vertex_handle ve2 = cdt.insert(Point_2(1.6, 0.0)); ve2->info().label = Str_label::CORNER;
		const Vertex_handle vf2 = cdt.insert(Point_2(2.2, 0.0)); vf2->info().label = Str_label::LINEAR;
		const Vertex_handle vg2 = cdt.insert(Point_2(2.8, 0.0)); vg2->info().label = Str_label::CORNER;
		const Vertex_handle vh2 = cdt.insert(Point_2(2.8, 0.2)); vh2->info().label = Str_label::CLUTTER;
		const Vertex_handle vi2 = cdt.insert(Point_2(2.8, 0.6)); vi2->info().label = Str_label::CORNER;
		const Vertex_handle vj2 = cdt.insert(Point_2(2.4, 0.6)); vj2->info().label = Str_label::CORNER;
		const Vertex_handle vk2 = cdt.insert(Point_2(2.4, 1.0)); vk2->info().label = Str_label::CORNER;
		const Vertex_handle vl2 = cdt.insert(Point_2(2.2, 1.0)); vl2->info().label = Str_label::LINEAR;
		const Vertex_handle vm2 = cdt.insert(Point_2(2.0, 1.0)); vm2->info().label = Str_label::CORNER;
		const Vertex_handle vn2 = cdt.insert(Point_2(2.0, 0.6)); vn2->info().label = Str_label::CORNER;
		const Vertex_handle vo2 = cdt.insert(Point_2(1.6, 0.6)); vo2->info().label = Str_label::CORNER;
		const Vertex_handle vp2 = cdt.insert(Point_2(1.6, 0.4)); vp2->info().label = Str_label::LINEAR;

		// cdt.insert_constraint(vg2, vh2); // unconstrained edge
		// cdt.insert_constraint(vh2, vi2); // unconstrained edge
		// cdt.insert_constraint(vj2, vk2); // unconstrained edge

		// constrained edges
		cdt.insert_constraint(ve2, vf2);
		cdt.insert_constraint(vf2, vg2);
		cdt.insert_constraint(vi2, vj2);
		cdt.insert_constraint(vk2, vl2);
		cdt.insert_constraint(vl2, vm2);
		cdt.insert_constraint(vm2, vn2);
		cdt.insert_constraint(vn2, vo2);
		cdt.insert_constraint(vo2, vp2);
		cdt.insert_constraint(vp2, ve2);


		// Second building - interior.
		const Vertex_handle vq3 = cdt.insert(Point_2(1.8, 0.2)); vq3->info().label = Str_label::CORNER;
		const Vertex_handle vr3 = cdt.insert(Point_2(2.6, 0.2)); vr3->info().label = Str_label::CORNER;
		const Vertex_handle vs3 = cdt.insert(Point_2(2.6, 0.4)); vs3->info().label = Str_label::CORNER;
		const Vertex_handle vt3 = cdt.insert(Point_2(2.2, 0.4)); vt3->info().label = Str_label::LINEAR;
		const Vertex_handle vu3 = cdt.insert(Point_2(1.8, 0.4)); vu3->info().label = Str_label::CORNER;

		cdt.insert_constraint(vq3, vr3);
		cdt.insert_constraint(vr3, vs3);
		cdt.insert_constraint(vs3, vt3);
		cdt.insert_constraint(vt3, vu3);
		cdt.insert_constraint(vu3, vq3);

		
		// Face and building information.
		const CGAL::Color r(169, 0, 0);
		const CGAL::Color g(0, 169, 0);
		const CGAL::Color b(0, 0, 169);

		const CGAL::Color grey(169, 169, 169);		

		Face_iterator fh = cdt.finite_faces_begin();

		fh->info().in = 1.0; fh->info().bu =  0; fh->info().bu_color = g; buildings[0].faces.push_back(static_cast<Face_handle>(fh)); buildings[0].color = g; ++fh; 
		fh->info().in = 1.0; fh->info().bu =  0; fh->info().bu_color = g; buildings[0].faces.push_back(static_cast<Face_handle>(fh)); buildings[0].color = g; ++fh;
		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;
		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;
		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;

		fh->info().in = 0.0; fh->info().bu = -1; fh->info().bu_color = grey; ++fh; 														 						   // outside

		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;
		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;
		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;

		fh->info().in = 0.0; fh->info().bu = -1; fh->info().bu_color = grey; ++fh; 														 						   // outside

		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;

		fh->info().in = 0.0; fh->info().bu = -1; fh->info().bu_color = grey; ++fh; 														 						   // outside

		fh->info().in = 1.0; fh->info().bu =  1; fh->info().bu_color = b; buildings[1].faces.push_back(static_cast<Face_handle>(fh)); buildings[1].color = b; ++fh;

		fh->info().in = 0.0; fh->info().bu = -1; fh->info().bu_color = grey; ++fh; 														 						   // outside

		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;

		fh->info().in = 0.0; fh->info().bu = -1; fh->info().bu_color = grey; ++fh; 														 						   // outside
		fh->info().in = 0.0; fh->info().bu = -1; fh->info().bu_color = grey; ++fh; 														 						   // outside

		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;
		fh->info().in = 1.0; fh->info().bu =  1; fh->info().bu_color = b; buildings[1].faces.push_back(static_cast<Face_handle>(fh)); buildings[1].color = b; ++fh;
		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;
		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;
		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;
		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;
		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;
		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;
		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;
		fh->info().in = 1.0; fh->info().bu =  1; fh->info().bu_color = b; buildings[1].faces.push_back(static_cast<Face_handle>(fh)); buildings[1].color = b; ++fh;
		fh->info().in = 1.0; fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; ++fh;

		// Log log; 
		// log.save_cdt_ply(cdt, "tmp" + std::string(PS) + "cdt_outliner", "bu");
		// log.save_buildings_info(cdt, buildings, "tmp" + std::string(PS) + "buildings_outliner_before");
	}
};

TEST_F(LOD_BuildingOutlinerTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_BuildingOutlinerTest, VerifiesOrientedMethod) {

	lodBuildingOutliner.set_boundary_type(CGAL::LOD::Building_boundary_type::ORIENTED);
	lodBuildingOutliner.find_boundaries(cdt, buildings);

	// Extra.	
	Log log;
	log.save_buildings_info(cdt, buildings, "tmp" + std::string(PS) + "buildings_oriented_outliner_after");

	ASSERT_THAT(static_cast<int>(buildings.size()), Eq(3));
	Boundaries &b = buildings[0].boundaries;

	// Main.
	assert(!b.empty());
	ASSERT_THAT(b[0].size(), Eq(bound[0].size()));

	ASSERT_THAT(b[0][0], Eq(bound[0][0]));
	ASSERT_THAT(b[0][1], Eq(bound[0][1]));
	ASSERT_THAT(b[0][2], Eq(bound[0][2]));
	ASSERT_THAT(b[0][3], Eq(bound[0][3]));
}

TEST_F(LOD_BuildingOutlinerTest, VerifiesUnorientedMethod) {

	lodBuildingOutliner.save_info(false);
	lodBuildingOutliner.set_boundary_type(CGAL::LOD::Building_boundary_type::UNORIENTED);

	lodBuildingOutliner.find_boundaries(cdt, buildings);
	
	// Extra.
	Log log;
	log.save_buildings_info(cdt, buildings, "tmp" + std::string(PS) + "buildings_unoriented_outliner_after");

	ASSERT_THAT(static_cast<int>(buildings.size()), Eq(3));
	Boundaries &b = buildings[0].boundaries;

	// Main.
	assert(!b.empty());
	ASSERT_THAT(static_cast<int>(b[0].size()), Eq(8));
}

TEST_F(LOD_BuildingOutlinerTest, VerifiesBuildingNeighbours) {

	lodBuildingOutliner.save_info(false);
	lodBuildingOutliner.find_boundaries(cdt, buildings);

	// Extra.
	ASSERT_THAT(static_cast<int>(buildings.size()), Eq(3));

	// Main.
	ASSERT_THAT(static_cast<int>(buildings[0].neighbours.size()), Eq(0));
	
	ASSERT_THAT(*(buildings[1].neighbours.begin()), Eq(2));
	ASSERT_THAT(*(buildings[2].neighbours.begin()), Eq(1));
}