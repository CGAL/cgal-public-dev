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
#include <CGAL/Point_set_3.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Buildings/Level_of_detail_building_roof_fitter_2.h>

using namespace testing;

class LOD_BuildingRoofFitterTest: public Test {

public:
	using FT 	 = double;
	using Traits = CGAL::Simple_cartesian<FT>;
	
	using Point_2 = Traits::Point_2;
	using Point_3 = Traits::Point_3;
	using Plane_3 = Traits::Plane_3;

	using Container = CGAL::Point_set_3<Point_3>;
	using Str_label = CGAL::LOD::Structured_label;

	using My_vertex_info = CGAL::LOD::My_vertex_info<Str_label>; 
	using My_face_info   = CGAL::LOD::My_face_info<FT>;

	using VB 		   = CGAL::Triangulation_vertex_base_with_info_2<My_vertex_info, Traits>;
	using FB_with_info = CGAL::Triangulation_face_base_with_info_2<My_face_info, Traits>;
	using FB 		   = CGAL::Constrained_triangulation_face_base_2<Traits, FB_with_info>;

	using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
	using CDT = CGAL::Constrained_Delaunay_triangulation_2<Traits, TDS>;

	using LodMinHeightFitter = CGAL::LOD::Level_of_detail_min_height_fitter<Traits>;
	using LodAvgHeightFitter = CGAL::LOD::Level_of_detail_avg_height_fitter<Traits>;
	using LodMaxHeightFitter = CGAL::LOD::Level_of_detail_max_height_fitter<Traits>;
	
	using LodBuildingRoofFitterMin = CGAL::LOD::Level_of_detail_building_roof_fitter_2<Traits, CDT, Container, LodMinHeightFitter>;
	using LodBuildingRoofFitterAvg = CGAL::LOD::Level_of_detail_building_roof_fitter_2<Traits, CDT, Container, LodAvgHeightFitter>;
	using LodBuildingRoofFitterMax = CGAL::LOD::Level_of_detail_building_roof_fitter_2<Traits, CDT, Container, LodMaxHeightFitter>;

	using Vertex_handle = CDT::Vertex_handle;
	using Face_handle   = CDT::Face_handle;
	using Face_iterator = CDT::Finite_faces_iterator;

	using Log = CGAL::LOD::Mylog;

	using Building  = CGAL::LOD::Building<FT, Vertex_handle, Face_handle>;
	using Buildings = std::map<int, Building>;

	using Label     = int; 
	using Label_map = Container:: template Property_map<Label>; 

	using Point_iterator = Container::iterator;

	using Point_index = typename Container::Index;
	using Face_points_map = std::map<Face_handle, std::vector<Point_index> >;

	CDT cdt; Buildings buildings; Container input; Face_points_map fp_map; Plane_3 ground_fitted_plane;
	
	LodBuildingRoofFitterMin lodBuildingRoofFitterMin;
	LodBuildingRoofFitterAvg lodBuildingRoofFitterAvg;
	LodBuildingRoofFitterMax lodBuildingRoofFitterMax;

	LOD_BuildingRoofFitterTest() { 
		create_data();
	}

	void create_data() {

		cdt.clear();
		buildings.clear();
		input.clear();
		fp_map.clear();

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

		
		// Input container with 3D points + see below.
		const Label ground = 0;
		const Label facade = 1;
		const Label roof   = 2;

		Label_map labels;
		set_labels_property(input, labels);


		// Face, input, and building information.
		const CGAL::Color r(169, 0, 0);
		const CGAL::Color g(0, 169, 0);
		const CGAL::Color b(0, 0, 169);

		const CGAL::Color grey(169, 169, 169);		
		Face_iterator fh = cdt.finite_faces_begin();

		// Face 0.
		fh->info().bu =  0; fh->info().bu_color = g; buildings[0].faces.push_back(static_cast<Face_handle>(fh)); buildings[0].color = g; 
		Point_iterator 
		pit = input.insert(Point_3(0.80, 0.80, 1.0)); labels[*pit] =   roof; fp_map[static_cast<Face_handle>(fh)].push_back(*pit); 
		pit = input.insert(Point_3(0.60, 0.60, 9.0)); labels[*pit] = ground; fp_map[static_cast<Face_handle>(fh)].push_back(*pit); 			   												   // outlier
		++fh;

		// Face 1.
		fh->info().bu =  0; fh->info().bu_color = g; buildings[0].faces.push_back(static_cast<Face_handle>(fh)); buildings[0].color = g; 
		pit = input.insert(Point_3(0.20, 0.40, 1.0)); labels[*pit] = roof; fp_map[static_cast<Face_handle>(fh)].push_back(*pit);
		++fh;

		// Face 2.
		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; 
		pit = input.insert(Point_3(1.93, 0.07, 1.0)); labels[*pit] = roof; fp_map[static_cast<Face_handle>(fh)].push_back(*pit);
		++fh;

		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh;
		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh;
		fh->info().bu = -1; fh->info().bu_color = grey; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh; 														 						   // outside
		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh;

		// Face 7.
		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; 
		pit = input.insert(Point_3(2.74, 0.33, 1.1)); labels[*pit] = roof; fp_map[static_cast<Face_handle>(fh)].push_back(*pit);
		++fh;

		// Face 8.
		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; 
		pit = input.insert(Point_3(2.35, 0.85, 1.0)); labels[*pit] = roof; fp_map[static_cast<Face_handle>(fh)].push_back(*pit);
		++fh;

		fh->info().bu = -1; fh->info().bu_color = grey; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh; 														 						   // outside
		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh;
		fh->info().bu = -1; fh->info().bu_color = grey; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh; 														 						   // outside

		// Face 12.
		fh->info().bu =  1; fh->info().bu_color = b; buildings[1].faces.push_back(static_cast<Face_handle>(fh)); buildings[1].color = b; 
		pit = input.insert(Point_3(2.00, 0.36, 1.9)); labels[*pit] = roof; fp_map[static_cast<Face_handle>(fh)].push_back(*pit);
		++fh;

		fh->info().bu = -1; fh->info().bu_color = grey; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh; 														 						   // outside
		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh;
		fh->info().bu = -1; fh->info().bu_color = grey; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh; 														 						   // outside
		fh->info().bu = -1; fh->info().bu_color = grey; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh; 														 						   // outside
		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh;

		// Face 18.
		fh->info().bu =  1; fh->info().bu_color = b; buildings[1].faces.push_back(static_cast<Face_handle>(fh)); buildings[1].color = b; 
		pit = input.insert(Point_3(2.40, 0.36, 2.1)); labels[*pit] = roof; fp_map[static_cast<Face_handle>(fh)].push_back(*pit);
		++fh;

		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh;
		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh;		
		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh;

		// Face 22.
		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; 
		pit = input.insert(Point_3(2.13, 0.55, 0.9)); labels[*pit] = roof; fp_map[static_cast<Face_handle>(fh)].push_back(*pit);
		++fh;

		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh;

		// Face 24.
		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; 
		pit = input.insert(Point_3(1.80, 0.60, 1.0)); labels[*pit] = facade; fp_map[static_cast<Face_handle>(fh)].push_back(*pit); 			    											    // outlier
		++fh;

		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh;

		// Face 26.
		fh->info().bu =  1; fh->info().bu_color = b; buildings[1].faces.push_back(static_cast<Face_handle>(fh)); buildings[1].color = b; 
		pit = input.insert(Point_3(2.12, 0.29, 2.0)); labels[*pit] = roof; fp_map[static_cast<Face_handle>(fh)].push_back(*pit);
		++fh;

		fh->info().bu =  2; fh->info().bu_color = r; buildings[2].faces.push_back(static_cast<Face_handle>(fh)); buildings[2].color = r; fp_map[static_cast<Face_handle>(fh)].clear(); ++fh;


		// Create ground plane.
		ground_fitted_plane = Plane_3(Point_3(0.0, 0.0, 0.0), Point_3(1.0, 0.0, 0.0), Point_3(0.0, 1.0, 0.0));

		/*
		Log log; 

		log.save_cdt_ply(cdt, "tmp" + std::string(PS) + "cdt_roof_fitter", "bu");
		log.save_buildings_info(cdt, buildings, "tmp" + std::string(PS) + "buildings_roof_fitter");
		
		std::cout << std::endl;
		for (Point_iterator it = input.begin(); it != input.end(); ++it)
			std::cout << input.point(*it) << ": " << labels[*it] << std::endl;
		std::cout << std::endl;

		std::cout << std::endl;
		for (typename Face_points_map::const_iterator it = fp_map.begin(); it != fp_map.end(); ++it) {
			for (size_t i = 0; i < (*it).second.size(); ++i)
				std::cout << (*it).first->info().bu << " with points " << input.point((*it).second[i]) << "; ";
			std::cout << std::endl;
		}
		std::cout << std::endl; */
	}

	void set_labels_property(Container &input, Label_map &labels) {
		
		auto success = false;
		boost::tie(labels, success)  = input. template add_property_map<Label>("label", -1);
		assert(success);
	}
};

TEST_F(LOD_BuildingRoofFitterTest, Compiles) {
   
	// Empty test.
}

TEST_F(LOD_BuildingRoofFitterTest, VerifiesMinHeightsOfThreeBuildings) {

	lodBuildingRoofFitterMin.fit_roof_heights(cdt, input, fp_map, ground_fitted_plane, buildings);

	ASSERT_THAT(buildings[0].height, Eq(1.0));
	ASSERT_THAT(buildings[1].height, Eq(1.9));
	ASSERT_THAT(buildings[2].height, Eq(0.9));
}

TEST_F(LOD_BuildingRoofFitterTest, VerifiesAvgHeightsOfThreeBuildings) {

	lodBuildingRoofFitterAvg.fit_roof_heights(cdt, input, fp_map, ground_fitted_plane, buildings);

	ASSERT_THAT(buildings[0].height, Eq(1.0));
	ASSERT_THAT(buildings[1].height, Eq(2.0));
	ASSERT_THAT(buildings[2].height, Eq(1.0));
}

TEST_F(LOD_BuildingRoofFitterTest, VerifiesMaxHeightsOfThreeBuildings) {

	lodBuildingRoofFitterMax.fit_roof_heights(cdt, input, fp_map, ground_fitted_plane, buildings);

	ASSERT_THAT(buildings[0].height, Eq(1.0));
	ASSERT_THAT(buildings[1].height, Eq(2.1));
	ASSERT_THAT(buildings[2].height, Eq(1.1));
}