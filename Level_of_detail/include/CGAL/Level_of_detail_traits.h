#ifndef CGAL_LEVEL_OF_DETAIL_TRAITS_H
#define CGAL_LEVEL_OF_DETAIL_TRAITS_H

// STL includes.
#include <map>
#include <vector>

// CGAL includes.
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

// New CGAL includes.
#include <CGAL/Loader/Level_of_detail_loader.h>
#include <CGAL/Loader/Level_of_detail_loader_stub.h>
#include <CGAL/Loader/Level_of_detail_loader_eth.h>
#include <CGAL/Preprocessor/Level_of_detail_preprocessor.h>
#include <CGAL/Selector/Level_of_detail_selector.h>
#include <CGAL/Selector/Level_of_detail_selection_strategy.h>
#include <CGAL/Regularizer/Level_of_detail_vertical_regularizer.h>
#include <CGAL/Regularizer/Level_of_detail_line_regularizer_jean_philippe.h>
#include <CGAL/Regularizer/Segment_regularizer_2/Level_of_detail_segment_regularizer_2.h>
#include <CGAL/Regularizer/Level_of_detail_polygonizer_jean_philippe.h>
#include <CGAL/Projector/Level_of_detail_projector.h>
#include <CGAL/Utils/Level_of_detail_utils.h>
#include <CGAL/Structuring_2/Level_of_detail_structuring_2.h>
#include <CGAL/Visibility_2/Level_of_detail_visibility_2.h>
#include <CGAL/Graphcut/Level_of_detail_graphcut.h>
#include <CGAL/Reconstruction/Level_of_detail_reconstruction.h>
#include <CGAL/Buildings/Level_of_detail_buildings.h>
#include <CGAL/Clutter/Level_of_detail_grid_simplify.h>
#include <CGAL/Clutter/Level_of_detail_thinning.h>
#include <CGAL/Clutter/Level_of_detail_clutter_filtering.h>
#include <CGAL/Clutter/Level_of_detail_clutter_processor.h>
#include <CGAL/Region_growing/Level_of_detail_region_growing.h>
#include <CGAL/Tools/Level_of_detail_complexity.h>
#include <CGAL/Tools/Level_of_detail_distortion.h>
#include <CGAL/Tools/Level_of_detail_coverage.h>
#include <CGAL/Tools/Level_of_detail_parameters.h>
#include <CGAL/Tools/Level_of_detail_parameters_estimator.h>
#include <CGAL/Container/Level_of_detail_container.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class OutputContainer>
		struct Level_of_detail_traits {

			typedef KernelTraits 	Kernel;
			typedef OutputContainer Container_3D;

			typedef CGAL::LOD::Level_of_detail_loader_eth<Kernel, Container_3D>   Loader;
			typedef CGAL::LOD::Level_of_detail_preprocessor<Kernel, Container_3D> Preprocessor;

			typedef CGAL::LOD::Level_of_detail_clutter<Kernel, Container_3D> 		   Clutter_strategy;
			typedef CGAL::LOD::Level_of_detail_ground<Kernel, Container_3D> 		   Ground_strategy;
			typedef CGAL::LOD::Level_of_detail_building_boundary<Kernel, Container_3D> Building_boundary_strategy;
			typedef CGAL::LOD::Level_of_detail_building_interior<Kernel, Container_3D> Building_interior_strategy;

			typedef CGAL::LOD::Level_of_detail_selector<Kernel, Clutter_strategy> 		    Clutter_selector;
			typedef CGAL::LOD::Level_of_detail_selector<Kernel, Ground_strategy> 		    Ground_selector;
			typedef CGAL::LOD::Level_of_detail_selector<Kernel, Building_boundary_strategy> Building_boundary_selector;
			typedef CGAL::LOD::Level_of_detail_selector<Kernel, Building_interior_strategy> Building_interior_selector;

			typedef std::map<int, std::vector<int> >          				   		   			  		Planes;
			typedef std::map<int, typename Kernel::Point_2>           	    				  		 	Projected_points;
			typedef CGAL::LOD::Level_of_detail_simple_projector<Kernel, Container_3D, Projected_points> Ground_projector;
			
			typedef CGAL::LOD::Level_of_detail_vertical_regularizer<Kernel, Container_3D, Planes> Vertical_regularizer;
			typedef CGAL::LOD::Level_of_detail_segment_regularizer_2<Kernel> 					  Line_regularizer;

			typedef CGAL::LOD::My_vertex_info<Structured_label>  My_vertex_info; 
	    	typedef CGAL::LOD::My_face_info<typename Kernel::FT> My_face_info;

			typedef CGAL::Triangulation_vertex_base_with_info_2<My_vertex_info, Kernel> VB;
	 	    typedef CGAL::Triangulation_face_base_with_info_2<My_face_info, Kernel>     FB_with_info;
			typedef CGAL::Constrained_triangulation_face_base_2<Kernel, FB_with_info>   FB;

			   typedef CGAL::Exact_predicates_tag                                TAG;
			// typedef CGAL::Exact_intersections_tag                             TAG;

			typedef CGAL::Triangulation_data_structure_2<VB, FB> 			  	 TDS;
			typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, TAG> CDT;
			
			// typedef CGAL::Constrained_triangulation_plus_2<DDT> CDT;

			typedef CGAL::LOD::Level_of_detail_utils<Kernel, Container_3D, CDT> 				   				 Utils;
			typedef CGAL::LOD::Level_of_detail_clutter_processor<Kernel, Planes, Projected_points, Container_3D> Clutter_processor;
			typedef CGAL::LOD::Level_of_detail_region_growing<Kernel, Planes, Projected_points, Container_3D>    Region_growing;

			typedef CGAL::LOD::Level_of_detail_grid_simplify<Kernel, Planes, Projected_points> 			Grid_simplifier;
			typedef CGAL::LOD::Level_of_detail_thinning<Kernel, Planes, Projected_points, Container_3D> Thinning;
			typedef CGAL::LOD::Level_of_detail_clutter_filtering<Kernel, Planes, Projected_points> 		Clutter_filtering;

			typedef int Label; 
			typedef std::vector< std::pair<typename Kernel::Point_2, Label> > Container_2D;

			typedef CGAL::LOD::Level_of_detail_structuring_2<Kernel> Structuring_2;
			typedef CGAL::LOD::Level_of_detail_graphcut<Kernel, CDT> Graph_cut;

			typedef CGAL::LOD::Level_of_detail_visibility_from_classification_2<Kernel, Container_2D, CDT> Visibility_2;
			
			typedef CGAL::LOD::Building<typename Kernel::FT, typename CDT::Vertex_handle, typename CDT::Face_handle> Building;
			typedef std::map<int, Building> 																	     Buildings;

			typedef CGAL::Polyhedron_3<Kernel> 											    Mesh;
			typedef CGAL::LOD::Level_of_detail_reconstruction<Kernel, CDT, Buildings, Mesh> Lods;

			typedef typename Lods::Mesh_facet_colors Mesh_facet_colors;

			typedef CGAL::LOD::Level_of_detail_building_splitter_2<Kernel, CDT> Building_splitter;
			typedef CGAL::LOD::Level_of_detail_building_outliner_2<Kernel, CDT> Building_outliner;

			typedef CGAL::LOD::Level_of_detail_min_height_fitter<Kernel> Min_height_fitter;
			typedef CGAL::LOD::Level_of_detail_avg_height_fitter<Kernel> Avg_height_fitter;
			typedef CGAL::LOD::Level_of_detail_max_height_fitter<Kernel> Max_height_fitter;
	
			typedef CGAL::LOD::Level_of_detail_building_roof_fitter_2<Kernel, CDT, Container_3D, Min_height_fitter> Building_min_roof_fitter;
			typedef CGAL::LOD::Level_of_detail_building_roof_fitter_2<Kernel, CDT, Container_3D, Avg_height_fitter> Building_avg_roof_fitter;
			typedef CGAL::LOD::Level_of_detail_building_roof_fitter_2<Kernel, CDT, Container_3D, Max_height_fitter> Building_max_roof_fitter;

			typedef CGAL::LOD::Level_of_detail_parameters<typename Kernel::FT> Level_of_detail_parameters;
			typedef typename Level_of_detail_parameters::Input_parameters 	   Parameters;

			typedef CGAL::LOD::Level_of_detail_parameters_estimator<Kernel, Container_3D, Parameters> Parameters_estimator;

			typedef CGAL::LOD::Level_of_detail_complexity<Kernel, Container_3D, Lods> Lod_complexity;
			typedef CGAL::LOD::Level_of_detail_distortion<Kernel, Container_3D, Lods> Lod_distortion;
			typedef CGAL::LOD::Level_of_detail_coverage<Kernel, Container_3D, Lods>   Lod_coverage;

			typedef CGAL::LOD::Level_of_detail_container<Kernel> 								     Lod_data_structure;
			typedef CGAL::LOD::Level_of_detail_polygonizer_jean_philippe<Kernel, Lod_data_structure> Polygonizer;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_TRAITS_H