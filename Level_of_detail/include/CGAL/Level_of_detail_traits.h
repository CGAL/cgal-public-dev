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

// New CGAL includes.
#include <CGAL/Loader/Level_of_detail_loader_stub.h>
#include <CGAL/Preprocessor/Level_of_detail_preprocessor.h>
#include <CGAL/Selector/Level_of_detail_selector.h>
#include <CGAL/Selector/Level_of_detail_selection_strategy.h>
#include <CGAL/Regularizer/Level_of_detail_regularizer.h>
#include <CGAL/Projector/Level_of_detail_projector.h>
#include <CGAL/Structuring_2/Level_of_detail_structuring_2.h>
#include <CGAL/Visibility_2/Level_of_detail_visibility_2.h>
#include <CGAL/Lod_0/Level_of_detail_reconstruction_0.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class OutputContainer>
		struct Level_of_detail_traits {

			typedef KernelTraits 	Kernel;
			typedef OutputContainer Container_3D;

			typedef Level_of_detail_loader_stub<Kernel, Container_3D>  Loader;
			typedef Level_of_detail_preprocessor<Kernel, Container_3D> Preprocessor;

			typedef Level_of_detail_clutter<Kernel, Container_3D> 			Clutter_strategy;
			typedef Level_of_detail_ground<Kernel, Container_3D> 			Ground_strategy;
			typedef Level_of_detail_building_boundary<Kernel, Container_3D> Building_boundary_strategy;
			typedef Level_of_detail_building_interior<Kernel, Container_3D> Building_interior_strategy;

			typedef Level_of_detail_selector<Kernel, Clutter_strategy> 		     Clutter_selector;
			typedef Level_of_detail_selector<Kernel, Ground_strategy> 		     Ground_selector;
			typedef Level_of_detail_selector<Kernel, Building_boundary_strategy> Building_boundary_selector;
			typedef Level_of_detail_selector<Kernel, Building_interior_strategy> Building_interior_selector;

			typedef std::map<int, std::vector<int> >          				   		   Planes;
			typedef Level_of_detail_vertical_regularizer<Kernel, Container_3D, Planes> Vertical_regularizer;
			
			typedef std::map<int, typename Kernel::Point_2>           	    				  		 Projected_points;
			typedef Level_of_detail_simple_projector<Kernel, Container_3D, Planes, Projected_points> Ground_projector;

			typedef CGAL::LOD::My_vertex_info<Structured_label>  My_vertex_info; 
	    	typedef CGAL::LOD::My_face_info<typename Kernel::FT> My_face_info;

			typedef CGAL::Triangulation_vertex_base_with_info_2<My_vertex_info, Kernel> VB;
	 	    typedef CGAL::Triangulation_face_base_with_info_2<My_face_info, Kernel>     FB_with_info;
			typedef CGAL::Constrained_triangulation_face_base_2<Kernel, FB_with_info>   FB;

			typedef CGAL::Triangulation_data_structure_2<VB, FB> 			TDS;
			typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS> CDT;

			typedef int Label; 
			typedef std::vector< std::pair<typename Kernel::Point_2, Label> > Container_2D;

			typedef Level_of_detail_structuring_2<Kernel>                                    	Structuring_2;
			typedef Level_of_detail_visibility_from_classification_2<Kernel, Container_2D, CDT> Visibility_2;
			
			typedef Level_of_detail_reconstruction_0<Kernel, CDT> Lod_0;
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_TRAITS_H