#ifndef CGAL_LEVEL_OF_DETAIL_H
#define CGAL_LEVEL_OF_DETAIL_H

// STL includes.
#include <memory>
#include <iostream>
#include <algorithm>

// LOD includes.
#include <CGAL/Level_of_detail/Parameters.h>
#include <CGAL/Level_of_detail/Enumerations.h>
#include <CGAL/Level_of_detail/Data_structure.h>

#include <CGAL/Level_of_detail/internal/Property_maps/Point_from_reference_property_map_2.h>
#include <CGAL/Level_of_detail/internal/Property_maps/Dereference_property_map.h>
#include <CGAL/Level_of_detail/internal/Property_maps/Estimated_normal_property_map_2.h>
#include <CGAL/Level_of_detail/internal/Property_maps/Segment_from_region_property_map_2.h>
#include <CGAL/Level_of_detail/internal/Property_maps/Building_height_property_map.h>
#include <CGAL/Level_of_detail/internal/Property_maps/Building_with_segment_constraints_property_map.h>

#include <CGAL/Level_of_detail/internal/Filtering/Alpha_shapes_filtering.h>
#include <CGAL/Level_of_detail/internal/Filtering/Grid_based_filtering.h>

#include <CGAL/Level_of_detail/internal/Data/Kd_tree_with_data_creator.h>

#include <CGAL/Level_of_detail/internal/Estimations/Tree_based_lines_estimator.h>

#include <CGAL/Level_of_detail/internal/Triangulations/Constrained_triangulation_creator.h>

#include <CGAL/Level_of_detail/internal/Shape_detection/Points_based_region_growing_2.h>

#include <CGAL/Level_of_detail/internal/Regularization/Segment_regularizer_parameters.h>
#include <CGAL/Level_of_detail/internal/Regularization/Segment_regularizer_2.h>

#include <CGAL/Level_of_detail/internal/Partitioning/Kinetic_based_partitioning_2.h>

#include <CGAL/Level_of_detail/internal/Visibility/Facet_visibility_estimator.h>
#include <CGAL/Level_of_detail/internal/Visibility/Visibility_consistency.h>

#include <CGAL/Level_of_detail/internal/Buildings/Buildings_creator.h>
#include <CGAL/Level_of_detail/internal/Buildings/Buildings_outliner.h>

#include <CGAL/Level_of_detail/internal/utils.h>

#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/linear_least_squares_fitting_3.h>

namespace CGAL {

	namespace Level_of_detail {

    /*!
      \ingroup PkgLevelOfDetail

      \brief Construct levels of detail from an input point cloud.

      \tparam GeomTraits model of \cgal Kernel.
      \tparam PointRange model of `ConstRange`. The value type of its
      iterator is the key type of `PointMap`.
      \tparam PointMap model of `ReadablePropertyMap` whose key
      type is the value type of the iterator of `PointRange` and value type
      is `GeomTraits::Point_3`.
    */
		template<class GeomTraits, class PointRange, class PointMap>
		class Level_of_detail {

		public:

      /// \cond SKIP_IN_MANUAL
			using Kernel      = GeomTraits;
			using Input_range = PointRange;
			using Point_map   = PointMap;
      /// \endcond

			typedef typename GeomTraits::FT FT;
			typedef typename GeomTraits::Point_3 Point_3;
			typedef Parameters<FT> Parameters;
      
      /// \cond SKIP_IN_MANUAL
			using Point_2 = typename Kernel::Point_2;
			using Plane_3 = typename Kernel::Plane_3;
      using Triangle_3 = typename Kernel::Triangle_3;


			using Data_structure = Data_structure<Kernel, Input_range, Point_map>;

			using Point_identifier  = typename Data_structure::Point_identifier;
			using Point_identifiers = typename Data_structure::Point_identifiers;
			using Triangulation = typename Data_structure::Triangulation;

			using Point_map_2 = Point_from_reference_property_map_2<Point_identifier, Point_2, Point_map>;
			using Point_map_3 = Dereference_property_map<Point_identifier, Point_map>;
			
			using Mesh 	= typename Data_structure::Mesh;
			using Lod_0 = typename Data_structure::Lod_0;
			using Lod_1 = typename Data_structure::Lod_1;

      /// \endcond

      /// \name Constructor
      /// @{

      /*!
        \brief Initializes data structures for Level Of Detail computation.
      */
			Level_of_detail(const PointRange& point_range, const PointMap& point_map, const Parameters& parameters) :
			m_data_structure(point_range, point_map),
			m_parameters(parameters),
			m_point_map_2(m_data_structure.point_map()),
			m_point_map_3(m_data_structure.point_map()) { 
				
				CGAL_assertion(input_range.size() != 0);
			}

      /// @}
      
      /// \name Complete Generation
      /// @{

      /*!
        \brief Generates LOD0 and LOD1.

        This method computes everything needed for LOD0 and LOD1. It
        is equivalent to calling these methods in the following order:

         - `split_semantic_data()` with `semantic_element_map` as a parameter
         - `fit_ground_plane()`
         - `extract_building_boundaries()`
         - `simplify_building_boundaries()`
         - `detect_lines()`
         - `regularize_segments()`
         - `create_partitioning()`
         - `compute_visibility()` with `visibility_map` as a parameter
         - `create_triangulation()`
         - `find_buildings()`
         - `find_building_walls()`
         - `fit_flat_building_roofs()`
         - `compute_triangulation_vertices_heights()`
      */
			template<class SemanticElementMap, class VisibilityMap>
			void build (SemanticElementMap semantic_element_map, VisibilityMap visibility_map) {

				if (m_parameters.verbose()) std::cout << std::endl << "... building LOD data ..." << std::endl << std::endl;

				split_semantic_data(semantic_element_map);
				
				fit_ground_plane();
				
				extract_building_boundaries();
				
				simplify_building_boundaries();

				detect_lines();

				regularize_segments();

				create_partitioning();

				compute_visibility(visibility_map);

				create_triangulation();

				find_buildings();

				find_building_walls();

				fit_flat_building_roofs();
        
        compute_triangulation_vertices_heights();
			}

      /// @}

      /// \cond SKIP_IN_MANUAL

			template<class Lod>
			void get_lod(Lod &lod) {
				
				if (m_parameters.verbose()) std::cout << "* constructing " << lod.name() << std::endl;
				lod.reconstruct(m_data_structure.buildings(), m_data_structure.ground_bounding_box());
			}

      /// \endcond

      /// \name Step by Step Generation
      /// @{

			template<class SemanticElementMap>
			void split_semantic_data (SemanticElementMap semantic_element_map) {
				if (m_parameters.verbose()) std::cout << "* splitting semantic data" << std::endl;

				// In this step, we split only ground, building interior, and building boundaries.
        m_data_structure.ground_points().clear();
				m_data_structure.building_boundary_points().clear();
				m_data_structure.building_interior_points().clear();
        m_data_structure.vegetation_points().clear();

				for (typename Input_range::const_iterator point
          = m_data_structure.input_range().begin(); point != m_data_structure.input_range().end(); ++point)
        {
          const Semantic_label label = get (semantic_element_map, *point);
					switch (label) {

						case Semantic_label::GROUND:
							m_data_structure.ground_points().push_back(point);
							break;

						case Semantic_label::BUILDING_BOUNDARY:
							m_data_structure.building_boundary_points().push_back(point);
							break;

						case Semantic_label::BUILDING_INTERIOR:
							m_data_structure.building_interior_points().push_back(point);
							break;

						case Semantic_label::VEGETATION:
							m_data_structure.vegetation_points().push_back(point);
							break;

						default:
							break;
					}
				}
			}

			void fit_ground_plane() {

        using Local_kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
        using Local_FT 		= typename Local_kernel::FT;
        using Local_point_3 = typename Local_kernel::Point_3;
        using Local_plane_3 = typename Local_kernel::Plane_3;

				if (m_parameters.verbose()) std::cout << "* fitting ground plane" << std::endl;

				// Here, we fit a plane to all ground points.
				CGAL_precondition(m_data_structure.ground_points().size() > 2);
        
				std::size_t i = 0;
				std::vector<Local_point_3> points(m_data_structure.ground_points().size());
				
				for (typename Point_identifiers::const_iterator ce_it = m_data_structure.ground_points().begin();
             ce_it != m_data_structure.ground_points().end(); ++ce_it, ++i)
        {
					const Point_3 &point = get(m_point_map_3, *ce_it);

					const Local_FT x = static_cast<Local_FT>(CGAL::to_double(point.x()));
					const Local_FT y = static_cast<Local_FT>(CGAL::to_double(point.y()));
					const Local_FT z = static_cast<Local_FT>(CGAL::to_double(point.z()));

					points[i] = Local_point_3(x, y, z);
				}

        Local_point_3 centroid;
				Local_plane_3 fitted_plane;

				CGAL::linear_least_squares_fitting_3 (points.begin(), points.end(),
                                              fitted_plane, centroid, CGAL::Dimension_tag<0>(), Local_kernel(),
                                              CGAL::Eigen_diagonalize_traits<Local_FT, 3>());

        m_data_structure.ground_plane() =
          Plane_3(static_cast<FT>(fitted_plane.a()), static_cast<FT>(fitted_plane.b()),
                  static_cast<FT>(fitted_plane.c()), static_cast<FT>(fitted_plane.d()));
        
				internal::compute_bounding_box_3(m_data_structure.ground_points(), m_point_map_3,
                                         m_data_structure.ground_plane(), m_data_structure.ground_bounding_box());

				auto  it = m_data_structure.ground_bounding_box().begin();
				const Point_3 &p1 = *it; ++it;
				const Point_3 &p2 = *it; ++it;
				const Point_3 &p3 = *it; ++it;

				m_data_structure.ground_plane() = Plane_3(p1, p2, p3);
//				m_data_structure.ground_points().clear();
			}

			void extract_building_boundaries() {
				if (m_parameters.verbose()) std::cout << "* extracting building boundary points" << std::endl;

				// In this step, we apply alpha shapes to extract only building boundary points.
				CGAL_precondition(m_data_structure.building_boundary_points().size() > 2 ||
                          m_data_structure.building_interior_points().size() > 2);
				
				m_data_structure.filtered_building_boundary_points().clear();
				const Alpha_shapes_filtering<Kernel, Point_identifier>
          alpha_shapes_filtering(m_parameters.alpha_shape_size());
				
				if (m_data_structure.building_boundary_points().size() > 2)
					alpha_shapes_filtering.add_points(m_data_structure.building_boundary_points(), m_point_map_2,
                                            m_data_structure.filtered_building_boundary_points());

				if (m_data_structure.building_interior_points().size() > 2)
					alpha_shapes_filtering.add_points(m_data_structure.building_interior_points(), m_point_map_2,
                                            m_data_structure.filtered_building_boundary_points());

				m_data_structure.building_boundary_points().clear();

        if (m_parameters.verbose())
          std::cout << " -> " << m_data_structure.filtered_building_boundary_points().size()
                    << " boundary point(s) extracted" << std::endl;
			}

			void simplify_building_boundaries() {
				if (m_parameters.verbose()) std::cout << "* simplifying building boundary points";

				// Here, we apply grid-based simplification to all building boundary points.
				if (m_parameters.no_simplification()) {
					
					m_data_structure.simplified_building_boundary_points() = m_data_structure.filtered_building_boundary_points();
					if (m_parameters.verbose()) std::cout << " - skipped";

				} else {

					const Grid_based_filtering<Kernel, Point_identifier>
            grid_based_filtering(m_parameters.grid_cell_width());
          
					grid_based_filtering.apply(m_data_structure.filtered_building_boundary_points(),
                                     m_point_map_2, m_data_structure.simplified_building_boundary_points());	
				}

				if (m_parameters.verbose()) std::cout << std::endl;
				m_data_structure.filtered_building_boundary_points().clear();
        
        if (m_parameters.verbose())
          std::cout << " -> " << m_data_structure.simplified_building_boundary_points().size()
                    << " boundary point(s) remaining" << std::endl;
			}

			void detect_lines() {
        
        using Points_tree_2 		       = Kd_tree_with_data_creator<Kernel, Point_identifier, Point_identifiers, Point_map_2>;
        using Tree_based_lines_estimator_2 = Tree_based_lines_estimator<Kernel, Point_identifiers, Point_map_2, Points_tree_2>;

        using Region_growing_2_normal_map = Estimated_normal_property_map_2<Kernel, Tree_based_lines_estimator_2>;
        using Region_growing_2       	  = Points_based_region_growing_2<Kernel, Points_tree_2>;

				if (m_parameters.verbose()) std::cout << "* detecting lines along building boundaries" << std::endl;

				// In this step, we apply region growing to detect 2D lines that form building boundaries.
				const Points_tree_2 points_tree_2(
            m_data_structure.simplified_building_boundary_points(),
            m_point_map_2, 
            m_parameters.region_growing_2_cluster_epsilon());

				const Tree_based_lines_estimator_2 lines_estimator_2(
            m_data_structure.simplified_building_boundary_points(), 
            m_point_map_2, 
            points_tree_2);

        const typename Tree_based_lines_estimator_2::Sorter sorter = lines_estimator_2.sorter();

				std::stable_sort(m_data_structure.simplified_building_boundary_points().begin(), m_data_structure.simplified_building_boundary_points().end(), sorter);
				
				const Region_growing_2_normal_map normal_map(
            m_data_structure.simplified_building_boundary_points(), 
            m_point_map_2, 
            lines_estimator_2.lines_2());

				Region_growing_2 region_growing_2(
					m_parameters.region_growing_2_epsilon(), 
					m_parameters.region_growing_2_cluster_epsilon(), 
					m_parameters.region_growing_2_normal_threshold(),
					m_parameters.region_growing_2_min_points(),
					points_tree_2);

				region_growing_2.detect(m_data_structure.simplified_building_boundary_points(), m_point_map_2, normal_map, m_data_structure.detected_2d_regions());
				m_data_structure.simplified_building_boundary_points().clear();
        
        if (m_parameters.verbose())
          std::cout << " -> " << m_data_structure.detected_2d_regions().size()
                    << " line(s) detected" << std::endl;
			}

			void regularize_segments() {
        
				if (m_parameters.verbose()) std::cout << "* regularizing segments detected along building boundaries";

				// Here, we regularize segments that form building boundaries wrt to angles and ordinates.
				const Segment_from_region_property_map_2<Kernel, Point_identifiers, Point_map_2,
                                                 typename Data_structure::Regularized_segments>
          segment_from_region_map_2(m_data_structure.detected_2d_regions(), m_point_map_2);

				Segment_regularizer_parameters<FT> segment_regularizer_parameters;				
				segment_regularizer_parameters.max_angle_in_degrees() 	  = m_parameters.segment_regularizer_2_max_angle_in_degrees();
				segment_regularizer_parameters.max_difference_in_meters() = m_parameters.segment_regularizer_2_max_difference_in_meters();

				if (m_parameters.no_regularization()) {
					if (m_parameters.verbose()) std::cout << " - skipped ";

					segment_regularizer_parameters.optimize_angles() 	= false;
					segment_regularizer_parameters.optimize_ordinates() = false;
				}
				
				Segment_regularizer_2<Kernel> segment_regularizer_2(segment_regularizer_parameters);
				segment_regularizer_2.regularize(m_data_structure.detected_2d_regions(), segment_from_region_map_2, m_data_structure.regularized_segments());

				if (m_parameters.verbose()) std::cout << std::endl;
				m_data_structure.detected_2d_regions().clear();
        
        if (m_parameters.verbose())
          std::cout << " -> " << m_data_structure.regularized_segments().size()
                    << " line(s) after regularization" << std::endl;
			}

			void create_partitioning() {

				if (m_parameters.verbose()) std::cout << "* computing partitioning" << std::endl;

				// In this step, we create a 2D partitioning of the domain.
				const Kinetic_based_partitioning_2<Kernel, typename Data_structure::Partition_face_2>
          kinetic_based_partitioning_2(
            m_parameters.kinetic_partitioning_2_num_intersections(),
            m_parameters.kinetic_partitioning_2_min_face_width());
					
				kinetic_based_partitioning_2.compute(m_data_structure.regularized_segments(), m_data_structure.regularized_segment_map(), m_data_structure.partition_faces_2());
        
        if (m_parameters.verbose())
          std::cout << " -> " << m_data_structure.partition_faces_2().size()
                    << " partition face(s) created" << std::endl;
			}

			template<class VisibilityMap>
			void compute_visibility(VisibilityMap visibility_map) {
				if (m_parameters.verbose()) std::cout << "* computing visibility" << std::endl;

        Facet_visibility_estimator<Kernel, Input_range, Point_map, VisibilityMap>
          visibility_estimator (m_data_structure.input_range(), m_data_structure.point_map(), visibility_map);

				// Here, we try to guess which of the partitioning faces is inside or outside the building.
        typename Data_structure::Partition_faces_2& facets_range
          = m_data_structure.partition_faces_2();
        for (typename Data_structure::Partition_faces_2::iterator
               facet = facets_range.begin(); facet != facets_range.end(); ++facet)
					visibility_estimator.estimate_visibility (*facet);

			}

			void create_triangulation() {
        
				if (m_parameters.verbose()) std::cout << "* creating triangulation" << std::endl;

				// In this step, we build constrained Delaunay triangulation.
				Identity_property_map<Point_2> partition_point_map;

				const Constrained_triangulation_creator<Kernel, Triangulation>
          constrained_triangulation_creator;
				constrained_triangulation_creator.make_triangulation_with_info(
          m_data_structure.ground_bounding_box(), m_parameters.scale() * 100.,
					m_data_structure.partition_faces_2(), 
					partition_point_map,
					m_data_structure.triangulation());

				if (!m_parameters.no_consistent_visibility()) {

					const Visibility_consistency<Triangulation>
            triangulation_visibility_consistency;
					triangulation_visibility_consistency.make_consistent(m_data_structure.triangulation());
				}

				m_data_structure.partition_faces_2().clear();
			}

			void find_buildings() {
        
				if (m_parameters.verbose()) std::cout << "* searching for buildings" << std::endl;

				// Here, we search for sets of triangles that form buildings.
				const Building_with_segment_constraints_property_map<Kernel, Triangulation>
          face_to_building_map(
            m_data_structure.triangulation(),
            m_data_structure.regularized_segments(),
            m_data_structure.regularized_segment_map(),
            m_parameters.segment_constraints_threshold());
				
        Triangulation& triangulation = m_data_structure.triangulation();
				for (typename Triangulation::Finite_faces_iterator
               face = triangulation.finite_faces_begin(); face != triangulation.finite_faces_end(); ++face)
        {
          typename Triangulation::Face_handle face_handle
            = static_cast<typename Triangulation::Face_handle>(face);
					put(face_to_building_map, face_handle);
				}
        
				const Buildings_creator<Kernel, typename Data_structure::Building>
          buildings_creator(m_parameters.min_num_building_floor_faces());
				buildings_creator.create(m_data_structure.triangulation(), m_data_structure.buildings());

				m_data_structure.regularized_segments().clear();
        if (m_parameters.verbose())
          std::cout << " -> " << m_data_structure.buildings().size()
                    << " buildings(s) found" << std::endl;
			}

			void find_building_walls() {
        
				if (m_parameters.verbose()) std::cout << "* searching for building walls" << std::endl;

				// In this step, we search for sets of segments that form building walls.
				const Buildings_outliner<Kernel, typename Data_structure::Building> buildings_outliner;
				buildings_outliner.find_walls(m_data_structure.triangulation(), m_data_structure.buildings());
			}

			void fit_flat_building_roofs() {
        
				if (m_parameters.verbose()) std::cout << "* fitting flat building roofs" << std::endl;

				// Here, we fit flat roofs to all buildings with the average (see parameters) building height.
        for (typename Triangulation::Finite_faces_iterator
               tf_it = m_data_structure.triangulation().finite_faces_begin();
             tf_it != m_data_structure.triangulation().finite_faces_end(); ++tf_it)
          tf_it->info().elements().clear();

        CGAL_precondition(m_data_structure.building_interior_points().size() > 0);

        for (typename Data_structure::Point_identifiers::const_iterator ce_it
               = m_data_structure.building_interior_points().begin();
             ce_it != m_data_structure.building_interior_points().end(); ++ce_it)
        {
          const Point_2 &query = get (m_point_map_2, *ce_it);

          typename Triangulation::Locate_type locate_type;
          int locate_index_stub;
          const typename Triangulation::Face_handle face_handle
            = m_data_structure.triangulation().locate(query, locate_type, locate_index_stub);

					if (locate_type == Triangulation::FACE ||
              locate_type == Triangulation::EDGE ||
              locate_type == Triangulation::VERTEX)
						face_handle->info().add_element(*ce_it);
        }


				const Building_height_property_map<Kernel, Triangulation,
                                           typename Data_structure::Building>
          lod_building_height_map(
            m_data_structure.triangulation(), 
            m_point_map_3,
            m_data_structure.ground_plane(),
            m_parameters.flat_roof_type());

        typename Data_structure::Buildings& buildings = m_data_structure.buildings();
				for (typename Data_structure::Buildings::iterator
               bu_it = buildings.begin(); bu_it != buildings.end(); ++bu_it)
          put (lod_building_height_map, *bu_it);
        
				m_data_structure.building_interior_points().clear();
			}
      
			void compute_triangulation_vertices_heights() {
        
				if (m_parameters.verbose()) std::cout << "* computing triangulation vertices heights" << std::endl;

        // First pass: init everything to ground plane
				for (typename Triangulation::Finite_faces_iterator
               face = m_data_structure.triangulation().finite_faces_begin();
             face != m_data_structure.triangulation().finite_faces_end(); ++face)
          for (std::size_t j = 0; j < 3; ++ j)
            face->info().height(j) = internal::position_on_plane (m_data_structure.ground_plane(),
                                                                  face->vertex(j)->point()).z();

        // Second pass: init building heights
        for (typename Data_structure::Buildings::iterator it = m_data_structure.buildings().begin();
             it != m_data_structure.buildings().end(); ++ it)
          for (typename Data_structure::Building::Floor_face_handles::iterator fit = it->floor_face_handles().begin();
               fit != it->floor_face_handles().end(); ++ fit)
            for (std::size_t j = 0; j < 3; ++ j)
              (*fit)->info().height(j) = it->height() + m_data_structure.ground_bounding_box().begin()->z();

        // Third pass: compute ground real heights
        using Points_tree_2 = Kd_tree_with_data_creator<Kernel, Point_identifier, Point_identifiers, Point_map_2>;
        
				const Points_tree_2 points_tree_2(
            m_data_structure.ground_points(),
            m_point_map_2, 
            int(6));

				for (typename Triangulation::Finite_faces_iterator
               face = m_data_structure.triangulation().finite_faces_begin();
             face != m_data_structure.triangulation().finite_faces_end(); ++face)
          if (face->info().visibility_label() == Visibility_label::OUTSIDE)
            for (std::size_t j = 0; j < 3; ++ j)
            {
              const Point_2& p2 = face->vertex(j)->point();
              typename Points_tree_2::Neighbours neighbors;
              points_tree_2.search_knn_2 (p2, neighbors);
              double h_mean = 0.;
              for (std::size_t i = 0; i < neighbors.size(); ++ i)
                h_mean += get (m_point_map_3, neighbors[i].second).z();
              face->info().height(j) = h_mean / neighbors.size();
            }

        // Fourth pass, refine ground
        std::vector<std::pair<FT, Point_identifier> > out_of_tolerance;
        typename Triangulation::Face_handle hint;
        for (typename Point_identifiers::const_iterator ce_it = m_data_structure.ground_points().begin();
             ce_it != m_data_structure.ground_points().end(); ++ce_it)
        {
					const Point_2& point_2 = get(m_point_map_2, *ce_it);
          const Point_3& point_3 = get(m_point_map_3, *ce_it);
          
          hint = m_data_structure.triangulation().locate (point_2, hint);
          if (hint->info().visibility_label() != Visibility_label::OUTSIDE)
            continue;

          Triangle_3 triangle = internal::triangle_3<Triangle_3>(hint);

          FT sq_dist = CGAL::squared_distance (point_3, triangle);
          if (sq_dist > m_parameters.scale() * m_parameters.scale())
            out_of_tolerance.push_back (std::make_pair (sq_dist, *ce_it));
				}

        std::sort (out_of_tolerance.begin(), out_of_tolerance.end());

        std::size_t previous_out_side = out_of_tolerance.size();
        while (!out_of_tolerance.empty())
        {
          typename Triangulation::Vertex_handle
            v = m_data_structure.triangulation().insert (get(m_point_map_2, out_of_tolerance.back().second));

          typename Triangulation::Face_circulator circ = m_data_structure.triangulation().incident_faces(v);
          typename Triangulation::Face_circulator start = circ;
          do
          {
            if (circ->info().visibility_label() == Visibility_label::OUTSIDE)
            {
              for (std::size_t j = 0; j < 3; ++ j)
              {
                const Point_2& p2 = circ->vertex(j)->point();
                typename Points_tree_2::Neighbours neighbors;
                points_tree_2.search_knn_2 (p2, neighbors);
                double h_mean = 0.;
                for (std::size_t i = 0; i < neighbors.size(); ++ i)
                  h_mean += get (m_point_map_3, neighbors[i].second).z();
                circ->info().height(j) = h_mean / neighbors.size();
              }
            }
            ++ circ;
          }
          while (circ != start);

          std::vector<std::pair<FT, Point_identifier> > new_out_of_tolerance;

          for (std::size_t i = 0; i < out_of_tolerance.size() - 1; ++ i)
          {
            const Point_2& point_2 = get(m_point_map_2, out_of_tolerance[i].second);
            const Point_3& point_3 = get(m_point_map_3, out_of_tolerance[i].second);
          
            hint = m_data_structure.triangulation().locate (point_2, hint);
            if (hint->info().visibility_label() != Visibility_label::OUTSIDE)
              continue;

            Triangle_3 triangle = internal::triangle_3<Triangle_3>(hint);

            FT sq_dist = CGAL::squared_distance (point_3, triangle);
            if (sq_dist > m_parameters.scale() * m_parameters.scale())
              new_out_of_tolerance.push_back (std::make_pair (sq_dist, out_of_tolerance[i].second));
          }
          new_out_of_tolerance.swap (out_of_tolerance);

          if (previous_out_side == out_of_tolerance.size())
            break;
          previous_out_side = out_of_tolerance.size();

          std::sort (out_of_tolerance.begin(), out_of_tolerance.end());
        }
      }

      /// @}

      /// \name Output
      /// @{

      template <typename Polygon>
      std::size_t
      output_lod0_to_polygon_soup (std::vector<Point_3>& vertices,
                                   std::vector<Polygon>& polygons) const
      {
        std::vector<typename Triangulation::Face_handle> ground_faces;
        std::vector<typename Triangulation::Face_handle> roof_faces;
        std::vector<typename Triangulation::Face_handle> vegetation_faces;

        internal::segment_semantic_faces (m_data_structure.triangulation(),
                                          ground_faces, roof_faces,
                                          vegetation_faces);

        internal::Indexer<Point_2> indexer;

        std::size_t out = 0;
        
        for (std::size_t i = 0; i < ground_faces.size(); ++ i)
        {
          polygons.push_back (std::vector<std::size_t>());

          for (std::size_t j = 0; j < 3; ++ j)
          {
            std::size_t idx = indexer(ground_faces[i]->vertex(j)->point());
            if (idx == vertices.size())
              vertices.push_back (internal::position_on_plane (m_data_structure.ground_plane(),
                                                               ground_faces[i]->vertex(j)->point()));
            polygons.back().push_back (idx);
          }
        }

        out = polygons.size();

        for (std::size_t i = 0; i < roof_faces.size(); ++ i)
        {
          polygons.push_back (std::vector<std::size_t>());

          for (std::size_t j = 0; j < 3; ++ j)
          {
            std::size_t idx = indexer(roof_faces[i]->vertex(j)->point());
            if (idx == vertices.size())
              vertices.push_back (internal::position_on_plane (m_data_structure.ground_plane(),
                                                               roof_faces[i]->vertex(j)->point()));
            polygons.back().push_back (idx);
          }
        }

        return out;
      }

      template <typename Polygon>
      std::pair<std::size_t, std::size_t>
      output_lod1_to_polygon_soup (std::vector<Point_3>& vertices,
                                   std::vector<Polygon>& polygons) const
      {
        std::vector<typename Triangulation::Face_handle> ground_faces;
        std::vector<typename Triangulation::Face_handle> roof_faces;
        std::vector<typename Triangulation::Face_handle> vegetation_faces;

        internal::segment_semantic_faces (m_data_structure.triangulation(),
                                          ground_faces, roof_faces,
                                          vegetation_faces);

        internal::Indexer<Point_3> indexer;

        std::pair<std::size_t, std::size_t> out;
        
        for (std::size_t i = 0; i < ground_faces.size(); ++ i)
        {
          polygons.push_back (std::vector<std::size_t>());

          for (std::size_t j = 0; j < 3; ++ j)
          {
            std::size_t idx = indexer(internal::point_3<Point_3>(ground_faces[i], j));
            if (idx == vertices.size())
              vertices.push_back (internal::point_3<Point_3> (ground_faces[i], j));
            polygons.back().push_back (idx);
          }
        }

        out.first = polygons.size();

        for (std::size_t i = 0; i < roof_faces.size(); ++ i)
        {
          polygons.push_back (std::vector<std::size_t>());

          for (std::size_t j = 0; j < 3; ++ j)
          {
            std::size_t idx = indexer(internal::point_3<Point_3>(roof_faces[i], j));
            if (idx == vertices.size())
              vertices.push_back (internal::point_3<Point_3>(roof_faces[i], j));
            polygons.back().push_back (idx);
          }
        }

        out.second = polygons.size();

        // Get wall faces
				for (typename Triangulation::Finite_edges_iterator
               e = m_data_structure.triangulation().finite_edges_begin();
             e != m_data_structure.triangulation().finite_edges_end(); ++ e)
        {
          typename Triangulation::Face_handle f0 = e->first;
          typename Triangulation::Face_handle f1 = e->first->neighbor(e->second);

          if (m_data_structure.triangulation().is_infinite(f0) ||
              m_data_structure.triangulation().is_infinite(f1))
            continue;

          typename Triangulation::Vertex_handle va = e->first->vertex((e->second + 1)%3);
          typename Triangulation::Vertex_handle vb = e->first->vertex((e->second + 2)%3);
            
          std::vector<Point_3> points; points.reserve(4);

          Point_3 p0a = internal::point_3<Point_3>(f0, f0->index(va));
          points.push_back (p0a);
          Point_3 p1a = internal::point_3<Point_3>(f1, f1->index(va));
          if (p1a != p0a) points.push_back (p1a);
          Point_3 p1b = internal::point_3<Point_3>(f1, f1->index(vb));
          points.push_back (p1b);
          Point_3 p0b = internal::point_3<Point_3>(f0, f0->index(vb));
          if (p0b != p1b) points.push_back (p0b);
          
          if (points.size() > 2) // facet if > than 2 vertices
          {
            polygons.push_back (std::vector<std::size_t>());
            
            for (std::size_t j = 0; j < points.size(); ++ j)
            {
              std::size_t idx = indexer(points[j]);
              if (idx == vertices.size())
                vertices.push_back (points[j]);
              polygons.back().push_back (idx);
            }
          }
        }
        
        return out;
      }

      /// @}

      /// \cond SKIP_IN_MANUAL

			//////////////////////////////////
			// Functions to be not documented!

			inline const Data_structure& get_internal_data_structure() const {
				return m_data_structure;
			}
      inline const Point_map_2& point_map_2() const {
        return m_point_map_2;
      }

      /// \endcond

		private:
			Data_structure m_data_structure;
			const Parameters &m_parameters;

			const Point_map_2 m_point_map_2;
			const Point_map_3 m_point_map_3;
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_H
