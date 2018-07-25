#ifndef CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H
#define CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H

// STL indludes.
#include <list>
#include <vector>

// CGAL includes.
#include <CGAL/Polygon_2.h>
#include <CGAL/property_map.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>
#include <CGAL/Level_of_detail/internal/Buildings/Building.h>
#include <CGAL/Level_of_detail/internal/Reconstruction/Lod_0.h>
#include <CGAL/Level_of_detail/internal/Reconstruction/Lod_1.h>
#include <CGAL/Level_of_detail/internal/Partitioning/Partition_element.h>
#include <CGAL/Level_of_detail/internal/Triangulations/Triangulation_face_info.h>
#include <CGAL/Level_of_detail/internal/Triangulations/Triangulation_vertex_info.h>

#include <CGAL/Iterator_range.h>
#include <boost/iterator/filter_iterator.hpp>


namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

		template<typename GeomTraits, typename InputRange, typename PointMap,
             typename SemanticMap, typename VisibilityMap>
		struct Data_structure {

		public:			
            using Kernel      = GeomTraits;
      			using Input_range = InputRange;
            using Point_map   = PointMap;

            typedef typename InputRange::const_iterator iterator;
            typedef typename iterator::value_type Range_type;

            struct Filter_points_by_label
            {
              Semantic_label label;
              SemanticMap semantic_map;
              Filter_points_by_label () { }
              Filter_points_by_label (Semantic_label label, SemanticMap semantic_map)
                : label (label), semantic_map (semantic_map) { }

              bool operator() (const typename SemanticMap::key_type& k) const
              {
                Semantic_label l = get (semantic_map, k);
                return (get (semantic_map, k) == label); // Keep elements with good label
              }
            };

            typedef boost::filter_iterator<Filter_points_by_label, iterator> Filtered_iterator;
            typedef Iterator_range<Filtered_iterator> Filtered_range;

            using FT        = typename Kernel::FT;
            using Point_3   = typename Kernel::Point_3;
            using Point_2   = typename Kernel::Point_2;
            using Plane_3   = typename Kernel::Plane_3;
            using Segment_2 = typename Kernel::Segment_2;

            using Polygon_2 = CGAL::Polygon_2<Kernel>;
            using Polygon_3 = std::vector<Point_3>;

            using Detected_regions = std::vector<std::vector<std::size_t> >;
            
            using Regularized_segments = std::vector<Segment_2>;

            using Partition_face_2  = LOD::Partition_element<Kernel, Polygon_2>;
            using Partition_faces_2 = std::vector<Partition_face_2>;

            using Face_info   = LOD::Triangulation_face_info<Filtered_iterator, FT>;
            using Vertex_info = LOD::Triangulation_vertex_info;

            using VB           = CGAL::Triangulation_vertex_base_with_info_2<Vertex_info, Kernel>;
            using FB_with_info = CGAL::Triangulation_face_base_with_info_2<Face_info, Kernel>;
            using FB           = CGAL::Constrained_triangulation_face_base_2<Kernel, FB_with_info>;

            using TAG = CGAL::Exact_predicates_tag;
            using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;

            using Triangulation             = CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, TAG>;
            using Triangulation_face_handle = typename Triangulation::Face_handle;

            using Building  = LOD::Building<Kernel, Triangulation_face_handle>;
            using Buildings = std::vector<Building>;

            using Mesh = CGAL::Polyhedron_3<Kernel>;

            using Lod_0 = LOD::Lod_0<Kernel, Building, Mesh>;
            using Lod_1 = LOD::Lod_1<Kernel, Building, Mesh>;

            Data_structure (const Input_range& input_range, Point_map point_map,
                            SemanticMap semantic_map, VisibilityMap visibility_map)
              : m_input_range(input_range),
                m_point_map(point_map),
                m_semantic_map (semantic_map),
                m_visibility_map (visibility_map)
            { }

            // Input.
            inline const Input_range& input_range() const {
                return m_input_range;
            }

            inline const Point_map& point_map() const {
                return m_point_map;
            }

            inline SemanticMap& semantic_map() {
                return m_semantic_map;
            }

            inline VisibilityMap& visibility_map() {
                return m_visibility_map;
            }

            // Points.
            inline Filtered_range ground_points() const {
              return make_range (boost::make_filter_iterator
                                 (Filter_points_by_label(Semantic_label::GROUND, m_semantic_map),
                                  m_input_range.begin(), m_input_range.end()),
                                 boost::make_filter_iterator
                                 (Filter_points_by_label(Semantic_label::GROUND, m_semantic_map),
                                  m_input_range.end(), m_input_range.end()));
            }

            inline Filtered_range building_boundary_points() const {
              return make_range (boost::make_filter_iterator
                                 (Filter_points_by_label(Semantic_label::BUILDING_BOUNDARY, m_semantic_map),
                                  m_input_range.begin(), m_input_range.end()),
                                 boost::make_filter_iterator
                                 (Filter_points_by_label(Semantic_label::BUILDING_BOUNDARY, m_semantic_map),
                                  m_input_range.end(), m_input_range.end()));
            }

            inline Filtered_range building_interior_points() const {
              return make_range (boost::make_filter_iterator
                                 (Filter_points_by_label(Semantic_label::BUILDING_INTERIOR, m_semantic_map),
                                  m_input_range.begin(), m_input_range.end()),
                                 boost::make_filter_iterator
                                 (Filter_points_by_label(Semantic_label::BUILDING_INTERIOR, m_semantic_map),
                                  m_input_range.end(), m_input_range.end()));
            }

            inline Filtered_range vegetation_points() const {
              return make_range (boost::make_filter_iterator
                                 (Filter_points_by_label(Semantic_label::VEGETATION, m_semantic_map),
                                  m_input_range.begin(), m_input_range.end()),
                                 boost::make_filter_iterator
                                 (Filter_points_by_label(Semantic_label::VEGETATION, m_semantic_map),
                                  m_input_range.end(), m_input_range.end()));
            }

            inline std::vector<Point_2>& filtered_building_boundary_points() {
                return m_filtered_building_boundary_points;
            }

            inline const std::vector<Point_2>& filtered_building_boundary_points() const {
                return m_filtered_building_boundary_points;
            }

            inline std::vector<Point_2>& simplified_building_boundary_points() {
                return m_simplified_building_boundary_points;
            }

            inline const std::vector<Point_2>& simplified_building_boundary_points() const {
                return m_simplified_building_boundary_points;
            }

            // Ground.
            inline Plane_3& ground_plane() {
                return m_ground_plane;
            }

            inline const Plane_3& ground_plane() const {
                return m_ground_plane;
            }

            inline Polygon_3& ground_bounding_box() {
                return m_ground_bounding_box;
            }

            inline const Polygon_3& ground_bounding_box() const {
                return m_ground_bounding_box;
            }

            // Regions.
            inline Detected_regions& detected_2d_regions() {
                return m_detected_2d_regions;
            }

            inline const Detected_regions& detected_2d_regions() const {
                return m_detected_2d_regions;
            }

            // Regularized segments.
            inline Regularized_segments& regularized_segments() {
                return m_regularized_segments;
            }

            inline const Regularized_segments& regularized_segments() const {
                return m_regularized_segments;
            }

            // Partitioning.
            inline Partition_faces_2& partition_faces_2() {
                return m_partition_faces_2;
            }

            inline const Partition_faces_2& partition_faces_2() const {
                return m_partition_faces_2;
            }

            // Triangulation.
            inline Triangulation& triangulation() {
                return m_triangulation;
            }

            inline const Triangulation& triangulation() const {
                return m_triangulation;
            }

            // Buildings.
            inline Buildings& buildings() {
                return m_buildings;
            }

            inline const Buildings& buildings() const {
                return m_buildings;
            }

        private:
            const Input_range &m_input_range;
            Point_map   m_point_map;
            SemanticMap m_semantic_map;
            VisibilityMap m_visibility_map;

            Plane_3   m_ground_plane;
            Polygon_3 m_ground_bounding_box;

            std::vector<Point_2> m_filtered_building_boundary_points;
            std::vector<Point_2> m_simplified_building_boundary_points;

            Detected_regions        m_detected_2d_regions;
            Regularized_segments    m_regularized_segments;

            Partition_faces_2 m_partition_faces_2;
            Triangulation     m_triangulation;

            Buildings m_buildings;
        };
    
    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_DATA_STRUCTURE_H
