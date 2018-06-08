#ifndef CGAL_LEVEL_OF_DETAIL_PARAMETERS_H
#define CGAL_LEVEL_OF_DETAIL_PARAMETERS_H

// STL include.
#include <string>

namespace CGAL {

	namespace Level_of_detail {

		template<typename NumberType>
		struct Parameters {

		public:
			using FT = NumberType;

			Parameters() :
			m_path_to_input("default_path"),
			m_silent(false),
			m_verbose(false),
			m_no_simplification(false),
			m_no_regularization(false),
			m_scale(FT(5)),
			m_epsilon(FT(16) / FT(5)),
			m_region_growing_2_normal_threshold(FT(7) / FT(10)),
			m_region_growing_2_min_points(10),
			m_segment_regularizer_2_max_angle_in_degrees(FT(25)),
			m_kinetic_partitioning_2_num_intersections(2) { 
				
				update_dependent();
			}


			//////////////////////////////////
			// Functions to be not documented:

			inline std::string& path_to_input() {
				return m_path_to_input;
			}

			inline const std::string& path_to_input() const {
				return m_path_to_input;
			}

			inline bool& silent() {
				return m_silent;
			}

			inline const bool& silent() const {
				return m_silent;
			}


			//////////////////////////////////
			// Functions to be documented:


			// Flags.
			inline bool& verbose() {
				return m_verbose;
			}

			inline const bool& verbose() const {
				return m_verbose;
			}

			inline bool& no_simplification() {
				return m_no_simplification;
			}

			inline const bool& no_simplification() const {
				return m_no_simplification;
			}

			inline bool& no_regularization() {
				return m_no_regularization;
			}

			inline const bool& no_regularization() const {
				return m_no_regularization;
			}


			// The main parameters.
			inline FT& scale() {
				return m_scale;
			}

			inline const FT& scale() const {
				return m_scale;
			}

			inline FT& epsilon() {
				return m_epsilon;
			}

			inline const FT& epsilon() const {
				return m_epsilon;
			}


			// Filtering.
			inline FT& alpha_shape_size() {
				return m_alpha_shape_size;
			}

			inline const FT& alpha_shape_size() const {
				return m_alpha_shape_size;
			}

			inline FT& grid_cell_width() {
				return m_grid_cell_width;
			}

			inline const FT& grid_cell_width() const {
				return m_grid_cell_width;
			}


			// Points based region growing in 2D.
			inline FT& region_growing_2_normal_threshold() {
				return m_region_growing_2_normal_threshold;
			}

			inline const FT& region_growing_2_normal_threshold() const {
				return m_region_growing_2_normal_threshold;
			}

			inline size_t& region_growing_2_min_points() {
				return m_region_growing_2_min_points;
			}

			inline const size_t& region_growing_2_min_points() const {
				return m_region_growing_2_min_points;
			}

			inline FT& region_growing_2_epsilon() {
				return m_region_growing_2_epsilon;
			}

			inline const FT& region_growing_2_epsilon() const {
				return m_region_growing_2_epsilon;
			}

			inline FT& region_growing_2_cluster_epsilon() {
				return m_region_growing_2_cluster_epsilon;
			}

			inline const FT& region_growing_2_cluster_epsilon() const {
				return m_region_growing_2_cluster_epsilon;
			}


			// Segment regularizer in 2D.
			inline FT& segment_regularizer_2_max_angle_in_degrees() {
				return m_segment_regularizer_2_max_angle_in_degrees;
			}

			inline const FT& segment_regularizer_2_max_angle_in_degrees() const {
				return m_segment_regularizer_2_max_angle_in_degrees;
			}

			inline FT& segment_regularizer_2_max_difference_in_meters() {
				return m_segment_regularizer_2_max_difference_in_meters;
			}

			inline const FT& segment_regularizer_2_max_difference_in_meters() const {
				return m_segment_regularizer_2_max_difference_in_meters;
			}


			// Kinetic based partitioning in 2D.
			inline size_t& kinetic_partitioning_2_num_intersections() {
				return m_kinetic_partitioning_2_num_intersections;
			}

			inline const size_t& kinetic_partitioning_2_num_intersections() const {
				return m_kinetic_partitioning_2_num_intersections;
			}

			inline FT& kinetic_partitioning_2_min_face_width() {
				return m_kinetic_partitioning_2_min_face_width;
			}

			inline const FT& kinetic_partitioning_2_min_face_width() const {
				return m_kinetic_partitioning_2_min_face_width;
			}


			// Automatically defined parameters.
			void update_dependent() {
				m_alpha_shape_size = m_scale;
				m_grid_cell_width  = FT(26) * m_scale / FT(100);

				m_region_growing_2_epsilon 		   = m_epsilon;
				m_region_growing_2_cluster_epsilon = FT(58) * m_scale / FT(100);

				m_segment_regularizer_2_max_difference_in_meters = m_scale / FT(4);

				m_kinetic_partitioning_2_min_face_width = m_scale / FT(2);
			}

		private:
			std::string m_path_to_input;
			
			bool m_silent;
			bool m_verbose;
			bool m_no_simplification;
			bool m_no_regularization;

			FT m_scale;
			FT m_epsilon;
			
			FT m_alpha_shape_size;
			FT m_grid_cell_width;

			FT 	   m_region_growing_2_normal_threshold;
			size_t m_region_growing_2_min_points;
			FT 	   m_region_growing_2_epsilon;
			FT 	   m_region_growing_2_cluster_epsilon;
			
			FT m_segment_regularizer_2_max_angle_in_degrees;
			FT m_segment_regularizer_2_max_difference_in_meters;

			size_t m_kinetic_partitioning_2_num_intersections;
			FT     m_kinetic_partitioning_2_min_face_width;
		};
	
	} // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_PARAMETERS_H