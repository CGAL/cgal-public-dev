#ifndef CGAL_LEVEL_OF_DETAIL_INTERIOR_BOUNDARY_EXTRACTOR_H
#define CGAL_LEVEL_OF_DETAIL_INTERIOR_BOUNDARY_EXTRACTOR_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\" 
#else 
#define PSR "/" 
#endif 

// STL includes.
#include <vector>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/constructions_d.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/property_map.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Simple_cartesian.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Projector/Level_of_detail_projector.h>
#include <CGAL/Utils/Level_of_detail_utils.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class InputIndices>
		class Level_of_detail_interior_boundary_extractor {
		
		public:
			typedef KernelTraits   Kernel;
			typedef InputContainer Container;
			typedef InputIndices   Indices;

			typedef typename Kernel::FT 	   FT;
			typedef typename Kernel::Point_2   Point_2;
			typedef typename Kernel::Point_3   Point_3;
			typedef typename Kernel::Plane_3   Plane_3;
			typedef typename Kernel::Segment_2 Segment_2;

			typedef CGAL::Alpha_shape_vertex_base_2<Kernel> Vb;
			typedef CGAL::Alpha_shape_face_base_2<Kernel>   Fb;

			typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
			typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>  Triangulation_2;
			typedef CGAL::Alpha_shape_2<Triangulation_2>         Alpha_shape_2;

			typedef typename Alpha_shape_2::Alpha_shape_vertices_iterator Alpha_vertex_iterator;
			typedef typename Alpha_shape_2::Alpha_shape_edges_iterator    Alpha_edge_iterator;

			typedef Level_of_detail_utils_simple<Kernel> Simple_utils;


			using Index            = int;
			using Const_iterator   = typename Container::const_iterator;
			using Index_map        = typename Container:: template Property_map<Index>;
			using Log 			   = CGAL::LOD::Mylog;
			using Projected_points = std::map<int, Point_2>;
			using Point_iterator   = typename Projected_points::const_iterator;


			Level_of_detail_interior_boundary_extractor() : m_alpha(-FT(1)), m_silent(false) { }

			void make_silent(const bool new_state) {
				m_silent = new_state;
			}

			void set_alpha(const FT new_value) {

				assert(new_value > FT(0));
				m_alpha = new_value;
			}


			int extract(const Container & /* input */, const Indices & /* mapping */, const Projected_points &projected, Indices &result) {

				// Some preconditions.
				assert(result.empty());
				int number_of_extracted_points = -1;

				// assert(input.number_of_points() != 0);
				// if (mapping.empty()) return number_of_extracted_points;

				if (projected.empty()) return number_of_extracted_points;


				// Set default parameters.
				set_default_parameters(projected);


				// Create temporary points. Optimize/remove this function!
				std::vector<Point_2> points;
				create_points(projected, points);


				// Extract boundary points using alpha shapes.
				std::vector<Point_2> extracted;
				number_of_extracted_points = apply_alpha_shape(points, extracted);


				// Map points to indices. Optimize/remove this function!
				map_extracted_points_to_indices(projected, extracted, result);


				// Save extracted points.
				if (!m_silent) {
					Log log;
					log.export_points("tmp" + std::string(PSR) + "extracted_boundary_points_2d", extracted);
				}


				// Return number of extracted points.
				assert(number_of_extracted_points <= static_cast<int>(projected.size()));
				return number_of_extracted_points;
			}

		private:
			FT m_alpha;
			Simple_utils m_simple_utils;
			bool m_silent;

			void set_default_parameters(const Projected_points &projected_points) {

				set_default_alpha(projected_points);
			}

			void set_default_alpha(const Projected_points &projected_points) {

				if (m_alpha > FT(0)) return;

				Point_2 bbmin, bbmax;
				m_simple_utils.compute_bounding_box_in_2d(bbmin, bbmax, projected_points);

				const FT diagonal = m_simple_utils.compute_2d_bounding_box_diagonal(bbmin, bbmax);
				m_alpha = diagonal / FT(10);

				assert(m_alpha > FT(0));
			}

			void create_points(const Projected_points &projected_points, std::vector<Point_2> &points) {
				assert(!projected_points.empty());

				points.clear();
				points.resize(projected_points.size());

				size_t count = 0;
				for (Point_iterator pit = projected_points.begin(); pit != projected_points.end(); ++pit, ++count) {

					const Point_2 &p = (*pit).second;
					points[count] = p;
				}
				assert(count == projected_points.size());
			}

			int apply_alpha_shape(const std::vector<Point_2> &points, std::vector<Point_2> &extracted) {

				assert(!points.empty());
				extracted.clear();

				assert(m_alpha > FT(0));
				Alpha_shape_2 alpha_shape(points.begin(), points.end(), m_alpha, Alpha_shape_2::GENERAL);

				const size_t number_of_extracted_points = std::distance(alpha_shape.alpha_shape_vertices_begin(), alpha_shape.alpha_shape_vertices_end());
				extracted.resize(number_of_extracted_points);

				size_t count = 0;
				for (Alpha_vertex_iterator ait = alpha_shape.alpha_shape_vertices_begin(); ait != alpha_shape.alpha_shape_vertices_end(); ++ait, ++count) {

					const Point_2 &p = (*ait)->point();
					extracted[count] = p;
				}

				assert(count == number_of_extracted_points);
				return number_of_extracted_points;
			}

			void map_extracted_points_to_indices(
			const Projected_points &projected, const std::vector<Point_2> &extracted, Indices &mapping) {

				assert(!projected.empty());
				const size_t num_extracted = extracted.size();

				mapping.clear();
				mapping.resize(num_extracted);

				for (size_t i = 0; i < num_extracted; ++i) {
					const Point_2 &p = extracted[i];

					bool found = false;
					for (Point_iterator pit = projected.begin(); pit != projected.end(); ++pit) {
						
						const Point_2 &q = (*pit).second;
						if (p == q) {

							found = true;
							mapping[i] = (*pit).first;

							break;
						}
					}
					assert(found);
				}
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_INTERIOR_BOUNDARY_EXTRACTOR_H