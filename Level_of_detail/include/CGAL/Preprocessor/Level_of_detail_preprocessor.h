#ifndef CGAL_LEVEL_OF_DETAIL_PREPROCESSOR_H
#define CGAL_LEVEL_OF_DETAIL_PREPROCESSOR_H

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


			Level_of_detail_interior_boundary_extractor() : m_alpha(-FT(1)) { }

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
				Log log;
				log.export_points("tmp/extracted_boundary_points", extracted);


				// Return number of extracted points.
				assert(number_of_extracted_points <= static_cast<int>(projected.size()));
				return number_of_extracted_points;
			}

		private:
			FT m_alpha;
			Simple_utils m_simple_utils;

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


		template<class KernelTraits, class InputContainer>
		class Level_of_detail_preprocessor {

		public:
			typedef KernelTraits   Traits;
			typedef InputContainer Container;

			typedef typename Traits::FT 	 FT;
			typedef typename Traits::Point_2 Point_2;
			typedef typename Traits::Point_3 Point_3;
			typedef typename Traits::Plane_3 Plane_3;

			typedef CGAL::Search_traits_3<Traits>                     Search_traits;
			typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> Neighbor_search;
			typedef CGAL::Fuzzy_sphere<Search_traits>                 Fuzzy_sphere;
			typedef typename Neighbor_search::Tree 					  Tree;

			typedef std::map<int, Point_2> 												  Projected_points;
			typedef Level_of_detail_simple_projector<Traits, Container, Projected_points> Ground_projector;

			using Index          = int;
			using Const_iterator = typename Container::const_iterator;
			using Index_map      = typename Container:: template Property_map<Index>;
			using Log 			 = CGAL::LOD::Mylog;

			Level_of_detail_preprocessor() : m_scale(FT(1)), m_alpha(-FT(1)) { }

			void set_scale(const FT new_value) {
				
				assert(new_value >= FT(0));
				m_scale = new_value;
			}

			void set_alpha(const FT new_value) {

				assert(new_value > FT(0));
				m_alpha = new_value;
			}

			template<class Planes>
			auto get_planes(const Container &input, Planes &planes) {

				auto number_of_planes = -1;
				create_indices(input);
				
				planes.clear();

				for (Const_iterator it = input.begin(); it != input.end(); ++it)
					if (m_indices[*it] >= 0) 
						planes[m_indices[*it]].push_back(*it);

				number_of_planes = planes.size();
				return number_of_planes;
			}

			// Later I can adapt boundary_clutter to some other data structure where I also take
			// into account the type of the clutter: undetected, cylinder, sphere and so on. The type is saved in property_map<Types>.
			template<class Indices, class Boundary_data>
			auto get_boundary_points(
				const Container &input, const Indices &boundary_mapping, const Indices &interior_mapping, const bool with_shape_detection, 
				Boundary_data &building_boundaries, Boundary_data &boundary_clutter) {

				boundary_clutter.clear();
				create_indices(input);

				add_interior_boundary_to_clutter(input, interior_mapping, boundary_clutter);
				return add_boundary_points(boundary_mapping, with_shape_detection, building_boundaries, boundary_clutter);
			}

			// This is very slow algorithm. Should be improved later.
			template<class Projected_points, class Point_sets>
			auto clean_projected_points(Projected_points &projected_points, Point_sets &point_sets) {

				if (projected_points.size() == 0) return 0;
				assert(!point_sets.empty());

				Projected_points cleaned_points;
				Point_sets cleaned_sets;

				const auto num_neighbours = 2;

				CGAL::Identity_property_map<Point_3> pmap;
				for (typename Point_sets::const_iterator it = point_sets.begin(); it != point_sets.end(); ++it) {
					const auto set_index = (*it).first;

					const size_t num_points = (*it).second.size();
					std::vector<Point_3> points(num_points);

					for (size_t i = 0; i < num_points; ++i) {
						const auto point_index = (*it).second[i];

						const Point_2 &p = projected_points.at(point_index);
						points[i] = Point_3(p.x(), p.y(), FT(0));
					}
						
					const FT average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(points.begin(), points.end(), pmap, num_neighbours, Traits());
					Tree tree(points.begin(), points.end());

					std::vector<Point_3> qs;
					for (size_t i = 0; i < num_points; ++i) {
						const auto point_index = (*it).second[i];
						
						const Point_3 &p = points[i];
						const Fuzzy_sphere sphere(p, m_scale * average_spacing);

						qs.clear();	
						tree.search(std::back_inserter(qs), sphere);

						if (qs.empty() || (qs.size() == 1 && qs[0] == p)) continue;

						cleaned_points[point_index] = Point_2(p.x(), p.y());
						cleaned_sets[set_index].push_back(point_index);
					}
				}

				const auto number_of_outliers = static_cast<int>(projected_points.size() - cleaned_points.size());

				projected_points = cleaned_points;
				point_sets = cleaned_sets;

				return number_of_outliers;
			}

		private:
			Index_map m_indices;

			FT m_scale;
			FT m_alpha;

			void create_indices(const Container &input) {
				boost::tie(m_indices, boost::tuples::ignore) = input. template property_map<Index>("index");
			}

			template<class Indices, class Boundary_data>
			auto add_boundary_points(const Indices &mapping, const bool with_shape_detection, Boundary_data &building_boundaries, Boundary_data &boundary_clutter) {

				auto number_of_boundaries = -1;
				building_boundaries.clear();

				if (mapping.empty()) return number_of_boundaries;

				for (size_t i = 0; i < mapping.size(); ++i) {
					if (m_indices[mapping[i]] >= 0 && with_shape_detection) building_boundaries[m_indices[mapping[i]]].push_back(mapping[i]);
					else boundary_clutter[0].push_back(mapping[i]);
				}

				number_of_boundaries = building_boundaries.size();
				return number_of_boundaries;
			}


			template<class Indices, class Boundary_data>
			void add_interior_boundary_to_clutter(const Container &input, const Indices &interior_mapping, Boundary_data &boundary_clutter) {
				if (interior_mapping.empty()) return;

				// (1) Project points onto the ground.
				Projected_points building_interior_projected;
				project_points_onto_ground<Indices, Boundary_data>(input, interior_mapping, building_interior_projected);


				// (2) Extract boundaries.
				Level_of_detail_interior_boundary_extractor<Traits, Container, Indices> extractor;

				Indices result;
				extractor.set_alpha(m_alpha);
				extractor.extract(input, interior_mapping, building_interior_projected, result);


				// (3) Add extracted points to the clutter.
				Boundary_data stub;
				add_boundary_points(result, false, stub, boundary_clutter);
			}

			template<class Indices, class Boundary_data>
			void project_points_onto_ground(const Container &input, const Indices &interior_mapping, Projected_points &building_interior_projected) {
				if (interior_mapping.empty()) return;

				Boundary_data building_interior, stub;
				add_boundary_points(interior_mapping, false, stub, building_interior);

				Plane_3 base_ground_plane(FT(0), FT(0), FT(1), FT(0));
				building_interior_projected.clear();
				
				Ground_projector projector; 
				projector.project(input, building_interior, base_ground_plane, building_interior_projected);

				Log log;
				log.export_projected_points_as_xyz("tmp/building_interior", building_interior_projected, "stub");
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_PREPROCESSOR_H