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

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer>
		class Level_of_detail_preprocessor {

		public:
			typedef KernelTraits   Traits;
			typedef InputContainer Container;

			typedef typename Traits::FT 	 FT;
			typedef typename Traits::Point_2 Point_2;
			typedef typename Traits::Point_3 Point_3;

			typedef CGAL::Search_traits_3<Traits>                     Search_traits;
			typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> Neighbor_search;
			typedef CGAL::Fuzzy_sphere<Search_traits>                 Fuzzy_sphere;
			typedef typename Neighbor_search::Tree 					  Tree;

			using Index          = int;
			using Const_iterator = typename Container::const_iterator;
			using Index_map      = typename Container:: template Property_map<Index>;

			Level_of_detail_preprocessor() : m_scale(FT(1)) { }

			void set_scale(const FT new_value) {
				m_scale = new_value;
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
			auto get_boundaries(const Container &input, const Indices &mapping, Boundary_data &building_boundaries, Boundary_data &boundary_clutter) {

				auto number_of_boundaries = -1;
				create_indices(input);

				building_boundaries.clear();
				boundary_clutter.clear();

				for (size_t i = 0; i < mapping.size(); ++i) {
					if (m_indices[mapping[i]] >= 0) building_boundaries[m_indices[mapping[i]]].push_back(mapping[i]);
					else boundary_clutter[0].push_back(mapping[i]);
				}

				number_of_boundaries = building_boundaries.size();
				return number_of_boundaries;
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

			void create_indices(const Container &input) {
				boost::tie(m_indices,  boost::tuples::ignore) = input. template property_map<Index>("index");
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_PREPROCESSOR_H