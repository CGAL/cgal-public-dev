#ifndef CGAL_LEVEL_OF_DETAIL_ESTIMATED_NORMAL_PROPERTY_MAP_2_H
#define CGAL_LEVEL_OF_DETAIL_ESTIMATED_NORMAL_PROPERTY_MAP_2_H

// STL includes.
#include <map>
#include <list>

// CGAL includes.
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>
#include <CGAL/Search_traits_2.h>

// LOD includes.
#include <CGAL/Level_of_detail/Tools/Fitters/Line_to_points_fitter.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

		template<typename KeyType, class InputKernel, class InputElements, class PointMap>
		class Estimated_normal_property_map_2 {

		public:
            using Kernel    = InputKernel;
            using Elements  = InputElements;
			using Point_map = PointMap;
            
            using FT        = typename Kernel::FT;
            using Point_2   = typename Kernel::Point_2;
            using Vector_2  = typename Kernel::Vector_2;
            using Line_2    = typename Kernel::Line_2;
            using ValueType = Vector_2;

            using Normals    = std::map<KeyType, Vector_2>;
            using Neighbours = std::list<Point_2>;

			using Search_traits_2 = CGAL::Search_traits_2<Kernel>;
			using Search_circle   = CGAL::Fuzzy_sphere<Search_traits_2>;
			using Search_tree     = CGAL::Kd_tree<Search_traits_2>;

            using Identity_point_map    = CGAL::Identity_property_map<Point_2>;
            using Line_to_points_fitter = LOD::Line_to_points_fitter<Kernel>;

			Estimated_normal_property_map_2(const Elements &elements, const Point_map &point_map, const FT local_search_radius) :
			m_elements(elements),
            m_point_map(point_map),
            m_local_search_radius(local_search_radius) { 

                estimate_normals();
            }

            const Normals &normals() const {
                return m_normals;
            }

		private:
            using key_type   = KeyType;
            using value_type = ValueType;
            using reference  = const value_type&;
            
            using Self = Estimated_normal_property_map_2<key_type, Kernel, Elements, Point_map>;

            friend reference get(const Self &self, const key_type &key) {
				return self.normals().at(key);
            }

            friend void put(const Self &, key_type &, const value_type &) { }
			
            const Elements  &m_elements;
            const Point_map &m_point_map;
            
            Normals  m_normals;
            const FT m_local_search_radius;

            void estimate_normals() {
                
                CGAL_precondition(m_elements.size() > 0);
                m_normals.clear();

				Search_tree tree;
                create_tree(tree);
                
                for (typename Elements::const_iterator element = m_elements.begin(); element != m_elements.end(); ++element) {
                    
                    const Point_2 &point = get(m_point_map, *element);
                    estimate_normal(tree, point, *element);
                }
            }

            void create_tree(Search_tree &tree) const {

                tree.clear();
                for (typename Elements::const_iterator element = m_elements.begin(); element != m_elements.end(); ++element) {
                    
                    const Point_2 &point = get(m_point_map, *element);
                    tree.insert(point);
                }
            }

            void estimate_normal(const Search_tree &tree, const Point_2 &query, const key_type &key) {

				Neighbours neighbours;
				Search_circle circle(query, m_local_search_radius);
				tree.search(std::back_inserter(neighbours), circle);

                Line_2 line;
                Identity_point_map identity_point_map;

                Line_to_points_fitter line_to_points_fitter;
                line_to_points_fitter.fit_line_2(neighbours, identity_point_map, line);
				
				const Vector_2 vector = line.to_vector();
				Vector_2 normal       = vector.perpendicular(CGAL::COUNTERCLOCKWISE);
				const FT length       = static_cast<FT>(CGAL::sqrt(CGAL::to_double(normal * normal)));

				CGAL_precondition(length != FT(0));
				normal /= length;
                m_normals[key] = normal;
            }
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ESTIMATED_NORMAL_PROPERTY_MAP_2_H