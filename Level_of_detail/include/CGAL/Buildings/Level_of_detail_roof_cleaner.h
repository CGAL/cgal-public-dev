#ifndef CGAL_LEVEL_OF_DETAIL_ROOF_CLEANER_H
#define CGAL_LEVEL_OF_DETAIL_ROOF_CLEANER_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\" 
#else 
#define PSR "/" 
#endif 

// STL includes.
#include <map>
#include <vector>
#include <cassert>
#include <algorithm>

// CGAL includes.
#include <CGAL/number_utils.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class InputCDT, class InputBuildings>
		class Level_of_detail_roof_cleaner {

        public:
            typedef KernelTraits   Kernel;
            typedef InputContainer Input;
            typedef InputCDT       CDT;
            typedef InputBuildings Buildings;

            using FT       = typename Kernel::FT;
            using Point_3  = typename Kernel::Point_3;

            using Vertex_handle   = typename CDT::Vertex_handle;
            using Face_handle     = typename CDT::Face_handle;

            using Building          = CGAL::LOD::Building<FT, Vertex_handle, Face_handle, Point_3>;
            using Building_iterator = typename Buildings::iterator;
            
            using Log = CGAL::LOD::Mylog;
            using Indices = std::vector<int>;

            using Shapes = typename Building::Shapes;
            using Shape_indices = typename Building::Indices;
            using Heights = std::map<size_t, FT>;

            Level_of_detail_roof_cleaner(const Input &input, const FT ground_height) : 
            m_input(input), 
            m_ground_height(ground_height), 
            m_silent(false),
            m_max_percentage(FT(80)),
            m_apply_size_criteria(true),
            m_apply_height_criteria(true) { }

            void clean_shapes(Buildings &buildings) const {
                
                if (buildings.size() == 0) return;
                for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {
                    
                    Building &building = (*bit).second;
                    if (building.is_valid) clean_building_roofs(building);
                }

                if (!m_silent) {
                    Log exporter; exporter.export_shapes_inside_buildings(buildings, m_input, "tmp" + std::string(PSR) + "inside_buildings_filtered_shapes");
                }
            }

            void make_silent(const bool new_state) {
                m_silent = new_state;
            }

        private:
            const Input &m_input;
            const FT m_ground_height;

            bool  m_silent;
            const FT m_max_percentage;

            const bool m_apply_size_criteria;
            const bool m_apply_height_criteria;

            class Size_comparator {
                
                public:
                    Size_comparator(const Shapes &shapes) : m_shapes(shapes) { }
                
                    bool operator() (const size_t i, const size_t j) const { 

                        assert(i >= 0 && i < m_shapes.size());
                        assert(j >= 0 && j < m_shapes.size());

                        return m_shapes[i].size() > m_shapes[j].size();
                    }

                private:
                    const Shapes &m_shapes;
            };

            class Height_comparator {
                
                public:
                    Height_comparator(const Heights &heights) : m_heights(heights) { }
                
                    bool operator() (const size_t i, const size_t j) const { 
                        return m_heights.at(i) < m_heights.at(j);
                    }

                private:
                    const Heights &m_heights;
            };

            void clean_building_roofs(Building &building) const {

                Shapes &shapes = building.shapes;
                const size_t num_shapes = shapes.size();

                Indices indices;
                set_default_indices(indices, num_shapes);

                if (m_apply_size_criteria)   apply_size_criteria(shapes, indices);
                if (m_apply_height_criteria) apply_height_criteria(shapes, indices);

                update_shapes(indices, shapes);
                if (shapes.size() == 0) building.is_valid = false;
            }

            void set_default_indices(Indices &indices, const size_t num_indices) const {
                
                indices.clear();
                indices.resize(num_indices);

                for (size_t i = 0; i < num_indices; ++i)
                    indices[i] = i;
            }

            void apply_size_criteria(const Shapes &shapes, Indices &indices) const {
                
                sort_indices_by_size(shapes, indices);
                remove_outliers(shapes, m_max_percentage, indices);
            }

            void sort_indices_by_size(const Shapes &shapes, Indices &indices) const {

                Size_comparator size_comparator(shapes);
                std::sort(indices.begin(), indices.end(), size_comparator);
            }

            void remove_outliers(const Shapes &shapes, const FT percentage, Indices &indices) const {
                
                const size_t num_total_points = get_total_number_of_points(shapes, indices);
                const FT scale = percentage / FT(100);

                const size_t num_points_to_keep = static_cast<size_t>(std::ceil(CGAL::to_double(scale * static_cast<FT>(num_total_points))));
                
                size_t curr_num_points = 0;
                for (size_t i = 0; i < indices.size(); ++i) {
                    
                    const size_t index = indices[i];
                    curr_num_points += shapes[index].size();

                    if (curr_num_points >= num_points_to_keep) {
                        
                        indices.erase(indices.begin() + i + 1, indices.end());
                        break;
                    }
                }
            }

            size_t get_total_number_of_points(const Shapes &shapes, const Indices &indices) const {

                size_t num_total_points = 0;
                for (size_t i = 0; i < indices.size(); ++i)
                    num_total_points += shapes[indices[i]].size();

                return num_total_points;
            }

            void apply_height_criteria(const Shapes &shapes, Indices &indices) const {

                Heights heights;
                compute_heights(shapes, indices, heights);

                sort_indices_by_height(heights, indices);
                remove_outliers(shapes, m_max_percentage - FT(20), indices);
            }

            void compute_heights(const Shapes &shapes, const Indices &indices, Heights &heights) const {

                heights.clear();
                for (size_t i = 0; i < indices.size(); ++i)
                    heights[indices[i]] = compute_height(shapes[indices[i]]);
            }

            FT compute_height(const Shape_indices &shape_indices) const {
                return compute_min_height(shape_indices);
            }

            FT compute_min_height(const Shape_indices &shape_indices) const {

                const FT ground_z = m_ground_height;
                FT min_height = FT(100000000000000);

                for (size_t i = 0; i < shape_indices.size(); ++i) {
                    const Point_3 &p = m_input.point(shape_indices[i]);

                    const FT height = p.z() - ground_z;
                    min_height = CGAL::min(min_height, height);
                }

                return min_height;
            }

            void sort_indices_by_height(const Heights &heights, Indices &indices) const {

                Height_comparator height_comparator(heights);
                std::sort(indices.begin(), indices.end(), height_comparator);
            }

            void update_shapes(const Indices &indices, Shapes &shapes) const {

                Shapes new_shapes(indices.size());
                for (size_t i = 0; i < indices.size(); ++i)
                    if (shapes[indices[i]].size() > 2) 
                        new_shapes[i] = shapes[indices[i]];
                shapes = new_shapes;
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_ROOF_CLEANER_H