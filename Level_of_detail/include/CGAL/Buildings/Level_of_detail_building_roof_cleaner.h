#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_CLEANER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_CLEANER_H

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
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/number_utils.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class InputCDT, class InputBuildings>
		class Level_of_detail_building_roof_cleaner {

        public:
            typedef KernelTraits   Kernel;
            typedef InputContainer Input;
            typedef InputCDT       CDT;
            typedef InputBuildings Buildings;

			typename Kernel::Compute_squared_length_3 		  squared_length;
			typename Kernel::Compute_scalar_product_3 		  dot_product;
			typename Kernel::Construct_cross_product_vector_3 cross_product;

            using FT       = typename Kernel::FT;
            using Point_3  = typename Kernel::Point_3;
            using Vector_3 = typename Kernel::Vector_3;
            using Plane_3  = typename Kernel::Plane_3;

            using Vertex_handle = typename CDT::Vertex_handle;
            using Face_handle   = typename CDT::Face_handle;

            using Building          = CGAL::LOD::Building<Kernel, Vertex_handle, Face_handle>;
            using Building_iterator = typename Buildings::iterator;
            
            using Log     = CGAL::LOD::Mylog;
            using Indices = std::vector<int>;
            using Points  = std::vector<Point_3>;

            using Shapes        = typename Building::Shapes;
            using Shape_indices = typename Building::Indices;
            using Heights       = std::map<size_t, FT>;

            using Search_traits   = CGAL::Search_traits_3<Kernel>;
			using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Search_traits>;
			using Fuzzy_sphere    = CGAL::Fuzzy_sphere<Search_traits>;
			using Tree            = typename Neighbor_search::Tree;

            Level_of_detail_building_roof_cleaner(const Input &input, const FT ground_height) : 
            m_input(input), 
            m_ground_height(ground_height), 
            m_silent(false),
            m_scale_upper_bound(-FT(1)),
            m_max_percentage(FT(80)),
            m_angle_threshold(FT(5)),
            m_apply_size_criteria(true),
            m_apply_height_criteria(false),
            m_apply_vertical_criteria(true),
            m_apply_scale_based_criteria(true) { }

            void clean_shapes(Buildings &buildings) const {
                
                if (buildings.size() == 0) return;
                for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {
                    
                    Building &building = (*bit).second;
                    if (building.is_valid) clean_building_roofs(building);
                }

                if (!m_silent) {
                    Log exporter; exporter.export_shapes_inside_buildings(buildings, m_input, "tmp" + std::string(PSR) + "lod_2" + std::string(PSR) + "filtered_roof_shapes");
                }
            }

            void make_silent(const bool new_state) {
                m_silent = new_state;
            }

            void set_scale_upper_bound(const FT new_value) {
                
                assert(new_value > FT(0));
                m_scale_upper_bound = new_value;
            }

            void set_max_percentage(const FT new_value) {
                
                assert(new_value >= FT(0) && new_value <= FT(100));
                m_max_percentage = new_value;
            }

        private:
            const Input &m_input;
            const FT m_ground_height;

            bool m_silent;

            FT m_scale_upper_bound;
            FT m_max_percentage;

            const FT m_angle_threshold;
            
            const bool m_apply_size_criteria;
            const bool m_apply_height_criteria;
            const bool m_apply_vertical_criteria;
            const bool m_apply_scale_based_criteria;

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

                if (m_apply_size_criteria)        apply_size_criteria(shapes, indices);
                if (m_apply_scale_based_criteria) apply_scale_based_criteria(shapes, indices);
                if (m_apply_vertical_criteria)    apply_vertical_criteria(shapes, indices);
                if (m_apply_height_criteria)      apply_height_criteria(shapes, indices);
                

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
                
                assert(m_max_percentage >= FT(0) && m_max_percentage <= FT(100));

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

                assert(m_max_percentage >= FT(0) && m_max_percentage <= FT(100));

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
            
            void apply_vertical_criteria(const Shapes &shapes, Indices &indices) const {

                Indices new_indices;
                for (size_t i = 0; i < indices.size(); ++i) {

                    const size_t index = indices[i];
                    if (!is_vertical_shape(shapes[index])) new_indices.push_back(index);
                }
                indices = new_indices;
            }

            bool is_vertical_shape(const Shape_indices &shape_indices) const {

				Vector_3 shape_normal;
				set_shape_normal(shape_indices, shape_normal);

				Vector_3 ground_normal;
				set_ground_normal(ground_normal);

                const FT angle      = compute_angle(shape_normal, ground_normal);
                const FT angle_diff = FT(90) - CGAL::abs(angle);

                if (angle_diff < m_angle_threshold) return true;
                return false;
            }

            void set_shape_normal(const Shape_indices &shape_indices, Vector_3 &m) const {
				
                FT x = FT(0), y = FT(0), z = FT(0);
                for (size_t i = 0; i < shape_indices.size(); ++i) {
                    
                    const auto index = shape_indices[i];
                    const Vector_3 &normal = m_input.normal(index);

                    x += normal.x();
                    y += normal.y();
                    z += normal.z();
                }

                x /= static_cast<FT>(shape_indices.size());
                y /= static_cast<FT>(shape_indices.size());
                z /= static_cast<FT>(shape_indices.size());

                m = Vector_3(x, y, z);
			}

			void set_ground_normal(Vector_3 &n) const {
				
                const Plane_3 ground = Plane_3(FT(0), FT(0), FT(1), FT(0));				
				n = ground.orthogonal_vector();
			}

            FT compute_angle(const Vector_3 &m, const Vector_3 &n) const {
				
				const auto cross = cross_product(m, n);
				const FT length  = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_length(cross))));
				const FT dot     = dot_product(m, n);

				const FT angle_rad = static_cast<FT>(std::atan2(CGAL::to_double(length), CGAL::to_double(dot)));
                const FT angle_deg = angle_rad * FT(180) / static_cast<FT>(CGAL_PI);
                
                return angle_deg;
			}

            void apply_scale_based_criteria(const Shapes &shapes, Indices &indices) const {

                Indices new_indices;
                for (size_t i = 0; i < indices.size(); ++i) {

                    const size_t index = indices[i];
                    if (!is_within_scale_bounds(shapes[index])) new_indices.push_back(index);
                }
                indices = new_indices;
            }

            bool is_within_scale_bounds(const Shape_indices &shape_indices) const {
                
                assert(m_scale_upper_bound > FT(0));

                Point_3 barycentre;
                compute_barycentre(shape_indices, barycentre);

                Points points;
                set_points(shape_indices, points);

                Tree tree(points.begin(), points.end());
                const Fuzzy_sphere sphere(barycentre, m_scale_upper_bound);

                Points result;
                tree.search(std::back_inserter(result), sphere);

                if (result.size() == points.size()) return true;
                return false;
            }

            void compute_barycentre(const Shape_indices &shape_indices, Point_3 &barycentre) const {
                
                FT x = FT(0), y = FT(0), z = FT(0);
                for (size_t i = 0; i < shape_indices.size(); ++i) {
                    const Point_3 &p = m_input.point(shape_indices[i]);

                    x += p.x();
                    y += p.y();
                    z += p.z();
                }

                x /= static_cast<FT>(shape_indices.size());
                y /= static_cast<FT>(shape_indices.size());
                z /= static_cast<FT>(shape_indices.size());

                barycentre = Point_3(x, y, z);
            }

            void set_points(const Shape_indices &shape_indices, Points &points) const {

                points.clear();
                points.resize(shape_indices.size());

                for (size_t i = 0; i < shape_indices.size(); ++i) 
                    points[i] = m_input.point(shape_indices[i]);
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

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_CLEANER_H