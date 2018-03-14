#ifndef CGAL_LEVEL_OF_DETAIL_REGION_GROWING_3_H
#define CGAL_LEVEL_OF_DETAIL_REGION_GROWING_3_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\"
#else 
#define PSR "/" 
#endif 

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>
#include <unordered_set>

// CGAL includes.
#include <CGAL/Point_set_3.h>
#include <CGAL/number_utils.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Region_growing/Region_growing_simon_3/Plane.h>
#include <CGAL/Region_growing/Region_growing_simon_3/Shape_base.h>
#include <CGAL/Region_growing/Region_growing_simon_3/Region_growing.h>
#include <CGAL/Region_growing/Region_growing_simon_3/Shape_detection_traits.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class InputCDT, class InputBuildings>
		class Level_of_detail_region_growing_3 {

		public:
            typedef KernelTraits   Kernel;
            typedef InputContainer Input;
            typedef InputCDT       CDT;
            typedef InputBuildings Buildings;

            using Index   = int;
			using Indices = std::vector<Index>;

            using FT       = typename Kernel::FT;
            using Point_3  = typename Kernel::Point_3;
            using Vector_3 = typename Kernel::Vector_3;

            using Vertex_handle   = typename CDT::Vertex_handle;
            using Face_handle     = typename CDT::Face_handle;

            using Building          = CGAL::LOD::Building<FT, Vertex_handle, Face_handle, Point_3>;
            using Building_iterator = typename Buildings::iterator;

            using Log = CGAL::LOD::Mylog;

            using Point_set  = CGAL::Point_set_3<Point_3>;
            using Point_map  = typename Point_set::Point_map;
            using Vector_map = typename Point_set::Vector_map;

            using Traits         = CGAL::Shape_detection_simon_3::Shape_detection_traits<Kernel, Point_set, Point_map, Vector_map>;
            using Region_growing = CGAL::Shape_detection_simon_3::Region_growing<Traits>;
            using Plane          = CGAL::Shape_detection_simon_3::Plane<Traits>;

            using Parameters     = typename CGAL::Shape_detection_simon_3::Region_growing<Traits>::Parameters;
            using Shape_iterator = typename Region_growing::Shape_range::iterator;
            using Shape          = typename Region_growing::Shape;
            using Index_iterator = std::vector<size_t>::const_iterator;

            Level_of_detail_region_growing_3(const Input &input) : 
            m_input(input), 
            m_epsilon(-FT(1)),
			m_cluster_epsilon(-FT(1)),
			m_normal_threshold(-FT(1)),
			m_min_points(0),
            m_silent(false) { }

            void detect(Buildings &buildings) const {
                
                assert(buildings.size() > 0);
                for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {
                    
                    Building &building = (*bit).second;
                    grow_regions(building);
                }

                if (!m_silent) {
                    Log exporter; exporter.export_shapes_inside_buildings(buildings, m_input, "tmp" + std::string(PSR) + "shapes_inside_buildings");
                }
            }
            
            void make_silent(const bool new_state) {
                m_silent = new_state;
            }

            void set_epsilon(const FT new_value) {

				assert(new_value > FT(0));
				m_epsilon = new_value;
			}

			void set_cluster_epsilon(const FT new_value) {

				assert(new_value > FT(0));
				m_cluster_epsilon = new_value;
			}

			void set_normal_threshold(const FT new_value) {

				assert(new_value > FT(0) && new_value < FT(1));
				m_normal_threshold = new_value;
			}

			void set_minimum_shape_points(const size_t new_value) {

				assert(new_value > 0);
				m_min_points = new_value;
			}

        private:
            const Input &m_input;
            bool  m_silent;

            FT m_epsilon;
            FT m_cluster_epsilon;
            FT m_normal_threshold;
            size_t m_min_points;

            void grow_regions(Building &building) const {

                const Indices &indices = building.interior_indices;
                if (indices.size() == 0) return;

                Point_set points;
                set_points(indices, points);

                apply_3d_region_growing(indices, points, building);
            }

            void set_points(const Indices &indices, Point_set &points) const {
                points.clear();
                points.add_normal_map();

                for (size_t i = 0; i < indices.size(); ++i) {
                    const Index index = indices[i];

                    const Point_3  &point  = m_input.point(index);
                    const Vector_3 &normal = m_input.normal(index);

                    points.insert(point, normal);
                }
            }

            void apply_3d_region_growing(const Indices &indices, Point_set &points, Building &building) const {

                Region_growing region_growing;

                region_growing.set_input(points, points.point_map(), points.normal_map());
                region_growing. template add_shape_factory<Plane>();

                assert(m_epsilon          > FT(0));
                assert(m_cluster_epsilon  > FT(0));
                assert(m_normal_threshold > FT(0) && m_normal_threshold < FT(1));
                assert(m_min_points       > 0);

                Parameters parameters;
                set_parameters(parameters);

                region_growing.detect(parameters);
                set_shapes_to_building(region_growing, indices, building);
            }

            void set_parameters(Parameters &parameters) const {

                parameters.epsilon          = m_epsilon / FT(4);
                parameters.cluster_epsilon  = m_cluster_epsilon;
                parameters.normal_threshold = m_normal_threshold;
                parameters.min_points       = m_min_points * FT(6);
            }

            void set_shapes_to_building(const Region_growing &region_growing, const Indices &indices, Building &building) const {
                
                building.clear_shapes();
                const size_t number_of_shapes = static_cast<size_t>(region_growing.shapes().end() - region_growing.shapes().begin());

                building.shapes.resize(number_of_shapes);
                Shape_iterator sit = region_growing.shapes().begin();

                size_t count = 0;
                while (sit != region_growing.shapes().end()) {
                    boost::shared_ptr<Shape> shape = *sit;

                    const size_t number_of_indices = static_cast<size_t>(shape->indices_of_assigned_points().end() - shape->indices_of_assigned_points().begin());
                    building.shapes[count].resize(number_of_indices);

                    Index_iterator index_it = shape->indices_of_assigned_points().begin(); size_t i = 0;
                    while (index_it != shape->indices_of_assigned_points().end()) {
                        const auto index = *index_it;

                        assert(index >= 0);
                        assert(index < indices.size());

                        building.shapes[count][i] = indices[index];
                           
                        ++index_it;
                        ++i;
                    }
                    ++sit;
                    ++count;
                }

                building.clear_interior_indices();
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_REGION_GROWING_3_H