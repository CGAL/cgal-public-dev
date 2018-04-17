#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_ESTIMATOR_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Simple_cartesian.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Buildings/Level_of_detail_diagonalize_traits.h>

namespace CGAL {

	namespace LOD {

		// Main class.
		template<class KernelTraits, class ContainerInput, class BuildingsInput, class EstimationStrategy>
		class Level_of_detail_building_roof_estimator {
            
        public:
            typedef KernelTraits       Kernel;
            typedef ContainerInput     Input;
            typedef BuildingsInput     Buildings;
            typedef EstimationStrategy Strategy;

            using FT      = typename Kernel::FT;
            using Point_3 = typename Kernel::Point_3;
            using Plane_3 = typename Kernel::Plane_3;

            using Local_kernel       = CGAL::Simple_cartesian<double>;
            using Diagonalize_traits = CGAL::LOD::Eigen_diagonalize_traits_lod<double, 3>;

			using Point_3ft    = Local_kernel::Point_3;
			using Plane_3ft    = Local_kernel::Plane_3;

            using Building          = typename Strategy::Building;
            using Building_iterator = typename Buildings::iterator;

            using Index   = int;
			using Indices = std::vector<Index>;

            Level_of_detail_building_roof_estimator(const Input &input) : 
            m_input(input), m_strategy(input), m_alpha(-FT(1)) { }

            void estimate(Buildings &buildings) {
                assert(buildings.size() > 0);

                std::cerr << "This class is turned off due to its very constrained application!" << std::endl;
                exit(1);

				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {

                    auto &building = bit->second;
					process_building(building);
                }
            }

            void set_alpha(const FT new_value) {
                assert(new_value > FT(0));
                m_alpha = new_value;
            }

            bool is_face_based() const {
				return m_strategy.is_face_based();
			}

        private:
            const Input &m_input;
            Strategy  m_strategy;
            FT        m_alpha;

            void process_building(Building &building) {

                const auto &shapes = building.shapes;
                if (shapes.size() == 0) return;
					
				building.clear_roofs();
				for (size_t i = 0; i < shapes.size(); ++i) {
					
                    const Indices &indices = shapes[i];
                    process_roof(indices, building);
                }
            }

            void process_roof(const Indices &indices, Building &building) {
                if (indices.size() < 3) return;

                Plane_3 plane;
                fit_plane_to_roof_points(indices, plane);

                std::vector<Point_3> points;
                
                const FT points_height = project_points_onto_plane(indices, plane, points);
                const FT height_difference = get_translation(points_height, building);

                translate_points(height_difference, points);

                m_strategy.set_alpha(m_alpha);
                m_strategy.estimate_roof(points, plane, building);
            }

            void fit_plane_to_roof_points(const Indices &indices, Plane_3 &plane) const {
                
                assert(indices.size() > 2);
                double bx = 0.0, by = 0.0, bz = 0.0;

                std::vector<Point_3ft> points(indices.size());
				for (size_t i = 0; i < indices.size(); ++i) {

					const Point_3 &p = m_input.point(indices[i]);

					const double x = CGAL::to_double(p.x());
					const double y = CGAL::to_double(p.y());
					const double z = CGAL::to_double(p.z());

					points[i] = Point_3ft(x, y, z);

                    bx += x;
                    by += y;
                    bz += z;
				}

                bx /= static_cast<double>(indices.size());
                by /= static_cast<double>(indices.size());
                bz /= static_cast<double>(indices.size());

				Plane_3ft tmp_plane;
                Point_3ft centroid = Point_3ft(bx, by, bz);

				CGAL::linear_least_squares_fitting_3(points.begin(), points.end(), tmp_plane, centroid, CGAL::Dimension_tag<0>(), Local_kernel(), Diagonalize_traits());
				plane = Plane_3(static_cast<FT>(tmp_plane.a()), static_cast<FT>(tmp_plane.b()), static_cast<FT>(tmp_plane.c()), static_cast<FT>(tmp_plane.d()));
            }

            FT project_points_onto_plane(const Indices &indices, const Plane_3 &plane, std::vector<Point_3> &points) const {
                assert(indices.size() > 2);

                points.clear();
                points.resize(indices.size());

                FT points_height = FT(100000000000000);
                for (size_t i = 0; i < indices.size(); ++i) {			
					
                    const Point_3 &p = m_input.point(indices[i]);
					points[i] = plane.projection(p);

                    points_height = CGAL::min(points_height, points[i].z());
                }

                return points_height;
            }

            FT get_translation(const FT points_height, const Building &building) const {
                
                const FT building_height = building.height;
                return points_height - building_height;
            }

            void translate_points(const FT height_difference, std::vector<Point_3> &points) const {

				for (size_t i = 0; i < points.size(); ++i) {

					Point_3 &p = points[i];
					const FT z = p.z() - height_difference;

					points[i] = Point_3(p.x(), p.y(), z);
				}
			}
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_ESTIMATOR_H