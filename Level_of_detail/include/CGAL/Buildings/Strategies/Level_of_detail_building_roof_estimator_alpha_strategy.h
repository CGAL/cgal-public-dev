#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_ESTIMATOR_ALPHA_STRATEGY_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_ESTIMATOR_ALPHA_STRATEGY_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		// Main class.
		template<class KernelTraits, class ContainerInput, class BuildingInput>
		class Level_of_detail_building_roof_estimator_alpha_strategy {
            
        public:
            typedef KernelTraits   Kernel;
            typedef ContainerInput Input;
            typedef BuildingInput  Building;

            using FT       = typename Kernel::FT;
            using Point_2  = typename Kernel::Point_2;
			using Point_3  = typename Kernel::Point_3;
            using Plane_3  = typename Kernel::Plane_3;
            using Points   = std::vector<Point_3>;

            using Roof = typename Building::Roof;

			using Vbwi 			= CGAL::Triangulation_vertex_base_with_info_2<int, Kernel>;
			using Vb   			= CGAL::Alpha_shape_vertex_base_2<Kernel, Vbwi>;
			using Fb   			= CGAL::Alpha_shape_face_base_2<Kernel>;
			using Tds 		 	= CGAL::Triangulation_data_structure_2<Vb, Fb>;
			using Dt   			= CGAL::Delaunay_triangulation_2<Kernel, Tds>;
			using Alpha_shape_2 = CGAL::Alpha_shape_2<Dt>;

			using Face_iterator = typename Alpha_shape_2::Finite_faces_iterator;
			using Vertex_handle = typename Dt::Vertex_handle;

            Level_of_detail_building_roof_estimator_alpha_strategy(const Input &input) : m_input(input), m_alpha(-FT(1)) { }

            void estimate_roof(const Points &roof_points, const Plane_3 &, Building &building) const {                
				if (roof_points.size() < 3) return;

				// Set triangulation.
				Dt dt;
				set_delaunay_triangulation(roof_points, dt);

				// Apply alpha shapes
				assert(m_alpha > FT(0));
				Alpha_shape_2 alpha_shape(dt, m_alpha, Alpha_shape_2::GENERAL);

				// Get roof elements.
				Points elements;
				get_roof_elements(alpha_shape, roof_points, elements);

				// Set estimated roof.
				Roof roof;
                roof.boundary = elements;
                building.roofs.push_back(roof);
            }

			void set_alpha(const FT new_value) { 
				assert(new_value > FT(0));
				m_alpha = new_value;
			}

			bool is_face_based() const {
				return true;
			}

			std::string name() const {
				return "alpha";
			}

        private:
            const Input &m_input;
			FT m_alpha;

			void set_delaunay_triangulation(const Points &roof_points, Dt &dt) const {
				assert(roof_points.size() > 2);

				for (size_t i = 0; i < roof_points.size(); ++i) {
					const Point_3 &p = roof_points[i];

					Vertex_handle vh = dt.insert(Point_2(p.x(), p.y()));
					vh->info() = static_cast<int>(i);
				}
			}

			void get_roof_elements(const Alpha_shape_2 &alpha_shape, const Points &roof_points, Points &elements) const {
				assert(elements.size() == 0);

				for (Face_iterator fit = alpha_shape.finite_faces_begin(); fit != alpha_shape.finite_faces_end(); ++fit) {
					if (alpha_shape.classify(fit) != Alpha_shape_2::INTERIOR) continue;

					for (size_t i = 0; i < 3; ++ i) {
						
						const Point_2 &p   = fit->vertex(i)->point();
						const size_t index = fit->vertex(i)->info();

						elements.push_back(Point_3(p.x(), p.y(), roof_points[index].z()));
					}
				}
			}
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_ESTIMATOR_ALPHA_STRATEGY_H