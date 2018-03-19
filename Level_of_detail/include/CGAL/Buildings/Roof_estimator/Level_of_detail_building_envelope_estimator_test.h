#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_ESTIMATOR_TEST_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_ESTIMATOR_TEST_H

// STL includes.
#include <map>
#include <list>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Cartesian.h>
#include <CGAL/envelope_3.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Env_plane_traits_3.h>
#include <CGAL/Env_triangle_traits_3.h>
#include <CGAL/Env_surface_data_traits_3.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		// Main class.
		template<class KernelTraits, class ContainerInput, class BuildingsInput>
		class Level_of_detail_building_envelope_estimator_test {
            
        public:
            typedef KernelTraits       Kernel;
            typedef ContainerInput     Input;
            typedef BuildingsInput     Buildings;

            typedef CGAL::Exact_rational         Number_type;
            typedef CGAL::Cartesian<Number_type> CKernel;
            typedef typename CKernel::Point_3    CPoint_3;
            typedef typename CKernel::Plane_3    CPlane_3;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Point_3 = typename Kernel::Point_3;
            using Plane_3 = typename Kernel::Plane_3;

            using Local_kernel = CGAL::Simple_cartesian<double>;
			using Point_3ft    = Local_kernel::Point_3;
			using Plane_3ft    = Local_kernel::Plane_3;

            using Building_iterator = typename Buildings::iterator;

            using Index   = int;
			using Indices = std::vector<Index>;

            using Envelope_traits  = CGAL::Env_plane_traits_3<CKernel>;
            using Surface          = typename Envelope_traits::Surface_3;
            using Envelope_diagram = CGAL::Envelope_diagram_2<Envelope_traits>;

            using Planes = std::list<Surface>;
			using Shapes = std::vector<Indices>;

            using Face_iterator       = typename Envelope_diagram::Face_const_iterator;
            using Halfedge_circulator = typename Envelope_diagram::Ccb_halfedge_const_circulator;
            using Surface_iterator    = typename Envelope_diagram::Surface_const_iterator;

            using Points = std::vector<Point_3>;

            Level_of_detail_building_envelope_estimator_test(const Input &input) : 
            m_input(input), m_add_example(true) { }

            void estimate(Buildings &buildings) {
                assert(buildings.size() > 0);

				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {

                    auto &building = bit->second;
					process_building(building);
                }

                if (m_add_example) print_example();
            }

            void set_alpha(const FT) { }

            bool is_face_based() const {
				return false;
			}

        private:
            const Input &m_input;
            const bool m_add_example;

            void print_example() const {

                typedef CGAL::Exact_rational                              Number_type_ex;
                typedef CGAL::Cartesian<Number_type_ex>                   Kernel_ex;
                typedef CGAL::Env_triangle_traits_3<Kernel_ex>            Traits_3;
                typedef typename Kernel_ex::Point_3                       Point_3_ex;
                typedef typename Traits_3::Surface_3                      Triangle_3;
                typedef CGAL::Env_surface_data_traits_3<Traits_3, char>   Data_traits_3;
                typedef typename Data_traits_3::Surface_3                 Data_triangle_3;
                typedef CGAL::Envelope_diagram_2<Data_traits_3>           Envelope_diagram_2;

                std::cout << std::endl << "Example: " << std::endl << std::endl;

                // Construct the input triangles, makred A and B.
                std::list<Data_triangle_3> triangles;
                triangles.push_back(Data_triangle_3(Triangle_3(Point_3_ex(0, 0, 0),
                                                               Point_3_ex(0, 6, 0),
                                                               Point_3_ex(5, 3, 4)), 'A'));

                triangles.push_back(Data_triangle_3(Triangle_3(Point_3_ex(6, 0, 0),
                                                               Point_3_ex(6, 6, 0),
                                                               Point_3_ex(1, 3, 4)), 'B'));
                
                // Compute and print the minimization diagram.
                Envelope_diagram_2 diag;
                CGAL::lower_envelope_3(triangles.begin(), triangles.end(), diag);

                Envelope_diagram_2::Face_const_iterator            fit;
                Envelope_diagram_2::Ccb_halfedge_const_circulator  ccb;
                Envelope_diagram_2::Surface_const_iterator         sit;

                for (fit = diag.faces_begin(); fit != diag.faces_end(); ++fit) {
                    
                    // Print the face boundary.
                    if (fit->is_unbounded()) {
                        std::cout << "[Unbounded face]";
                    } else {
                    
                        // Print the vertices along the outer boundary of the face.
                        ccb = fit->outer_ccb(); std::cout << "[Face]  ";
                        do {
                            std::cout << '(' << CGAL::to_double(ccb->target()->point().x()) << " " << CGAL::to_double(ccb->target()->point().y()) << ")  ";
                            ++ccb;
                        } while (ccb != fit->outer_ccb());
                    }

                    // Print the labels of the triangles that induce the envelope on this face.
                    std::cout << "-->  " << fit->number_of_surfaces() << " triangles:";
                    for (sit = fit->surfaces_begin(); sit != fit->surfaces_end(); ++sit) std::cout << ' ' << sit->data();
                    
                    std::cout << std::endl << std::endl;
                }
            }

            template<class Building>
            void process_building(Building &building) const {

                Planes planes;

                add_roofs(building, planes);
                add_walls(building, planes);
                
                CPoint_3 minb, maxb;
                compute_bounding_box(building, minb, maxb);
                add_box(minb, maxb, planes);

                Envelope_diagram min_diag;
                CGAL::lower_envelope_3(planes.begin(), planes.end(), min_diag);

                get_back_result(min_diag, building);
            }

            template<class Building>
            void add_roofs(const Building &building, Planes &planes) const {
                
                const auto &shapes = building.shapes;
                if (shapes.size() == 0) return;

                CPlane_3 plane;
				for (size_t i = 0; i < shapes.size(); ++i) {
                    
                    const Indices &indices = shapes[i];
                    if (indices.size() < 3) continue;    

                    fit_plane_to_roof_points(indices, plane);
                    planes.push_back(Surface(plane));
                }
            }

            void fit_plane_to_roof_points(const Indices &indices, CPlane_3 &plane) const {
                assert(indices.size() > 2);

                std::vector<Point_3ft> points(indices.size());
				for (size_t i = 0; i < indices.size(); ++i) {

					const Point_3 &p = m_input.point(indices[i]);

					const double x = CGAL::to_double(p.x());
					const double y = CGAL::to_double(p.y());
					const double z = CGAL::to_double(p.z());

					points[i] = Point_3ft(x, y, z);
				}

				Plane_3ft tmp_plane;
				CGAL::linear_least_squares_fitting_3(points.begin(), points.end(), tmp_plane, CGAL::Dimension_tag<0>());
				plane = CPlane_3(tmp_plane.a(), tmp_plane.b(), tmp_plane.c(), tmp_plane.d());
            }

            template<class Building>
            void add_walls(const Building &building, Planes &planes) const {

                const auto &boundary = building.boundaries[0];
                if (boundary.size() == 0) return;

                const size_t num_vertices = boundary.size();
				for (size_t i = 0; i < num_vertices; i += 2) {
					
					const size_t ip = i + 1;
					add_wall_plane(boundary[i], boundary[ip], FT(0), FT(10), planes);
				}
            }

            template<class Vertex_handle>
            void add_wall_plane(const Vertex_handle &v1, const Vertex_handle &v2, const FT height_base, const FT height_roof, Planes &planes) const {

                const Point_2 &va = v1->point();
				const Point_2 &vb = v2->point();

				const Point_3 a = Point_3(va.x(), va.y(), height_base);
				const Point_3 b = Point_3(vb.x(), vb.y(), height_base);
				const Point_3 c = Point_3(vb.x(), vb.y(), height_roof);

                planes.push_back(Surface(CPlane_3(CPoint_3(a.x(), a.y(), a.z()), CPoint_3(b.x(), b.y(), b.z()), CPoint_3(c.x(), c.y(), c.z()))));
            }

            template<class Building>
            void get_back_result(const Envelope_diagram &diag, Building &building) const {
                building.clear_roofs();

                Points boundary;
                for (Face_iterator fit = diag.faces_begin(); fit != diag.faces_end(); ++fit) {
                    boundary.clear();

                    if (fit->is_unbounded()) continue;

                    Halfedge_circulator ccb = fit->outer_ccb();
                    do {

                        const double x = CGAL::to_double(ccb->target()->point().x());
                        const double y = CGAL::to_double(ccb->target()->point().y());

                        if (x > 10000 || y > 10000 || x < -10000 || y < -10000) {
                         
                            boundary.clear();
                            break;
                        }

                        const Point_2 p = Point_2(static_cast<FT>(x), static_cast<FT>(y));
                        boundary.push_back(Point_3(p.x(), p.y(), FT(100)));

                        ++ccb;
                    } while (ccb != fit->outer_ccb());

                    using Roof = typename Building::Roof;
                    Roof  roof;

                    roof.boundary = boundary;
                    building.roofs.push_back(roof);
                }
            }

            template<class Building>
            void compute_bounding_box(const Building &building, CPoint_3 &minb, CPoint_3 &maxb) const {
                
                const auto &boundary = building.boundaries[0];
                if (boundary.size() == 0) {
                 
                    minb = CPoint_3( 0,  0,  0);
                    maxb = CPoint_3(10, 10, 10);

                    return;
                }

                const FT big_value = FT(100000000000000);
                FT minx =  big_value, miny =  big_value;
                FT maxx = -big_value, maxy = -big_value;

                const size_t num_vertices = boundary.size();
				for (size_t i = 0; i < num_vertices; ++i) {

                    minx = CGAL::min(minx, boundary[i]->point().x());
                    miny = CGAL::min(miny, boundary[i]->point().y());

                    maxx = CGAL::max(maxx, boundary[i]->point().x());
                    maxy = CGAL::max(maxy, boundary[i]->point().y());
                }

                minb = CPoint_3(minx, miny, 0);
                maxb = CPoint_3(maxx, maxy, building.height + 20);
            }

            void add_box(const CPoint_3 &minb, const CPoint_3 &maxb, Planes &planes) const {

                const CPoint_3 p1 = CPoint_3(minb.x(), minb.y(), 0);
                const CPoint_3 p2 = CPoint_3(maxb.x(), minb.y(), 0);
                const CPoint_3 p3 = CPoint_3(maxb.x(), maxb.y(), 0);
                const CPoint_3 p4 = CPoint_3(minb.x(), maxb.y(), 0);

                const CPoint_3 p5 = CPoint_3(minb.x(), minb.y(), 50);
                const CPoint_3 p6 = CPoint_3(maxb.x(), minb.y(), 50);
                const CPoint_3 p7 = CPoint_3(maxb.x(), maxb.y(), 50);
                const CPoint_3 p8 = CPoint_3(minb.x(), maxb.y(), 50);

                planes.push_back(CPlane_3(p1, p2, p3));
                planes.push_back(CPlane_3(p5, p6, p7));
                
                planes.push_back(CPlane_3(p1, p2, p6));
                planes.push_back(CPlane_3(p3, p4, p8));

                planes.push_back(CPlane_3(p2, p3, p7));
                planes.push_back(CPlane_3(p4, p1, p5));
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_ESTIMATOR_TEST_H