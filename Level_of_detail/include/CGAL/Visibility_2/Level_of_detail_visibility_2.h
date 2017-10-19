#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_2_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_2_H

// STL includes.
#include <utility>
#include <cassert>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/Barycentric_coordinates_2.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <CGAL/Random.h>
#include <CGAL/number_utils.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/constructions_d.h>
#include <CGAL/utils.h>
#include <CGAL/IO/Color.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class ContainerInput, class CDTInput>
		class Level_of_detail_visibility_2 {

		public:
			typedef KernelTraits 		Kernel;
			typedef typename Kernel::FT FT;

			typedef ContainerInput Container;
			typedef CDTInput 	   CDT;
			
			virtual void save_info(const bool) = 0;
			virtual void show_progress(const bool) = 0;

			virtual int compute(const Container &, CDT &) = 0;

			virtual ~Level_of_detail_visibility_2() { }
		};

		// This is the ray shooting based classification algorithm.
		template<class KernelTraits, class ContainerInput, class CDTInput>
		class Level_of_detail_visibility_ray_shooting_2 : public Level_of_detail_visibility_2<KernelTraits, ContainerInput, CDTInput> {

		public:
			// Typedefs.
			typedef Level_of_detail_visibility_2<KernelTraits, ContainerInput, CDTInput> Base;
			
			typedef typename Base::Kernel 	   Kernel;
			typedef typename Base::FT     	   FT;
			typedef typename Kernel::Point_2   Point_2;
			typedef typename Kernel::Vector_2  Vector_2;
			typedef typename Kernel::Ray_2     Ray_2;
			typedef typename Kernel::Segment_2 Segment_2;

			typedef typename Base::Container Container;			
			typedef typename Base::CDT 		 CDT;

			typedef typename CDT::Vertex_handle 	   Vertex_handle;
			typedef typename CDT::Face_handle   	   Face_handle;
			typedef typename CDT::Line_face_circulator Line_face_circulator;


			// Extra.
			using Log     = CGAL::LOD::Mylog;
			using Bbox    = std::vector<Point_2>;
			using Samples = std::vector<Point_2>;
			using Ray     = Segment_2;
			using Rays    = std::vector<Ray>;

			using Cdt_vertex_iterator = typename CDT::Finite_vertices_iterator;
			using Cdt_face_iterator   = typename CDT::Finite_faces_iterator;

			using Sample_iterator = typename Samples::const_iterator;
			using Ray_iterator    = typename Rays::const_iterator;


			Level_of_detail_visibility_ray_shooting_2() : 
			m_save_info(true), 
			m_show_progress(true),
			m_num_samples(3),
			m_num_rays_per_side(2) { }

			void save_info(const bool new_state) override {
				m_save_info = new_state;
			}

			void show_progress(const bool new_state) override {
				m_show_progress = new_state;
			}

			void set_number_of_samples(const size_t new_value) {
				m_num_samples = new_value;
			}

			void set_number_of_rays_per_side(const size_t new_value) {
				m_num_rays_per_side = new_value;
			}

			int compute(const Container &, CDT &cdt) override {

				Log log;

				// (1) Compute bounding box.
				std::vector<Point_2> bbox;
				compute_bounding_box(cdt, log, bbox);

				if (m_save_info) log.out << "(1) Bounding box is computed!" << std::endl << std::endl;

				// (2) Compute visibility.
				compute_visibility(cdt, log, bbox);

				if (m_save_info) {
					
					log.out << "\n(2) Visibility is computed!" << std::endl;
					log.save("tmp/visibility_ray_shooting");
				}

				return static_cast<int>(cdt.number_of_faces());
			}

		private:
			bool m_save_info;
			bool m_show_progress;
			size_t m_num_samples;
			size_t m_num_rays_per_side;

			CGAL::Random m_rand;

			void compute_bounding_box(const CDT &cdt, Log &log, Bbox &bbox) {

				bbox.clear();
				bbox.resize(4);

				const FT big_value = FT(1000000); // change it

				FT minx =  big_value, miny =  big_value;
				FT maxx = -big_value, maxy = -big_value;

				for (typename CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
					const Point_2 &p = vit->point();

					const FT x = p.x();
					const FT y = p.y();

					minx = CGAL::min(minx, x);
					miny = CGAL::min(miny, y);

					maxx = CGAL::max(maxx, x);
					maxy = CGAL::max(maxy, y);
				}

				const FT min_length = CGAL::min(CGAL::abs(maxx - minx), CGAL::abs(maxy - miny));
				const FT threshold  = min_length / FT(4); 

				bbox[0] = Point_2(minx - threshold, miny - threshold);
				bbox[1] = Point_2(maxx + threshold, miny - threshold);
				bbox[2] = Point_2(maxx + threshold, maxy + threshold);
				bbox[3] = Point_2(minx - threshold, maxy + threshold);

				if (m_save_info) {
					
					log.out << "Bounding box: ";
					log.skip_line();

					log.out << bbox[0] << "; ";
					log.out << bbox[1] << "; ";
					log.out << bbox[2] << "; ";
					log.out << bbox[3] << "; ";

					log.skip_line();
					log.skip_line();
				}
			}

			void compute_visibility(CDT &cdt, Log &log, const Bbox &bbox) {

				// Can be removed -->
				size_t num_stages = 1;
				size_t stage      = 1;

				if (m_show_progress) {

					num_stages = 100;
					stage = cdt.number_of_faces() / num_stages;

					std::cout << std::endl;
					std::cout << "Progress: " << std::endl;
				}

				size_t progress = 0; // <--

				size_t face_index = 0;
				for (Cdt_face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); 
					++fit, 
					++face_index, 
					++progress) {

					if (m_show_progress) { // can be removed
						
						if (stage != 0 && progress % stage == 0)
							std::cout << (progress * 100) / cdt.number_of_faces() << "% " << std::endl; 
					}

					if (m_save_info) log.out << "\nFace: " << face_index << std::endl;

					const Face_handle fh = static_cast<Face_handle>(fit);

					fh->info().in       = estimate_face_visibility(log, cdt, fh, bbox, face_index);
					fh->info().in_color = get_color(fh->info().in);
				}

				if (m_show_progress) { // can be removed

					std::cout << "100%" << std::endl;
					std::cout << std::endl;
				}
			}

			FT estimate_face_visibility(Log &log, const CDT &cdt, const Face_handle &fh, const Bbox &bbox, const size_t /* face_index */ ) {

				// (a) Generate samples.
				Samples samples(m_num_samples);
				generate_samples(cdt, fh, samples);

				/* // Debugging.
				Log saver; saver.save_triangle_with_points_eps(
					cdt.triangle(fh).vertex(0), 
					cdt.triangle(fh).vertex(1), 
					cdt.triangle(fh).vertex(2), samples, "tmp/triangle_" + std::to_string(face_index)); */

				if (m_save_info) log.out << "(a) Samples are generated." << std::endl;


				// (b) Handle all samples.
				const FT face_visibility = estimate_face_visibility_from_samples(log, cdt, fh, bbox, samples);
				if (m_save_info) log.out << "\n(b) Face visibility is estimated from the given samples: " << face_visibility << std::endl << std::endl;


				return face_visibility;
			}

			FT estimate_face_visibility_from_samples(Log &log, const CDT &cdt, const Face_handle &fh, const Bbox &bbox, const Samples &samples) {

				FT face_visibility = FT(0);
				for (Sample_iterator sit = samples.begin(); sit != samples.end(); ++sit) {

					if (m_save_info) log.out << "\nNew sample: " << std::endl;

					const Point_2 &sample = *sit;
					face_visibility += estimate_visibility_from_sample(log, sample, cdt, fh, bbox);
				}

				return face_visibility / static_cast<FT>(m_num_samples);
			}

			FT estimate_visibility_from_sample(Log &log, const Point_2 &sample, const CDT &cdt, const Face_handle &fh, const Bbox &bbox) {

				assert(bbox.size() == 4);
				FT sample_visibility = FT(0);


				// Handle bottom side of the bounding box.
				if (m_save_info) log.out << "\nBottom side: " << std::endl;
				sample_visibility += get_side_visibility(log, sample, cdt, fh, bbox[0], bbox[1]);


				// Handle right side of the bounding box.
				if (m_save_info) log.out << "\nRight side: "  << std::endl;
				sample_visibility += get_side_visibility(log, sample, cdt, fh, bbox[1], bbox[2]);


				// Handle top side of the bounding box.
				if (m_save_info) log.out << "\nTop side: "    << std::endl;
				sample_visibility += get_side_visibility(log, sample, cdt, fh, bbox[2], bbox[3]);


				// Handle left side of the bounding box.
				if (m_save_info) log.out << "\nLeft side: "   << std::endl;
				sample_visibility += get_side_visibility(log, sample, cdt, fh, bbox[3], bbox[0]);


				if (m_save_info) log.out << std::endl;
				return sample_visibility / static_cast<FT>(4);
			}

			FT get_side_visibility(Log &log, const Point_2 &sample, const CDT &cdt, const Face_handle &fh, const Point_2 &a, const Point_2 &b) {

				if (m_save_info) log.out << "side: (" << a << ") -- (" << b << ")" << std::endl;


				// Generate rays.
				Rays rays(m_num_rays_per_side);
				generate_side_rays(sample, a, b, rays);

				if (m_save_info) log.out << "rays are generated" << std::endl;


				// Handle all rays.
				return get_visibility_from_rays(log, rays, cdt, fh) / static_cast<FT>(m_num_rays_per_side);
			}

			void generate_side_rays(const Point_2 &source, const Point_2 &a, const Point_2 &b, Rays &rays) {

				assert(m_num_rays_per_side > 0);
				assert(rays.size() == m_num_rays_per_side);

				for (size_t i = 0; i < m_num_rays_per_side; ++i) {

					const FT s = static_cast<FT>(m_rand.uniform_real(0.0, 1.0));
					const FT t = FT(1) - s;

					const FT x = s * a.x() + t * b.x();
					const FT y = s * a.y() + t * b.y();

					const Point_2 target = Point_2(x, y);

					rays[i] = Ray(source, target);
				}
			}

			FT get_visibility_from_rays(Log &log, const Rays &rays, const CDT &cdt, const Face_handle &fh) {

				if (m_save_info) log.out << std::endl;

				FT rays_visibility = FT(0);
				for (Ray_iterator rit = rays.begin(); rit != rays.end(); ++rit) {
					const Ray &ray = *rit;

					if (m_save_info) log.out << "New ray: " << ray;
					rays_visibility += get_visibility_from_ray(log, ray, cdt, fh);
				}

				if (m_save_info) log.out << std::endl;
				return rays_visibility;
			}

			FT get_visibility_from_ray(Log &log, const Ray &ray, const CDT &cdt, const Face_handle &fh) {

				const Point_2 &p = ray.source();
				const Point_2 &q = ray.target();

				const Line_face_circulator lfc = cdt.line_walk(p, q);

				assert(lfc != NULL);
				const FT ray_visibility = traverse_ray_faces(log, lfc, cdt, fh);

				if (m_save_info) log.out << " ray visibility: " << ray_visibility << std::endl;
				return ray_visibility;
			}

			FT traverse_ray_faces(Log &log, const Line_face_circulator &lfc, const CDT &cdt, const Face_handle &) {

				if (m_save_info) log.out << std::endl << " first face: " << cdt.triangle(lfc) << std::endl;

				// to be implemented!

				return FT(1) / FT(2);
			}

			void generate_samples(const CDT &cdt, const Face_handle &fh, Samples &samples) {				
				generate_samples_uniform_0(cdt, fh, samples);
			}

			void generate_samples_uniform_0(const CDT &cdt, const Face_handle &fh, Samples &samples) {

				assert(!samples.empty() && samples.size() == m_num_samples);

				const Point_2 &a = cdt.triangle(fh).vertex(0);
				const Point_2 &b = cdt.triangle(fh).vertex(1);
				const Point_2 &c = cdt.triangle(fh).vertex(2);

				Vector_2 ab = b - a;
				Vector_2 ac = c - a;

				for (size_t i = 0; i < m_num_samples; ++i) {

					FT s = static_cast<FT>(m_rand.uniform_real(0.0, 1.0));
					FT t = static_cast<FT>(m_rand.uniform_real(0.0, 1.0));

					if (s + t > FT(1)) {

						s = FT(1) - s;
						t = FT(1) - t;
					}
					samples[i] = a + s * ab + t * ac;
				}
			}

			CGAL::Color get_color(const FT visibility) {

				const FT half  = FT(1) / FT(2);
				
				FT scale_in  = FT(1) - visibility;
				FT scale_out = visibility;

				const FT scale_thres = FT(2) / FT(10);

				if (scale_in  < scale_thres) scale_in  = scale_thres;
				if (scale_out < scale_thres) scale_out = scale_thres;

				if (visibility > half)      return CGAL::Color(51 * scale_in, 255 * scale_in, 51 * scale_in);  	 // INSIDE
				else if (visibility < half) return CGAL::Color(255 * scale_out, 51 * scale_out, 51 * scale_out); // OUTSIDE
									  
				return CGAL::Color(255, 204, 0); // UNKNOWN
			}
		};

		// This class works only with the xy aligned ground plane that is Plane(0, 0, 1, 0).
		template<class KernelTraits, class ContainerInput, class CDTInput>
		class Level_of_detail_visibility_from_classification_2 : public Level_of_detail_visibility_2<KernelTraits, ContainerInput, CDTInput> {

		public:
			typedef Level_of_detail_visibility_2<KernelTraits, ContainerInput, CDTInput> Base;

			typedef typename Base::Kernel 	  Kernel;
			typedef typename Base::FT     	  FT;
			typedef typename Kernel::Point_2  Point_2;
			typedef typename Kernel::Vector_2 Vector_2;

			typedef typename Base::Container  Container;			
			typedef typename Base::CDT 		  CDT;
			
			typedef typename CDT::Vertex_handle Vertex_handle;
			typedef typename CDT::Face_handle   Face_handle;

			typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<Kernel> 							 Triangle_coordinates;
			typedef CGAL::Barycentric_coordinates::Mean_value_2<Kernel> 									 Mean_value;
			typedef CGAL::Barycentric_coordinates::Generalized_barycentric_coordinates_2<Mean_value, Kernel> Mean_value_coordinates;

			typedef int Label;
			typedef typename std::pair<Point_2, Label> 							Point_with_label;
			typedef typename CGAL::First_of_pair_property_map<Point_with_label> Point_map;

			typedef CGAL::Search_traits_2<Kernel>                       					  Search_traits_2;
			typedef CGAL::Search_traits_adapter<Point_with_label, Point_map, Search_traits_2> Search_traits;
			typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> 						  Neighbor_search;
			typedef typename Neighbor_search::Tree 											  Tree;
			typedef CGAL::Fuzzy_sphere<Search_traits>                    					  Fuzzy_circle;
			typedef CGAL::Kd_tree<Search_traits>					                          Fuzzy_tree;

			typedef CGAL::Delaunay_triangulation_2<Kernel> Delaunay_triangulation;
			typedef CGAL::Interpolation_traits_2<Kernel>   Interpolation_traits;
			
			typedef std::map<Point_2, FT, typename Kernel::Less_xy_2> Function_type;
			typedef CGAL::Data_access<Function_type> 		 		  Value_access;

			using Visibility = std::map<Face_handle, std::vector<Visibility_label> >;
			using Log = CGAL::LOD::Mylog;

			using Point_iterator = typename Container::const_iterator;
			using Face_iterator  = typename CDT::Finite_faces_iterator;

			enum class Radius_type { MIN, MAX };

			Level_of_detail_visibility_from_classification_2() : 
			m_approach(Visibility_approach::FACE_BASED), 
			m_method(Visibility_method::FACE_BASED_COUNT),
			m_sampler(Visibility_sampler::UNIFORM_0),
			m_radius_type(Radius_type::MIN),
			m_num_samples(200),
			m_k(3), 
			m_save_info(true), 
			m_show_progress(true),
			m_norm_threshold(FT(1000)) { }

			void set_number_of_samples(const size_t new_value) {
				m_num_samples = new_value;
			}

			void set_number_of_neighbours(const size_t new_value) {
				assert(!"This value is valid only for barycentric method, which currently does not work!");
				m_k = new_value;
			}

			void set_norm_threshold(const FT new_value) {
				m_norm_threshold = new_value;
			}

			void set_approach(const Visibility_approach new_approach) {
				m_approach = new_approach;
			}

			void set_method(const Visibility_method new_method) {
				m_method = new_method;
			}

			void save_info(const bool new_state) override {
				m_save_info = new_state;
			}

			void show_progress(const bool new_state) override {
				m_show_progress = new_state;
			}

			int compute(const Container &input, CDT &cdt) override {

				const FT half = FT(1) / FT(2);
				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
					
					fit->info().in       = half;
					fit->info().in_color = get_color(fit->info().in);
				}

				switch(m_approach) {

					case Visibility_approach::POINT_BASED:
					compute_point_based_visibility(input, cdt);
					break;

					case Visibility_approach::FACE_BASED:
					compute_face_based_approach(input, cdt);
					break;

					default:
					assert(!"Wrong approach!");
					break;
				}

				// Remove later.
				if (m_save_info) {

					Log log;
					log.out << "Visibility labels: " << std::endl;

					int count = 0;
					for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit, ++count) {

						const FT result = fit->info().in;
						std::string labelName = "default";

						if (result >  half) labelName = "IN";
						if (result <  half) labelName = "OUT";
						if (result == half) labelName = "UNKNOWN";

						log.out << "face index: " << count << " with label: " << labelName << " and visibility: " << result << std::endl;
					}
					log.save("tmp/visibility_classification");
				}

				return static_cast<int>(cdt.number_of_faces());
			}

		private:
			Visibility_approach m_approach;
			Visibility_method   m_method;

			const Visibility_sampler m_sampler;
			const Radius_type 		 m_radius_type;

			size_t m_num_samples;
			size_t m_k; 		  // change it to autodetection later!
			bool m_save_info;
			bool m_show_progress;
			FT m_norm_threshold;

			CGAL::Random m_rand;

			void compute_point_based_visibility(const Container &input, CDT &cdt) {

				Visibility visibility;
				for (Point_iterator it = input.begin(); it != input.end(); ++it) {

					const Point_2 &p = (*it).first;
					typename CDT::Locate_type locate_type;
					
					int locate_index = -1;
					const Face_handle face_handle = cdt.locate(p, locate_type, locate_index);

					if (locate_type == CDT::VERTEX ||
						locate_type == CDT::EDGE /*||
						on_the_border(cdt, fh, p) */) {

						// Improve this part of the code if possible!
						set_unknown(face_handle, visibility);
						continue;
					}

					// assert(locate_type != CDT::OUTSIDE_CONVEX_HULL);
					// assert(locate_type != CDT::OUTSIDE_AFFINE_HULL);

					if (locate_type == CDT::OUTSIDE_CONVEX_HULL ||
						locate_type == CDT::OUTSIDE_AFFINE_HULL) continue;

					assert(locate_type == CDT::FACE);
					switch(m_method) {

						case Visibility_method::POINT_BASED_CLASSIFICATION: {

							const Label point_label = (*it).second;
							estimate_with_classification(face_handle, point_label, visibility);
							break;
						}

						default:
						assert(!"Wrong visibility method!");
						break;
					}
				}
				postprocess(visibility);
			}

			void estimate_with_classification(const Face_handle face_handle, const Label point_label, Visibility &visibility) {

				const Label ground     = 0;
				const Label facade     = 1;
				const Label roof       = 2;
				const Label vegetation = 3; 

				switch (point_label) {

					case ground:
					set_outside(face_handle, visibility);
					break;

					case facade:
					set_unknown(face_handle, visibility);
					break;

					case roof:
					set_inside(face_handle, visibility);
					break;

					case vegetation:
					set_outside(face_handle, visibility);
					break;

					default: {
						// std::cout << "WARNING! CLASSIFICATION LABEL IS MISSING!" << std::endl; 
						// assert(!"Classification label is missing!");
						
						set_unknown(face_handle, visibility);
						break;
					}
				}
			}

			void compute_face_based_approach(const Container &input, CDT &cdt) {

				switch(m_method) {

					case Visibility_method::FACE_BASED_BARYCENTRIC:
						assert(!"Barycentric method does not work!");
						compute_face_based_approach_barycentric(input, cdt);
						break;

					case Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS:
						compute_face_based_approach_natural_neighbours(input, cdt);
						break;

					case Visibility_method::FACE_BASED_AFFINE:
						compute_face_based_approach_affine(input, cdt);
						break;

					case Visibility_method::FACE_BASED_COUNT:
						compute_face_based_approach_count(input, cdt);
						break;

					default:
						assert(!"Wrong visibility method!");
						break;
				}
			}

			void compute_face_based_approach_affine(const Container &, CDT &) {

				assert(!"This approach is not yet implemented!");
				// to be implemented!
			}

			void compute_face_based_approach_natural_neighbours(const Container &input, CDT &cdt) {

				Delaunay_triangulation dt;
				Function_type function_values;
				set_delunay_and_function_values(input, dt, function_values);

				// Can be removed -->
				size_t num_stages = 1;
				size_t stage      = 1;

				if (m_show_progress) {

					num_stages = 100;
					stage = cdt.number_of_faces() / num_stages;

					std::cout << std::endl;
					std::cout << "Progress: " << std::endl;
				}

				size_t progress = 0; // <--
				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit, ++progress) {
					
					if (m_show_progress) { // can be removed
						
						if (stage != 0 && progress % stage == 0)
							std::cout << (progress * 100) / cdt.number_of_faces() << "% " << std::endl; 
					}

					fit->info().in       = compute_interpolated_value(cdt, fit, dt, function_values);
					fit->info().in_color = get_color(fit->info().in);
				}

				if (m_show_progress) { // can be removed

					std::cout << "100%" << std::endl;
					std::cout << std::endl;
				}
			}

			CGAL::Color get_color(const FT visibility) {

				const FT half  = FT(1) / FT(2);
				
				FT scale_in  = FT(1) - visibility;
				FT scale_out = visibility;

				const FT scale_thres = FT(2) / FT(10);

				if (scale_in  < scale_thres) scale_in  = scale_thres;
				if (scale_out < scale_thres) scale_out = scale_thres;

				if (visibility > half)      return CGAL::Color(51 * scale_in, 255 * scale_in, 51 * scale_in);  	 // INSIDE
				else if (visibility < half) return CGAL::Color(255 * scale_out, 51 * scale_out, 51 * scale_out); // OUTSIDE
									  
				return CGAL::Color(255, 204, 0); // UNKNOWN
			}

			void set_delunay_and_function_values(const Container &input, Delaunay_triangulation &dt, Function_type &function_values) {

				const Label ground     = 0;
				const Label facade     = 1;
				const Label roof       = 2;
				const Label vegetation = 3; 

				const FT half = FT(1) / FT(2);
				for (size_t i = 0; i < input.size(); ++i) {

					const Point_2 &p  = input[i].first;
					const Label label = input[i].second;

					FT inside = half;
					switch(label) {

						case ground:
						inside = FT(0);
						break;

						case facade:
						inside = half;
						break;

						case roof:
						inside = FT(1);
						break;

						case vegetation:
						inside = FT(0);
						break;

						default:
						assert(!"Wrong label!");
						break;
					}

					dt.insert(p);
					function_values.insert(std::make_pair(p, inside));
				}
			}

			FT compute_interpolated_value(const CDT &cdt, const Face_iterator fh, const Delaunay_triangulation &dt, const Function_type &function_values) {

				std::vector<Point_2> samples(m_num_samples);
				generate_samples(cdt, fh, samples);

				// Log log; log.save_triangle_with_points_eps(cdt.triangle(fh).vertex(0), cdt.triangle(fh).vertex(1), cdt.triangle(fh).vertex(2), samples); // debugging info

				const FT half    = FT(1) / FT(2);
				FT result        = FT(0);
				size_t full_size = 0;

				for (size_t i = 0; i < samples.size(); ++i) {

					const Point_2 &query = samples[i];
					std::vector<std::pair<Point_2, FT> > coords;

					// May bug for some samples, gives division by zero assertion, for basic data set, probably because not enough natural neighbours can be found!
					const auto triple = CGAL::natural_neighbor_coordinates_2(dt, query, std::back_inserter(coords));

					const bool success = triple.third;
					const FT norm 	   = triple.second;

					// If a sample point gives an invalid result, skip it.
					if (!success) 			   continue;
					if (is_invalid_norm(norm)) continue; 

					assert(norm > FT(0));
					const FT intp = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(function_values));

					full_size += coords.size();
					result += intp;
				}
				assert(samples.size() != 0);

				result /= FT(samples.size());
				if (full_size == 0) result = half;				

				return result;
			}

			bool is_invalid_norm(const FT norm) {
				return (!std::isfinite(norm) || norm <= FT(0) || norm > m_norm_threshold);
			}

			void compute_face_based_approach_count(const Container &input, CDT &cdt) {		

				Fuzzy_tree tree(input.begin(), input.end());
				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {

					const Point_2 &a = cdt.triangle(fit).vertex(0);
					const Point_2 &b = cdt.triangle(fit).vertex(1);
					const Point_2 &c = cdt.triangle(fit).vertex(2);

					const Fuzzy_circle circle = compute_circle(a, b, c);

					Container pwl;
					tree.search(std::back_inserter(pwl), circle);

					fit->info().in       = compute_estimator(pwl);
					fit->info().in_color = get_color(fit->info().in);
				}
			}

			FT compute_estimator(const Container &pwl) {

				const Label ground     = 0;
				const Label facade     = 1;
				const Label roof       = 2;
				const Label vegetation = 3; 

				FT half = FT(1) / FT(2);
				
				FT inside  = FT(0);
				FT outside = FT(0);

				if (pwl.empty()) return half;
				for (size_t i = 0; i < pwl.size(); ++i) {

					const Label label = pwl[i].second;
					switch(label) {						// THIS CODE MAY GIVE COUNTERINTUITIVE RESULTS IF LABELS ARE DIFFERENT!

						case ground:
							outside += FT(1);
							break;

						case facade: {

							inside  += half;
							outside += half;

							break;
						}

						case roof:
							inside  += FT(1);
							break;

						case vegetation:
							outside += FT(1);
							break;

						default:
							assert(!"Wrong label!");
							break;
					}
				}

				const FT sum = FT(inside + outside); 
				if (sum == FT(0)) return half;

				return inside / sum;
			}

			// Here I can probably improve by taking the outscribed circle.
			Fuzzy_circle compute_circle(const Point_2 &a, const Point_2 &b, const Point_2 &c) {

				const Point_2 centre = compute_circle_centre(a, b, c);
				const FT radius = compute_circle_radius(a, b, c);

				return Fuzzy_circle(centre, radius);
			}

			Point_2 compute_circle_centre(const Point_2 &a, const Point_2 &b, const Point_2 &c) {

				const FT third = FT(1) / FT(3);
				return Point_2(third * a.x() + third * b.x() + third * c.x(), third * a.y() + third * b.y() + third * c.y());
			}

			FT compute_circle_radius(const Point_2 &a, const Point_2 &b, const Point_2 &c) {

				FT result = FT(0);
				switch(m_radius_type) {

					case Radius_type::MIN: {
					
					    FT min_length = CGAL::squared_distance(a, b);

						min_length = CGAL::min(min_length, CGAL::squared_distance(b, c));
						min_length = CGAL::min(min_length, CGAL::squared_distance(c, a));

						const FT threshold = min_length / FT(3);
						result = CGAL::sqrt(min_length) - threshold;

						break;
					}

					case Radius_type::MAX: {
					
						FT max_length = CGAL::squared_distance(a, b);

						max_length = CGAL::max(max_length, CGAL::squared_distance(b, c));
						max_length = CGAL::max(max_length, CGAL::squared_distance(c, a));

						const FT threshold = max_length / FT(5);
						result = CGAL::sqrt(max_length) + threshold;

						break;
					}

					default:
						assert(!"Wrong radius type!");
						break;
				}
				
				const FT threshold = FT(0); // (FT(1) / FT(5)) * max_length;
				return result + threshold;
			}

			void compute_face_based_approach_barycentric(const Container &input, CDT &cdt) {

				// Create a tree.
				Tree tree(input.begin(), input.end());

				// Iterate over all faces.
				std::vector<Point_2> samples(m_num_samples);
				std::vector<FT> visibility(m_num_samples);

				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {

					generate_samples(cdt, fit, samples);
					compute_visibility(tree, samples, visibility);
				}
			}

			void compute_visibility(const Tree &tree, const std::vector<Point_2> &samples, std::vector<FT> &visibility) {

				assert(visibility.size() == samples.size());
				for (size_t i = 0; i < samples.size(); ++i) {

					const Point_2 &p = samples[i];
					Neighbor_search search(tree, p, m_k);

					std::vector<FT> bc;
					compute_barycentric_coordinates(search, p, bc);
				}
			}

			// For the moment it works only with m_k = 3! For the bigger number of nearest neighbours, we need to use triangulation, because
			// the found nearest neighbours do not necesserary form an oriented polygon!
			// It can also give negative values.
			void compute_barycentric_coordinates(const Neighbor_search &search, const Point_2 &query, std::vector<FT> &bc) {

				const size_t num_found_neighbours = static_cast<size_t>(std::distance(search.begin(), search.end()));
				
				std::vector<Point_2> points(num_found_neighbours);
				bc.reserve(num_found_neighbours);

				size_t i = 0;
				for (typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it, ++i) points[i] = it->first.first;

				Mean_value_coordinates mean_value_coordinates(points.begin(), points.end());
				mean_value_coordinates(query, std::back_inserter(bc));

				assert(bc.size() == num_found_neighbours);
			}

			void generate_samples(const CDT &cdt, const Face_iterator fit, std::vector<Point_2> &samples) {
				
				switch(m_sampler) {

					case Visibility_sampler::UNIFORM_0:
						generate_samples_uniform_0(cdt, fit, samples);
						break;

					case Visibility_sampler::UNIFORM_1:
						generate_samples_uniform_1(cdt, fit, samples);
						break;

					default:
						assert(!"Wrong sampler!");
						break;
				}
			}

			void generate_samples_uniform_0(const CDT &cdt, const Face_iterator fh, std::vector<Point_2> &samples) {
				
				assert(!samples.empty() && samples.size() == m_num_samples);

				const Point_2 &a = cdt.triangle(fh).vertex(0);
				const Point_2 &b = cdt.triangle(fh).vertex(1);
				const Point_2 &c = cdt.triangle(fh).vertex(2);

				Vector_2 ab = b - a;
				Vector_2 ac = c - a;

				for (size_t i = 0; i < m_num_samples; ++i) {

					FT s = static_cast<FT>(m_rand.uniform_real(0.0, 1.0));
					FT t = static_cast<FT>(m_rand.uniform_real(0.0, 1.0));

					if (s + t > FT(1)) {

						s = FT(1) - s;
						t = FT(1) - t;
					}
					samples[i] = a + s * ab + t * ac;
				}

				// Log log;
				// log.save_triangle_with_points_eps(a, b, c, samples, "tmp/triangle_0");
			}

			void generate_samples_uniform_1(const CDT &cdt, const Face_iterator fh, std::vector<Point_2> &samples) {

				assert(!samples.empty() && samples.size() == m_num_samples);

				const Point_2 &a = cdt.triangle(fh).vertex(0);
				const Point_2 &b = cdt.triangle(fh).vertex(1);
				const Point_2 &c = cdt.triangle(fh).vertex(2);

				for (size_t i = 0; i < m_num_samples; ++i) {

					FT s = static_cast<FT>(m_rand.uniform_real(0.0, 1.0));
					FT t = static_cast<FT>(m_rand.uniform_real(0.0, 1.0));

					FT u =  FT(1)      - CGAL::sqrt(t);
					FT v = (FT(1) - s) * CGAL::sqrt(t);
					FT w =  		s  * CGAL::sqrt(t);

					samples[i] = Point_2(u * a.x() + v * b.x() + w * c.x(), u * a.y() + v * b.y() + w * c.y());
				}

				// Log log;
				// log.save_triangle_with_points_eps(a, b, c, samples, "tmp/triangle_1");
			}

			void set_inside(const Face_handle face_handle, Visibility &visibility) {
				visibility[face_handle].push_back(Visibility_label::IN);
			}

			void set_outside(const Face_handle face_handle, Visibility &visibility) {
				visibility[face_handle].push_back(Visibility_label::OUT);
			}

			void set_unknown(const Face_handle face_handle, Visibility &visibility) {
				visibility[face_handle].push_back(Visibility_label::UNKNOWN);
			}

			void postprocess(Visibility &visibility) {

				for (typename Visibility::const_iterator it = visibility.begin(); it != visibility.end(); ++it) {
					
					FT inside  = FT(0);
					FT outside = FT(0);

					const FT half = FT(1) / FT(2);
					for (size_t i = 0; i < (*it).second.size(); ++i) {

						switch ((*it).second[i]) {

							case Visibility_label::IN:
								inside  += FT(1);
								break;

							case Visibility_label::OUT:
								outside += FT(1);
								break;

							case Visibility_label::UNKNOWN: {

								inside  += half;
								outside += half;

								break;
							}

							default:
								assert(!"Should never get here!");
								break;
						}
					}

					const FT sum = FT(inside + outside); 
					if (sum == FT(0)) {

						(*it).first->info().in       = half;
						(*it).first->info().in_color = get_color((*it).first->info().in);
						
						continue;	
					}

					(*it).first->info().in 		 = inside / sum;
					(*it).first->info().in_color = get_color((*it).first->info().in);
				}
			}

			bool on_the_border(const CDT &cdt, const Face_handle fh, const Point_2 &p) {

				const Point_2 &a = cdt.triangle(fh).vertex(0);
				const Point_2 &b = cdt.triangle(fh).vertex(1);
				const Point_2 &c = cdt.triangle(fh).vertex(2);

				Triangle_coordinates tric(a, b, c);
				std::vector<FT> bc;
				tric(p, bc);

				const FT tol = 1.0e-6;

				if (bc[0] < tol || bc[1] < tol || bc[2] < tol) return true;
				return false;
			}
		};		
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_2_H