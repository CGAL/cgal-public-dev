#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_CLASSIFICATION_2_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_CLASSIFICATION_2_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\" 
#else 
#define PS "/" 
#endif 

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
#include <CGAL/Utils/Level_of_detail_utils.h>
#include <CGAL/Visibility_2/Level_of_detail_visibility_base_2.h>

// Eigen includes.
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

namespace CGAL {

	namespace LOD {

		// This class works only with the xy aligned ground plane that is Plane(0, 0, 1, 0).
		template<class KernelTraits, class ContainerInput, class CDTInput>
		class Level_of_detail_visibility_from_classification_2 : public Level_of_detail_visibility_base_2<KernelTraits, ContainerInput, CDTInput> {

		public:
			// Typedefs.
			typedef Level_of_detail_visibility_base_2<KernelTraits, ContainerInput, CDTInput> Base;

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


			// Extra.
			using Visibility = std::map<Face_handle, std::vector<Visibility_label> >;
			using Log = CGAL::LOD::Mylog;

			using Point_iterator = typename Container::const_iterator;
			using Face_iterator  = typename CDT::Finite_faces_iterator;

			using Samples = std::vector<Point_2>;
			using Sample_generator = Uniform_sample_generator<Kernel>;

			enum class Radius_type { MIN, MAX };


			// Eigen typedefs.
        	typedef Eigen::VectorXd VectorXd;
        	typedef Eigen::MatrixXd MatrixXd;


			// Constructor.
			Level_of_detail_visibility_from_classification_2() : 
			m_approach(Visibility_approach::POINT_BASED), 
			m_method(Visibility_method::POINT_BASED_CLASSIFICATION),
			m_sampler(Visibility_sampler::UNIFORM_SUBDIVISION),
			m_radius_type(Radius_type::MAX),
			m_num_samples(3),
			m_k(6), 
			m_save_info(true), 
			m_show_progress(true),
			m_norm_threshold(FT(1000)),
			m_without_bc(true) { }


			// Public functions.
			void save_info(const bool new_state) override {
				m_save_info = new_state;
			}

			void set_approach(const Visibility_approach new_approach) override {
				m_approach = new_approach;
			}

			void set_method(const Visibility_method new_method) override {
				m_method = new_method;
			}

			void set_number_of_samples(const size_t new_value) override {
				assert(new_value >= 0);
				m_num_samples = new_value;
			}

			void show_progress(const bool new_state) override {
				m_show_progress = new_state;
			}

			void set_norm_threshold(const FT new_value) override {
				assert(new_value >= FT(0));
				m_norm_threshold = new_value;
			}

			void set_number_of_neighbours(const size_t new_value) override {
				assert(new_value >= 0);
				m_k = new_value;
			}

			void set_sampler_type(const Visibility_sampler new_sampler) override {
				m_sampler = new_sampler;
			}

			void set_number_of_rays_per_side(const size_t) override {
				return;
			}

			void set_small_edge_threshold(const FT) override {
				return;
			}

			std::string name() override {
				return "classification";
			}


			int compute(const Container &input, CDT &cdt) override {

				m_sample_generator.set_number_of_samples(m_num_samples);

				const FT half = FT(1) / FT(2);
				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
					
					fit->info().in       = half;
					fit->info().in_color = this->get_color(fit->info().in);
				}

				switch (m_approach) {

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
					log.save("tmp" + std::string(PS) + "visibility_classification");
				}

				this->global_postprocess(cdt);
				return static_cast<int>(cdt.number_of_faces());
			}

		private:
			Visibility_approach m_approach;
			Visibility_method   m_method;

			Visibility_sampler m_sampler;
			const Radius_type  m_radius_type;

			size_t    m_num_samples;
			size_t              m_k; // change it to the autodetection later!
			bool        m_save_info;
			bool    m_show_progress;
			FT     m_norm_threshold;

			Sample_generator m_sample_generator;
			const bool m_without_bc;

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
					switch (m_method) {

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

				switch (m_method) {

					case Visibility_method::FACE_BASED_BARYCENTRIC:
						assert(!"Does not work!");
						compute_face_based_approach_barycentric(input, cdt);
						break;

					case Visibility_method::FACE_BASED_NATURAL_NEIGHBOURS:
						compute_face_based_approach_natural_neighbours(input, cdt);
						break;

					case Visibility_method::FACE_BASED_COUNT:
						assert(!"Does not work!");
						compute_face_based_approach_count(input, cdt);
						break;

					default:
						assert(!"Wrong visibility method!");
						break;
				}
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

					fit->info().in       = compute_interpolated_value_natural_neighbours(cdt, fit, dt, function_values);
					fit->info().in_color = this->get_color(fit->info().in);
				}

				if (m_show_progress) { // can be removed

					std::cout << "100%" << std::endl;
					std::cout << std::endl;
				}
			}

			void set_delunay_and_function_values(const Container &input, Delaunay_triangulation &dt, Function_type &function_values) {

				const Label ground     =  0;
				const Label facade     =  1;
				const Label roof       =  2;
				const Label vegetation =  3; 
				const Label fix 	   = -1; // remove later!

				const FT half = FT(1) / FT(2);
				for (size_t i = 0; i < input.size(); ++i) {

					const Point_2 &p  = input[i].first;
					Label label = input[i].second;

					FT inside = half;
					switch (label) {

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

						case fix:
							inside = half;
							break;

						default:
							assert(!"Wrong label!");
							break;
					}

					dt.insert(p);
					function_values.insert(std::make_pair(p, inside));
				}
			}

			FT compute_interpolated_value_natural_neighbours(const CDT &cdt, const Face_iterator fh, const Delaunay_triangulation &dt, const Function_type &function_values) {

				Samples samples;
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

				if (result > FT(1) && result < FT(1) + FT(1) / FT(1000000)) result = FT(1);
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
					fit->info().in_color = this->get_color(fit->info().in);
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

					// THIS CODE MAY GIVE COUNTERINTUITIVE RESULTS IF LABELS ARE DIFFERENT!
					const Label label = pwl[i].second;
					switch (label) {						

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
				switch (m_radius_type) {

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
				std::vector<Point_2> samples;

				size_t num_stages = 1;
				size_t stage      = 1;

				if (m_show_progress) {

					num_stages = 100;
					stage = cdt.number_of_faces() / num_stages;

					std::cout << std::endl;
					std::cout << "Progress: " << std::endl;
				}

				size_t progress = 0;
				for (Face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit, ++progress) {
					
					if (m_show_progress) {
						
						if (stage != 0 && progress % stage == 0)
							std::cout << (progress * 100) / cdt.number_of_faces() << "% " << std::endl; 
					}

					generate_samples(cdt, fit, samples);

					fit->info().in       = compute_visibility_barycentric(tree, samples);
					fit->info().in_color = this->get_color(fit->info().in);
				}

				if (m_show_progress) {

					std::cout << "100%" << std::endl;
					std::cout << std::endl;
				}
			}

			FT compute_visibility_barycentric(const Tree &tree, const std::vector<Point_2> &samples) {

				assert(m_k >= 3);
				assert(!samples.empty());
 
				FT result = FT(0);
				for (size_t i = 0; i < samples.size(); ++i) {

					const Point_2 &query = samples[i];
					Neighbor_search search(tree, query, m_k);

					const FT intp = interpolate_sample_with_affine_coordinates(search, query);
					result += intp;
				}

				result /= static_cast<FT>(samples.size());				
				return result;
			}

			FT interpolate_sample_with_affine_coordinates(const Neighbor_search &search, const Point_2 &query) {

				const size_t num_found_neighbours = static_cast<size_t>(std::distance(search.begin(), search.end()));
				if (num_found_neighbours == 0) return FT(1) / FT(2);
				
				std::vector<Point_2> points(num_found_neighbours);
				std::vector<FT>      function_values(num_found_neighbours);

				size_t i = 0; 
				for (typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it, ++i) {
					assert(i < points.size());

					points[i]          = it->first.first;
					function_values[i] = get_function_value(it->first.second);
				}
				
				std::vector<FT> bc;
				compute_affine_coordinates(points, query, bc);

				return interpolate_with_barycentric_coordinates(function_values, bc);
			}

			FT get_function_value(const Label label) {

				const Label ground     = 0;
				const Label facade     = 1;
				const Label roof       = 2;
				const Label vegetation = 3; 

				switch (label) {						

					case ground:
						return FT(0);

					case facade:
						return FT(1) / FT(2);

					case roof:
						return FT(1);

					case vegetation:
						return FT(0);

					default:
						assert(!"Wrong label!");
						return FT(0);
				}
			}

			void compute_affine_coordinates(const std::vector<Point_2> &points, const Point_2 &query, std::vector<FT> &bc) {

				bc.clear();
				if (m_without_bc) return;

				assert(!points.empty());
				const size_t n = points.size();
				bc.resize(n, FT(0));

				// Compute the barycentre.
				Point_2 c;
				compute_barycentre(points, c);

				// Set the matrices.
            	MatrixXd V(2, n);
            	for (size_t i = 0; i < n; ++i) {

                	const FT x = points[i].x() - c.x();
                	const FT y = points[i].y() - c.y();

                	V(0, i) = x;
                	V(1, i) = y;
            	}

            	const MatrixXd Vs  = V.adjoint();
            	const MatrixXd mat = V * Vs;
            	const MatrixXd inv = mat.inverse();

            	VectorXd vert(2);
            	Point_2 diff;

            	const FT num_points = static_cast<FT>(n);
            	for (size_t i = 0; i < n; ++i) {
                	
            		const FT x = query.x() - c.x();
            		const FT y = query.y() - c.y();

                	vert(0) = V(0, i);
                	vert(1) = V(1, i);

                	VectorXd res = inv * vert;
    	            bc[i] = x * res(0) + y * res(1) + FT(1) / num_points;
	            }
			}

			void compute_barycentre(const std::vector<Point_2> &points, Point_2 &c) {

				FT x = FT(0); 
				FT y = FT(0);

				for (size_t i = 0; i < points.size(); ++i) {

					x += points[i].x();
					y += points[i].y();
				}

				const FT num_points = static_cast<FT>(points.size());

				x /= num_points;
				y /= num_points;

				c = Point_2(x, y);
			}

			FT interpolate_with_barycentric_coordinates(const std::vector<FT> &function_values, const std::vector<FT> &bc) {

				if (bc.empty() && m_without_bc) return apply_averaging(function_values);

				assert(!bc.empty());
				assert(!function_values.empty());
				assert(function_values.size() == bc.size());

				FT result = FT(0);
				for (size_t i = 0; i < function_values.size(); ++i) {
					
					assert(function_values[i] >= FT(0) && function_values[i] <= FT(1));
					result += function_values[i] * bc[i];
				}

				assert(result >= FT(0) && result <= FT(1));
				if (result < FT(0) || result > FT(1)) {
				
					// std::cerr << "The value is truncated!" << std::endl;
					return FT(1) / FT(2);
				}
				return result;
			}

			FT apply_averaging(const std::vector<FT> &function_values) {

				const FT half = FT(1) / FT(2);
				if (function_values.empty()) return half;

				FT result = FT(0);
				for (size_t i = 0; i < function_values.size(); ++i) result += function_values[i];
				result /= static_cast<FT>(function_values.size());

				return result;
			}

			void generate_samples(const CDT &cdt, const Face_iterator fit, std::vector<Point_2> &samples) {
				
				switch (m_sampler) {

					case Visibility_sampler::BARYCENTRE:
						generate_barycentre(cdt, fit, samples);
						break;

					case Visibility_sampler::RANDOM_UNIFORM_0:
						generate_samples_random_uniform_0(cdt, fit, samples);
						break;

					case Visibility_sampler::RANDOM_UNIFORM_1:
						generate_samples_random_uniform_1(cdt, fit, samples);
						break;

					case Visibility_sampler::UNIFORM_SUBDIVISION:
						generate_samples_uniform_subdivision(cdt, fit, samples);
						break;

					default:
						assert(!"Wrong sampler!");
						break;
				}
			}

			void generate_barycentre(const CDT &cdt, const Face_iterator &fh, Samples &samples) {

				const Point_2 &a = cdt.triangle(fh).vertex(0);
				const Point_2 &b = cdt.triangle(fh).vertex(1);
				const Point_2 &c = cdt.triangle(fh).vertex(2);

				m_sample_generator.set_number_of_samples(0);
				m_sample_generator.create_uniform_subdivision_samples(a, b, c, samples);
			}

			void generate_samples_random_uniform_0(const CDT &cdt, const Face_iterator &fh, Samples &samples) {

				const Point_2 &a = cdt.triangle(fh).vertex(0);
				const Point_2 &b = cdt.triangle(fh).vertex(1);
				const Point_2 &c = cdt.triangle(fh).vertex(2);

				m_sample_generator.create_random_uniform_samples_0(a, b, c, samples);

				// Log log;
				// log.save_triangle_with_points_eps(a, b, c, samples, "tmp" + std::string(PS) + "triangle_0");
			}

			void generate_samples_random_uniform_1(const CDT &cdt, const Face_iterator &fh, Samples &samples) {

				const Point_2 &a = cdt.triangle(fh).vertex(0);
				const Point_2 &b = cdt.triangle(fh).vertex(1);
				const Point_2 &c = cdt.triangle(fh).vertex(2);

				m_sample_generator.create_random_uniform_samples_1(a, b, c, samples);

				// Log log;
				// log.save_triangle_with_points_eps(a, b, c, samples, "tmp" + std::string(PS) + "triangle_1");
			}

			void generate_samples_uniform_subdivision(const CDT &cdt, const Face_iterator &fh, Samples &samples) {

				const Point_2 &a = cdt.triangle(fh).vertex(0);
				const Point_2 &b = cdt.triangle(fh).vertex(1);
				const Point_2 &c = cdt.triangle(fh).vertex(2);

				m_sample_generator.create_uniform_subdivision_samples(a, b, c, samples);
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
						(*it).first->info().in_color = this->get_color((*it).first->info().in);
						
						continue;	
					}

					(*it).first->info().in 		 = inside / sum;
					(*it).first->info().in_color = this->get_color((*it).first->info().in);
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

#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_FROM_CLASSIFICATION_2_H