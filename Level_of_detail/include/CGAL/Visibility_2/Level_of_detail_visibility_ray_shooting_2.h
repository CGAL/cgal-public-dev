#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_RAY_SHOOTING_2_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_RAY_SHOOTING_2_H

#if defined(WIN32) || defined(_WIN32) 
#define PS "\\"
#define PN "\r\n"
#else 
#define PS "/" 
#define PN "\n"
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

		// This is the ray shooting based classification algorithm.
		template<class KernelTraits, class ContainerInput, class CDTInput>
		class Level_of_detail_visibility_ray_shooting_2 : public Level_of_detail_visibility_base_2<KernelTraits, ContainerInput, CDTInput> {

		public:
			// Typedefs.
			typedef Level_of_detail_visibility_base_2<KernelTraits, ContainerInput, CDTInput> Base;
			
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
			typedef typename CDT::Edge 				   Edge;


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

			using Sample_generator = CGAL::LOD::Uniform_sample_generator<Kernel>;


			// Constructor.
			Level_of_detail_visibility_ray_shooting_2() : 
			m_save_info(true), 
			m_show_progress(true),
			m_num_samples(1),
			m_num_rays_per_side(100),
			m_small_edge_thres(FT(0)),
			m_sampler(Visibility_sampler::BARYCENTRE) { }


			// Public functions.
			void save_info(const bool new_state) override {
				m_save_info = new_state;
			}

			void set_approach(const Visibility_approach) override {
				return;
			}

			void set_method(const Visibility_method) override {
				return;
			}

			void set_number_of_samples(const size_t new_value) override {
				assert(new_value >= 0);
				m_num_samples = new_value;
			}

			void show_progress(const bool new_state) override {
				m_show_progress = new_state;
			}

			void set_norm_threshold(const FT) override {
				return;
			}

			void set_number_of_neighbours(const size_t) override {
				return;
			}

			void set_sampler_type(const Visibility_sampler new_sampler) override {
				m_sampler = new_sampler;
			}

			void set_number_of_rays_per_side(const size_t new_value) override {
				assert(new_value > 0);
				m_num_rays_per_side = new_value;
			}

			// NOTE: new_value can be negative, too!
			void set_small_edge_threshold(const FT new_value) override {
				m_small_edge_thres = new_value;
			}

			std::string name() override {
				return "ray shooting";
			}


			int compute(const Container &, CDT &cdt) override {

				m_sample_generator.set_number_of_samples(m_num_samples);
				m_sample_generator.set_number_of_rays(m_num_rays_per_side);	

				Log log;

				// (1) Compute bounding box.
				std::vector<Point_2> bbox;
				compute_bounding_box(cdt, log, bbox);

				if (m_save_info) log.out << "(1) Bounding box is computed!" << std::endl << std::endl;

				// (2) Compute visibility.
				compute_visibility(cdt, log, bbox);

				if (m_save_info) {
					
					log.out << "" + std::string(PN) + "(2) Visibility is computed!" << std::endl;
					log.save("tmp" + std::string(PS) + "visibility_ray_shooting");
				}

				this->global_postprocess(cdt);
				return static_cast<int>(cdt.number_of_faces());
			}

		private:
			bool   m_save_info;
			bool   m_show_progress;
			size_t m_num_samples;
			size_t m_num_rays_per_side;
			FT     m_small_edge_thres;

			Visibility_sampler m_sampler;
			Sample_generator   m_sample_generator;

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

					if (m_save_info) log.out << "" + std::string(PN) + "Face: " << face_index << std::endl;

					const Face_handle fh = static_cast<Face_handle>(fit);

					fh->info().in       = estimate_face_visibility(log, cdt, fh, bbox, face_index);
					fh->info().in_color = this->get_color(fh->info().in);
				}

				if (m_show_progress) { // can be removed

					std::cout << "100%" << std::endl;
					std::cout << std::endl;
				}
			}

			FT estimate_face_visibility(Log &log, const CDT &cdt, const Face_handle &fh, const Bbox &bbox, const size_t /* face_index */ ) {

				// (a) Generate samples.
				Samples samples;
				generate_samples(cdt, fh, samples);

				/* // Debugging.
				Log saver; saver.save_triangle_with_points_eps(
					cdt.triangle(fh).vertex(0), 
					cdt.triangle(fh).vertex(1), 
					cdt.triangle(fh).vertex(2), samples, "tmp" + std::string(PS) + "triangle_" + std::to_string(face_index)); */

				if (m_save_info) log.out << "(a) Samples are generated." << std::endl;


				// (b) Handle all samples.
				const FT face_visibility = estimate_face_visibility_from_samples(log, cdt, fh, bbox, samples);
				if (m_save_info) log.out << "" + std::string(PN) + "(b) Face visibility is estimated from the given samples: " << face_visibility << std::endl << std::endl;


				return face_visibility;
			}

			FT estimate_face_visibility_from_samples(Log &log, const CDT &cdt, const Face_handle &fh, const Bbox &bbox, const Samples &samples) {

				FT face_visibility = FT(0);
				for (Sample_iterator sit = samples.begin(); sit != samples.end(); ++sit) {

					if (m_save_info) log.out << "" + std::string(PN) + "New sample: " << std::endl;

					const Point_2 &sample = *sit;
					face_visibility += estimate_visibility_from_sample(log, sample, cdt, fh, bbox);
				}

				assert(!samples.empty());
				return face_visibility / static_cast<FT>(samples.size());
			}

			FT estimate_visibility_from_sample(Log &log, const Point_2 &sample, const CDT &cdt, const Face_handle &fh, const Bbox &bbox) {

				assert(bbox.size() == 4);
				FT sample_visibility = FT(0);


				// Handle bottom side of the bounding box.
				if (m_save_info) log.out << "" + std::string(PN) + "Bottom side: " << std::endl;
				sample_visibility += get_side_visibility(log, sample, cdt, fh, bbox[0], bbox[1]);


				// Handle right side of the bounding box.
				if (m_save_info) log.out << "" + std::string(PN) + "Right side: "  << std::endl;
				sample_visibility += get_side_visibility(log, sample, cdt, fh, bbox[1], bbox[2]);


				// Handle top side of the bounding box.
				if (m_save_info) log.out << "" + std::string(PN) + "Top side: "    << std::endl;
				sample_visibility += get_side_visibility(log, sample, cdt, fh, bbox[2], bbox[3]);


				// Handle left side of the bounding box.
				if (m_save_info) log.out << "" + std::string(PN) + "Left side: "   << std::endl;
				sample_visibility += get_side_visibility(log, sample, cdt, fh, bbox[3], bbox[0]);


				if (m_save_info) log.out << std::endl;
				return sample_visibility / static_cast<FT>(4);
			}

			FT get_side_visibility(Log &log, const Point_2 &sample, const CDT &cdt, const Face_handle &fh, const Point_2 &a, const Point_2 &b) {

				if (m_save_info) log.out << "side: (" << a << ") -- (" << b << ")" << std::endl;


				// Generate rays.
				Rays rays;
				generate_side_rays(sample, a, b, rays);

				if (m_save_info) log.out << "rays are generated" << std::endl;


				// Handle all rays.
				assert(!rays.empty());
				return get_visibility_from_rays(log, rays, cdt, fh) / static_cast<FT>(rays.size());
			}

			void generate_side_rays(const Point_2 &source, const Point_2 &a, const Point_2 &b, Rays &rays) {
				
				if (m_sampler == Visibility_sampler::RANDOM_UNIFORM_0 ||
					m_sampler == Visibility_sampler::RANDOM_UNIFORM_1) {

					m_sample_generator.create_random_uniform_rays(source, a, b, rays);
				} else m_sample_generator.create_uniform_rays(source, a, b, rays);
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

				Line_face_circulator lfc;
				if (cdt.oriented_side(fh, p) == CGAL::ON_NEGATIVE_SIDE) {
				
					if (cdt.oriented_side(fh, q) == CGAL::ON_NEGATIVE_SIDE) return FT(1) / FT(2);
					lfc = cdt.line_walk(q, p, fh);

				} else lfc = cdt.line_walk(p, q, fh);

				const FT max_small_edge = compute_max_small_edge(lfc, cdt);
				const FT ray_visibility = traverse_ray_faces(log, lfc, cdt, fh, max_small_edge);

				if (m_save_info) log.out << "max small edge: " << max_small_edge << "; ray visibility: " << ray_visibility << std::endl << std::endl;
				return ray_visibility;
			}

			FT compute_max_small_edge(const Line_face_circulator &lfc, const CDT &cdt) {

				Line_face_circulator tmp_lfc = lfc;

				Edge edge;
				const Face_handle end = tmp_lfc;
				
				FT sum_edge_length = FT(0);
				size_t num_edges   = 0;

				do {
					assert(is_valid_traversal_face(tmp_lfc, cdt));

					const Face_handle curr = static_cast<Face_handle>(tmp_lfc);
					const Face_handle next = static_cast<Face_handle>(++tmp_lfc);

					if (adjacent_by_edge(curr, next, edge)) {

						sum_edge_length += compute_edge_length(edge, cdt);
						++num_edges;
					}

				} while (is_valid_traversal_face(tmp_lfc, cdt) && tmp_lfc != end);

				return sum_edge_length / static_cast<FT>(num_edges) - m_small_edge_thres;
			}

			FT compute_edge_length(const Edge &edge, const CDT &cdt) {

				const Segment_2 &segment = cdt.segment(edge);
				return CGAL::sqrt(segment.squared_length());
			}

			FT traverse_ray_faces(Log &log, Line_face_circulator &lfc, const CDT &cdt, const Face_handle &fh, const FT max_small_edge) {

				if (m_save_info) log.out << std::endl;
				const Face_handle end = lfc;

				assert(lfc != NULL);
				assert(lfc == fh);
				
				if (lfc != fh) {
					std::cerr << "Error: lfc != fh, traverse_ray_faces function ray shooting!" << std::endl;
					exit(EXIT_FAILURE);
				}

				size_t tmp_sign_changes, sign_changes = 0;
				do {
					assert(is_valid_traversal_face(lfc, cdt));
					if (m_save_info) log.out << "f: " << cdt.triangle(lfc) << "; ";

					tmp_sign_changes = check_face_and_its_neighbour(log, lfc, cdt, max_small_edge);
					if (m_save_info) log.out << "sign changed: " << tmp_sign_changes << " times; " << std::endl;

					sign_changes += tmp_sign_changes;

				} while (is_valid_traversal_face(lfc, cdt) && lfc != end);
				assert(tmp_sign_changes >= 0);

				const FT ray_visibility = compute_ray_visibility(sign_changes);
				return ray_visibility;
			}

			size_t check_face_and_its_neighbour(Log &log, Line_face_circulator &lfc, const CDT &cdt, const FT max_small_edge) {

				const Face_handle curr = static_cast<Face_handle>(lfc);
				const Face_handle next = static_cast<Face_handle>(++lfc);

				return get_sign_changes_from_two_faces(log, curr, next, cdt, max_small_edge);
			}

			size_t get_sign_changes_from_two_faces(Log &log, const Face_handle &curr, const Face_handle &next, const CDT &cdt, const FT max_small_edge) {

				Edge edge;
				if (!adjacent_by_edge(curr, next, edge)) return get_sign_changes_from_two_faces_adjacent_by_vertex();

				if (m_save_info) log.out << "edge: " << cdt.segment(edge) << "; ";
				return get_sign_changes_from_two_faces_ajacent_by_edge(edge, next, cdt, max_small_edge);
			}

			bool adjacent_by_edge(const Face_handle &curr, const Face_handle &next, Edge &edge) {

				const size_t num_equal_vertices = number_of_equal_vertices(curr, next);
				if (num_equal_vertices == 1) return false;

				edge = std::make_pair(curr, curr->index(next));
				return true;
			}

			size_t number_of_equal_vertices(const Face_handle &curr, const Face_handle &next) {
				
				size_t num_equal_vertices = 0;
				for (size_t i = 0; i < 3; ++i) {	

					for (size_t j = 0; j < 3; ++j) {
						if (curr->vertex(i) == next->vertex(j)) {
							
							++num_equal_vertices;
							break;
						}
					}
				}

				assert(num_equal_vertices != 0 && num_equal_vertices != 3);
				return num_equal_vertices;
			}

			size_t get_sign_changes_from_two_faces_adjacent_by_vertex() {

				// For the moment we do not handle this case. Since I use statistics, it should not affect the final result very much.
				// Moreover it is very unlikely that this case will happen often!

				std::cerr << "WARNING: Two faces are adjacent at vertex!" << std::endl;
				return 0;
			}

			size_t get_sign_changes_from_two_faces_ajacent_by_edge(const Edge &edge, const Face_handle &next, const CDT &cdt, const FT /* max_small_edge */ ) {

				// (1) Structured points.

				// The constrained edge with the infinite face changes the sign twice.
				// It basically means that this edge is boundary edge of the building, which
				// is also the boundary edge of the convex hull of CDT.
				if (cdt.is_infinite(next) && cdt.is_constrained(edge)) return 2;


				// Constrained edge changes the sign.
				if (cdt.is_constrained(edge)) return 1;


				// ------------------------------------


				// (2) No more constrained edges. Handle clutter here!

				// Similar to above, but here we estimate if this edge is good enough
				// to change the sign or not.
				// if (cdt.is_infinite(next) && is_valid_edge_for_sign_change(edge, cdt, max_small_edge)) return 2;


				// Infinite boundary face changes the sign.
				if (cdt.is_infinite(next)) return 1;


				// Estimate validity of the edge to change the sign.
				// if (is_valid_edge_for_sign_change(edge, cdt, max_small_edge)) return 1;


				// Otherwise, the sign does not change.
				return 0;
			}

			bool is_valid_edge_for_sign_change(const Edge &edge, const CDT &cdt, const FT max_small_edge) {
				
				const Segment_2 &segment = cdt.segment(edge);
				const FT edge_length = CGAL::sqrt(segment.squared_length());

				return edge_length < max_small_edge;
			}

			FT compute_ray_visibility(const size_t sign_changes) {
				
				if (sign_changes % 2 == 0) return FT(1);
				return FT(0);
			}

			bool is_valid_traversal_face(const Line_face_circulator &lfc, const CDT &cdt) {
				return !cdt.is_infinite(lfc);
			}

			void generate_samples(const CDT &cdt, const Face_handle &fh, Samples &samples) {				

				switch (m_sampler) {

					case Visibility_sampler::BARYCENTRE:
						generate_barycentre(cdt, fh, samples);
						break;

					case Visibility_sampler::RANDOM_UNIFORM_0:
						generate_samples_random_uniform_0(cdt, fh, samples);
						break;

					case Visibility_sampler::RANDOM_UNIFORM_1:
						generate_samples_random_uniform_1(cdt, fh, samples);
						break;

					case Visibility_sampler::UNIFORM_SUBDIVISION:
						generate_samples_uniform_subdivision(cdt, fh, samples);
						break;

					default:
						assert(!"Wrong visibility sampler");
						break;
				}
			}

			void generate_barycentre(const CDT &cdt, const Face_handle &fh, Samples &samples) {

				/*
				const Point_2 &a = cdt.triangle(fh).vertex(0);
				const Point_2 &b = cdt.triangle(fh).vertex(1);
				const Point_2 &c = cdt.triangle(fh).vertex(2); */

				m_sample_generator.set_number_of_samples(0);
				m_sample_generator.create_uniform_subdivision_samples(cdt.triangle(fh).vertex(0), cdt.triangle(fh).vertex(1), cdt.triangle(fh).vertex(2), samples);
			}

			void generate_samples_random_uniform_0(const CDT &cdt, const Face_handle &fh, Samples &samples) {

				/*
				const Point_2 &a = cdt.triangle(fh).vertex(0);
				const Point_2 &b = cdt.triangle(fh).vertex(1);
				const Point_2 &c = cdt.triangle(fh).vertex(2); */

				m_sample_generator.create_random_uniform_samples_0(cdt.triangle(fh).vertex(0), cdt.triangle(fh).vertex(1), cdt.triangle(fh).vertex(2), samples);
			}

			void generate_samples_random_uniform_1(const CDT &cdt, const Face_handle &fh, Samples &samples) {

				/*
				const Point_2 &a = cdt.triangle(fh).vertex(0);
				const Point_2 &b = cdt.triangle(fh).vertex(1);
				const Point_2 &c = cdt.triangle(fh).vertex(2); */

				m_sample_generator.create_random_uniform_samples_1(cdt.triangle(fh).vertex(0), cdt.triangle(fh).vertex(1), cdt.triangle(fh).vertex(2), samples);
			}

			void generate_samples_uniform_subdivision(const CDT &cdt, const Face_handle &fh, Samples &samples) {

				/*
				const Point_2 &a = cdt.triangle(fh).vertex(0);
				const Point_2 &b = cdt.triangle(fh).vertex(1);
				const Point_2 &c = cdt.triangle(fh).vertex(2); */

				m_sample_generator.create_uniform_subdivision_samples(cdt.triangle(fh).vertex(0), cdt.triangle(fh).vertex(1), cdt.triangle(fh).vertex(2), samples);
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_RAY_SHOOTING_2_H