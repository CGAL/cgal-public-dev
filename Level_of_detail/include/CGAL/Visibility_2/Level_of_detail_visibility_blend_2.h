#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_BLEND_2_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_BLEND_2_H

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
#include <CGAL/Visibility_2/Level_of_detail_visibility_ray_shooting_2.h>
#include <CGAL/Visibility_2/Level_of_detail_visibility_from_classification_2.h>

// Eigen includes.
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

namespace CGAL {

	namespace LOD {

		// This is a visibility class that blends both ray shooting and classification visibility estimation algorithms below.
		template<class KernelTraits, class ContainerInput, class CDTInput>
		class Level_of_detail_visibility_blend_2 : public Level_of_detail_visibility_base_2<KernelTraits, ContainerInput, CDTInput> {

		public:
			// Typedefs.
			typedef KernelTraits 		Kernel;
			typedef typename Kernel::FT FT;

			typedef ContainerInput Container;
			typedef CDTInput 	   CDT;

			using LOD_classification = Level_of_detail_visibility_from_classification_2<Kernel, Container, CDT>;
			using LOD_ray_shooting   = Level_of_detail_visibility_ray_shooting_2<Kernel, Container, CDT>;


			// Extra.
			using Face_handle 		  = typename CDT::Face_handle;
			using Cdt_vertex_iterator = typename CDT::Finite_vertices_iterator;
			using Cdt_face_iterator   = typename CDT::Finite_faces_iterator;


			// Constructor.
			Level_of_detail_visibility_blend_2() : 
			m_save_info(true), 
			m_show_progress(true),
			m_num_samples(1),
			m_num_rays_per_side(100),
			m_k(6),
			m_norm_threshold(FT(1000)),
			m_small_edge_thres(FT(0)),
			m_sampler(Visibility_sampler::BARYCENTRE),
			m_approach(Visibility_approach::POINT_BASED),
			m_method(Visibility_method::POINT_BASED_CLASSIFICATION) { }


			// Public functions.
			void save_info(const bool new_state) override {
				m_save_info = new_state;
			}

			void show_progress(const bool new_state) override {
				m_show_progress = new_state;
			}

			void set_number_of_samples(const size_t new_value) override {
				assert(new_value >= 0);
				m_num_samples = new_value;
			}

			void set_number_of_rays_per_side(const size_t new_value) override {
				assert(new_value > 0);
				m_num_rays_per_side = new_value;
			}

			void set_number_of_neighbours(const size_t new_value) override {
				assert(new_value >= 0);
				m_k = new_value;
			}

			void set_norm_threshold(const FT new_value) override {
				assert(new_value >= FT(0));
				m_norm_threshold = new_value;
			}

			void set_small_edge_threshold(const FT new_value) override {
				m_small_edge_thres = new_value;
			}

			void set_sampler_type(const Visibility_sampler new_sampler) override {
				m_sampler = new_sampler;
			}

			void set_approach(const Visibility_approach new_approach) override {
				m_approach = new_approach;
			}

			void set_method(const Visibility_method new_method) override {
				m_method = new_method;
			}

			std::string name() override {
				return "blend";
			}


			int compute(const Container &input, CDT &cdt) override {

				set_classification_lod();
				set_ray_shooting_lod();

				CDT tmp = cdt;

				m_lod_classification.compute(input, tmp);
				m_lod_ray_shooting.compute(input, cdt);

				update_results(tmp, cdt);

				this->global_postprocess(cdt);
				return static_cast<int>(cdt.number_of_faces());
			}

		private:
			bool m_save_info;
			bool m_show_progress;

			size_t m_num_samples;
			size_t m_num_rays_per_side;
			size_t m_k;

			FT m_norm_threshold;
			FT m_small_edge_thres;

			Visibility_sampler  m_sampler;
			Visibility_approach m_approach;
			Visibility_method   m_method;

			LOD_classification m_lod_classification;
			LOD_ray_shooting   m_lod_ray_shooting;

			void set_classification_lod() {

				m_lod_classification.save_info(m_save_info);
				m_lod_classification.show_progress(m_show_progress);
				m_lod_classification.set_number_of_samples(m_num_samples);
				m_lod_classification.set_number_of_rays_per_side(m_num_rays_per_side);
				m_lod_classification.set_number_of_neighbours(m_k);
				m_lod_classification.set_norm_threshold(m_norm_threshold);
				m_lod_classification.set_small_edge_threshold(m_small_edge_thres);
				m_lod_classification.set_sampler_type(m_sampler);
				m_lod_classification.set_approach(m_approach);
				m_lod_classification.set_method(m_method);
			}

			void set_ray_shooting_lod() {

				m_lod_ray_shooting.save_info(m_save_info);
				m_lod_ray_shooting.show_progress(m_show_progress);
				m_lod_ray_shooting.set_number_of_samples(m_num_samples);
				m_lod_ray_shooting.set_number_of_rays_per_side(m_num_rays_per_side);
				m_lod_ray_shooting.set_number_of_neighbours(m_k);
				m_lod_ray_shooting.set_norm_threshold(m_norm_threshold);
				m_lod_ray_shooting.set_small_edge_threshold(m_small_edge_thres);
				m_lod_ray_shooting.set_sampler_type(m_sampler);
				m_lod_ray_shooting.set_approach(m_approach);
				m_lod_ray_shooting.set_method(m_method);
			}

			void update_results(const CDT &tmp, CDT &cdt) {
				assert(tmp.number_of_faces() == cdt.number_of_faces());

				Cdt_face_iterator fit_tmp = tmp.finite_faces_begin();
				for (Cdt_face_iterator fit_cdt = cdt.finite_faces_begin(); fit_cdt != cdt.finite_faces_end(); ++fit_tmp, ++fit_cdt) { 

					const Face_handle fh_tmp = static_cast<Face_handle>(fit_tmp);
					const Face_handle fh_cdt = static_cast<Face_handle>(fit_cdt);
					
					fh_cdt->info().in       = estimate_face_visibility(fh_tmp, fh_cdt);
					fh_cdt->info().in_color = this->get_color(fh_cdt->info().in);
				}
			}

			FT estimate_face_visibility(const Face_handle &fh_tmp, const Face_handle &fh_cdt) {

				const FT half = FT(1) / FT(2);

				const FT in_tmp = fh_tmp->info().in;
				const FT in_cdt = fh_cdt->info().in;

				// Simple averaging.

				// return (in_tmp + in_cdt) / FT(2);

				// More complex version.

				if (in_tmp >= half && in_cdt >= half) return (in_tmp + in_cdt) / FT(2);
				if (in_tmp <= half && in_cdt <= half) return (in_tmp + in_cdt) / FT(2);

				assert(in_tmp != half && in_cdt != half);

				if (in_tmp > half && in_cdt < half) {

					const FT val_tmp = FT(1) - in_tmp;
					const FT val_cdt = in_cdt;

					return compare_two_vals(val_tmp, val_cdt, in_tmp, in_cdt);
				}

				if (in_tmp < half && in_cdt > half) {
					
					const FT val_tmp = in_tmp;
					const FT val_cdt = FT(1) - in_cdt;

					return compare_two_vals(val_tmp, val_cdt, in_tmp, in_cdt);
				}

				assert(!"Cannot be here!");
				return half;
			}

			FT compare_two_vals(const FT val_tmp, const FT val_cdt, const FT in_tmp, const FT in_cdt) {

				if (val_tmp > val_cdt) return in_cdt;
				else if (val_tmp < val_cdt) return in_tmp;
				
				return (in_tmp + in_cdt) / FT(2);
			}
		};	
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_BLEND_2_H