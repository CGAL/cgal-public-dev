#ifndef CGAL_LEVEL_OF_DETAIL_VISIBILITY_BASE_2_H
#define CGAL_LEVEL_OF_DETAIL_VISIBILITY_BASE_2_H

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

// Eigen includes.
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class ContainerInput, class CDTInput>
		class Level_of_detail_visibility_base_2 {

		public:
			typedef KernelTraits 		Kernel;

			typedef typename Kernel::FT 	  FT;
			typedef typename Kernel::Point_2  Point_2;
			typedef typename Kernel::Vector_2 Vector_2;

			typedef ContainerInput Container;
			typedef CDTInput 	   CDT;

			typedef CGAL::Color 			  VColor;
			typedef typename CDT::Face_handle Face_handle;

			typename Kernel::Compute_scalar_product_2 dot_product;
			
			// Public functions.
			virtual void save_info(const bool) = 0;
			virtual void show_progress(const bool) = 0;
			virtual void set_approach(const Visibility_approach) = 0;
			virtual void set_method(const Visibility_method) = 0;
			virtual void set_number_of_samples(const size_t) = 0;
			virtual void set_norm_threshold(const FT) = 0;
			virtual void set_number_of_neighbours(const size_t) = 0;
			virtual void set_sampler_type(const Visibility_sampler) = 0;
			virtual void set_number_of_rays_per_side(const size_t) = 0;
			virtual void set_small_edge_threshold(const FT new_value) = 0;
			virtual std::string name() = 0;

			virtual int compute(const Container &, CDT &) = 0;

			Level_of_detail_visibility_base_2() : m_angle_eps(FT(0)) { }
			virtual ~Level_of_detail_visibility_base_2() { }

			void set_angle_eps(const FT new_value) {

				assert(new_value >= FT(0));
				m_angle_eps = new_value;
			}

		protected:
			VColor get_color(const FT visibility) {

				assert(visibility >= FT(0) && visibility <= FT(1));
				const FT half  = FT(1) / FT(2);
				
				FT scale_in  = FT(1) - visibility;
				FT scale_out = visibility;

				const FT scale_thres = FT(2) / FT(10);

				if (scale_in  < scale_thres) scale_in  = scale_thres;
				if (scale_out < scale_thres) scale_out = scale_thres;

				if (visibility > half)      return VColor(51 * scale_in, 255 * scale_in, 51 * scale_in);  	 // INSIDE
				else if (visibility < half) return VColor(255 * scale_out, 51 * scale_out, 51 * scale_out); // OUTSIDE
									  
				return VColor(255, 204, 0); // UNKNOWN
			}

			void global_postprocess(CDT &cdt) {

				using Cdt_face_iterator = typename CDT::Finite_faces_iterator;
				for (Cdt_face_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
					
					const Face_handle fh = static_cast<Face_handle>(fit);
					if (!is_valid_face(cdt, fh)) {
						
						fh->info().in = FT(0);
						fh->info().in_color = get_color(fh->info().in);
					}
				}
			}

		private:
			FT m_angle_eps;

			bool is_valid_face(const CDT &cdt, const Face_handle &fh) {

				const auto &triangle = cdt.triangle(fh);
				for (size_t i = 0; i < 3; ++i) {

					const size_t im = (i + 2) % 3;
					const size_t ip = (i + 1) % 3;

					const Point_2 &a = triangle.vertex(im);
					const Point_2 &b = triangle.vertex(i);
					const Point_2 &c = triangle.vertex(ip);

					Vector_2 ba = Vector_2(b, a);
					Vector_2 bc = Vector_2(b, c);

					ba /= CGAL::sqrt(ba.squared_length());
					bc /= CGAL::sqrt(bc.squared_length());

					const FT prod = dot_product(bc, ba);

					if (CGAL::abs(FT(1) + prod) < m_angle_eps) return false;
				}

				return true;
			}

			FT compute_face_area(const CDT &cdt, const Face_handle &fh) {

				const auto &triangle = cdt.triangle(fh);
				return triangle.area();
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_VISIBILITY_BASE_2_H