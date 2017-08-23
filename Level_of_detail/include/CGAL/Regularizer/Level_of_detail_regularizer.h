#ifndef CGAL_LEVEL_OF_DETAIL_REGULARIZER_H
#define CGAL_LEVEL_OF_DETAIL_REGULARIZER_H

// STL includes.
#include <iostream>
#include <cmath>
#include <vector>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/number_utils.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class Planes>
		class Level_of_detail_regularizer {

		public:
			virtual int regularize(Planes&, InputContainer &) const = 0;
		};

		template<class KernelTraits, class InputContainer, class PlanesContainer>
		class Level_of_detail_vertical_regularizer : public Level_of_detail_regularizer<KernelTraits, InputContainer, PlanesContainer> {

		private:
			enum class Regularization_status { REJECT, REGULARIZE, NOACTION };
			const static auto M_NORMAL_INDEX = 0;

		public:
			typedef KernelTraits    Traits;
			typedef InputContainer  Container;
			typedef PlanesContainer Planes;

			typedef typename Traits::FT FT;

			typedef typename Traits::Plane_3  Plane;
			typedef typename Traits::Vector_3 Vector;
			typedef typename Traits::Vector_3 Normal;
			typedef typename Traits::Point_3  Point;

			typename Traits::Compute_squared_length_3 		  squared_length;
			typename Traits::Compute_scalar_product_3 		  dot_product;
			typename Traits::Construct_cross_product_vector_3 cross_product;

			using Const_iterator = typename Planes::const_iterator;
			using Plane_map = typename Container:: template Property_map<Plane>;

			// Here I assume that all points in the plane have the same normals and these normals
			// are the same as the associated plane's normal.
			int regularize(Planes &planes, Container &input) const override { 

				auto number_of_regularized_planes = 0;

				std::vector<bool> rejected(planes.size(), false);
				bool update_planes = false;

				int count = 0;
				for (Const_iterator it = planes.begin(); it != planes.end(); ++it, ++count) {
					
					Normal m;
					set_plane_normal(input, it, m);

					Normal n;
					set_ground_normal(n);

					// Compute rotation angle and rotation axis.
					FT rad_angle; Vector axis;
					const Regularization_status status = check_status(m, n, rad_angle, axis);

					switch (status) {
						
						case Regularization_status::REJECT: {
							
							rejected[count] = true;	
							update_planes   = true;

							break;
						}

						case Regularization_status::REGULARIZE: {
						
							rotate_plane(it, rad_angle, axis, input);
							++number_of_regularized_planes;

							break;
						}

						case Regularization_status::NOACTION:
							break;

						default:
							break;
					}
				}

				if (update_planes) update_rejected_planes(rejected, planes);
				return number_of_regularized_planes;
			}

		private:
			void set_plane_normal(const Container &input, const Const_iterator &it, Normal &plane_normal) const {
				plane_normal = input.normal((*it).second[M_NORMAL_INDEX]);
			} 

			void set_ground_normal(Normal &ground_normal) const {
				ground_normal = Normal(0, 0, 1);
			}

			Regularization_status check_status(const Normal &m, const Normal &n, FT &rad_angle, Vector &axis) const {

				const auto cross  = cross_product(m, n);
				const auto length = CGAL::sqrt(squared_length(cross));
				const auto dot    = dot_product(m, n);

				rad_angle = std::atan2(length, dot) - (FT(M_PI) / FT(2));
				axis = cross / length;

				const FT deg_angle = CGAL::abs(rad_angle * FT(180) / FT(M_PI));

				if (deg_angle < FT(1)) return Regularization_status::NOACTION;
				else if (deg_angle < FT(10)) return Regularization_status::REGULARIZE;
				
				return Regularization_status::REJECT;
			}

			void rotate_plane(const Const_iterator &it, const FT &rad_angle, const Vector &axis, Container &input) const {

				const double angle = static_cast<double>(rad_angle);

				const auto c = static_cast<FT>(std::cos(angle));
				const auto s = static_cast<FT>(std::sin(angle));

				const auto C = FT(1) - c;

				Plane_map planes;
				boost::tie(planes, boost::tuples::ignore) = input. template property_map<Plane>("plane");
				
				rotate_element(c, s, C, axis, input.normal((*it).second[M_NORMAL_INDEX]));
				for (size_t i = 0; i < (*it).second.size(); ++i) {
					
					// Update point positions and normals.
					const auto index = (*it).second[i];
					rotate_element(c, s, C, axis, input.point(index));
					input.normal(index) = input.normal((*it).second[M_NORMAL_INDEX]); // Here we enforce the same normals for all points!

					// Update associated planes. Most-likely can be removed later.
					planes[index] = Plane(input.point(index), input.normal(index));
				}
			}

			template<class Element>
			void rotate_element(const FT c, const FT s, const FT C, const Vector &axis, Element &e) const {

				const auto x = axis.x();
				const auto y = axis.y();
				const auto z = axis.z();

				e = Element((x * x * C + c)     * e.x() + (x * y * C - z * s) * e.y() + (x * z * C + y * s) * e.z(),
					  		(y * x * C + z * s) * e.x() + (y * y * C + c)     * e.y() + (y * z * C - x * s) * e.z(),
					  		(z * x * C - y * s) * e.x() + (z * y * C + x * s) * e.y() + (z * z * C + c)     * e.z());
			}

			void update_rejected_planes(const std::vector<bool> &rejected, Planes &planes) const {

				Planes new_planes;
				assert(planes.size() == rejected.size());

				size_t i = 0;
				for (Const_iterator it = planes.begin(); it != planes.end(); ++i, ++it)
					if (!rejected[i]) new_planes[(*it).first] = (*it).second;

				planes = new_planes;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_REGULARIZER_H