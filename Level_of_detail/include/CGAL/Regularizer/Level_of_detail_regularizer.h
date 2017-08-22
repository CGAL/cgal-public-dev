#ifndef CGAL_LEVEL_OF_DETAIL_REGULARIZER_H
#define CGAL_LEVEL_OF_DETAIL_REGULARIZER_H

#include <iostream>
#include <cmath>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer>
		class Level_of_detail_regularizer {

		public:
			typedef KernelTraits    		  Traits;
			typedef typename Traits::Plane_3  Plane;
			typedef InputContainer  		  Container;
			typedef typename Traits::Vector_3 Vector;
			typedef typename Traits::Point_3  Point;

			typename Traits::Compute_squared_length_3 squared_length;
			typename Traits::Compute_scalar_product_3 dot_product;

			template<class Planes>
			auto regularize(Container &, const Planes &) { 

				const auto number_of_regularized_planes = -1;
				// using Const_iterator = typename Planes::const_iterator;

				// Const_iterator it = planes.end();
				// --it;
				// std::cout << input.point((*it).second[1]) << std::endl;

				Vector m = Vector(-0.36, 0.0, -0.03);
				Vector n(0, 0, 1);

				Point p1(1.0, 0.2, 0.9);
				Point p2(1.05, 0.5, 0.3);
				Point p3(1.0, 0.8, 0.9);

				Point newp1, newp2, newp3;

				get_point(m, n, p1, newp1);
				get_point(m, n, p2, newp2);
				get_point(m, n, p3, newp3);

				std::cout << newp1 << std::endl;
				std::cout << newp2 << std::endl;
				std::cout << newp3 << std::endl;

				return number_of_regularized_planes;
			}

		private:
			void get_point(const Vector &m, const Vector &n, const Point &p, Point &newp) {

				const auto cross  = CGAL::cross_product(m, n);
				const auto length = std::sqrt(squared_length(cross));
				const auto dot    = dot_product(m, n);

				// std::cout << std::atan2(length, dot) * 180.0 / M_PI << std::endl; 

				auto angle = -(M_PI / 2.0) + std::atan2(length, dot);

				// const auto cosa = cos(angle);
				// const auto cosa = dot / (std::sqrt(squared_length(m))*std::sqrt(squared_length(n)));

				const auto axis = cross / length;

				const auto c = std::cos(angle);
				const auto s = std::sin(angle);
				const auto C = 1.0 - c;

				const auto x = axis.x();
				const auto y = axis.y();
				const auto z = axis.z();

				newp = 
				Point((x*x*C+c)*p.x()   + (x*y*C-z*s)*p.y() + (x*z*C+y*s)*p.z(),
					  (y*x*C+z*s)*p.x() + (y*y*C+c)*p.y()   + (y*z*C-x*s)*p.z(),
					  (z*x*C-y*s)*p.x() + (z*y*C+x*s)*p.y() + (z*z*C+c)*p.z());

				// std::cout << newp << std::endl;
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_REGULARIZER_H