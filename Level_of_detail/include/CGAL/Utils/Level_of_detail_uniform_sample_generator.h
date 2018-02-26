#ifndef CGAL_LEVEL_OF_DETAIL_UNIFORM_SAMPLE_GENERATOR_H
#define CGAL_LEVEL_OF_DETAIL_UNIFORM_SAMPLE_GENERATOR_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\" 
#else 
#define PSR "/" 
#endif 

// STL includes.
#include <map>
#include <cassert>
#include <vector>
#include <memory>
#include <string>

// CGAL includes.
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/utils.h>
#include <CGAL/Random.h>
#include <CGAL/number_utils.h>
#include <CGAL/constructions_d.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/property_map.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>

// Boost includes.
#include <boost/tuple/tuple.hpp>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits>
		class Uniform_sample_generator {

		public:
			typedef KernelTraits  			   Kernel;
			typedef typename Kernel::FT 	   FT;
			typedef typename Kernel::Point_2   Point_2;
			typedef typename Kernel::Vector_2  Vector_2;
			typedef typename Kernel::Segment_2 Ray;

			Uniform_sample_generator() : m_num_samples(3), m_num_rays(2) { }

			void set_number_of_samples(const size_t new_value) {
				
				assert(new_value >= 0);
				m_num_samples = new_value;
			}

			void set_number_of_rays(const size_t new_value) {

				assert(new_value > 0);
				m_num_rays = new_value;
			}

			template<class Samples>
			void create_uniform_subdivision_samples(const Point_2 &a, const Point_2 &b, const Point_2 &c, Samples &samples) {

				assert(m_num_samples < 10); // here m_num_samples means the number of subdivision steps

				samples.clear();
				size_t flood_count = 0;

				flood_triangles(a, b, c, samples, flood_count);
			}

			template<class Samples>
			void create_random_uniform_samples_0(const Point_2 &a, const Point_2 &b, const Point_2 &c, Samples &samples) {

				assert(m_num_samples > 0);

				samples.clear();
				samples.resize(m_num_samples);

				assert(!samples.empty() && samples.size() == m_num_samples);

				const Vector_2 ab = b - a;
				const Vector_2 ac = c - a;

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

			template<class Samples>
			void create_random_uniform_samples_1(const Point_2 &a, const Point_2 &b, const Point_2 &c, Samples &samples) {

				assert(m_num_samples > 0);

				samples.clear();
				samples.resize(m_num_samples);

				assert(!samples.empty() && samples.size() == m_num_samples);

				for (size_t i = 0; i < m_num_samples; ++i) {

					FT s = static_cast<FT>(m_rand.uniform_real(0.0, 1.0));
					FT t = static_cast<FT>(m_rand.uniform_real(0.0, 1.0));

					FT u =  FT(1)      - static_cast<FT>(CGAL::sqrt(CGAL::to_double(t)));
					FT v = (FT(1) - s) * static_cast<FT>(CGAL::sqrt(CGAL::to_double(t)));
					FT w =  		s  * static_cast<FT>(CGAL::sqrt(CGAL::to_double(t)));

					samples[i] = Point_2(u * a.x() + v * b.x() + w * c.x(), u * a.y() + v * b.y() + w * c.y());
				}
			}

			template<class Rays>
			void create_random_uniform_rays(const Point_2 &source, const Point_2 &a, const Point_2 &b, Rays &rays) {

				assert(m_num_rays > 0);

				rays.clear();
				rays.resize(m_num_rays);

				for (size_t i = 0; i < m_num_rays; ++i) {

					const FT s = static_cast<FT>(m_rand.uniform_real(0.0, 1.0));
					const FT t = FT(1) - s;

					const FT x = s * a.x() + t * b.x();
					const FT y = s * a.y() + t * b.y();

					const Point_2 target = Point_2(x, y);
					rays[i] = Ray(source, target);
				}
			}

			template<class Rays>
			void create_uniform_rays(const Point_2 &source, const Point_2 &a, const Point_2 &b, Rays &rays) {

				assert(m_num_rays > 0);

				rays.clear();
				rays.resize(m_num_rays);

				const FT num_rays = static_cast<FT>(m_num_rays + 1);

				size_t count = 0;
				for (size_t i = 1; i < num_rays; ++i) {

					const FT t = static_cast<FT>(i) / num_rays;
					const FT s = FT(1) - t;

					const FT x = s * a.x() + t * b.x();
					const FT y = s * a.y() + t * b.y();

					const Point_2 target = Point_2(x, y);

					assert(count < rays.size());
					rays[count++] = Ray(source, target);
				}
			}

		private:
			size_t m_num_samples;
			size_t m_num_rays;
			
			CGAL::Random m_rand;

			template<class Samples>
			void flood_triangles(const Point_2 &a, const Point_2 &b, const Point_2 &c, Samples &samples, size_t flood_count) {

				if (flood_count >= m_num_samples) {
					
					// Insert new sample.
					const FT third = FT(1) / FT(3);

					const FT x = a.x() + b.x() + c.x();
					const FT y = a.y() + b.y() + c.y();

					const Point_2 new_sample = Point_2(x * third, y * third);
					samples.push_back(new_sample);

					return;
				}

				++flood_count;

				// Subdivide.
				std::vector<Point_2> tri(3);
				tri[0] = a; tri[1] = b; tri[2] = c;

				std::vector<Point_2> mids(3);
				const FT half = FT(1) / FT(2);

				for (size_t i = 0; i < 3; ++i) {
					const size_t ip = (i + 1) % 3;					

					const FT x = tri[i].x() + tri[ip].x();
					const FT y = tri[i].y() + tri[ip].y();

					mids[i] = Point_2(half * x, half * y);
				}

				for (size_t i = 0; i < 3; ++i) {
					
					const size_t im = (i + 2) % 3;
					flood_triangles(mids[im], tri[i], mids[i], samples, flood_count);
				}

				flood_triangles(mids[0], mids[1], mids[2], samples, flood_count);
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_UNIFORM_SAMPLE_GENERATOR_H