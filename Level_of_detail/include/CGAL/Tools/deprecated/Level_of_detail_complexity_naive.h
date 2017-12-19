#ifndef CGAL_LEVEL_OF_DETAIL_COMPLEXITY_NAIVE_H
#define CGAL_LEVEL_OF_DETAIL_COMPLEXITY_NAIVE_H

// STL includes.
#include <string>
#include <vector>
#include <cassert>
#include <map>

// New CGAL includes.
#include <CGAL/Selector/Level_of_detail_selector.h>
#include <CGAL/Selector/Level_of_detail_selection_strategy.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class LodReconstruction>
		class Level_of_detail_complexity_naive {

		public:
			typedef KernelTraits   	  Kernel;
			typedef InputContainer    Container;
			typedef LodReconstruction LODS;

			typedef typename Kernel::FT FT;

			using Index   = int;
			using Indices = std::vector<Index>;

			typedef CGAL::LOD::Level_of_detail_building_boundary<Kernel, Container> Boundary_point_counter_strategy;
			typedef CGAL::LOD::Level_of_detail_building_interior<Kernel, Container> Interior_point_counter_strategy;

			typedef CGAL::LOD::Level_of_detail_selector<Kernel, Boundary_point_counter_strategy> Boundary_point_counter;
			typedef CGAL::LOD::Level_of_detail_selector<Kernel, Interior_point_counter_strategy> Interior_point_counter;

			enum class Normalization_type { POINT_BASED, FIXED_ONE };
			enum class Estimation_method  { UNIFORM, SPLIT };


			Level_of_detail_complexity_naive(const Container &input, const LODS &lods) 
			: m_input(input), m_lods(lods), m_complexity(-FT(1)), m_normalize(false),
			m_normalization_type(Normalization_type::POINT_BASED),
			m_estimation_method(Estimation_method::SPLIT) { }


			void estimate() {

				switch (m_estimation_method) {

					case Estimation_method::SPLIT:
						estimate_split();
						break;

					case Estimation_method::UNIFORM:
						estimate_uniform();
						break;

					default:
						assert(!"Wrong estimation method!");
						break;
				}
			}

			FT get() const {
				return m_complexity;
			}


		private:
			const Container &m_input;
			const LODS 		&m_lods;
			
			FT m_complexity;
			bool m_normalize;

			Normalization_type m_normalization_type;
			Estimation_method  m_estimation_method;


			void estimate_uniform() {

				const FT num_elements 		  = get_number_of_elements();
				const FT normalization_factor = get_normalization_factor();

				m_complexity = num_elements / normalization_factor;
			}

			void estimate_split() {

				FT num_roofs = get_number_of_roofs();
				if (m_normalize) num_roofs /= get_roofs_normalization_factor();

				FT num_walls = get_number_of_walls();
				if (m_normalize) num_walls /= get_walls_normalization_factor();

				m_complexity = num_roofs + num_walls;
			}

			FT get_number_of_roofs() const {

				const int num_roofs = m_lods.get_number_of_roofs();
				return static_cast<FT>(num_roofs);
			}

			FT get_number_of_walls() const {
				
				const int num_walls = m_lods.get_number_of_walls();
				return static_cast<FT>(num_walls);
			}

			FT get_number_of_elements() const {
				return get_number_of_roofs() + get_number_of_walls();
			}

			FT get_normalization_factor() const {

				switch (m_normalization_type) {
					case Normalization_type::POINT_BASED:
						return get_point_based_normalization_factor();

					case Normalization_type::FIXED_ONE:
						return get_fixed_one_normalization_factor();

					default:
						assert(!"Wrong normalization type!");
						return -FT(1);
				}
			}

			FT get_fixed_one_normalization_factor() const {
				return FT(1);
			}

			FT get_roofs_normalization_factor() const {			
				
				Indices indices;
				Interior_point_counter counter;
				counter.select_elements(m_input, std::back_inserter(indices));
	
				const int num_roofs_points = indices.size();
				return static_cast<FT>(num_roofs_points);
			}

			FT get_walls_normalization_factor() const {			
				
				Indices indices;
				Boundary_point_counter counter;
				counter.select_elements(m_input, std::back_inserter(indices));
	
				const int num_walls_points = indices.size();
				return static_cast<FT>(num_walls_points);
			}

			FT get_point_based_normalization_factor() const {
				return get_roofs_normalization_factor() + get_walls_normalization_factor();
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_NAIVE_H