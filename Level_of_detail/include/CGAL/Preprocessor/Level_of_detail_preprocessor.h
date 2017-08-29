#ifndef CGAL_LEVEL_OF_DETAIL_PREPROCESSOR_H
#define CGAL_LEVEL_OF_DETAIL_PREPROCESSOR_H

// Boost includes.
#include <boost/tuple/tuple.hpp>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer>
		class Level_of_detail_preprocessor {

		public:
			typedef KernelTraits   Traits;
			typedef InputContainer Container;

			using Index          = int;
			using Const_iterator = typename Container::const_iterator;
			using Index_map      = typename Container:: template Property_map<Index>;

			template<class Planes>
			auto get_planes(const Container &input, Planes &planes) {

				auto number_of_planes = -1;
				create_indices(input);
				
				planes.clear();

				for (Const_iterator it = input.begin(); it != input.end(); ++it)
					if (m_indices[*it] >= 0) 
						planes[m_indices[*it]].push_back(*it);

				number_of_planes = planes.size();
				return number_of_planes;
			}

			template<class Indices, class Planes>
			auto get_planes(const Container &input, const Indices &mapping, Planes &planes) {

				auto number_of_planes = -1;
				create_indices(input);

				planes.clear();

				for (size_t i = 0; i < mapping.size(); ++i)
					if (m_indices[mapping[i]] >= 0)
						planes[m_indices[mapping[i]]].push_back(mapping[i]);

				number_of_planes = planes.size();
				return number_of_planes;
			}

		private:
			Index_map m_indices;

			void create_indices(const Container &input) {
				boost::tie(m_indices,  boost::tuples::ignore) = input. template property_map<Index>("index");
			}
		};
	}
}

#endif // CGAL_LEVEL_OF_DETAIL_PREPROCESSOR_H