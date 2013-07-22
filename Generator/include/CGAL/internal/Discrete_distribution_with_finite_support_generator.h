#ifndef _DISCRETE_DISTRIBUTION_WITH_FINITE_SUPPORT_GENERATOR_H_
#define _DISCRETE_DISTRIBUTION_WITH_FINITE_SUPPORT_GENERATOR_H_
#include <iostream>
#include <vector>
#include <algorithm>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <iterator>
#include <cstdlib>

namespace CGAL { namespace internal {
template<typename Random_generator_with_weight>
class Discrete_distribution_with_finite_support_generator {
	private:
		std::vector<Random_generator_with_weight> container;
		std::vector<double> presums;
	public:
		typedef std::vector<Random_generator_with_weight> Container;
		typedef typename Random_generator_with_weight::result_type result_type;
		Discrete_distribution_with_finite_support_generator(Container &input) {
			int N = input.size();
			typename Container::iterator el_begin = input.begin();
			typename Container::iterator el_end = input.end();
			container.reserve(N);
			presums.reserve(N);
			typename std::vector<Random_generator_with_weight>::iterator it = el_begin;
			int i = 0;
			for (; it != el_end; it++) {
				container[i] = Random_generator_with_weight(*it);
				i++;
			}

			for (i = 0; i < N; i++) {
				presums[i] = (i == 0 ? container[i].getWeight() :
						container[i].getWeight() + presums[i-1]);
			}
		}

		void generate() {
			int N = presums.size();
			typename Container::iterator el_begin = container.begin();
			typename Container::iterator el_end = container.end();
			CGAL::Random rand;
			double tmp_presum = rand.get_double(0, presums[N-1]);
//			cout << tmp << '\n';
			//typename vector<CGAL::internal::Weighted_random_element<PointGeneratorFunctor> >::iterator SampleIterator = upper_bound(container.begin(), container.end(), tmp);
			typename std::vector<double>::iterator SampleIterator =
				upper_bound(presums.begin(), presums.end(),
						tmp_presum);

			int SampleIndex = SampleIterator - presums.begin();
			std::cout << "The picked Element is: " << SampleIndex <<
				std::endl;

			Random_generator_with_weight SampleElement = *(el_begin + SampleIndex);
			std::cout << "Coords of the picked element: " <<
				SampleElement.getWeight() << std::endl;

			result_type p = container[SampleIndex].getRand();
			std::cout << "The generated point is " << p.x() << " " <<
				p.y() << " " << p.z() << '\n';
		}
};
};
};
#endif //_DISCRETE_DISTRIBUTION_WITH_FINITE_SUPPORT_GENERATOR_H_

