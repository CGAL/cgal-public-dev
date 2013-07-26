#ifndef _DISCRETE_DISTRIBUTION_WITH_FINITE_SUPPORT_GENERATOR_H_
#define _DISCRETE_DISTRIBUTION_WITH_FINITE_SUPPORT_GENERATOR_H_
#include <iostream>
#include <vector>
#include <CGAL/Random.h>

//#define VERBOSE

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
			const int N = input.size();
			typename Container::iterator el_begin = input.begin();
			typename Container::iterator el_end = input.end();
			container.reserve(N);
			presums.reserve(N);
			typename std::vector<Random_generator_with_weight>::iterator it = el_begin;
			for (; it != el_end; it++) {
				container.push_back(Random_generator_with_weight(*it));
			}

			for (int i = 0; i < N; i++) {
				presums.push_back(i == 0 ? container[i].getWeight() :
						container[i].getWeight() + presums[i-1]);
			}
		}

		result_type generate(CGAL::Random &rand) {
			const int N = presums.size();
			typename Container::iterator el_begin = container.begin();
			typename Container::iterator el_end = container.end();
			double tmp_presum = rand.get_double(0, presums[N-1]);
			typename std::vector<double>::iterator SampleIterator =
				upper_bound(presums.begin(), presums.end(),
						tmp_presum);

			int SampleIndex = SampleIterator - presums.begin();
#ifdef VERBOSE
			std::cout << "The picked Element is: " << SampleIndex <<
				std::endl;
#endif
			Random_generator_with_weight SampleElement = *(el_begin + SampleIndex);
#ifdef VERBOSE
			std::cout << "Weight of the picked element: " <<
				SampleElement.getWeight() << std::endl;
#endif

			result_type p = container[SampleIndex].getRand();
#ifdef VERBOSE
			std::cout << "The generated point is " << p.x() << " " <<
				p.y() << " " << p.z() << std::endl;
#endif
			return p;
		}
};
};
};
#endif //_DISCRETE_DISTRIBUTION_WITH_FINITE_SUPPORT_GENERATOR_H_

