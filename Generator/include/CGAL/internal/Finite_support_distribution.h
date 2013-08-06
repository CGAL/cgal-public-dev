#ifndef _FINITE_SUPPORT_DISTRIBUTION_H_
#define _FINITE_SUPPORT_DISTRIBUTION_H_
#include <iostream>
#include <vector>
#include <CGAL/Random.h>

//#define VERBOSE

namespace CGAL { namespace internal {
template<typename Weighted_random_generator>
class Finite_support_distribution {
	private:
		std::vector<Weighted_random_generator> container;
		std::vector<double> presums;
	public:
		typedef std::vector<Weighted_random_generator> Container;
		typedef typename Weighted_random_generator::result_type result_type;
		Finite_support_distribution() {}
		Finite_support_distribution(Container &input) {
			const int N = input.size();
			typename Container::iterator el_begin = input.begin();
			typename Container::iterator el_end = input.end();
			container.reserve(N);
			presums.reserve(N);
			typename std::vector<Weighted_random_generator>::iterator it = el_begin;
			for (; it != el_end; it++) {
				container.push_back(Weighted_random_generator(*it));
			}

			for (int i = 0; i < N; i++) {
				presums.push_back(i == 0 ? container[i].getWeight() :
						container[i].getWeight() + presums[i-1]);
			}
		}

		Finite_support_distribution(const Finite_support_distribution &in) {
			container = Container(in.container);
			presums = std::vector<double>(in.presums);
			/*
			const int N = input.size();
			typename Container::iterator el_begin = input.begin();
			typename Container::iterator el_end = input.end();
			container.reserve(N);
			presums.reserve(N);
			typename std::vector<Weighted_random_generator>::iterator it = el_begin;
			for (; it != el_end; it++) {
				container.push_back(Weighted_random_generator(*it));
			}

			for (int i = 0; i < N; i++) {
				presums.push_back(i == 0 ? container[i].getWeight() :
						container[i].getWeight() + presums[i-1]);
			}
			*/
		}

		result_type generate(CGAL::Random &rand) {
			const int N = presums.size();
			typename Container::iterator el_begin = container.begin();
			typename Container::iterator el_end = container.end();
			double tmp_presum = rand.get_double(0, presums[N-1]);
#ifdef VERBOSE
			std::cout << "Random double: " << tmp_presum << std::endl;
#endif
			typename std::vector<double>::iterator SampleIterator =
				upper_bound(presums.begin(), presums.end(),
						tmp_presum);

			int SampleIndex = SampleIterator - presums.begin();
#ifdef VERBOSE
			std::cout << "The picked Element is: " << SampleIndex <<
				std::endl;
#endif
			Weighted_random_generator SampleElement = *(el_begin + SampleIndex);
#ifdef VERBOSE
			std::cout << "Weight of the picked element: " <<
				SampleElement.getWeight() << std::endl;
#endif

			result_type p = container[SampleIndex].getRand();
#ifdef VERBOSE
//			std::cout << "The generated point is " << p.x() << " " <<
//				p.y() << " " << p.z() << std::endl;
#endif
			return p;
		}

		Finite_support_distribution& operator=(const
				Finite_support_distribution &in) {
			this->container = Container(in.container);
			this->presums = std::vector<double>(in.presums);
			return *this;
		}
};
};
};
#endif //FINITE_SUPPORT_DISTRIBUTION_H_

