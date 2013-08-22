#ifndef _FINITE_SUPPORT_DISTRIBUTION_H_
#define _FINITE_SUPPORT_DISTRIBUTION_H_
#include <iostream>
#include <vector>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

//#define VERBOSE
CGAL::Timer t;

namespace CGAL { namespace internal {
template<typename Weighted_random_generator>
class Finite_support_distribution {
	private:
		Weighted_random_generator* container;
		double* presums;
		int size;
	public:
		typedef Weighted_random_generator* Container;
		typedef typename Weighted_random_generator::result_type result_type;
		Finite_support_distribution() {
			presums = NULL;
			container = NULL;
		}
		Finite_support_distribution(Container input, int size) {
			t.start();
			this->size = size;
			container = new Weighted_random_generator[size];
			presums = new double[size];
			for (int i = 0; i < size; i++) {
				container[i] = input[i];
			}
			for (int i = 0; i < size; i++) {
				presums[i] = (i == 0 ? container[i].getWeight() :
						container[i].getWeight() + presums[i-1]);
			}
		}

		Finite_support_distribution(const Finite_support_distribution &in) {
			size = in.size;
			container = new Weighted_random_generator[size];
			presums = new double[size];
			memmove(container, in.container, size * sizeof(in.container[0]));
			memmove(presums, in.presums, size * sizeof(in.presums[0]));
		}

		~Finite_support_distribution() {
			delete[] container;
			delete[] presums;
		}

		result_type generate(CGAL::Random &rand) {
			double tmp_presum = rand.get_double(0, presums[size-1]);
#ifdef VERBOSE
			std::cout << "Random double: " << tmp_presum << std::endl;
#endif
			double* SampleIterator = std::upper_bound<double*,
				double>(presums,
					presums+size-1, tmp_presum);

			int SampleIndex = SampleIterator - presums;
#ifdef VERBOSE
			std::cout << "The picked Element is: " << SampleIndex <<
				std::endl;
#endif
			Weighted_random_generator SampleElement = container[SampleIndex];
#ifdef VERBOSE
			std::cout << "Weight of the picked element: " <<
				SampleElement.getWeight() << std::endl;
#endif
			result_type p;
			p = container[SampleIndex].getRand();
#ifdef VERBOSE
//			std::cout << "The generated point is " << p.x() << " " <<
//				p.y() << " " << p.z() << std::endl;
#endif
			return p;
		}

		Finite_support_distribution& operator=(const
				Finite_support_distribution &in) {
			this->size = in.size;

			if (presums != NULL) delete[] presums;
			if (container != NULL) delete[] container;
			container = new Weighted_random_generator[size];
			presums = new double[size];
			memmove(container, in.container, size * sizeof(in.container[0]));
			memmove(presums, in.presums, size * sizeof(in.presums[0]));
			return *this;
		}
};
};
};
#endif //FINITE_SUPPORT_DISTRIBUTION_H_

