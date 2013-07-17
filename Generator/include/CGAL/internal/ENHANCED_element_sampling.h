#ifndef _ENHANCED_ELEMENT_SAMPLING_H_
#define _ENHANCED_ELEMENT_SAMPLING_H_
#include <iostream>
#include <vector>
#include <algorithm>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <iterator>
#include <cstdlib>
#include <time.h>
#include <CGAL/Timer.h>
#include <CGAL/internal/weighted_random_element.h>

using namespace std;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;

namespace CGAL { namespace internal {
template<typename Element_RandomAccessIterator, typename VolumeElementFunctor,
	typename PointGeneratorFunctor>
class EnhancedElementSampling {
	private: vector<CGAL::internal::Weighted_random_element<PointGeneratorFunctor> > container;
	public:
		EnhancedElementSampling(int N, Element_RandomAccessIterator el_begin,
				Element_RandomAccessIterator el_end) {
			container.reserve(N);
			Element_RandomAccessIterator it = el_begin;
			int i = 0;
			for (; it != el_end; it++) {
				VolumeElementFunctor volElem(*it);
				double weight = volElem();
				double presum = (i == 0 ? weight : weight +
						container[i-1].getPresum());
				PointGeneratorFunctor randGen(*it);
				CGAL::internal::Weighted_random_element<PointGeneratorFunctor>
					aux(randGen, presum);
				container[i] = aux;
				i++;
			}


			CGAL::Random rand;
			double tmp_presum = rand.get_double(0,
					container[N-1].getPresum());
			CGAL::internal::Weighted_random_element<PointGeneratorFunctor>
				tmp(tmp_presum);
//			cout << tmp << '\n';
			typename vector<CGAL::internal::Weighted_random_element<PointGeneratorFunctor> >::iterator SampleIterator = upper_bound(container.begin(), container.end(), tmp);

			int SampleIndex = SampleIterator - container.begin();
			cout << "The picked Element is: " << SampleIndex << '\n';

			Element_RandomAccessIterator SampleElement = el_begin + SampleIndex;
			cout << "Coords of the picked element: " <<
				SampleElement[0] << " " << SampleElement[1] <<
				" " << SampleElement[2] << " " <<
				SampleElement[3] << '\n';

			Point p = container[SampleIndex].getRand()(6);
			cout << "The generated point is " << p.x() << " " <<
				p.y() << " " << p.z() << '\n';
		}
};
};
};
#endif //_ENHANCED_ELEMENT_SAMPLING_H_
