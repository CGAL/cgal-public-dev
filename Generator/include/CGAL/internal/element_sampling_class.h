#ifndef _ELEMENT_SAMPLING_H_
#define _ELEMENT_SAMPLING_H_
#include <iostream>
#include <vector>
#include <algorithm>
#include <CGAL/Cartesian_d.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Random.h>
#include <CGAL/algorithm.h>
#include <iterator>
#include <cstdlib>
#include <time.h>
#include <CGAL/Timer.h>
#include "weighted_random_element.h"

using namespace std;

namespace CGAL { namespace internal {
template<typename Element_RandomAccessIterator, typename VolumeElementFunctor>
class ElementSampling {
	private:
	public:
		ElementSampling(int N, Element_RandomAccessIterator el_begin,
				Element_RandomAccessIterator el_end) {
			CGAL::Random rand;
			CGAL::Timer tmp;
			double *volumes, *sums;
			int i = 0;
			volumes = new double[N];
			sums = new double[N];
			Element_RandomAccessIterator it = el_begin;
			tmp.start();
			for (; it != el_end; it++) {
				VolumeElementFunctor volElem(*it);
				volumes[i] = volElem();

				cout << i << " " << volumes[i] << '\n';
				i++;
			}
			tmp.stop();
			cout << "Total time: " << tmp.time() << '\n';

			sums[0] = volumes[0];
			for (i = 1; i < N; i++) {
				sums[i] = volumes[i] + sums[i - 1];
				cout << i << "-th sum " << sums[i] << '\n';
			}

			double aux = rand.get_double(0, sums[N-1]);
			cout << aux << '\n';

			double* SampleIterator = upper_bound(sums, sums+N, aux);
			int SampleIndex = SampleIterator - sums;
			cout << "The picked Element is: " << SampleIndex << '\n';

			//Element SampleElement = *(el_begin + SampleIndex);

			//TODO: I haven't done yet the splitting into atomic elements of the
			//picked Element; So far I have only tested on examples where Element is
			//an atomic element itself



			delete[] volumes;
			delete[] sums;
		}
};
};
};
#endif //_ELEMENT_SAMPLING_H_
