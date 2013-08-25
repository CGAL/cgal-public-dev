#ifndef _FINITE_SUPPORT_DISTRIBUTION_H_
#define _FINITE_SUPPORT_DISTRIBUTION_H_
#include <iostream>
#include <vector>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

//#define VERBOSE
CGAL::Timer t;

namespace CGAL { namespace internal {
template<typename Weighted_random_generator, int N>
class Finite_support_distribution {
	private:
		Weighted_random_generator* container;
		double *presums;
		int **offsets;
		double **aux_array;
		int *length; // the number of elements in offsets[i]
		int size;
	public:
		typedef Weighted_random_generator* Container;
		typedef typename Weighted_random_generator::result_type result_type;
		Finite_support_distribution() {
			presums = NULL;
			container = NULL;
			offsets = NULL;
			aux_array = NULL;
			length = NULL;
		}
		Finite_support_distribution(Container input, int size) {
			int parts = N / sizeof(std::vector<Weighted_random_generator>*);
#ifdef VERBOSE
			std::cout << "N = " << N << std::endl;
			std::cout << "sizeof = " << sizeof(std::vector<Weighted_random_generator>*) << std::endl;
			std::cout << "Size = " << size << std::endl;
			std::cout << "Parts = " << parts << std::endl;
#endif
			this->size = size;
			container = new Weighted_random_generator[size];
			presums = new double[size];
			offsets = new int*[parts];
			aux_array = new double*[parts];
			length = new int[parts];
			for (int i = 0; i < size; i++) {
				container[i] = input[i];
			}
			for (int i = 0; i < size; i++) {
				presums[i] = (i == 0 ? container[i].getWeight() :
						container[i].getWeight() + presums[i-1]);
			}
			double step = presums[size-1] / parts;
#ifdef VERBOSE
			std::cout << "Step = " << step << std::endl;
#endif
			int j = 0;
			for (int i = 0; i < parts; i++) {
				offsets[i] = new int[10];
				aux_array[i] = new double[10];
				int k = 0;
				while (presums[j] < (i+1) * step) {
					aux_array[i][k] = presums[j];
					offsets[i][k] = j;
					j++;
					k++;
				}
				aux_array[i][k] = step * (i+1);
				offsets[i][k] = j;
				length[i] = k+1;
			}
		}

		Finite_support_distribution(const Finite_support_distribution &in) {
			size = in.size;
			int parts = N / sizeof(std::vector<Weighted_random_generator>*);
			container = new Weighted_random_generator[size];
			presums = new double[size];
			offsets = new int*[parts];
			aux_array = new double*[parts];
			length = new int[parts];
			memmove(container, in.container, size * sizeof(in.container[0]));
			memmove(presums, in.presums, size * sizeof(double));
			memmove(length, in.length, parts * sizeof(int));
			for (int i = 0; i < parts; i++) {
				offsets[i] = new int[length[i]];
				memmove(offsets[i], in.offsets[i], length[i] *
						sizeof(in.offsets[i][0]));
				aux_array[i] = new double[length[i]];
				memmove(aux_array[i], in.aux_array[i], length[i] *
						sizeof(in.aux_array[i][0]));
			}
		}

		~Finite_support_distribution() {
			int parts = N / sizeof(std::vector<Weighted_random_generator>*);
			delete[] container;
			delete[] presums;
			for (int i = 0; i < parts; i++) {
				delete[] offsets[i];
				delete[] aux_array[i];
			}
			delete[] offsets;
			delete[] aux_array;
			delete[] length;
		}

		Finite_support_distribution& operator=(const Finite_support_distribution &in) {
			this->size = in.size;
			int parts = N / sizeof(std::vector<Weighted_random_generator>*);
			if (presums != NULL) delete[] presums;
			if (container != NULL) delete[] container;
			if (offsets != NULL) delete[] offsets;
			if (aux_array != NULL) delete[] aux_array;
			if (length != NULL) delete[] length;
			container = new Weighted_random_generator[size];
			presums = new double[size];
			offsets = new int*[parts];
			aux_array = new double*[parts];
			length = new int[parts];
			memmove(container, in.container, size * sizeof(in.container[0]));
			memmove(presums, in.presums, size * sizeof(double));
			memmove(length, in.length, parts * sizeof(int));
			for (int i = 0; i < parts; i++) {
				if (offsets[i] != NULL) delete[] offsets[i];
				offsets[i] = new int[length[i]];
				memmove(offsets[i], in.offsets[i], length[i] *
						sizeof(in.offsets[i][0]));
				if (aux_array[i] != NULL) delete[] aux_array[i];
				aux_array[i] = new double[length[i]];
				memmove(aux_array[i], in.aux_array[i], length[i] *
						sizeof(in.aux_array[i][0]));
			}
			return *this;
		}

		result_type generate(CGAL::Random &rand) {
			double tmp_presum = rand.get_double(0, presums[size-1]);
#ifdef VERBOSE
			std::cout << "Random double: " << tmp_presum << std::endl;
#endif
			int parts = N / sizeof(std::vector<Weighted_random_generator>*);
			double step = presums[size-1] / parts;
			int index = tmp_presum / step; // the index of the
						       // picked array
#ifdef VERBOSE
			std::cout<<"Index in the array: "<<index<<std::endl;
#endif
			double *SampleIterator = aux_array[index];
			if (length[index] > 1) {
				SampleIterator = std::upper_bound<double*,
					double>(aux_array[index],
						aux_array[index]+length[index]-1, tmp_presum);
			}
			int SampleIndex = SampleIterator - aux_array[index];
#ifdef VERBOSE
			std::cout << "Length of picked array " << length[index]
				<< std::endl;
			std::cout << "The picked Element is: " << SampleIndex <<
				std::endl;
#endif
			return container[offsets[index][SampleIndex]].getRand();
		}
};
};
};
#endif //FINITE_SUPPORT_DISTRIBUTION_H_

