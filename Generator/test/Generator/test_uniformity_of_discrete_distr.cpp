#include <CGAL/generators.h>
#include <CGAL/function_objects.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/Random.h>
#include <iostream>
#include <CGAL/internal/Finite_support_distribution.h>
#include <CGAL/internal/Weighted_random_generator.h>
#include <CGAL/point_generators_3.h>
#include <vector>
#include <algorithm>
#include <cassert>

//#define TEST_VERBOSE

//Probability 1 generator
template < class P >
class Probability_1_generator : public CGAL::Random_generator_base<P> {
	P _i;
	void generate_point();
public:
	typedef P result_type;
	typedef Probability_1_generator<P> This;
	Probability_1_generator() {}
	Probability_1_generator( const This& x,CGAL::Random& rnd = CGAL::default_random)
	: CGAL::Random_generator_base<P>( 1, rnd ),_i(x._i) {
		generate_point();
	}
	Probability_1_generator( const P& i, CGAL::Random& rnd = CGAL::default_random)
	: CGAL::Random_generator_base<P>( 1, rnd ),_i(i) {
		generate_point();
	}
	This operator=(This x) {
		_i = x._i;
		return *this;
	}
	This& operator++() {
		generate_point();
		return *this;
	}
	This operator++(int) {
		This tmp = *this;
		++(*this);
		return tmp;
	}
};

template<class P>
void Probability_1_generator<P>::generate_point() {
	this->d_item = _i;
}

bool inside_generation_range(const int lower, const int upper, int gen) {
	if (gen < lower || gen > upper) {
		return false;
	}
	return true;
}

bool all_posibilities_covered(int *begin, int size) {
	for (int i = 0; i < size; i++) {
		if (!begin[i]) return false;
	}
	return true;
}

bool is_uniform(int *begin, int size, int average) {
	for(int i = 0; i < size; i++) {
		if (1.20 * average < begin[i] || 0.80 * average > begin[i])
		{
			std::cout << average << " " << begin[i] << " " << i<<'\n';
			return false;
		}
	}
	return true;
}

//Random generator
typedef
CGAL::internal::Weighted_random_generator<Probability_1_generator<int> > GeneratorWithWeight;

int main()
{
// Testing for small values of N
	for (int cont = 1; cont < 10; cont++) {
		CGAL::Random rand;
		const int N = cont;
#ifdef TEST_VERBOSE
		std::cout << "N = " << N << std::endl;
#endif
		const int MIN_POINTS = 1000;
		const int MAX_POINTS = 1000000;
		const int number_points = rand.get_int(MIN_POINTS, MAX_POINTS);
#ifdef TEST_VERBOSE
		std::cout << "number_points = " << number_points << std::endl;
#endif
		for (int cont2 = 1; cont2 < 10; cont2++) {
			const int parts = cont2;
			GeneratorWithWeight *containing_structure = new GeneratorWithWeight[N];
			for (int i = 0; i < N; i++) {
				int aux = i;
				Probability_1_generator<int> randGen( aux );
				GeneratorWithWeight tmp = GeneratorWithWeight (randGen, (double) 1/N);
				containing_structure[i] = tmp;
			}

			CGAL::internal::Finite_support_distribution<GeneratorWithWeight
				> randomGen(containing_structure, N, parts);

			int *ret = (int *) calloc(N, sizeof(int));
			for (int i = 0; i < number_points; i++) {
				int index = randomGen.generate(rand);
				assert(inside_generation_range(0, N, index));
				ret[index]++;
			}

#ifdef TEST_VERBOSE
			for(int i = 0; i < N; i++) {
				std::cout << ret[i] << std::endl;
			}
#endif
	//		int avg = 0;
	//		for(int i = 0; i < N; i++) {
	//			avg += ret[i];
	//		}
	//		avg = avg / N;
	//		assert(is_uniform(ret, N, avg));

			assert(all_posibilities_covered(ret, N));
			free(ret);
		}
	}

// Testing for small values of parts and N
	for (int cont = 1; cont < 10; cont++) {
		CGAL::Random rand;
		const int N = cont;
#ifdef TEST_VERBOSE
		std::cout << "N = " << N << std::endl;
#endif
		const int MIN_POINTS = 1000;
		const int MAX_POINTS = 1000000;
		const int number_points = rand.get_int(MIN_POINTS, MAX_POINTS);
#ifdef TEST_VERBOSE
		std::cout << "number_points = " << number_points << std::endl;
#endif
		const int MIN_PARTS = 1;
		const int MAX_PARTS = 1<<20;
		const int parts = rand.get_int(MIN_PARTS, MAX_PARTS);
		GeneratorWithWeight *containing_structure = new GeneratorWithWeight[N];
		for (int i = 0; i < N; i++) {
			int aux = i;
			Probability_1_generator<int> randGen( aux );
			GeneratorWithWeight tmp = GeneratorWithWeight (randGen, (double) 1/N);
			containing_structure[i] = tmp;
		}

		CGAL::internal::Finite_support_distribution<GeneratorWithWeight
			> randomGen(containing_structure, N, parts);

		int *ret = (int *) calloc(N, sizeof(int));
		for (int i = 0; i < number_points; i++) {
			int index = randomGen.generate(rand);
			assert(inside_generation_range(0, N, index));
			ret[index]++;
		}

#ifdef TEST_VERBOSE
		for(int i = 0; i < N; i++) {
			std::cout << ret[i] << std::endl;
		}
#endif
//		int avg = 0;
//		for(int i = 0; i < N; i++) {
//			avg += ret[i];
//		}
//		avg = avg / N;
//		assert(is_uniform(ret, N, avg));

		assert(all_posibilities_covered(ret, N));
		free(ret);
	}

// Testing for small values of number_points
	for (int cont = 1; cont < 10; cont++) {
		CGAL::Random rand;
		const int MIN_N = 1;
		const int MAX_N = cont+1;
		const int N = rand.get_int(MIN_N, MAX_N);
#ifdef TEST_VERBOSE
		std::cout << "N = " << N << std::endl;
#endif
		const int number_points = cont;
#ifdef TEST_VERBOSE
		std::cout << "number_points = " << number_points << std::endl;
#endif
		const int MIN_PARTS = 1;
		const int MAX_PARTS = 1<<20;
		const int parts = rand.get_int(MIN_PARTS, MAX_PARTS);
		GeneratorWithWeight *containing_structure = new GeneratorWithWeight[N];
		for (int i = 0; i < N; i++) {
			int aux = i;
			Probability_1_generator<int> randGen( aux );
			GeneratorWithWeight tmp = GeneratorWithWeight (randGen, (double) 1/N);
			containing_structure[i] = tmp;
		}

		CGAL::internal::Finite_support_distribution<GeneratorWithWeight
			> randomGen(containing_structure, N, parts);

		int *ret = (int *) calloc(N, sizeof(int));
		for (int i = 0; i < number_points; i++) {
			int index = randomGen.generate(rand);
			assert(inside_generation_range(0, N, index));
			ret[index]++;
		}

#ifdef TEST_VERBOSE
		for(int i = 0; i < N; i++) {
			std::cout << ret[i] << std::endl;
		}
#endif
//		int avg = 0;
//		for(int i = 0; i < N; i++) {
//			avg += ret[i];
//		}
//		avg = avg / N;
//		assert(is_uniform(ret, N, avg));

		assert(all_posibilities_covered(ret, N));
		free(ret);
	}

	for (int cont = 0; cont < 10; cont++) {
		CGAL::Random rand;
		const int MIN_N = 90;
		const int MAX_N = 100;
		const int N = rand.get_int(MIN_N, MAX_N);
#ifdef TEST_VERBOSE
		std::cout << "N = " << N << std::endl;
#endif
		const int MIN_POINTS = 1000;
		const int MAX_POINTS = 1000000;
		const int number_points = rand.get_int(MIN_POINTS, MAX_POINTS);
#ifdef TEST_VERBOSE
		std::cout << "number_points = " << number_points << std::endl;
#endif
		const int MIN_PARTS = 1;
		const int MAX_PARTS = 1<<20;
		const int parts = rand.get_int(MIN_PARTS, MAX_PARTS);
		GeneratorWithWeight *containing_structure = new GeneratorWithWeight[N];
		for (int i = 0; i < N; i++) {
			int aux = i;
			Probability_1_generator<int> randGen( aux );
			GeneratorWithWeight tmp = GeneratorWithWeight (randGen, (double) 1/N);
			containing_structure[i] = tmp;
		}

		CGAL::internal::Finite_support_distribution<GeneratorWithWeight
			> randomGen(containing_structure, N, parts);

		int *ret = (int *) calloc(N, sizeof(int));
		for (int i = 0; i < number_points; i++) {
			int index = randomGen.generate(rand);
			assert(inside_generation_range(0, N, index));
			ret[index]++;
		}

#ifdef TEST_VERBOSE
		for(int i = 0; i < N; i++) {
			std::cout << ret[i] << std::endl;
		}
#endif
//		int avg = 0;
//		for(int i = 0; i < N; i++) {
//			avg += ret[i];
//		}
//		avg = avg / N;
//		assert(is_uniform(ret, N, avg));

		assert(all_posibilities_covered(ret, N));
		free(ret);
		
// Testing the copy-constructor of FSD
		CGAL::internal::Finite_support_distribution<GeneratorWithWeight
			> randomGen_copy(randomGen);

		ret = (int *) calloc(N, sizeof(int));
		for (int i = 0; i < number_points; i++) {
			int index = randomGen.generate(rand);
			assert(inside_generation_range(0, N, index));
			ret[index]++;
		}

#ifdef TEST_VERBOSE
		for(int i = 0; i < N; i++) {
			std::cout << ret[i] << std::endl;
		}
#endif
//		int avg = 0;
//		for(int i = 0; i < N; i++) {
//			avg += ret[i];
//		}
//		avg = avg / N;
//		assert(is_uniform(ret, N, avg));

		assert(all_posibilities_covered(ret, N));
		free(ret);
	}
	return 0;
}

