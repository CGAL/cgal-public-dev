#include <CGAL/generators.h>
#include <CGAL/function_objects.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/Random.h>
#include <iostream>
#include <CGAL/internal/Discrete_distribution_with_finite_support_generator.h>
#include <CGAL/internal/Random_generator_with_weight.h>
#include <CGAL/point_generators_3.h>
#include <vector>
#include <algorithm>

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

template<class P >
void Probability_1_generator<P>::generate_point() {
	this->d_item = _i;
}

//Random generator
typedef
CGAL::internal::Random_generator_with_weight<Probability_1_generator<int> > GeneratorWithWeight;

int main()
{
	CGAL::Random rand;
	const int MIN_N = 90;
	const int MAX_N = 100;
	const int N = rand.get_int(MIN_N, MAX_N);
	std::cout << "N = " << N << std::endl;
	const int MIN_POINTS = 1000;
	const int MAX_POINTS = 1000000;
	const int number_points = rand.get_int(MIN_POINTS, MAX_POINTS);
	std::cout << "number_points = " << number_points << std::endl;
	std::vector<GeneratorWithWeight> containing_structure;
	containing_structure.reserve(N);
	for (int i = 0; i < N; i++) {
		int aux = i;
		Probability_1_generator<int> randGen( aux );
		GeneratorWithWeight tmp = GeneratorWithWeight (randGen, (double) 1/N);
		containing_structure.push_back(tmp);
	}

	CGAL::internal::Discrete_distribution_with_finite_support_generator<GeneratorWithWeight
		> randomGen(containing_structure);

	int *ret = (int *) calloc(N, sizeof(int));
	for (int i = 0; i < number_points; i++) {
		int index = randomGen.generate(rand);
		ret[index]++;
	}
	for(int i = 0; i < N; i++) {
		std::cout << ret[i] << std::endl;
	}
	free(ret);
	return 0;
}

