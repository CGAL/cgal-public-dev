#include <CGAL/Random.h>
#include <iostream>
#include <CGAL/internal/ALTERED_Discrete_distribution_with_finite_support_generator.h>
#include <CGAL/internal/ALTERED_Random_generator_with_weight.h>
#include <CGAL/point_generators_3.h>

//Random generator
typedef CGAL::internal::ALTERED_Random_generator_with_weight GeneratorWithWeight;

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
		GeneratorWithWeight tmp = GeneratorWithWeight ((float) 1/N);
		containing_structure.push_back(tmp);
	}

	CGAL::internal::ALTERED_Discrete_distribution_with_finite_support_generator<GeneratorWithWeight
		> randomGen(containing_structure);

	int *ret = (int *) calloc(N, sizeof(int));
	for (int i = 0; i < number_points; i++) {
		ret[randomGen.generate(rand)]++;
	}
	for(int i = 0; i < N; i++) {
		std::cout << ret[i] << std::endl;
	}
	free(ret);
	return 0;
}

