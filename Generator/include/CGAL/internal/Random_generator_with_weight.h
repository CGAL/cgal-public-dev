#ifndef _RANDOM_GENERATOR_WITH_WEIGHT_H_
#define _RANDOM_GENERATOR_WITH_WEIGHT_H_
#include <vector>
#include <iostream>
#include <CGAL/algorithm.h>

namespace CGAL { namespace internal {
template <typename PointGeneratorClass>
class Random_generator_with_weight {
	private:
		PointGeneratorClass _rand;
		double _weight;
	public:
		typedef typename PointGeneratorClass::result_type result_type;
		Random_generator_with_weight() {}

		Random_generator_with_weight(const Random_generator_with_weight<PointGeneratorClass> &x) {
			this->_rand = PointGeneratorClass(x._rand);
			this->_weight = x._weight;
		}

		Random_generator_with_weight(const PointGeneratorClass &rand,
				const double weight) {
			_rand = PointGeneratorClass(rand);
			_weight = weight;
		}

		double getWeight() const {
			return _weight;
		}

		result_type getRand() const {
			std::vector<result_type> output;
			output.reserve(1);
			CGAL::cpp11::copy_n( _rand, 1, std::back_inserter(output));
			return output[0];
		}

		Random_generator_with_weight& operator=(const Random_generator_with_weight &x) {
			this->_rand = PointGeneratorClass(x._rand);
			this->_weight = x._weight;
			return *this;
		}
};
};
};
#endif //_RANDOM_GENERATOR_WITH_WEIGHT_H_

