#ifndef _ALTERED_RANDOM_GENERATOR_WITH_WEIGHT_H_
#define _ALTERED_RANDOM_GENERATOR_WITH_WEIGHT_H_
#include <vector>
#include <CGAL/algorithm.h>
#include <iostream>

namespace CGAL { namespace internal {
class ALTERED_Random_generator_with_weight {
	private:
		double _weight;
	public:
		ALTERED_Random_generator_with_weight() {}

		ALTERED_Random_generator_with_weight(const ALTERED_Random_generator_with_weight &x) {
			this->_weight = x._weight;
		}

		ALTERED_Random_generator_with_weight(double weight) {
			_weight = weight;
		}

		double getWeight() {
			return _weight;
		}

		ALTERED_Random_generator_with_weight& operator=(const ALTERED_Random_generator_with_weight &x) {
			this->_weight = x._weight;
			return *this;
		}
};
};
};
#endif //_ALTERED_RANDOM_GENERATOR_WITH_WEIGHT_H_

