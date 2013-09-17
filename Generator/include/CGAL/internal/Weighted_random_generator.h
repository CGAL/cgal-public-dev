#ifndef _WEIGHTED_RANDOM_GENERATOR_H_
#define _WEIGHTED_RANDOM_GENERATOR_H_
#include <vector>
#include <iostream>
#include <CGAL/algorithm.h>

namespace CGAL { namespace internal {
template <typename PointGeneratorClass>
class Weighted_random_generator {
	private:
		PointGeneratorClass _rand;
		double _weight;
	public:
		typedef typename PointGeneratorClass::result_type result_type;
		Weighted_random_generator() {}

		Weighted_random_generator(const
				Weighted_random_generator<PointGeneratorClass>
				&x) : _weight(x._weight) {
			CGAL_precondition(_weight >= 0);
			this->_rand = PointGeneratorClass(x._rand);
		}

		Weighted_random_generator(const PointGeneratorClass &rand,
				const double weight) : _weight(weight) {
			CGAL_precondition(_weight >= 0);
			_rand = PointGeneratorClass(rand);
		}

		double getWeight() const {
			return _weight;
		}

		result_type getRand() const {
			std::vector<result_type> output;
			output.reserve(1);
			CGAL::cpp11::copy_n( _rand, 1, std::back_inserter(output));
			result_type ret = output[0];
			output.clear();
			return ret;
		}

		Weighted_random_generator& operator=(const Weighted_random_generator &x) {
			this->_rand = PointGeneratorClass(x._rand);
			this->_weight = x._weight;
			return *this;
		}
};
};
};
#endif //_WEIGHTED_RANDOM_GENERATOR_H_
