#ifndef _WEIGHTED_RANDOM_ELEMENT_H_
#define _WEIGHTED_RANDOM_ELEMENT_H_
namespace CGAL { namespace internal {
template <typename MyRandom>
class Weighted_random_element {
	private:
		MyRandom _rand;
		double _presum;
		double _weight; //volume
	public:
		Weighted_random_element() {}

		Weighted_random_element(Weighted_random_element &x) {
			this->_weight = x._weight;
			this->_rand = x._rand;
			this->_presum = x._presum;
		}

		Weighted_random_element(MyRandom &rand, double presum, double
				weight) {
			_rand = rand;
			_presum = presum;
			_weight = weight;
		}

		double getPresum() {
			return _presum;
		}

		MyRandom getRand() {
			return _rand;
		}

		bool operator< (Weighted_random_element &rhs) {
			return this->_weight < rhs._weight;
		}

		Weighted_random_element& operator=(Weighted_random_element &x) {
			this->_weight = x._weight;
			this->_rand = x._rand;
			this->_presum = x._presum;
			return *this;
		}
};
};
};
#endif //_WEIGHTED_RANDOM_ELEMENT_H_
