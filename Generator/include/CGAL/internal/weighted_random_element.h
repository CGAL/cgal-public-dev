#ifndef _WEIGHTED_RANDOM_ELEMENT_H_
#define _WEIGHTED_RANDOM_ELEMENT_H_
namespace CGAL { namespace internal {
template <typename MyRandom>
class Weighted_random_element {
	private:
		MyRandom _rand;
		double _presum; //sum of volumes
	public:
		Weighted_random_element() {}

		Weighted_random_element(Weighted_random_element &x) {
			this->_rand = x._rand;
			this->_presum = x._presum;
		}

		Weighted_random_element(MyRandom &rand, double presum) {
			_rand = rand;
			_presum = presum;
		}

		double getPresum() {
			return _presum;
		}

		MyRandom getRand() {
			return _rand;
		}

		bool operator< (Weighted_random_element &rhs) {
			return this->_presum < rhs._presum;
		}

		Weighted_random_element& operator=(Weighted_random_element &x) {
			this->_rand = x._rand;
			this->_presum = x._presum;
			return *this;
		}
};
};
};
#endif //_WEIGHTED_RANDOM_ELEMENT_H_
