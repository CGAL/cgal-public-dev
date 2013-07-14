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

		bool operator< (Weighted_random_element &rhs) {
			return this->_weight < rhs._weight;
		}
};
};
};
#endif //_WEIGHTED_RANDOM_ELEMENT_H_
