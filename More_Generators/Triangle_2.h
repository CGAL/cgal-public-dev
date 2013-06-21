#ifndef RANDOM_POINTS_IN_TRIANGLE_2_H
#define RANDOM_POINTS_IN_TRIANGLE_2_H

#include <CGAL/generators.h>
#include <iterator>

namespace CGAL {
template < class P, class Creator = 
           Creator_uniform_2<typename Kernel_traits<P>::Kernel::RT,P> >
class Random_points_in_triangle_2 : public Random_generator_base<P> {
	P _p,_q,_r;
	void generate_point();
public:
	typedef Random_points_in_triangle_2<P> This;
	Random_points_in_triangle_2() {}
	Random_points_in_triangle_2( const P& p, const P& q, const P& r, Random& rnd = default_random)
		: Random_generator_base<P>( 1, rnd ),_p(p),_q(q),_r(r) {
    	generate_point();
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
	

template<class P, class Creator >
void Random_points_in_triangle_2<P, Creator>::generate_point() {
	typedef typename Creator::argument_type T;
	Creator creator;
	double a1 = this->_rnd.get_double(0,1);
	double a2 = this->_rnd.get_double(0,1);
	if(a1>a2) std::swap(a1,a2);
	double b1 = a1;
	double b2 = a2-a1;
	double b3 = 1.0-a2;
	this->d_item = creator(T(to_double(_p.x())*b1+to_double(_q.x())*b2+to_double(_r.x())*b3),
			       T(to_double(_p.y())*b1+to_double(_q.y())*b2+to_double(_r.y())*b3));
}

} //namespace CGAL
#endif //RANDOM_POINTS_IN_TRIANGLE_2_H
