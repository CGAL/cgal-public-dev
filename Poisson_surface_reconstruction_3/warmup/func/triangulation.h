#include <CGAL/basic.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include <CGAL/Vector_3.h>

template < typename Kernel,
typename Vb = CGAL::Triangulation_vertex_base_3<Kernel> >
class VB : public Vb{
public:
	typedef typename Kernel::FT FT;
	typedef typename Vb::Cell_handle Cell_handle;
	typedef typename Vb::Point Point;

	double m_value;

public:
	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
		typedef VB<Kernel, Vb2> Other;
	};

	VB()
		: Vb(), m_value(2.0) {
	}

	VB(const Point & p)
		: Vb(p), m_value(2.0) {
	}

	VB(const Point & p, Cell_handle c)
		: Vb(p,c), m_value(2.0) {
	}

	VB(Cell_handle c)
		: Vb(c), m_value(2.0) {
	}

	const double value() const {
		return m_value;
	}
	double& value() {
		return m_value;
	}

};

template <typename K, typename TDS>
class Min_triangulation_3D : public CGAL::Delaunay_triangulation_3<K, TDS> {
public:
	typedef typename K::FT FT;
	typedef Min_triangulation_3D < K, TDS > MT3;
	typedef typename MT3::Cell_handle Cell_handle;
	typedef typename MT3::Point Point;
	typedef typename K::Vector_3 Vector_3;

public:
	Min_triangulation_3D() {}
	~Min_triangulation_3D(){}

	void find_barycentric_coords(const Point& query, 
		Cell_handle ch, 
		FT& a, 
		FT& b, 
		FT& c, 
		FT& d){
		Point coords[4];
		for(int i = 0; i < 4; i++){
			coords[i] = ch->vertex(i)->point();
		}
		Vector_3 ap(coords[0], query);
		Vector_3 bp(coords[1], query);
		Vector_3 ab(coords[0], coords[1]);
		Vector_3 ac(coords[0], coords[2]);
		Vector_3 ad(coords[0], coords[3]);
		Vector_3 bc(coords[1], coords[2]);
		Vector_3 bd(coords[1], coords[3]);
		const double v = bp * CGAL::cross_product(bd, bc);
		const double w = ap * CGAL::cross_product(ac, ad);
		const double x = ap * CGAL::cross_product(ad, ab);
		const double y = ap * CGAL::cross_product(ab, ac);
		const double z = ab * CGAL::cross_product(ac, ad);

		// TODO: check case where z = 0.0
		a = v/z; b = w/z; c = x/z; d = y/z;
	}

	double compute_func_value(Point query){
		Cell_handle ch = this->locate(query);

		if(this->is_infinite(ch)) {
			return 2.0;
		}

		FT a,b,c,d;
		find_barycentric_coords(query,ch,a,b,c,d);
		return a * ch->vertex(0)->value() +
					 b * ch->vertex(1)->value() +
					 c * ch->vertex(2)->value() +
					 d * ch->vertex(3)->value();
	}

	double compute_func_value_BB(Point query){
		return 0;
	}
};
