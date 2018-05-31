#include <CGAL/basic.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include <CGAL/Vector_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Triangle_3.h>
#include"smooth.h"
#include <CGAL/Triangulation_cell_base_3.h>
#include<iterator>

template < typename Kernel,
typename Vb = CGAL::Triangulation_vertex_base_3<Kernel> >
class VB : public Vb{
public:
	typedef typename Kernel::FT FT;
	typedef typename Kernel::Vector_3 Vector_3;
	typedef typename Vb::Cell_handle Cell_handle;
	typedef typename Vb::Point Point;

	double m_value;
	Vector_3 m_grad;
public:
	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
		typedef VB<Kernel, Vb2> Other;
	};

	VB()
		: Vb(), m_value(2.0), m_grad(CGAL::NULL_VECTOR) {
	}

	VB(const Point & p)
		: Vb(p), m_value(2.0), m_grad(CGAL::NULL_VECTOR) {
	}

	VB(const Point & p, Cell_handle c)
		: Vb(p,c), m_value(2.0), m_grad(CGAL::NULL_VECTOR) {
	}

	VB(Cell_handle c)
		: Vb(c), m_value(2.0), m_grad(CGAL::NULL_VECTOR) {
	}

	const double f() const {
		return m_value;
	}
	double& f() {
		return m_value;
	}

	const Vector_3 df() const {
		return m_grad;
	}
	Vector_3& df() {
		return m_grad;
	}

};

template < typename Kernel,
typename Cb = CGAL::Triangulation_cell_base_3<Kernel> >
class CB : public Cb{
public:
	typedef typename Kernel::FT FT;
	typedef typename Cb::Cell_handle Cell_handle;
	typedef typename Cb::Vertex_handle Vertex_handle;
	typedef typename Cb::Point Point;
	typedef typename Kernel::Vector_3 Vector_3;
	typedef typename Kernel::Triangle_3 Triangle_3;
	typedef typename Kernel::Tetrahedron_3 Tetrahedron_3;
	Vector_3 m_grad;

public:
	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename Cb::template Rebind_TDS<TDS2>::Other Cb2;
		typedef CB<Kernel, Cb2> Other;
	};

	CB()
		: Cb(), m_grad(CGAL::NULL_VECTOR) {
	}

 	CB(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3)
	 : Cb(v0, v1, v2, v3) , m_grad(CGAL::NULL_VECTOR){
	 }

 	CB(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3, Cell_handle n0, Cell_handle n1,
	 	 Cell_handle n2, Cell_handle n3): Cb(v0, v1, v2, v3, n0, n1, n2, n3), m_grad(CGAL::NULL_VECTOR) {
	 }

 const Vector_3 df() const {
   	 return m_grad;
 	}

 	Vector_3& df() {
 		return m_grad;
 	}

	FT compute_volume() const{
		const Point& pa = this->vertex(0)->point();
		const Point& pb = this->vertex(1)->point();
		const Point& pc = this->vertex(2)->point();
		const Point& pd = this->vertex(3)->point();
		Tetrahedron_3 tet(pa, pb, pc, pd);
		const double v = tet.volume();
		return v;
	}

	Vector_3 unnormalized_ingoing_normal(const int index){
		const Point& p1 = this->vertex((index+1)%4)->point();
		const Point& p2 = this->vertex((index+2)%4)->point();
		const Point& p3 = this->vertex((index+3)%4)->point();
		Vector_3 cross = CGAL::cross_product(p2 - p1, p3 - p1);
		if(index%2 == 0)
			return -cross;
		else
			return cross;
	}

	void compute_grad(){
			int index = 0; Vector_3 grad = CGAL::NULL_VECTOR;
			FT findex = this->vertex(index)->f();
			for(int i = 0; i < 3; i++){ //face opposite each of i
				FT fi = this->vertex(i)->f();
				Vector_3 normal = this->unnormalized_ingoing_normal(i);
				grad = grad + (fi - findex ) * normal;
			}
			m_grad = grad;
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
	typedef typename K::Tetrahedron_3 Tetrahedron_3;
	typedef typename K::Triangle_3 Triangle_3;
	typedef typename MT3::Vertex_handle Vertex_handle;
	typedef typename MT3::Facet Facet;

public:
	Min_triangulation_3D() {}
	~Min_triangulation_3D(){}

	void read_xyz(const char* filename){
		std::ifstream ifile(filename);
	  int num_vertices;
	  FT x, y, z, f;
	  ifile >> num_vertices;
	  std::cout << num_vertices << std::endl;

	  for(int i = 0; i < num_vertices; i++){
	    ifile >> x >> y >> z >> f;
	    std::cout << x << " " << y << " " <<  z << " " << f << std::endl;

	    // insert to triangulation
	    Point p(x, y, z);
	    Vertex_handle v = this->insert(p);
	    v->f() = f;
	  }

	  ifile.close();
	}

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

	void compute_grad(Vertex_handle v){
		std::vector<Cell_handle> cells;
		this->incident_cells(v, std::back_inserter(cells));

		Vector_3 vec =CGAL::NULL_VECTOR;
		FT sum_volumes = 0.0;

		typename std::vector<Cell_handle>::iterator it;
		for(it = cells.begin(); it != cells.end(); it++){
			Cell_handle c = *it;
			if(this->is_infinite(c))
				continue;
			// use cell (get gradient df and compute cell volume)
			const FT volume = c->compute_volume();
			const Vector_3 df = c->df();
			vec = vec + df;
			sum_volumes += volume;
		}
		v->df() = vec / sum_volumes;
	}

	void compute_grad_per_cell(){
		for(auto it = this->all_cells_begin(); it != this->all_cells_end(); it++){
			it->compute_grad();
		}
	}

	void compute_grad_per_vertex(){
		for(auto it = this->all_vertices_begin(); it != this->all_vertices_end(); it++){
			compute_grad(it);
		}
	}

	double compute_func_value(Point query){
		Cell_handle ch = this->locate(query);

		if(this->is_infinite(ch)) {
			return 2.0;
		}

		FT a,b,c,d;
		find_barycentric_coords(query,ch,a,b,c,d);
		return a * ch->vertex(0)->f() +
					 b * ch->vertex(1)->f() +
					 c * ch->vertex(2)->f() +
					 d * ch->vertex(3)->f();
	}

	double compute_func_value_BB(Point query){
		Cell_handle ch = this->locate(query);

		if(this->is_infinite(ch)) {
			return 2.0;
		}

		double x[3*4];
		double f[4];
		double gradf[3*4];

		for (int i=0;i<4;++i){
			Point p = ch->vertex(i)->point();
			x[3*i] = p[0];
			x[3*i + 1] = p[1];
			x[3*i + 2] = p[2];
		}

	  for (int i=0;i<4;++i){
			f[i] = ch->vertex(i)->f();
		}

	  for (int i=0;i<4;++i){
			Vector_3 p = ch->vertex(i)->df();
			gradf[3*i] = p[0];
			gradf[3*i + 1] = p[1];
			gradf[3*i + 2] = p[2];
		}

		double b[20];
	  control_points(b,x,f,gradf);

		FT w[4]; //barycentric coordinates
		find_barycentric_coords(query, ch, w[0], w[1], w[2], w[3]);
		return eval_bernstein3(b,w);
	}
};
