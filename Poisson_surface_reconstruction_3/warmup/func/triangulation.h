#include <CGAL/basic.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include <CGAL/Vector_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Triangle_3.h>
#include "smooth.h"
#include <CGAL/Triangulation_cell_base_3.h>
#include <iterator>
#include <fstream>
#include "random.h"

template < typename Kernel,
typename Vb = CGAL::Triangulation_vertex_base_3<Kernel> >
class VB : public Vb{
public:
	typedef typename Kernel::FT FT;
	typedef typename Kernel::Vector_3 Vector;
	typedef typename Vb::Cell_handle Cell_handle;
	typedef typename Vb::Point Point;

	double m_value;
	Vector m_grad;
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

	const Vector df() const {
		return m_grad;
	}
	Vector& df() {
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
	typedef typename Kernel::Vector_3 Vector;
	typedef typename Kernel::Triangle_3 Triangle;
	typedef typename Kernel::Tetrahedron_3 Tetrahedron;
	Vector m_grad;

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

 const Vector df() const {
   	 return m_grad;
 	}

 	Vector& df() {
 		return m_grad;
 	}

	FT compute_volume() const{
		const Point& pa = this->vertex(0)->point();
		const Point& pb = this->vertex(1)->point();
		const Point& pc = this->vertex(2)->point();
		const Point& pd = this->vertex(3)->point();
		Tetrahedron tet(pa, pb, pc, pd);
		return std::fabs(tet.volume()); // abs of signed volume
	}

	Vector unnormalized_ingoing_normal(const int index)
	{
		const Point& p1 = this->vertex((index+1)%4)->point();
		const Point& p2 = this->vertex((index+2)%4)->point();
		const Point& p3 = this->vertex((index+3)%4)->point();
		Vector cross = CGAL::cross_product(p2 - p1, p3 - p1);

		// Triangle t(p1, p2, p3);
		// if(4.0 * t.squared_area() == cross * cross) std:: cout << "good area" << std::endl; //sanity check
		// else std::cout << "wrong area" << std::endl;
		// std:: cout << t.squared_area() << " " << cross*cross << std::endl;
		if(index%2 == 0)
			return -0.5 * cross; // area of triangle is 0.5 the cross product
		else
			return 0.5 * cross;
	}

	Vector compute_grad(const int ref)
	{
			FT fref = this->vertex(ref)->f();

			m_grad = CGAL::NULL_VECTOR;
			double volume = this->compute_volume();
			for(int i = 1; i <= 3; i++)
			{ //face opposite each of i
				const int other = (ref + i) % 4;
				FT fother = this->vertex(other)->f();
				const Vector normal = this->unnormalized_ingoing_normal(other) / (3.0 * volume);
				m_grad = m_grad + (fother - fref) * normal;
			}
			return m_grad;
	}

};

template <typename K, typename TDS>
class Min_triangulation_3D : public CGAL::Delaunay_triangulation_3<K, TDS> {
public:
	typedef typename K::FT FT;
	typedef Min_triangulation_3D < K, TDS > MT3;
	typedef typename MT3::Cell_handle Cell_handle;
	typedef typename MT3::Point Point;
	typedef typename K::Vector_3 Vector;
	typedef typename K::Tetrahedron_3 Tetrahedron;
	typedef typename K::Triangle_3 Triangle;
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
		const Point p = Point(x, y, z); //  +::random_vec<Vector>(1e-6);;
	    Vertex_handle v = this->insert(p);
	    v->f() = f;
	  }

	  ifile.close();
	}
	
	void read_xyzn(const char* filename){
	  std::ifstream ifile(filename);
	  int num_vertices;
	  FT x, y, z, dx, dy, dz, f;
	  //ifile >> num_vertices;
	  //std::cout << num_vertices << std::endl;

	  while(!ifile.eof()){
	    ifile >> x >> y >> z >> dx >> dy >> dz >> f;
	    std::cout << x << " " << y << " " <<  z << " " << dx << " " << dy << " " <<  dz << " " << f << std::endl;

	    // insert to triangulation
		const Point p = Point(x, y, z); //  +::random_vec<Vector>(1e-6);
		Vector grad(dx, dy, dz);
	    Vertex_handle v = this->insert(p);
	    v->f() = f;
	    v->df() = grad;
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
		Vector ap(coords[0], query);
		Vector bp(coords[1], query);
		Vector ab(coords[0], coords[1]);
		Vector ac(coords[0], coords[2]);
		Vector ad(coords[0], coords[3]);
		Vector bc(coords[1], coords[2]);
		Vector bd(coords[1], coords[3]);
		const double v = bp * CGAL::cross_product(bd, bc);
		const double w = ap * CGAL::cross_product(ac, ad);
		const double x = ap * CGAL::cross_product(ad, ab);
		const double y = ap * CGAL::cross_product(ab, ac);
		const double z = ab * CGAL::cross_product(ac, ad);

		if (v < 0.0) std::cout << "negative v" << std::endl;
		if (w < 0.0) std::cout << "negative w" << std::endl;
		if (x < 0.0) std::cout << "negative x" << std::endl;
		if (y < 0.0) std::cout << "negative y" << std::endl;
		if (z < 0.0) std::cout << "negative z" << std::endl;

		// TODO: check case where z = 0.0
		a = v/z; b = w/z; c = x/z; d = y/z;

		//const double sum = a + b + c + d;
		//std::cout << "sum: " << sum << std::endl;
	}

	void compute_grad(Vertex_handle v)
	{
		// get incident cells
		std::vector<Cell_handle> cells;
		this->incident_cells(v, std::back_inserter(cells));

		FT sum_volumes = 0.0;
		Vector sum_vec = CGAL::NULL_VECTOR;

		typename std::vector<Cell_handle>::iterator it;
		for(it = cells.begin(); it != cells.end(); it++)
		{
			Cell_handle c = *it;
			if(this->is_infinite(c))
				continue;

			// use cell (get gradient df and compute cell volume)
			const FT volume = c->compute_volume();
			sum_vec = sum_vec + volume * c->df();
			sum_volumes += volume;
		}

		if (sum_volumes != 0.0)
			v->df() = sum_vec / sum_volumes;
		else
			v->df() = CGAL::NULL_VECTOR;

		// DEBUG HARDCODED
		const Point& p = v->point();
		Vector vec = p - CGAL::ORIGIN;
		vec = vec / std::sqrt(vec*vec); // normalize
		v->df() = 5.0 * vec; // scale
	}

	void compute_grad_per_cell(){
		int i = 0;
		for(auto it = this->finite_cells_begin(); 
			it != this->finite_cells_end(); 
			it++, i++)
		{
			// std::cout << "CELL " << i << ":" << std::endl;
			//std::cout << it->vertex(0)->point() << " " << it->vertex(1)->point() << " " << it->vertex(2)->point() << " "<< it->vertex(3)->point() << " " << std::endl;

			it->compute_grad(0);

			 // DEBUGGING
			/*

			Vector vec0 = it->compute_grad(0);
			Vector vec1 = it->compute_grad(1);
			Vector vec2 = it->compute_grad(2);
			Vector vec3 = it->compute_grad(3);

			std::cout << "GRADIENT" << std::endl;
			std::cout << vec0 << std::endl;
			std::cout << vec1 << std::endl;
			std::cout << vec2 << std::endl;
			std::cout << vec3 << std::endl;
			*/
		}
	}

	void compute_grad_per_vertex(){
		for(auto it = this->finite_vertices_begin(); it != this->finite_vertices_end(); it++){
			compute_grad(it);
		}
	}

	double compute_func_value(Point query){
		Cell_handle ch = this->locate(query);

		if(this->is_infinite(ch)) {
			return 2.0;
		}

		FT a, b, c, d;
		find_barycentric_coords(query, ch, a, b, c, d);
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

		double f[4];
		double x[3*4];
		double gradf[3*4];

		for (int i = 0; i < 4; i++)
		{
			Vertex_handle v = ch->vertex(i);

			f[i] = v->f();

			const Point& p = v->point();
			x[3 * i] = p[0];
			x[3 * i + 1] = p[1];
			x[3 * i + 2] = p[2];

			Vector df = v->df(); // gradient per vertex
			gradf[3 * i] = df[0];
			gradf[3 * i + 1] = df[1];
			gradf[3 * i + 2] = df[2];
		}


 //barycentric coordinates
		double b[20];
	  control_points(b, x, f, gradf);

		double w[4]; //barycentric coordinates
		find_barycentric_coords(query, ch, w[0], w[1], w[2], w[3]);
		return eval_bernstein3(b, w);
	}

	void output_grads_to_off(){
		std::ofstream ofile("grad.off");
		std::cout << "Number of vertices (inside writing to off)" << this->number_of_vertices() << std::endl;
		ofile << "OFF" << std:: endl
			           << 7*(this->number_of_vertices()) << " "
			           << 2*(this->number_of_vertices()) << " 0" << std::endl;
		int i = 0;
		for(auto it = this->finite_vertices_begin(); it != this->finite_vertices_end(); it++, i++){
			//if (i == 0) continue;
			ofile << it->point()[0] << " " << it->point()[1] << " " << it->point()[2] << std::endl << it->point()[0] << " " << it->point()[1] << " " << it->point()[2] << std::endl;
			Vector grad = it->df();
			std::cout << "Point: " << it->point() <<  " Grad " << i << ": " << grad << std::endl;
			ofile << it->point()[0] + grad[0] << " " << it->point()[1] + grad[1] << " " << it->point()[2] + grad[2]<< std::endl << it->point()[0] + grad[0] << " " << it->point()[1] + grad[1] << " " << it->point()[2] + grad[2]<< std::endl;
			ofile << it->point()[0] + grad[0] + 0.01 << " " << it->point()[1] + grad[1] << " " << it->point()[2] + grad[2]<< std::endl;
			ofile << it->point()[0] + grad[0] - 0.01 << " " << it->point()[1] + grad[1] << " " << it->point()[2] + grad[2]<< std::endl;
			ofile << it->point()[0] + grad[0] << " " << it->point()[1] + grad[1] + 0.02 << " " << it->point()[2] + grad[2]<< std::endl;
		}
		i = 0;
		for(auto it = this->finite_vertices_begin(); it != this->finite_vertices_end(); it++, i++){
			ofile << "4 " << 7*i << " " << 7*i + 1 << " " << 7*i + 2 << " " << 7*i + 3 << std::endl;
			ofile << "3 " << 7*i + 4 << " " << 7*i + 5 << " " << 7*i + 6 << std::endl;
		}
		ofile.close();
	}
};
