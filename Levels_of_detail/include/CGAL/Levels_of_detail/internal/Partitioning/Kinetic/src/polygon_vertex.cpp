#include "../include/polygon_vertex.h"
#include "../include/polygon_edge.h"
#include "../include/intersection_line.h"
#include "../include/segment.h"
#include "../include/support_plane.h"
#include "../include/universe.h"
#include "../include/parameters.h"
#include "../include/event_queue.h"
#include "../include/stats.h"
// #include "sqrt_ft.h"

namespace Skippy 
{

	using CGAL::to_double;


	Polygon_Vertex::Polygon_Vertex(const int _id_plane)
		: Support_Plane_Object(_id_plane)
	{
		e1 = e2 = nullptr;
		constraints = std::list<Constraint>();
		polygon = nullptr;
	}



	Polygon_Vertex::~Polygon_Vertex()
	{
	
	}



	void Polygon_Vertex::add(Polygon_Edge* e)
	{
		if (e1 == nullptr) {
			e1 = e;
		} else {
			assert(e2 == nullptr);
			e2 = e;
		}
	}



	void Polygon_Vertex::remove(Polygon_Edge* e)
	{
		if (e1 == e) {
			e1 = nullptr;
		} else {
			assert(e2 != nullptr);
			e2 = nullptr;
		}
	}


	Polygon_Edge* Polygon_Vertex::other_edge(Polygon_Edge* e) const
	{
		assert(e1 == e || e2 == e);
		if (e1 == e) {
			return e2;
		} else {
			return e1;
		}
	}


	bool Polygon_Vertex::has_neighbor(Polygon_Vertex* v) const
	{
		return (e1 != nullptr && e1->other_vertex(this) == v)
			|| (e2 != nullptr && e2->other_vertex(this) == v);
	}


	bool Polygon_Vertex::unconstrained() const
	{
		return constraints.empty();
	}


	bool Polygon_Vertex::constrained() const
	{
		return !constraints.empty();
	}


	bool Polygon_Vertex::is_constrained_by(Intersection_Line* I) const
	{
		for (std::list<Constraint>::const_iterator it_l = constraints.begin(); it_l != constraints.end(); it_l++) {
			if (it_l->first == I) return true;
		}
		return false;
	}



	Constraint Polygon_Vertex::get_constraint() const
	{
		assert(!constraints.empty());
		return constraints.front();
	}




	Sign Polygon_Vertex::sign_of_constraint(Intersection_Line* I) const
	{
		assert(!constraints.empty());

		std::list<Constraint>::const_iterator it_l = constraints.begin();
		do {
			if (it_l->first == I) break;
		} while (++it_l != constraints.end());
		assert(it_l != constraints.end());

		return it_l->second;
	}



	bool Polygon_Vertex::is_constrained_neighbor(const std::list<Intersection_Line*> & I) const
	{
		for (std::list<Intersection_Line*>::const_iterator it_l = I.begin(); it_l != I.end(); it_l++) {
			if (is_constrained_by(*it_l)) {
				return true;
			}
		}

		return false;
	}



	Polygon* Polygon_Vertex::get_polygon()
	{
		return polygon;
	}


	Polygon_Vertex_R* Polygon_Vertex::to_r() 
	{
		return dynamic_cast<Polygon_Vertex_R*>(this);
	}


	Polygon_Vertex_S* Polygon_Vertex::to_s() 
	{
		return dynamic_cast<Polygon_Vertex_S*>(this);
	}
}