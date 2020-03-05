#include "../include/polygon_vertex.h"
#include "../include/universe.h"
#include "../include/support_plane.h"
#include "../include/parameters.h"


namespace Skippy {

	Polygon_Vertex_S::Polygon_Vertex_S(const int id_plane, const Constraint & C_1, const Constraint & C_2)
		: Polygon_Vertex(id_plane)
	{
		type = STILL_VERTEX;

		constraints.push_back(C_1);
		constraints.push_back(C_2);

		Support_Plane* SP = Universe::map_of_planes[id_plane];

		//SP->set_landmark(C_1.first, C_2.first, M);
		if (Universe::params->use_landmarks) {
			SP->get_landmark(C_1.first, C_2.first);
		}
	}



	Polygon_Vertex_S::~Polygon_Vertex_S()
	{
	}

	

	void Polygon_Vertex_S::set_polygon(Polygon* _polygon)
	{
		polygon = _polygon;
	}


	CGAL_Point_2 Polygon_Vertex_S::get_M() const
	{
		Constraint C_1 = constraints.front();
		Constraint C_2 = constraints.back();

		Support_Plane* SP = Universe::map_of_planes[id_plane];
		if (Universe::params->use_landmarks) {
			return SP->get_landmark(C_1.first, C_2.first);
		} else {
			return SP->get_intersection_point(C_1.first, C_2.first);
		}
	}


	/*const CGAL_Point_2 & Polygon_Vertex_S::get_M() const
	{
		Constraint C_1 = constraints.front();
		Constraint C_2 = constraints.back();

		Support_Plane* SP = Universe::map_of_planes[id_plane];

		//return SP->get_landmark(C_1.first, C_2.first);
		return SP->get_intersection_point(C_1.first, C_2.first);
	}*/

	
	CGAL_Point_2 Polygon_Vertex_S::pt(const FT & t) const
	{
		return get_M();
	}


	Constraint Polygon_Vertex_S::get_second_constraint() const
	{
		return constraints.back();
	}


	Constraint Polygon_Vertex_S::get_other_constraint(const Constraint & C) const
	{
		Constraint C_1 = constraints.front(), C_2 = constraints.back();
		return (C == C_1 ? C_2 : C_1);
	}



	bool Polygon_Vertex_S::represents_same_intersection(Polygon_Vertex_S* v1, Polygon_Vertex_S* v2) const
	{
		return (represents_same_intersection(v1) || represents_same_intersection(v2));
	}



	bool Polygon_Vertex_S::represents_same_intersection(Polygon_Vertex_S* v) const
	{
		// We assume that all vertices are stopped. Therefore they have two constraints.
		const Constraint C_1 = get_constraint(), C_2 = get_second_constraint();
		const Constraint C_v1 = v->get_constraint(), C_v2 = v->get_second_constraint();

		// We get the indices of the lines that correspond to such constraints.
		const Intersection_Line *I_1 = C_1.first, *I_2 = C_2.first;
		const Intersection_Line *I_v1 = C_v1.first, *I_v2 = C_v2.first;

		// If (C_1, C_2) and (C_v1, C_v2) represent the same lines, then this vertex and v represent the same intersection
		if ((I_1 == I_v1 && I_2 == I_v2) || (I_1 == I_v2 && I_2 == I_v1)) return true;

		// If we couldn't simultaneously match (C_1, C_2) and (C_v1, C_v2), 
		// then we get triplets of concurrent lines in which we find (C_1 and C_2)
		// If C_v1 or C_v2 is part of such a triplet, then we return true.
		std::list<Intersection_Line *> I_L;
		Universe::map_of_planes[v->id_plane]->get_concurrent_lines(I_1, I_2, I_L);

		for (Intersection_Line* L : I_L) {
			if (I_v1 == L || I_v2 == L) return true;
		}

		// The search failed
		return false;
	}


	Polygon_Vertex* Polygon_Vertex_S::to_base()
	{
		return dynamic_cast<Polygon_Vertex*>(this);
	}
}