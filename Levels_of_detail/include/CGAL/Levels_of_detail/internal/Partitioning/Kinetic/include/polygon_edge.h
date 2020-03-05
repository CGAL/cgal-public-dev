#pragma once
#include "support_plane_object.h"
#include "support_plane_objects.h"
#include "defs_cgal.h"


namespace Skippy
{
	class Polygon_Edge : public Support_Plane_Object 
	{
	public:
		Polygon_Edge(const int _id_plane, Polygon_Vertex* _v1, Polygon_Vertex* _v2);

		Polygon_Edge(const int _id_plane, Polygon_Vertex* _v1, Polygon_Vertex* _v2, const Constraint & C);

		Polygon_Edge(const int _id_plane, Polygon_Vertex* _v1, Polygon_Vertex* _v2, Intersection_Line* I, Sign epsilon);

		~Polygon_Edge();

		Polygon_Vertex* other_vertex(const Polygon_Vertex* v) const;

		bool is_constrained_by(Intersection_Line* I) const;

		// Different versions of intersection_pt_dir :
		// - we don't know if v1 or v2 already intersects I.
		void intersection_pt_dir(Intersection_Line* I, const FT & t, CGAL_Point_2 & M, CGAL_Vector_2 & dM) const;

		// - we know that v1 intersects I in V1_t.
		static void intersection_pt_dir(Intersection_Line* I, Polygon_Vertex* v1, Polygon_Vertex* v2, const FT & t, const CGAL_Point_2 & V1_t, CGAL_Point_2 & M, CGAL_Vector_2 & dM);

		Polygon_Vertex* intersection(Intersection_Line* I, Sign s, const FT & t, const int K) const;

		Polygon_Vertex* intersection(Intersection_Line* I, Sign s, const FT & t, const CGAL_Point_2 & M, const CGAL_Vector_2 & dM, const int K) const;

		// bool is_adjacent_to_segment() const;

		bool is_constrained() const;

		Constraint get_constraint() const;

		static bool two_edges_intersect_two_lines(const std::list<Event_Vertex_Line*> & E_ref, const std::list<Polygon_Vertex_R*> & V_ref,
			std::list<Event_Vertex_Line*> & E_1, std::list<Event_Vertex_Line*> & E, std::list<Event_Vertex_Line*> & E_2,
			Polygon_Vertex_R* & v1, Polygon_Vertex_R* & v, Polygon_Vertex_R* & v2);

	public:
		Polygon_Vertex* v1;
		Polygon_Vertex* v2;

	protected:
		Constraint constraint;
	};
}