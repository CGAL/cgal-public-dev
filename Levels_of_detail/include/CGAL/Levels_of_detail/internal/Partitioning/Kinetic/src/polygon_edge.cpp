#include "../include/polygon_edge.h"
#include "../include/polygon_vertex.h"
#include "../include/intersection_line.h"
#include "../include/segment.h"
#include "../include/support_plane.h"
#include "../include/universe.h"


namespace Skippy {

	Polygon_Edge::Polygon_Edge(const int _id_plane, Polygon_Vertex* _v1, Polygon_Vertex* _v2)
		: Support_Plane_Object(_id_plane),
		v1(_v1),
		v2(_v2)
	{
		v1->add(this);
		v2->add(this);

		constraint = Constraint();
	}



	Polygon_Edge::Polygon_Edge(const int _id_plane, Polygon_Vertex* _v1, Polygon_Vertex* _v2, const Constraint & C)
		: Polygon_Edge(_id_plane, _v1, _v2)
	{
		constraint = C;
	}



	Polygon_Edge::Polygon_Edge(const int _id_plane, Polygon_Vertex* _v1, Polygon_Vertex* _v2, Intersection_Line* I, Sign epsilon)
		: Polygon_Edge(_id_plane, _v1, _v2)
	{
		constraint = Constraint(I, epsilon);
	}



	Polygon_Edge::~Polygon_Edge()
	{
		v1->remove(this);
		v2->remove(this);
	}



	Polygon_Vertex* Polygon_Edge::other_vertex(const Polygon_Vertex* v) const
	{
		return (v1 != v ? v1 : v2);
	}


	bool Polygon_Edge::is_constrained_by(Intersection_Line* I) const
	{
		return (constraint.first != nullptr && constraint.first == I);
	}



	Polygon_Vertex* Polygon_Edge::intersection(Intersection_Line* I, Sign s, const FT & t, const int K) const
	{
		CGAL_Point_2 M;
		CGAL_Vector_2 dM;
		intersection_pt_dir(I, t, M, dM);

		return intersection(I, s, t, M, dM, K);
	}



	Polygon_Vertex* Polygon_Edge::intersection(Intersection_Line* I, Sign s, const FT & t, const CGAL_Point_2 & M, const CGAL_Vector_2 & dM, const int K) const
	{
		Polygon_Vertex* v = nullptr;
		Constraint C = Constraint(I, s);
		if (!is_constrained()) {
			v = new Polygon_Vertex_R(I->id_plane, t, M, dM, C, nullptr, nullptr, K, NO_SCHEDULE);
		} else {
			v = new Polygon_Vertex_S(I->id_plane, C, constraint);
		}

		return v;
	}



	void Polygon_Edge::intersection_pt_dir(Intersection_Line* I, const FT & t, CGAL_Point_2 & M, CGAL_Vector_2 & dM) const
	{
		// We assume that edge propagates within the frame of a homeomorphic transformation,
		// that's why the speed of a point that represents the intersection of an edge e and a line I is constant.
		// We suppose that, at the time when this function is called, the edge is colliding I.

		CGAL_Point_2 v1_t = v1->pt(t), v2_t = v2->pt(t);
		
		CGAL_Vector_2 v1_dM = (v1->type == RUNNING_VERTEX ? v1->to_r()->dM : CGAL_Vector_2(0, 0));
		CGAL_Vector_2 v2_dM = (v2->type == RUNNING_VERTEX ? v2->to_r()->dM : CGAL_Vector_2(0, 0));
		//assert(v1->type == RUNNING_VERTEX || v2->type == RUNNING_VERTEX);
		//assert(v1_dM.squared_length() > 0 || v2_dM.squared_length() > 0);

		const CGAL_Line_2 & L = I->line;
		CGAL_Segment_2 S_t(v1_t, v2_t);

		// Step 1.
		// Computes the intersection of L and S_t

		// For better efficiency, we do not compute a intersection of a line L and a segment S_t,
		// but the intersection of two lines : L and S_t.supporting_line()

		CGAL_Point_2 M_t;
		bool M_t_exists = false;

		CGAL_Line_2 S_t_line = S_t.supporting_line();

		CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object_1 = CGAL::intersection(L, S_t_line);
		if (object_1) {
			if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object_1)) {
				M_t = (*ptr);
				M_t_exists = true;
			}
		}

		if (!M_t_exists) {
			std::cout << "L = " << L << ", S_t_line = " << S_t_line << std::endl;
			throw std::logic_error("Error : unexpected behavior of CGAL::intersection.");
		}

		// Step 2.
		// Computes the support line of the edge at time t + 1, and its intersection with L.
		// This lets us compute the speed dM of the seeked intersection point, and finally the initial point M.

		CGAL_Point_2 v1_u = v1_t, v2_u = v2_t;
		if (v1->type == RUNNING_VERTEX) v1_u = v1_t + v1->to_r()->dM;
		if (v2->type == RUNNING_VERTEX) v2_u = v2_t + v2->to_r()->dM;

		bool half_factor = false;
		if (v1_u == v2_u) {
			v1_u = v1_t + v1->to_r()->dM / 2;
			v2_u = v2_t + v2->to_r()->dM / 2;
			half_factor = true;
		}

		CGAL_Line_2 S_u(v1_u, v2_u);

		CGAL_Point_2 M_u;
		bool M_u_exists = false;

		CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object_2 = CGAL::intersection(L, S_u);
		if (object_2) {
			if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object_2)) {
				M_u = (*ptr);
				M_u_exists = true;
			}
		}

		if (!M_u_exists) {
			std::cout.precision(15);
			std::cout << "v1_t = " << v1_t << ", v2_t = " << v2_t << std::endl;
			std::cout << "v1_dm = " << v1->to_r()->dM << ", v2_dm = " << v2->to_r()->dM << std::endl;
			std::cout << "v1_u = " << v1_u << ", v2_u = " << v2_u << std::endl;
			std::cout << "v1_type = " << v1->type << ", v2_type = " << v2->type << std::endl;

			std::cout << "L = " << L << ", S_u = " << S_u << std::endl;
			throw std::logic_error("Error : unexpected behavior of CGAL::intersection.");
		}

		dM = (half_factor ? 2 * (M_u - M_t) : M_u - M_t);
		M = M_t - t * dM;
	}



	void Polygon_Edge::intersection_pt_dir(Intersection_Line* I, Polygon_Vertex* v1, Polygon_Vertex* v2, const FT & t, const CGAL_Point_2 & v1_t, CGAL_Point_2 & M, CGAL_Vector_2 & dM)
	{
		const CGAL_Point_2 v2_t = v2->pt(t);
		const CGAL_Line_2 & L = I->line;

		// The first part of the previous function is now useless : M_t = v1_t.

		// We compute the support line of the edge (v1 v2) at time t + 1, and its intersection with L.
		// This lets us compute the speed dM of the seeked intersection point, and finally the initial point M.

		CGAL_Point_2 v1_u = v1_t, v2_u = v2_t;
		if (v1->type == RUNNING_VERTEX) v1_u = v1_t + v1->to_r()->dM;
		if (v2->type == RUNNING_VERTEX) v2_u = v2_t + v2->to_r()->dM;

		bool half_factor = false;
		if (v1_u == v2_u) {
			v1_u = v1_t + v1->to_r()->dM / 2;
			v2_u = v2_t + v2->to_r()->dM / 2;
			half_factor = true;
		}

		CGAL_Line_2 S_u(v1_u, v2_u);

		CGAL_Point_2 M_u;
		bool M_u_exists = false;

		CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object_2 = CGAL::intersection(L, S_u);
		if (object_2) {
			if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object_2)) {
				M_u = (*ptr);
				M_u_exists = true;
			}
		}

		//CGAL_Vector_2 v1_dM = (v1->type == RUNNING_VERTEX ? v1->to_r()->dM : CGAL_Vector_2(0, 0));
		//CGAL_Vector_2 v2_dM = (v2->type == RUNNING_VERTEX ? v2->to_r()->dM : CGAL_Vector_2(0, 0));
		//assert(v1->type == RUNNING_VERTEX || v2->type == RUNNING_VERTEX);
		//assert(v1_dM.squared_length() > 0 || v2_dM.squared_length() > 0);

		if (!M_u_exists) {
			throw std::logic_error("Error : unexpected behavior of CGAL::intersection.");
		}

		dM = (half_factor ? 2 * (M_u - v1_t) : M_u - v1_t);
		M = v1_t - t * dM;
	}



	/*bool Polygon_Edge::is_adjacent_to_segment() const
	{
		// We assume that vertices v1 and v2 are stopped
		// assert(!v1->is_active && !v2->is_active);
		Intersection_Line* I = constraint.first;
		assert(I != nullptr);

		// Gets a definition of the edge at t = + \infty, and the coordinates of its midpoint
		const CGAL_Point_2 & A = v1->get_M(), &B = v2->get_M();
		const CGAL_Point_2 M = CGAL::midpoint(A, B);

		for (std::list<Segment *>::iterator it_s = I->segments.begin(); it_s != I->segments.end(); it_s++) {
			if (Planar_Segment* pl_s = dynamic_cast<Planar_Segment*>(*it_s)) {
				// The segment is adjacent to an edge of the Support_Plane
				return true;
			}

			else if (Polygon_Segment* po_s = dynamic_cast<Polygon_Segment*>(*it_s)) {
				// We check if the point M computed before belongs to the segment
				const CGAL_Point_2 &S = po_s->origin(), &T = po_s->end();
				const CGAL_Segment_2 ST(S, T);
				//if (ST.has_on(M)) return true;

				FT x_s = S.x(), x_t = T.x(), x = M.x();
				if (x_s != x_t) {
					if ((x_s <= x && x <= x_t) || (x_t <= x && x <= x_s)) return true;
				} else {
					FT y_s = S.y(), y_t = T.y(), y = M.y();
					if ((y_s <= y && y <= y_t) || (y_t <= y && y <= y_s)) return true;
				}
			}
		}

		return false;
	}*/


	bool Polygon_Edge::is_constrained() const
	{
		assert(v1 != nullptr && v2 != nullptr);
		return constraint.first != nullptr;
	}


	Constraint Polygon_Edge::get_constraint() const
	{
		assert(v1 != nullptr && v2 != nullptr);
		return constraint;
	}


	bool Polygon_Edge::two_edges_intersect_two_lines(const std::list<Event_Vertex_Line*> & E_ref, const std::list<Polygon_Vertex_R*> & V_ref,
		std::list<Event_Vertex_Line*> & E_1, std::list<Event_Vertex_Line*> & E, std::list<Event_Vertex_Line*> & E_2,
		Polygon_Vertex_R* & v1, Polygon_Vertex_R* & v, Polygon_Vertex_R* & v2)
	{
		// In this function, provided a list of simultaneous events E_ref, and a list of adjacent vertices V_ref
		// we first identify three consecutive vertices (v1 v v2) such that e1 = (v1 v) and e2 = (v v2) intersect two lines.
	
		int pid = E_ref.front()->plane;
		Support_Plane* SP = Universe::map_of_planes[pid];

		std::list<Polygon_Vertex_R*>::const_iterator it_v = V_ref.begin();
		Polygon_Vertex_R* a = (*it_v); ++it_v;
		Polygon_Vertex_R* b = (*it_v); ++it_v;
		Polygon_Vertex_R* c = (*it_v); ++it_v;

		v1 = v = v2 = nullptr;
		if (a->has_neighbor(b) && a->has_neighbor(c)) {
			v = a; v1 = b; v2 = c;
		} else if (b->has_neighbor(a) && b->has_neighbor(c)) {
			v = b; v1 = a; v2 = c;
		} else if (c->has_neighbor(a) && c->has_neighbor(b)) {
			v = c; v1 = a; v2 = b;
		}

		if (v == nullptr) {
			// If a sequence of vertices couldn't be identified, then returns.
			return false;
		}

		// Now, we loop on the list of events E_ref
		// and decompose it into sub-lists depending on the vertex to which the events apply
		int id_1 = v1->id_object, id = v->id_object, id_2 = v2->id_object;
		for (std::list<Event_Vertex_Line*>::const_iterator it_e = E_ref.begin() ; it_e != E_ref.end() ; ++it_e) {
			int j = (*it_e)->intersectant;
			if (j == id_1) {
				E_1.push_back(*it_e);
			} else if (j == id) {
				E.push_back(*it_e);
			} else if (j == id_2) {
				E_2.push_back(*it_e);
			}
		}

		return true;
	}
}