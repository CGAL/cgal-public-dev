#include "../include/support_plane_objects.h"
#include "../include/polygon_vertex.h"
#include "../include/polygon_edge.h"
#include "../include/intersection_line.h"
#include "../include/polygon.h"
#include "../include/universe.h"
#include "../include/parameters.h"
#include "../include/polygon_node.h"



namespace Skippy 
{
	Polygon::Polygon(const int _seed, const std::list<Polygon_Vertex*> & _vertices, const std::list<Polygon_Edge*> & _edges)
	{
		running_vertices = 0;

		for (std::list<Polygon_Vertex *>::const_iterator it_v = _vertices.begin(); it_v != _vertices.end(); it_v++) {
			vertices.push_back(*it_v);
			(*it_v)->set_polygon(this);
		}

		for (std::list<Polygon_Edge *>::const_iterator it_e = _edges.begin(); it_e != _edges.end(); it_e++) {
			edges.push_back(*it_e);
		}

		seed = _seed;
		cell = nullptr;
	}



	Polygon::Polygon(const int id_plane, const int _seed, Polygon_Directions* D)
	{
		running_vertices = 0;

		size_t n = D->size();

		// Initializes a list of Polygon_Vertices

		for (size_t i = 0; i < n; i++) {
			const std::pair<CGAL_Point_2, CGAL_Vector_2> & D_i = D->get_vertex(i);
			Polygon_Vertex_R* v_i = new Polygon_Vertex_R(id_plane, 0, D_i.first, D_i.second, Universe::params->K, NO_SCHEDULE);

			vertices.push_back(v_i);
			v_i->set_polygon(this);
		}

		// Initializes the list of Polygon_Edges.

		std::list<Polygon_Vertex *>::iterator it_v1 = vertices.begin(), it_v2 = ++vertices.begin();
		while (it_v1 != vertices.end()) {
			Polygon_Vertex* v1 = (*it_v1);
			Polygon_Vertex* v2 = (*it_v2);

			Polygon_Edge* e = new Polygon_Edge(id_plane, v1, v2);
			edges.push_back(e);

			++it_v1;
			if (++it_v2 == vertices.end()) it_v2 = vertices.begin();
		}

		seed = _seed;
		cell = nullptr;
	}



	Polygon::Polygon(const int id_plane, const int _seed, Polygon_Directions* D, const CGAL_Vector_2 & OM, const FT & t)
	{
		running_vertices = 0;

		// Same algorithm as before.

		size_t n = D->size();

		for (size_t i = 0; i < n; i++) {
			const std::pair<CGAL_Point_2, CGAL_Vector_2> & D_i = D->get_vertex(i);
			Polygon_Vertex_R* v_i = new Polygon_Vertex_R(id_plane, t, D_i.first + OM, D_i.second, 0, NO_SCHEDULE);

			vertices.push_back(v_i);
			v_i->set_polygon(this);
		}

		std::list<Polygon_Vertex *>::iterator it_v1 = vertices.begin(), it_v2 = ++vertices.begin();
		while (it_v1 != vertices.end()) {
			Polygon_Vertex* v1 = (*it_v1);
			Polygon_Vertex* v2 = (*it_v2);

			Polygon_Edge* e = new Polygon_Edge(id_plane, v1, v2);
			edges.push_back(e);

			++it_v1;
			if (++it_v2 == vertices.end()) it_v2 = vertices.begin();
		}

		seed = _seed;
		cell = nullptr;
	}



	Polygon::~Polygon()
	{
		for (std::list<Polygon_Edge *>::iterator it_e = edges.begin(); it_e != edges.end(); it_e++) {
			delete (*it_e);
		}

		for (std::list<Polygon_Vertex *>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
			delete (*it_v);
		}

		edges.clear();
		vertices.clear();
	}



	void Polygon::clear()
	{
		edges.clear();
		vertices.clear();
	}



	bool Polygon::intersection_exists(Intersection_Line* I, const FT & t, std::map<int, Sign> & vertices_signs,
		std::list<Polygon_Vertex *> & vertices_plus, std::list<Polygon_Edge *> & edges_plus, Polygon_Edge* & edge_plus_front, Polygon_Edge* & edge_plus_back,
		std::list<Polygon_Vertex *> & vertices_minus, std::list<Polygon_Edge *> & edges_minus, Polygon_Edge* & edge_minus_front, Polygon_Edge* & edge_minus_back,
		std::list<Polygon_Vertex *> & vertices_zero)
	{
		// Computes the positives of the vertices relatively to the line I
		vertices_signs.clear();
		if (t == 0) {
			for (std::list<Polygon_Vertex *>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
				Polygon_Vertex* v = (*it_v);
				vertices_signs[v->id_object] = I->sign(v->get_M());
			}
		} else {
			for (std::list<Polygon_Vertex *>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
				Polygon_Vertex* v = (*it_v);
				vertices_signs[v->id_object] = I->sign(v->pt(t));
			}
		}

		// Exhibits two sequences, whose vertices lie strictly above or below the line I
		bool processed_plus = false, processed_minus = false;
		for (std::list<Polygon_Vertex *>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
			Polygon_Vertex* v = (*it_v);

			if (vertices_signs[v->id_object] == PLUS && !processed_plus) {
				exhibit_sequence(I, v, PLUS, vertices_signs, vertices_plus, edges_plus, edge_plus_front, edge_plus_back);
				processed_plus = true;

			} else if (vertices_signs[v->id_object] == MINUS && !processed_minus) {
				exhibit_sequence(I, v, MINUS, vertices_signs, vertices_minus, edges_minus, edge_minus_front, edge_minus_back);
				processed_minus = true;
			}

			if (processed_plus && processed_minus) break;
		}

		// Lists vertices which are on the line I
		for (std::list<Polygon_Vertex*>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
			Polygon_Vertex* v = (*it_v);
			if (vertices_signs[v->id_object] == ZERO) {
				vertices_zero.push_back(v);
			}
		}

		return (vertices_plus.size() > 0 && vertices_minus.size() > 0);
	}



	void Polygon::exhibit_sequence(Intersection_Line* I, Polygon_Vertex* v, Sign epsilon, std::map<int, Sign> & signs,
		std::list<Polygon_Vertex *> & sequence_vertices, std::list<Polygon_Edge *> & sequence_edges,
		Polygon_Edge* & excluded_edge_front, Polygon_Edge* & excluded_edge_back)
	{
		// We assume that v strictly falls above or below I
		// Let us initialize the sequence of vertices with it
		sequence_vertices.push_back(v);

		// v is connected to two edges
		// We loop on the vertices that delimit each of them
		Polygon_Edge *e1 = v->e1, *e2 = v->e2;
		Polygon_Vertex *v1 = (e1->v1 == v ? e1->v2 : e1->v1);
		while (signs[v1->id_object] == epsilon && v1 != v) {
			sequence_vertices.push_front(v1);
			sequence_edges.push_front(e1);
			e1 = (v1->e1 == e1 ? v1->e2 : v1->e1);
			v1 = (e1->v1 == v1 ? e1->v2 : e1->v1);
		}

		if (v1 == v) {
			// Case when we are back to the start
			// All the vertices of the polygon are on the same side of I
			// There is no need to loop on the second edge of v
			return;
		}


		Polygon_Vertex *v2 = (e2->v1 == v ? e2->v2 : e2->v1);
		while (signs[v2->id_object] == epsilon) {
			sequence_vertices.push_back(v2);
			sequence_edges.push_back(e2);
			e2 = (v2->e1 == e2 ? v2->e2 : v2->e1);
			v2 = (e2->v1 == v2 ? e2->v2 : e2->v1);
		}

		// Edges delimiting the sequence
		excluded_edge_front = e1;
		excluded_edge_back = e2;
	}



	// Events handlers
	// Case A


	void Polygon::split_unconstrained_vertex(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t,
		Polygon_Vertex_R* & v1, Polygon_Vertex_R* & v2, Constraint & C)
	{
		// Part 1.

		// We are going to replace it by two constrained vertices, v1 and v2, that move towards opposite directions.
		// From a connected sequence (a v b), we want to obtain (a v1 v2 b).
		Polygon_Edge *e_1 = v->e1, *e_2 = v->e2;
		Polygon_Vertex* a = (e_1->v1 == v ? e_1->v2 : e_1->v1);
		Polygon_Vertex* b = (e_2->v1 == v ? e_2->v2 : e_2->v1);

		// Computes the intersection of I and the two edges of v : initial virtual point, and speed vector

		CGAL_Point_2 M_1, M_2;
		CGAL_Vector_2 dM_1, dM_2;

		Polygon_Edge::intersection_pt_dir(I, v, a, t, V_t, M_1, dM_1);
		Polygon_Edge::intersection_pt_dir(I, v, b, t, V_t, M_2, dM_2);
		// e_1->intersection_pt_dir(I, t, M_1, dM_1);
		// e_2->intersection_pt_dir(I, t, M_2, dM_2);

		// We are going to create two vertices that will remain on the same side of the line as v(t - epsilon)
		// Normally, v->pt(0) is clearly one or another side of the line, except when it is on the line at t = 0

		C = Constraint(I, I->sign(v, 0));

		v1 = new Polygon_Vertex_R(v->id_plane, t, M_1, dM_1, C, nullptr, a->to_r(), 0, SCHEDULE);
		v2 = new Polygon_Vertex_R(v->id_plane, t, M_2, dM_2, C, nullptr, b->to_r(), 0, SCHEDULE);

		// Part 2.
		// Intersection of v1 and v2 into the polygon, removal of v.

		// Edges e1 and e2 are now deleted and replaced by new ones

		delete_object(e_1, e_2, edges, true);
		delete_object(v->to_base(), vertices, false);

		insert_object(v1->to_base(), v2->to_base(), vertices);

		Polygon_Edge* a_v1 = new Polygon_Edge(v->id_plane, a, v1);
		Polygon_Edge* v2_b = new Polygon_Edge(v->id_plane, v2, b);
		Polygon_Edge* v1_v2 = new Polygon_Edge(v->id_plane, v1, v2, C);

		insert_object(a_v1, v2_b, v1_v2, edges);

		// Note that v is NOT destroyed
	}



	Polygon* Polygon::build_as_unconstrained_vertex_propagates(Polygon_Vertex_R* v, const FT & t,
		const int polygon_seed, Polygon_Vertex_R* v1_ts, Polygon_Vertex_R* v2_ts, const Constraint & C_ts,
		Polygon_Vertex_R* & v1_os, Polygon_Vertex_R* & v2_os)
	{
		// This function is the sequel of the previous one
		// We assume that v is allowed to propagate on the other side of the line I, 
		// which is implicitly mentioned in the list of parameters, through C_ts.
		// We build two vertices v1_os, v2_os by reflection of v1_ts and v2_ts,
		// and initialize three edges.

		int pid = v->id_plane;
		const Constraint C_os = Constraint(C_ts.first, (C_ts.second == PLUS ? MINUS : PLUS));

		v1_os = new Polygon_Vertex_R(v1_ts, SCHEDULE);
		v2_os = new Polygon_Vertex_R(v2_ts, SCHEDULE);

		// v1_os = new Polygon_Vertex(pid, t, v1_ts->M, v1_ts->dM, C_os, nullptr, SCHEDULE);
		// v2_os = new Polygon_Vertex(pid, t, v2_ts->M, v2_ts->dM, C_os, nullptr, SCHEDULE);

		Polygon_Edge* v_v1 = new Polygon_Edge(pid, v, v1_os);
		Polygon_Edge* v_v2 = new Polygon_Edge(pid, v, v2_os);
		Polygon_Edge* v1_v2 = new Polygon_Edge(pid, v1_os, v2_os, C_os);

		std::list<Polygon_Vertex*> V;
		V.push_back(v); V.push_back(v1_os); V.push_back(v2_os);

		std::list<Polygon_Edge*> E;
		E.push_back(v_v1); E.push_back(v_v2); E.push_back(v1_v2);

		return new Polygon(polygon_seed, V, E);
	}



	// Handlers
	// Case B

	void Polygon::remove_vertices_colliding_on_intersection_line(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t, 
		const CGAL_Point_2 & V_t, Polygon_Vertex_R* & vc_ts, CGAL_Vector_2 & u_ref)
	{
		// Vertex v, which is unconstrained, intersects the line I.
		// It is assumed that one of this neighbor vertices is constrained by I, and therefore intersects v at time t.

		// Illustration :
		// 
		// I ------------------     I ----------------------
		//   vb ----- vc ->            vb ------- vc' ->
		//             \   ^                        \
		//              \ /     =>                   \
		//               v                            \
		//              /                             va
		//             va
		// 
		// We have a sequence of vertices (va v vc vb), and a sequence of edges (ea = va_v, ec = v_vc, eb = vc_vb)
		// which is going to become : (va vc' vb), (ea' = va_vc', eb' = vc'_vb).

		// Part 1.
		// Identification of the sequence

		Polygon_Vertex *v_a, *v_b, *v_c;
		Polygon_Edge *e_a, *e_b, *e_c;

		Polygon_Vertex* v1 = (v->e1->v1 == v ? v->e1->v2 : v->e1->v1);
		Polygon_Vertex* v2 = (v->e2->v1 == v ? v->e2->v2 : v->e2->v1);

		assert(v1->is_constrained_by(I) ^ v2->is_constrained_by(I));
		if (v1->is_constrained_by(I)) {
			v_c = v1; v_a = v2;
			e_c = v->e1; e_a = v->e2;
		} else {
			v_c = v2; v_a = v1;
			e_c = v->e2; e_a = v->e1;
		}
		e_b = (v_c->e1 == e_c ? v_c->e2 : v_c->e1);
		v_b = (e_b->v1 == v_c ? e_b->v2 : e_b->v1);

		// u_ref = v_c->pt(t + 1) - v->pt(t + 1)
		//       = v_c->pt(t) + v_c->dM - (v->pt(t) + v->dM)
		//       = v_c->dM - v->dM (since v_c(t) = v(t))
		u_ref = v_c->to_r()->dM - v->to_r()->dM;

		// Part 2.
		// Modification of the sequence

		// First we should create the new vertex v_c', which results from the intersection of I and e_a
		// and update the events involving the segments led by e_c
		CGAL_Point_2 M; CGAL_Vector_2 dM;
		Polygon_Edge::intersection_pt_dir(I, v, v_a, t, V_t, M, dM);
		// e_a->intersection_pt_dir(I, t, M, dM);

		const Constraint C = v_c->get_constraint();
		vc_ts = new Polygon_Vertex_R(v->id_plane, t, M, dM, C, nullptr, v_a->to_r(), 0, SCHEDULE);
		v_c->to_r()->transfer_segments(vc_ts->to_r(), t, V_t);

		// Deletes and inserts elements

		delete_object(e_a, e_b, e_c, edges, true);

		delete_object(v_c, vertices, true);
		delete_object(v->to_base(), vertices, false);

		insert_object(vc_ts->to_base(), vertices);

		Polygon_Edge* e_a_prime = new Polygon_Edge(v->id_plane, v_a, vc_ts);
		Polygon_Edge* e_b_prime = new Polygon_Edge(v->id_plane, vc_ts, v_b, C);

		insert_object(e_a_prime, e_b_prime, edges);

		// Warning : v is NOT destroyed
		// Since it may be transferred to the polygon on the other side of I
	}



	Polygon_Vertex_R* Polygon::get_tied_vertex(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t)
	{
		// See illustration below.
		// Unconstrained vertex v was previously on "this side" of line I, and we want to find the vertex vc in this polygon.
		// vc is a vertex constrained by I, which is at the same location as v at time t.

		for (std::list<Polygon_Vertex*>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
			if (Polygon_Vertex_R* v_it = (*it_v)->to_r()) {
				if (v_it->is_constrained_by(I) && ((v_it->pt(t) - V_t) == CGAL::NULL_VECTOR)) {
					return v_it;
				}
			}
		}

		return nullptr;
	}



	Polygon_Vertex_R* Polygon::get_tied_vertex(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t, const CGAL_Vector_2 & u)
	{
		// This function performs the same operation as the previous one :
		// it deals with the retrieval of a vertex vc constrained by I, 
		// such that v(t) = vc(t) and vc is located on the other side of I

		// However, here we perform an additional check :
		// we return vc if its non-constrained edge has integer r or vector u for reference vector

		Polygon_Vertex_R* v_c = get_tied_vertex(I, v, t, V_t);
		if (v_c == nullptr) return nullptr;

		Polygon_Edge* e = (v_c->e1->is_constrained_by(I) ? v_c->e2 : v_c->e1);
		CGAL_Vector_2 u_ref = e->v1->pt(t + 1) - e->v2->pt(t + 1);
		if (CGAL::determinant(u_ref, u) == 0) return v_c;

		return nullptr;
	}



	void Polygon::insert_unconstrained_vertex(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t,
		const CGAL_Point_2 & V_t, Polygon_Vertex_R* vc_ts, Polygon_Vertex_R* v_c)
	{
		// Unconstrained vertex v was previously on "this side" of line I.
		// The polygon to which it used to belong now includes a constrained vertex, vc_ts.
		// In this function, we aim at inserting v in polygon "this", which is located on the "other side" of I.

		// Illustration :
		//                            ... ---- va   ^
		//    ---- va                           \  /
		//          \                            v 
		//           \         =>                 \
		//            \                            \
		//   vb ----- vc ->          vb ---------- vc_os ->   [OS]
		// I -----------------     I -----------------------------
		//   v0 ----- vc_ts ->       v0 ---------- vc_ts ->   [TS]
		//             /                            /
		//            /                            /
		//          ...                          ...
		// 
		// We assume there exists a constrained vertex vc issued from another previous with I.
		// Necessarily there exists at least two other vertices va and vb, with vb constrained by I.
		// We are going to turn the sequence of vertices (va vc vb), and a sequence of edges (ea = va_vc, eb = vc_vb)
		// into (va v vc' vb) and (ea' = va_v, ec' = v_vc', eb' = vb_vc'.

		// Part 1.
		// Identification of the sequence

		Polygon_Vertex *v_a, *v_b;
		Polygon_Edge *e_a, *e_b;

		Polygon_Vertex *v1 = (v_c->e1->v1 == v_c ? v_c->e1->v2 : v_c->e1->v1), *v2 = (v_c->e2->v1 == v_c ? v_c->e2->v2 : v_c->e2->v1);
		assert(v1->is_constrained_by(I) ^ v2->is_constrained_by(I));
		if (v1->is_constrained_by(I)) {
			v_a = v2; e_a = v_c->e2;
			v_b = v1; e_b = v_c->e1;
		} else {
			v_a = v1; e_a = v_c->e1;
			v_b = v2; e_b = v_c->e2;
		}

		// Part 2.
		// Modification of the sequence

		const Constraint C_ts = vc_ts->get_constraint();
		const Constraint C_os = Constraint(C_ts.first, C_ts.second == PLUS ? MINUS : PLUS);

		Polygon_Vertex_R* vc_os = new Polygon_Vertex_R(vc_ts, SCHEDULE);
		//Polygon_Vertex* vc_os = new Polygon_Vertex(v->id_plane, t, vc_ts->M, vc_ts->dM, C_os, nullptr, SCHEDULE);

		v_c->to_r()->transfer_segments(vc_os->to_r(), t, V_t);

		delete_object(e_a, e_b, edges, true);
		delete_object(v_c->to_base(), vertices, true);

		insert_object(v->to_base() , vertices);
		insert_object(vc_os->to_base(), vertices);

		Polygon_Edge* e_a_prime = new Polygon_Edge(v->id_plane, v_a, v);
		Polygon_Edge* e_c_prime = new Polygon_Edge(v->id_plane, v, vc_os);

		Polygon_Edge* e_b_prime = new Polygon_Edge(v->id_plane, v_b, vc_os, C_os);
		insert_object(e_a_prime, e_c_prime, e_b_prime, edges);
	}



	// Handlers
	// Case C1


	void Polygon::remove_vertices_colliding_in_dead_end(Polygon_Vertex_R* v_1, Polygon_Vertex_R* v_2, const FT & t)
	{
		// v1 and v2 are constrained vertices on lines L1 and L2.
		// They meet at time t, as v1 intersects L2 and v2 intersects L1.
		// They should be merged as one still vertex.

		// Illustration :
		//
		//      / v1 --- v2 \               / v \
		//     / /         \ \             / / \ \
		//    / /           \ \     =>    / /   \ \
		//   / v0           v3 \         / v0   v3 \
		//  /                   \       /           \
		// L1                   L2     L1           L2

		// We should turn a sequence of vertices (v0, v1, v2, v3) and a sequence of edges (e_01 = v0_v1, e_12 = v1_v2, e_23 = v2_v3) 
		// into new sequences (v0, v, v3) and (e = v_v0, e' = v_v3).

		if (running_vertices == 2) {
			forget();
			return;
		}

		// Part 1. Identification of the sequence

		Polygon_Vertex *v_0, *v_3;
		Polygon_Edge *e_01, *e_12, *e_23;

		Polygon_Edge* e = v_1->e1;
		Polygon_Vertex* v_e = (e->v1 == v_1 ? e->v2 : e->v1);
		if (v_e != v_2) {
			e_01 = e;
			e_12 = v_1->e2;
			v_0 = v_e;
		} else {
			e_12 = e;
			e_01 = v_1->e2;
			v_0 = (e_01->v1 == v_1 ? e_01->v2 : e_01->v1);
		}

		e_23 = (v_2->e1 == e_12 ? v_2->e2 : v_2->e1);
		v_3 = (e_23->v1 == v_2 ? e_23->v2 : e_23->v1);

		// Part 2. Modifies the sequence

		const Constraint C_1 = v_1->get_constraint();
		const Constraint C_2 = v_2->get_constraint();

		Polygon_Vertex* v = new Polygon_Vertex_S(v_1->id_plane, C_1, C_2);

		delete_object(e_01, e_12, e_23, edges, true);
		delete_object(v_1->to_base(), v_2->to_base(), vertices, true);

		insert_object(v, vertices);

		Polygon_Edge* v_v0 = new Polygon_Edge(v->id_plane, v, v_0, C_1);
		Polygon_Edge* v_v3 = new Polygon_Edge(v->id_plane, v, v_3, C_2);

		insert_object(v_v0, v_v3, edges);

		/*if (running_vertices == 0) {
			forget();
		}*/
	}



	// Handlers
	// Case C2

	void Polygon::redirect_constrained_vertex(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t,
		Polygon_Vertex_S* & v_de, Polygon_Vertex_R* & v_ts, Constraint & C_ts)
	{
		// At time t, a constrained vertex v intersects line I.
		// Therefore we should create a redirected vertex that slides along I.

		// Illustration :
		//
		//           \   I_v                   \   I_v
		//            \  /         v2 ---- v_ts \  /
		//             \/                     \  \/
		// v2 ------ v /\     =>            v_de /\
		//          / /  \                   /  /  \
		//         / /    I                 /  /    I
		//        / /                      /  /
		//      v0 /                      v0 /

		// We turn a sequence of vertices (v0 v v2) into a sequence (v0 v_de v_ts v2)
		// and a sequence of edges (e_0, e_2) into a sequence (e_d0, e_td, e_t2)

		// Part 1. Identification of the sequence

		Polygon_Vertex *v_0, *v_2;
		Polygon_Edge *v_v0, *v_v2;

		Constraint C_v = v->get_constraint();
		Intersection_Line* I_v = C_v.first;

		if (v->e1->is_constrained_by(I_v)) {
			v_v0 = v->e1; v_v2 = v->e2;
		} else {
			v_v0 = v->e2; v_v2 = v->e1;
		}
		v_0 = (v_v0->v1 == v ? v_v0->v2 : v_v0->v1);
		v_2 = (v_v2->v1 == v ? v_v2->v2 : v_v2->v1);

		// Part 2. Modifies the sequence

		CGAL_Point_2 M; CGAL_Vector_2 dM;
		Polygon_Edge::intersection_pt_dir(I, v, v_2, t, V_t, M, dM);
		// v_v2->intersection_pt_dir(I, t, M, dM);

		C_ts = Constraint(I, I->sign(v, 0));

		v_de = new Polygon_Vertex_S(v->id_plane, C_v, C_ts);
		v_ts = new Polygon_Vertex_R(v->id_plane, t, M, dM, C_ts, I_v, v_2->to_r(), 0, SCHEDULE);

		delete_object(v_v0, v_v2, edges, true);
		delete_object(v->to_base(), vertices, false);

		insert_object(v_de->to_base(), v_ts->to_base(), vertices);

		Polygon_Edge* e_d0 = new Polygon_Edge(v->id_plane, v_0, v_de, C_v);
		Polygon_Edge* e_td = new Polygon_Edge(v->id_plane, v_ts, v_de, C_ts);
		Polygon_Edge* e_t2 = new Polygon_Edge(v->id_plane, v_ts, v_2);

		insert_object(e_d0, e_td, e_t2, edges);
	}



	Polygon* Polygon::build_as_constrained_vertex_propagates(Polygon_Vertex_R* v, const FT & t,
		const int seed, Polygon_Vertex_S* vd_ts, Polygon_Vertex_R* v_ts,
		const Constraint & C_os, const Constraint & C_v, Polygon_Vertex_R* & v_os)
	{
		// This is the sequel function of 'redirect_constrained_vertex'.
		// Vertex v, constrained by I_v, propagates beyond the line I.

		// We initialize a still vertex with the same constrained as vd_ts,
		// then a vertex v_os that is the reflection of v_ts,
		// and make a polygon with these three vertices.

		int pid = v->id_plane;

		Intersection_Line* I_v = C_v.first;

		v_os = new Polygon_Vertex_R(v_ts->to_r(), SCHEDULE);
		//v_os = new Polygon_Vertex(pid, t, v_ts->M, v_ts->dM, C_os, I_v, SCHEDULE);

		Polygon_Vertex* vd_os = new Polygon_Vertex_S(pid, C_v, C_os);

		Polygon_Edge* e_cos = new Polygon_Edge(pid, vd_os, v_os, C_os);
		Polygon_Edge* e_cv = new Polygon_Edge(pid, vd_os, v, C_v);
		Polygon_Edge* e = new Polygon_Edge(pid, v_os, v);

		std::list<Polygon_Vertex*> V;
		V.push_back(v); V.push_back(v_os); V.push_back(vd_os);

		std::list<Polygon_Edge*> E;
		E.push_back(e_cos); E.push_back(e_cv); E.push_back(e);

		return new Polygon(seed, V, E);
	}



	// Handlers
	// Case C2*


	void Polygon::redirect_constrained_vertex(Intersection_Line* I, Polygon_Vertex_R* v, Polygon_Edge* e, const FT & t,
		const std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & D, const std::vector<std::pair<Constraint, Constraint> > & C,
		const std::list<Intersection_Line*> & I_discarded, Polygon_Vertex_R* & v_ts)
	{
		// Same as before, except that we don't need to compute the directions and the constraints
		// of the new vertices : they are already known via vectors D and C

		// Step 1 : identifies the sequence

		Polygon_Vertex *v_0, *v_2;
		Polygon_Edge *v_v0, *v_v2;

		if (v->e1 == e) {
			v_v2 = v->e1; v_v0 = v->e2;
		} else {
			v_v0 = v->e1; v_v2 = v->e2;
		}
		v_0 = (v_v0->v1 == v ? v_v0->v2 : v_v0->v1);
		v_2 = (v_v2->v1 == v ? v_v2->v2 : v_v2->v1);

		// Step 2 : modifies the sequence

		//v->delimit_opposite_bidirectional_segments(I);

		delete_object(v_v0, v_v2, edges, true);
		delete_object(v->to_base(), vertices, false);

		const CGAL_Point_2 & M = D[1].first;
		const CGAL_Vector_2 & dM = D[1].second;
		const Constraint C_0 = C[0].first, C_1 = C[1].first;

		Polygon_Vertex_S* v_de = new Polygon_Vertex_S(v->id_plane, C_0, C_1);
		v_ts = new Polygon_Vertex_R(v->id_plane, t, M, dM, C_1, I_discarded, nullptr, 0, SCHEDULE);

		insert_object(v_de->to_base(), v_ts->to_base(), vertices);

		Polygon_Edge* e_d0 = new Polygon_Edge(v->id_plane, v_0, v_de, C_0);
		Polygon_Edge* e_td = new Polygon_Edge(v->id_plane, v_ts, v_de, C_1);
		Polygon_Edge* e_t2 = new Polygon_Edge(v->id_plane, v_ts, v_2);

		insert_object(e_d0, e_td, e_t2, edges);
	}



	Polygon* Polygon::build_as_constrained_vertex_propagates(Polygon_Vertex_R* v, const FT & t,
		const int seed, const int k,
		const std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & D, const std::vector<std::pair<Constraint, Constraint> > & C,
		const std::list<Intersection_Line*> & I_discarded, Polygon_Vertex_R* & v_k, Polygon_Vertex_R* & v_l)
	{
		const int n = int(C.size());
		const int pid = v->id_plane;

		// This function is called as a vertex constrained by I_0 intersects a multiple set of lines I = { I_1 .. I_n }.
		// Here, we build the polygon between lines I_k and I_{k + 1}.

		// First let us build the vertices

		const Constraint C_k = C[k].second;

		const CGAL_Point_2 & M_k = D[k].first;
		const CGAL_Vector_2 & dM_k = D[k].second;
		v_k = new Polygon_Vertex_R(pid, t, M_k, dM_k, C_k, I_discarded, nullptr, 0, SCHEDULE);

		const Constraint C_l = C[(k + 1) % n].first;
		v_l = nullptr;

		if (k + 1 == n) {
			v_l = v;
		} else {
			const CGAL_Point_2 & M_l = D[k + 1].first;
			const CGAL_Vector_2 & dM_l = D[k + 1].second;
			v_l = new Polygon_Vertex_R(pid, t, M_l, dM_l, C_l, I_discarded, nullptr, 0, SCHEDULE);
		}

		Polygon_Vertex* v_kl = new Polygon_Vertex_S(pid, C_k, C_l);

		// Then the edges

		Polygon_Edge* e_k = new Polygon_Edge(pid, v_k, v_kl, C_k);
		Polygon_Edge* e_l = new Polygon_Edge(pid, v_l, v_kl, C_l);
		Polygon_Edge* e_kl = new Polygon_Edge(pid, v_k, v_l);

		// Assembles everything into a polygon

		std::list<Polygon_Vertex*> V;
		V.push_back(v_k); V.push_back(v_l); V.push_back(v_kl);

		std::list<Polygon_Edge*> E;
		E.push_back(e_k); E.push_back(e_l); E.push_back(e_kl);

		return new Polygon(seed, V, E);
	}


	
	// Event handlers
	// Case D1


	void Polygon::remove_intersectant_edge(Intersection_Line* I, Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, const FT & t,
		const CGAL_Point_2 & V1_t, const CGAL_Point_2 & V2_t,
		Polygon_Edge* & e, Polygon_Vertex* & v1_ts, Polygon_Vertex* & v2_ts,
		Constraint & C_1, Constraint & C_2, Constraint & C_ts)
	{
		// At time t, the edge e = (v1 v2) intersects I.
		// We remove e from the polygon, and replace vertices v1 and v2,
		// by v1_ts and v2_ts which are constrained by I and take the same constraints as v1 and v2 if they exist.

		// Illustration :
		//
		//     |                              | ^                ^
		//     |                              | |               /
		//     |                              | v1 ---------- v2
		//     |                              | |              \
		//     |                              | v1_os -------- v2_os -->
		// I -------------------------    I --+-----------------------------
		//     |                              | v1_ts --------- v2_ts -->
		//     | ^                 ^          | |       e_ts      \
		//     | |                /           | |                  \
		//     | v1 ---------- v2             | |                   \
		//     | |      e       \             | | e1_ts       e2_ts  \
		//     | | e1         e2 \            | |                     \
		//     | |                \           | |                      \
		//     | v0                v3         | v0                     v3
		//

		// Step 1.
		// Identification.

		int id = v1->id_plane;

		Polygon_Edge *e1, *e2;
		Polygon_Vertex *v0, *v3;

		e = v1->e1;
		if ((e->v1 == v1 && e->v2 == v2) || (e->v1 == v2 && e->v2 == v1)) {
			e1 = v1->e2;
		} else {
			e = v1->e2;
			e1 = v1->e1;
		}
		e2 = (v2->e1 == e ? v2->e2 : v2->e1);
		v0 = (e1->v1 == v1 ? e1->v2 : e1->v1);
		v3 = (e2->v1 == v2 ? e2->v2 : e2->v1);

		C_1 = (v1->unconstrained() ? Constraint() : v1->get_constraint());
		C_2 = (v2->unconstrained() ? Constraint() : v2->get_constraint());

		//if (C_1.first != nullptr) v1->delimit_opposite_bidirectional_segments(I);
		//if (C_2.first != nullptr) v2->delimit_opposite_bidirectional_segments(I);

		// Step 2.
		// Generates v1_ts, v2_ts

		C_ts = Constraint(I, I->sign(v1, 0));

		if (C_1.first == nullptr) {
			CGAL_Point_2 M;
			CGAL_Vector_2 dM;
			Polygon_Edge::intersection_pt_dir(I, v1, v0, t, V1_t, M, dM);
			v1_ts = new Polygon_Vertex_R(id, t, M, dM, C_ts, I, v0->to_r(), 0, SCHEDULE);
		} else {
			v1_ts = new Polygon_Vertex_S(id, C_1, C_ts);
		}

		if (C_2.first == nullptr) {
			CGAL_Point_2 M;
			CGAL_Vector_2 dM;
			Polygon_Edge::intersection_pt_dir(I, v2, v3, t, V2_t, M, dM);
			v2_ts = new Polygon_Vertex_R(id, t, M, dM, C_ts, I, v3->to_r(), 0, SCHEDULE);
		} else {
			v2_ts = new Polygon_Vertex_S(id, C_2, C_ts);
		}

		// Step 3
		// Replaces

		delete_object(e, edges, false);
		delete_object(e1, e2, edges, true);

		delete_object(v1->to_base(), v2->to_base(), vertices, false);

		Polygon_Edge* e_ts = new Polygon_Edge(id, v1_ts, v2_ts, C_ts);
		Polygon_Edge *e1_ts, *e2_ts;

		if (C_1.first != nullptr) {
			e1_ts = new Polygon_Edge(id, v0, v1_ts, C_1);
		} else {
			e1_ts = new Polygon_Edge(id, v0, v1_ts);
		}

		if (C_2.first != nullptr) {
			e2_ts = new Polygon_Edge(id, v3, v2_ts, C_2);
		} else {
			e2_ts = new Polygon_Edge(id, v3, v2_ts);
		}

		insert_object(v1_ts, v2_ts, vertices);
		insert_object(e1_ts, e2_ts, e_ts, edges);
	}


	Polygon* Polygon::build_as_edge_propagates(Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, 
		Polygon_Vertex* v1_ts, Polygon_Vertex* v2_ts, Polygon_Edge* e, const FT & t,
		const int seed, const Constraint & C_ts, Polygon_Vertex* & v1_os, Polygon_Vertex* & v2_os)
	{
		int id = v1->id_plane;

		// Builds vertices v1_os and v2_os, cf. scheme above

		const Constraint C_1 = (v1->unconstrained() ? Constraint() : v1->get_constraint());
		const Constraint C_2 = (v2->unconstrained() ? Constraint() : v2->get_constraint());

		Intersection_Line* I = C_ts.first;
		const Constraint C_os = Constraint(I, C_ts.second == PLUS ? MINUS : PLUS);

		if (C_1.first == nullptr) {
			v1_os = new Polygon_Vertex_R(v1_ts->to_r(), SCHEDULE);
		} else {
			v1_os = new Polygon_Vertex_S(id, C_1, C_os);
		}

		if (C_2.first == nullptr) {
			v2_os = new Polygon_Vertex_R(v2_ts->to_r(), SCHEDULE);
		} else {
			v2_os = new Polygon_Vertex_S(id, C_2, C_os);
		}

		// Builds edges

		Polygon_Edge *e1, *e2;
		if (v1->unconstrained()) {
			e1 = new Polygon_Edge(id, v1, v1_os);
		} else {
			e1 = new Polygon_Edge(id, v1, v1_os, C_1);
		}

		if (v2->unconstrained()) {
			e2 = new Polygon_Edge(id, v2, v2_os);
		} else {
			e2 = new Polygon_Edge(id, v2, v2_os, C_2);
		}

		Polygon_Edge* e_os = new Polygon_Edge(id, v1_os, v2_os, C_os);

		// Builds polygon

		std::list<Polygon_Vertex*> V;
		V.push_back(v1); V.push_back(v2); V.push_back(v1_os); V.push_back(v2_os);

		std::list<Polygon_Edge*> E;
		E.push_back(e); E.push_back(e1); E.push_back(e2); E.push_back(e_os);

		return new Polygon(seed, V, E);
	}




	// Event handlers
	// Case D2


	void Polygon::remove_intersectant_edges(Intersection_Line* I_1, Intersection_Line* I_2, 
		Polygon_Vertex_R* v1, Polygon_Vertex_R* v, Polygon_Vertex_R* v2, const FT & t, 
		const CGAL_Point_2 & V1_t, const CGAL_Point_2 & V_t, const CGAL_Point_2 & V2_t,
		Polygon_Vertex* & v1_ts, Polygon_Vertex_S* & v_de, Polygon_Vertex* & v2_ts)
	{
		// At time t, two edges e1 = (v1 v) and e2 = (v v2) intersect two lines I_1 and I_2.

		// We do the following updates :
		// 1) v is unconstrained and gets replaced by a still vertex.
		// 2) if v1 is unconstrained, it is replaced by a vertex constrained by I_1, obtained by intersection of I_0 and (v0 v1).
		//    otherwise, it is replaced by a still vertex, constrained by I_1 and I_v1.
		// 3) same operation for v2.

		// 
		//                  I_1            I_v2                         I_1               I_v2
		//                  /               |                           /                  |
		//                 /                |                          /                   |
		//                /                 |                         /                    |
		// I_2 --------- + ---------------- + --   -->  I_2 -------- + ------------------- + --
		//              / ^             ^   |                       /  v_de -------- v2_ts |
		//             /   \            |   |                      /  /      e2_ts      |  |
		//            /     v -------- v2   |                     /  / e1_ts            |  |
		//           /     /    e2      |   |                    /  /                   |  |
		//          /     / e1       e3 |   |                   /  /             e3_ts  |  |
		//         /     /              |   |                  / v1_ts                  |  |
		//        /  <- v1             v3   |                 /  /  \                   |  |
		//       /       \                  |                /  v    \  e0_ts          v3  |
		//     ...        \ e0                              ...       \                    |
		//                 \                                           \
		//                 v0                                          v0
		//

		// Step 1.
		// Identification of the sequence of vertices and edges.

		int pid = v->id_plane;

		Polygon_Vertex *v0, *v3;
		Polygon_Edge *e0, *e1, *e2, *e3;

		if (v->e1->other_vertex(v) == v1) {
			e1 = v->e1; e2 = v->e2;
		} else {
			e1 = v->e2; e2 = v->e1;
		}

		e0 = (v1->e1 != e1 ? v1->e1 : v1->e2);
		e3 = (v2->e1 != e2 ? v2->e1 : v2->e2);

		v0 = e0->other_vertex(v1);
		v3 = e3->other_vertex(v2);

		Constraint C_v1 = (v1->unconstrained() ? Constraint() : v1->get_constraint());
		Constraint C_v2 = (v2->unconstrained() ? Constraint() : v2->get_constraint());

		// Step 2.
		// Creates new vertices.

		Sign eps_1 = I_1->sign(v, 0), eps_2 = I_2->sign(v, 0);
		Constraint C_1 = std::make_pair(I_1, eps_1);
		Constraint C_2 = std::make_pair(I_2, eps_2);

		v_de = new Polygon_Vertex_S(pid, C_1, C_2);

		if (C_v1.first == nullptr) {
			CGAL_Point_2 M_1;
			CGAL_Vector_2 dM_1;
			Polygon_Edge::intersection_pt_dir(I_1, v1, v0, t, V1_t, M_1, dM_1);
			v1_ts = new Polygon_Vertex_R(pid, t, M_1, dM_1, C_1, I_1, v0->to_r(), 0, SCHEDULE);
		} else {
			v1_ts = new Polygon_Vertex_S(pid, C_1, C_v1);
		}

		if (C_v2.first == nullptr) {
			CGAL_Point_2 M_2;
			CGAL_Vector_2 dM_2;
			Polygon_Edge::intersection_pt_dir(I_2, v2, v3, t, V2_t, M_2, dM_2);
			v2_ts = new Polygon_Vertex_R(pid, t, M_2, dM_2, C_2, I_2, v3->to_r(), 0, SCHEDULE);
		} else {
			v2_ts = new Polygon_Vertex_S(pid, C_2, C_v2);
		}

		// Step 3.
		// Inserts vertices.

		delete_object(e0, e1, e2, e3, edges, true);
		delete_object(v1->to_base(), v->to_base(), v2->to_base(), vertices, false);

		Polygon_Edge *e0_ts, *e1_ts, *e2_ts, *e3_ts;
		if (C_v1.first == nullptr) {
			e0_ts = new Polygon_Edge(pid, v0, v1_ts);
		} else {
			e0_ts = new Polygon_Edge(pid, v0, v1_ts, C_v1);
		}

		if (C_v2.first == nullptr) {
			e3_ts = new Polygon_Edge(pid, v2_ts, v3);
		} else {
			e3_ts = new Polygon_Edge(pid, v2_ts, v3, C_v2);
		}

		e1_ts = new Polygon_Edge(pid, v1_ts, v_de, C_1);
		e2_ts = new Polygon_Edge(pid, v_de, v2_ts, C_2);

		insert_object(v1_ts, v_de->to_base(), v2_ts, vertices);
		insert_object(e0_ts, e1_ts, e2_ts, e3_ts, edges);
	}


	/*bool Polygon::is_adjacent_to(Polygon* P, Polygon_Edge* & e, Polygon_Edge* & f)
	{
		// Two polygons P and Q are adjacent if and only if there have one edge in common
		// To determine if one edge is shared by this polygon and P, we need to find evidence that its delimiting
		// vertices are both in this object and P. To this end we interpret vertices as intersections of Intersection_Lines.

		// The common edge to this polygon and P is : e in 'this', f in P.

		// Part 1.
		// Finds shared edges e and f, or returns false

		e = nullptr, f = nullptr;

		for (std::list<Polygon_Edge*>::iterator it_e = edges.begin(); it_e != edges.end(); it_e++) {
			Constraint C_e = (*it_e)->get_constraint();

			for (std::list<Polygon_Edge*>::iterator it_f = P->edges.begin(); it_f != P->edges.end(); it_f++) {
				Constraint C_f = (*it_f)->get_constraint();

				if (C_e.first == C_f.first && C_e.second != C_f.second) {
					// Edges e of 'this', and f of polygon P, stand on opposite sides of the same intersection line
					// All polygons are convex, so if 'this' and P are adjacent, 
					// e is necessarily common to P and reciprocally
					e = (*it_e);
					f = (*it_f);
					goto next;
				}
			}
		}

		if (e == nullptr && f == nullptr) return false;

	next:
		// Part 2.
		// Matches vertices of e to vertices of f
		Polygon_Vertex_S *e_v1 = e->v1->to_s(), *e_v2 = e->v2->to_s(), *f_v1 = f->v1->to_s(), *f_v2 = f->v2->to_s();
		if (!(e_v1->represents_same_intersection(f_v1, f_v2) && e_v2->represents_same_intersection(f_v1, f_v2))) {
			return false;
		}

		// Part 3.
		// Determines the existence of a segment separating e from f

		Intersection_Line* I_e = e->get_constraint().first;
		if (I_e->exists_segment_adjacent_to_edge(e)) {
			return false;
		}

		return true;
	}*/



	CGAL_Point_2 Polygon::get_barycenter(const FT & t) const
	{
		FT x = 0, y = 0;
		int n = int(vertices.size());

		for (std::list<Polygon_Vertex*>::const_iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
			CGAL_Point_2 V_t = (*it_v)->pt(t);
			x += V_t.x(), y += V_t.y();
		}
		x /= n, y /= n;

		return CGAL_Point_2(x, y);
	}


	void Polygon::get_active_constrained_vertices(Intersection_Line* I_1, Intersection_Line* I_2, Polygon_Vertex_R* & v_1, Polygon_Vertex_R* & v_2)
	{
		// We assume that P has exactly one vertex constrained by I_1, and exactly one vertex constrained by I_2

		v_1 = nullptr;
		v_2 = nullptr;

		for (std::list<Polygon_Vertex*>::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
			Polygon_Vertex_R* v = (*it_v)->to_r();
			if (v != nullptr && !v->unconstrained()) {
				Intersection_Line* I = v->get_constraint().first;
				if (I == I_1) {
					assert(v_1 == nullptr);
					v_1 = v;
				} else if (I == I_2) {
					assert(v_2 == nullptr);
					v_2 = v;
				}
			}
		}

		assert(v_1 != nullptr && v_2 != nullptr);
	}


	void Polygon::shift_vertices(Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, const CGAL_Point_2 & M, const FT & t)
	{
		// This function is called during the algorithm of lateral propagation,
		// as we hierarchically divide a set of vertices which the same shape as the original polygon.
		// In this set of vertices, positions of all vertices are a bit shifted, ahead of time :
		// otherwise we would always get 0 when computing the signs of all vertices relatively to lines I and I_v.
		// That's what we correct here : at time t, v(t) = M

		v1->set_M(M - t * v1->dM);
		v2->set_M(M - t * v2->dM);
	}



	void Polygon::shift_remaining_vertices_and_schedule_events(Polygon* P_1, Polygon* P_2, Intersection_Line* I_1, Intersection_Line* I_2, const CGAL_Point_2 & M, const FT & t)
	{
		std::list<Polygon*> P(1, P_1);
		if (P_2 != nullptr) P.push_back(P_2);

		std::list<Intersection_Line*> discarded_lines;
		discarded_lines.push_back(I_1);
		discarded_lines.push_back(I_2);

		shift_remaining_vertices_and_schedule_events(P, discarded_lines, M, t);
	}



	void Polygon::shift_remaining_vertices_and_schedule_events(const std::list<Polygon*> P, std::list<Intersection_Line*> I, const CGAL_Point_2 & M, const FT & t)
	{
		// Step 1.
		// Gets all vertices of all elements of P and sorts them by index

		size_t n = 0;
		for (std::list<Polygon*>::const_iterator it_p = P.begin(); it_p != P.end(); it_p++) n += (*it_p)->vertices.size();

		std::vector<Polygon_Vertex_R*> V;
		V.reserve(n);

		for (std::list<Polygon*>::const_iterator it_p = P.begin(); it_p != P.end(); it_p++) {
			Polygon* P = (*it_p);
			for (std::list<Polygon_Vertex*>::iterator it_v = P->vertices.begin(); it_v != P->vertices.end(); it_v++) {
				Polygon_Vertex_R* v = (*it_v)->to_r();
				if (v != nullptr) V.push_back(v);
			}
		}

		V.shrink_to_fit();

		struct _Vertices_ID_Comparator {
			bool operator() (const Polygon_Vertex* vl, const Polygon_Vertex* vr) {
				return vl->id_object < vr->id_object;
			}
		} Vertices_ID_Comparator;

		std::sort(V.begin(), V.end(), Vertices_ID_Comparator);

		for (size_t i = 0; i < V.size(); i++) {
			Polygon_Vertex_R* v = V[i];

			// We only have to correct the locations of all unconstrained vertices
			// indeed, the locations of all inactive vertices are already equal to vertex M
			// and the locations of constrained vertices have been corrected before,
			// by a direct call to function shift_vertex
			if (v->unconstrained()) {
				v->set_M(M - t * v->dM);
			}

			// Schedule events but doesn't sort them
			if (v->is_independent()) {
				v->schedule_events(I, nullptr);
			} else {
				Polygon_Vertex_R* v_ts = v->get_paired_vertex();
				v->schedule_events(v_ts);
			}
		}
	}


	void Polygon::forget()
	{
		// This function is called where this polygon has no more active vertices,
		// or two vertices that are going to intersect.

		// There is no need to keep all the data related to vertices and edges in memory.
		// Actually, we are only interested in the sequence of contours that delimits the sequence :
		// vertices can be obtained later, via the 'landmarks' attribute in the class Support_Plane.

		if (cell->contours_are_empty()) {

			std::list<Constraint> circular_constraints;
			Polygon_Edge *e_init = edges.front();
			Polygon_Edge *e_prev = nullptr, *e_curr = e_init;
			do {
				Constraint C = e_curr->get_constraint();
				if (C.first != nullptr) circular_constraints.push_back(C);

				Polygon_Edge* e_next = e_curr->v1->other_edge(e_curr);
				if (e_next == e_prev) {
					e_next = e_curr->v2->other_edge(e_curr);
				}
				e_prev = e_curr;
				e_curr = e_next;
			} while (e_curr != e_init);

			cell->set_contours(circular_constraints);
		}

		// Cleans the polygon and frees the memory we no longer need.

		for (std::list<Polygon_Edge*>::iterator it_e = edges.begin() ; it_e != edges.end() ; ++it_e) {
			delete (*it_e);
		}
		edges.clear();

		for (std::list<Polygon_Vertex*>::iterator it_v = vertices.begin() ; it_v != vertices.end() ; ++it_v) {
			delete (*it_v);
		}
		vertices.clear();
	}
}