#include "../include/support_plane.h"
#include "../include/intersection_line.h"
#include "../include/polygon_vertex.h"
#include "../include/polygon_edge.h"
#include "../include/universe.h"
#include "../include/parameters.h"


namespace Skippy 
{
	using CGAL::to_double;

	//
	// HANDLERS CASE A
	//


	//#define PRINT_EVENTS


	void Support_Plane::unconstrained_vertex_intersects_line(const std::list<Event_Vertex_Line*> & E, Polygon_Vertex_R* v, const CGAL_Point_2 & V_t)
	{
		int E_size = int(E.size());
		if (E_size != 1) {
			throw std::logic_error("Error : unconstrained vertex intersects several lines : not implemented yet.");
		}

#ifdef PRINT_EVENTS
		append(E.front(), "A");
#endif

		Event_Vertex_Line* e = E.front();
		Intersection_Line* I = get_line_by_identifier(e->intersected);
		const FT t = e->t_intersectant;
		const CGAL_Point_3 W_t = backproject(V_t);

		// Deletes events
		v->delete_events(E);
		v->decrement_queued_events(E_size);

		// Part 1.
		// Splits v, which intersects I, into two vertices moving towards opposite directions.
		// Vertices and constraints are suffixed by 'ts' for 'this side' of I, and 'os' which means 'other side' here

		Polygon_Vertex_R *v1_ts, *v1_os, *v2_ts, *v2_os;
		Constraint C_ts;

		Polygon* P_ts = v->get_polygon();
		P_ts->split_unconstrained_vertex(I, v, t, V_t, v1_ts, v2_ts, C_ts);

		// Part 2.
		// If v keeps propagating beyond I, then we initialize another polygon with v, v1_os and v2_os, 
		// otherwise, it is deleted.

		Signature S_os = polygon_set->get_adjacent_polygons_signature(P_ts, I);

		if (propagation_continues_outside_intersections(I, v, false, V_t, W_t, t, S_os, P_ts->seed)) {

			Polygon* P_os = Polygon::build_as_unconstrained_vertex_propagates(v, t, P_ts->seed, v1_ts, v2_ts, C_ts, v1_os, v2_os);
			polygon_set->insert(S_os, P_os);

			// Two vertices on both sides of I lead four segments
			init_4_uplets_of_bidirectional_segments_of_null_length(I, V_t, W_t, v1_ts, v1_os, v2_ts, v2_os, t);

			// If necessary, computes more events for v.
			v->reschedule_events(); ++v->crossed_lines;

		} else {
			init_2_uplets_of_bidirectional_segments_of_null_length(I, V_t, W_t, v1_ts, v2_ts, t);
			v->stop(t);
			delete v;
		}
	}



	//
	// HANDLERS CASE B
	//



	void Support_Plane::unconstrained_vertex_intersects_line_kth_time(const std::list<Event_Vertex_Line*> & E, Polygon_Vertex_R* v, const CGAL_Point_2 & V_t)
	{
		int E_size = int(E.size());
		if (E_size != 1) {
			throw std::logic_error("Error : unconstrained vertex intersects several lines : not implemented yet.");
		}

#ifdef PRINT_EVENTS
		append(E.front(), "B");
#endif

		Event_Vertex_Line* e = E.front();
		Intersection_Line* I = get_line_by_identifier(e->intersected);
		const FT t = e->t_intersectant;

		// Deletes events
		v->delete_events(E);
		v->decrement_queued_events(E_size);

		// As v intersects I, it also intersects a constrained vertex along I, which is part of the same polygon as it.
		// We call a function that applies the following actions on the polygon : 
		// a) Update of the list of vertices and edges, in particular, creation of a new constrained vertex ;
		// b) Rescheduling of events involving the segment obtained by intersection of P and I.

		Polygon_Vertex_R* vc_ts;
		Polygon* P_ts = v->get_polygon();
		CGAL_Vector_2 u_ref;

		P_ts->remove_vertices_colliding_on_intersection_line(I, v, t, V_t, vc_ts, u_ref);

		// Even if a segment is intersected by v, the vertex is transferred to the polygon located on the other side of I.
		// Indeed, if we don't transfer the vertex, then there is a risk to create edges of a polygon that won't be adjacent to any segment
		// Only one possible exception : v intersects a planar segment

		if (!I->is_border) {

			Polygon_Vertex_R* vc = nullptr;
			Polygon* P_os = polygon_set->get_adjacent_polygon(P_ts, I);
			if (P_os != nullptr) vc = P_os->get_tied_vertex(I, v, t, V_t, u_ref);

			if (vc != nullptr) {
				P_os->insert_unconstrained_vertex(I, v, t, V_t, vc_ts, vc);

				// If necessary, computes more events for v.
				v->reschedule_events(); ++v->crossed_lines;

			} else {
				v->stop(t);
				delete v;
			}

		} else {
			v->stop(t);
			delete v;
		}
	}



	//
	// HANDLERS CASE C1
	//



	void Support_Plane::constrained_vertices_intersect(const std::list<Event_Vertex_Line*> & E_1, const std::list<Event_Vertex_Line*> & E_2, 
		Polygon_Vertex_R* v1, Polygon_Vertex_R* v2)
	{
		// If a constrained vertex v1 intersects another constrained vertex at the intersection of n lines, n > 2, 
		// then we make this n-uplet of concurrent lines known by the algorithm

#ifdef PRINT_EVENTS
		append(E.front(), "C1");
#endif

		int E_size_1 = int(E_1.size());
		int E_size_2 = int(E_2.size());

		if (E_size_1 > 1) {
			std::vector<Intersection_Line*> I;

			I.push_back(v1->get_constraint().first);
			for (std::list<Event_Vertex_Line*>::const_iterator it_e = E_1.begin(); it_e != E_1.end(); ++it_e) {
				I.push_back(get_line_by_identifier((*it_e)->intersected));
			}
			insert_triplets_of_concurrent_lines(I);
		}

		Event_Vertex_Line* e = E_1.front();
		const FT t = e->t_intersectant;

		// Deletes event
		v1->delete_events(E_1);
		v2->delete_events(E_2);
		v1->decrement_queued_events(E_size_1);
		v2->decrement_queued_events(E_size_2);

		// Stops segments of v1 and v2
		Polygon* P_ts = v1->get_polygon();
		int seed = P_ts->seed;

		Intersection_Line* I_1 = v1->get_constraint().first;
		Intersection_Line* I_2 = v2->get_constraint().first;
		v1->stop(t);
		v2->stop(t);
		v1->stop_segments(I_2, seed, t);
		v2->stop_segments(I_1, seed, t);

		// Merge vertices
		P_ts->remove_vertices_colliding_in_dead_end(v1, v2, t);
	}



	//
	// HANDLERS CASE C2
	//



	void Support_Plane::constrained_vertex_intersects_line(const std::list<Event_Vertex_Line*> & E, Polygon_Vertex_R* v, const CGAL_Point_2 & V_t)
	{
		int E_size = int(E.size());
		v->decrement_queued_events(E_size);

		if (E.size() == 1) {
			Event_Vertex_Line* e = E.front();
			constrained_vertex_intersects_single_line(e, v, V_t);
		} else {
			// throw std::logic_error("Error : constrained vertex intersects several lines : disabled case -> recompile.");
			constrained_vertex_intersects_multiple_lines(E, v, V_t);
		}
	}



	void Support_Plane::constrained_vertex_intersects_single_line(Event_Vertex_Line* e, Polygon_Vertex_R* v, const CGAL_Point_2 & V_t)
	{
#ifdef PRINT_EVENTS
		append(e, "C2");
#endif

		Intersection_Line* I_v = v->get_constraint().first;
		Intersection_Line* I = get_line_by_identifier(e->intersected);
		const FT t = e->t_intersectant;

		v->delete_event(e);

		CGAL_Point_3 W_t = CGAL::ORIGIN;
		if (Universe::params->stopping_condition >= 1) W_t = backproject(V_t);

		/*
		// At time t, we get the following situation :
		//             (H12)
		//  \             \                                          \ _     ^   \  ^
		//    \            \                                             \ _  \   \  \
		// [P]  \    (Q1)   \           (Q2)                       [P]        v_ts \ v_os
		//        \          \                                                  \   \  \  
		// ------- v -->      \               [C_v_ts]            ------------- v_de \  K ----- v ->      [C_v_ts]
		// (H14) ------------- M -------------- (H23)  I_v   -->  ------------------- M -------------------------- I_v
		//                      \             [C_v_os]                   <- vb --- K'' \  K' -- vl ->     [C_v_os]
		//                       \                                             \     \  \  \   /
		//     (Q4)               \         (Q3)                                  vl_ts  \ vl_os
		//                         \                                                  \   \  \
		//                          \                                                  v   \  v
		//                         (H34)                                                    \
		//                      [C_ts] I [C_os]                                       [C_ts] I [C_os]   

		// We observe four quarters of plane Q1 = (C_ts, C_v_ts), Q2 = (C_os, C_v_ts), Q3 = (C_os, C_v_os), and Q4 = (C_ts, C_v_os).
		// The vertex v arrives from a polygon P, located in Q1. 

		// The actions to perform are the following :
		// (Q1) : we create a still vertex in M, and a vertex constrained by I by ricochet of v on I.
		// Then if v can cross I :
		// (Q2) : inserts a node in the tree, corresponding to a new polygon that will contain v.
		// Then if there is no segment containing M in the halfline [C_os, ...) on I_v :
		// (Q3) : inserts a node in the tree, to prevent a hole to appear.
		// Finally, if there is no segment containing M in the halfline [C_v_os, ...) on I :
		// (Q4) : inserts a node in the tree, to prevent a hole to appear. */

		Polygon* P_1 = v->get_polygon();
		const int seed = P_1->seed;

		Constraint C_v_ts = v->get_constraint();
		Constraint C_v_os = Constraint(C_v_ts.first, (C_v_ts.second == PLUS ? MINUS : PLUS));
		Constraint C_ts, C_os;

		// Quarter of plane (Q1).
		// We create two vertices v_de and v_ts. v_de is at the intersection of I and I_v (dead-end),
		// whereas vertex v_ts that corresponds to the reflection of v on I. Both are inserted into polygon P.

		Polygon_Vertex_S *v_de;
		Polygon_Vertex_R *v_ts, *v_os, *vl, *vl_os, *vl_ts, *vb;
		P_1->redirect_constrained_vertex(I, v, t, V_t, v_de, v_ts, C_ts);
		C_os = Constraint(C_ts.first, (C_ts.second == PLUS ? MINUS : PLUS));
		const CGAL_Point_2 & M = v_de->get_M();

		// Quarter of plane (Q2).
		// We assert that v is allowed to cross I, by absence of segment along its way.
		// Now, we test the existence of a polygon on the other side of I. 
		// We define two booleans for answering this question, which say that I or I_v should be used to split the node. 
		// Indeed, P_ts and P_os, the polygon obtained as v crosses I are not necessarily two direct children of a same node :
		// it depends on the other of the Intersection_Lines as we go down the tree.

		Signature S_2 = polygon_set->get_adjacent_polygons_signature(P_1, I);
		if (propagation_continues_at_intersection(I, v, false, V_t, W_t, C_v_ts, t, S_2, seed)) {

			v->extend_segments(I);

			Polygon* P_2 = Polygon::build_as_constrained_vertex_propagates(v, t, seed, v_de, v_ts, C_os, C_v_ts, v_os);
			polygon_set->insert(S_2, P_2);

			// Two segments should be created on (H12), led by v_ts and v_os
			init_2_uplets_of_unidirectional_segments_of_null_length(I, V_t, v_ts, v_os, C_v_ts, t);

			// Quarter of plane (Q3).
			// We test if M is contained by any segment on line I_v, on the side C_os.
			// If not, and if there is no polygon in the quarter of plane (Q3), 
			// then we create one and initialize couples of segments on halfline (H23).
			// Otherwise, we only initialize a segment on halfline (H23).

			Signature S_3 = polygon_set->get_adjacent_polygons_signature(P_1, I, I_v);
			Signature S_4 = polygon_set->get_adjacent_polygons_signature(P_1, I_v);

			if (propagation_continues_at_intersection(I_v, nullptr, false, V_t, W_t, C_os, t, S_3, seed)) {

				// If the polygon is allowed to propagate laterally, then we initialize a local Polygon_Tree
				// with all vertices centered in M, that we hierarchically divide using I_v and I.
				// We then obtain the desired leaves of the tree, that contains the polygons in areas (Q3) and (Q4).

				CGAL_Vector_2 OM = M - polygon_directions[seed]->get_barycenter();
				Polygon_Tree* T = new Polygon_Tree(v->id_plane, seed, polygon_directions[seed], OM, t);

				T->split(I_v, t, 0, false);
				Polygon_Tree* T_34 = (C_v_os.second == PLUS ? T->subtree_plus : T->subtree_minus);

				T_34->split(I, t, 0, false);
				Polygon_Tree* T_3 = (C_os.second == PLUS ? T_34->subtree_plus : T_34->subtree_minus);
				Polygon_Tree* T_4 = (C_ts.second == PLUS ? T_34->subtree_plus : T_34->subtree_minus);

				// Removes the polygon P_3 from node T_3 and initializes a new cell
				Polygon *P_3 = T_3->remove_reference_to_polygon(), *P_4 = nullptr;
				polygon_set->insert(S_3, P_3);

				// Gets the vertex constrained by I_v and adjacent to v to initialize a couple of segments
				P_3->get_active_constrained_vertices(I_v, I, vl, vl_os);
				P_3->shift_vertices(vl, vl_os, M, t);

				// We cannot initialize a 2-uplet of segments representing v and vl, 
				// since they don't necessarily propagate with the same speed
				// init_1_uplets_of_unidirectional_segments(I_v, v, C_os, t);
				init_1_uplets_of_unidirectional_segments_of_null_length(I_v, V_t, vl, C_os, t);

				// Quarter of plane (Q4).
				// We test if M is contained by any segment in line I, on the side C_v_os.
				// If not, and if there is no polygon in the quarter of plane (Q4), we create one, 
				// and set two segments on (H34), and one segment on (H14).
				// Otherwise, we only initialize one segment on that halfline.

				if (propagation_continues_at_intersection(I, nullptr, false, V_t, W_t, C_v_os, t, S_4, seed)) {

					P_4 = T_4->remove_reference_to_polygon();
					polygon_set->insert(S_4, P_4);

					P_4->get_active_constrained_vertices(I_v, I, vb, vl_ts);
					if (vl_ts->id_object < vl_os->id_object) {
						Polygon_Vertex_R::set_paired_vertices(vl_ts, vl_os);
					} else {
						Polygon_Vertex_R::set_paired_vertices(vl_os, vl_ts);
					}
					P_4->shift_vertices(vb, vl_ts, M, t);

					init_2_uplets_of_unidirectional_segments_of_null_length(I, V_t, vl_os, vl_ts, C_v_os, t);
					init_1_uplets_of_unidirectional_segments_of_null_length(I_v, V_t, vb, C_ts, t);

				} else {
					init_1_uplets_of_unidirectional_segments_of_null_length(I, V_t, vl_os, C_v_os, t); // (H34)
				}

				// Half-plane (Q3) and (Q4)
				// We can now schedule events of all vertices included in P3 and P4.
				Polygon::shift_remaining_vertices_and_schedule_events(P_3, P_4, I_v, I, M, t);
				delete T;

			}

			// If necessary, computes more events for v.
			v->reschedule_events();
			++v->crossed_lines;

		} else {
			// v can't cross I.
			// We stop the segment carried by this vertex, and initialize one segment on halfline (H12).
			v->stop(t);
			v->stop_segments(I, seed, t);
			init_1_uplets_of_unidirectional_segments_of_null_length(I, V_t, v_ts, C_v_ts, t);
			delete v;
		}
	}



	void Support_Plane::constrained_vertex_intersects_multiple_lines(const std::list<Event_Vertex_Line*> & E, Polygon_Vertex_R* v, const CGAL_Point_2 & V_t)
	{
#ifdef PRINT_EVENTS
		append(E.front(), "C2 *");
#endif

		const FT t = E.front()->t_intersectant;
		CGAL_Point_3 W_t = CGAL::ORIGIN;
		if (Universe::params->stopping_condition >= 1) W_t = backproject(V_t);

		// Part 1.
		// At time t, constrained v intersects a set of lines I = {I_1, I_2, I_n}.
		// We first need to sort this line by angle made with the line I_0, to determine in what order polygons should be added to the support plane.

		std::vector<Intersection_Line*> I;
		std::vector<std::pair<Constraint, Constraint> > C;

		get_sorted_lines(v, V_t, E, I);
		get_sorted_constraints(v, I, C);
		insert_triplets_of_concurrent_lines(I);

		v->delete_events(E);

		// Part 2.
		// Considering the half-plane to which v belongs (i.e. C_0_ts), 
		// we iteratively build triangles that correspond to each fraction of plane,
		// except for the first fraction of plane, that contains the polygon including v.
		// If v was able to cross all the lines of I, and if there is no segment on I[0] beyond V_t,
		// then the polygon should be propagated laterally between all fractions of planes.

		bool propagate_laterally = constrained_vertex_intersects_multiple_lines_ts(v, t, V_t, W_t, I, C);

		// Part 3.
		// Propagates polygon laterally.

		if (propagate_laterally) {
			constrained_vertex_intersects_multiple_lines_os(v, t, V_t, W_t, I, C);
		}
}



	bool Support_Plane::constrained_vertex_intersects_multiple_lines_ts(Polygon_Vertex_R* v, const FT t, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t,
		const std::vector<Intersection_Line*> & I, const std::vector<std::pair<Constraint, Constraint> > & C)
	{
		int n = int(I.size());
		const Constraint & C_0_ts = C[0].first, C_0_os = C[0].second;

		std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > D;
		Polygon_Edge* e = (v->e1->is_constrained_by(I[0]) ? v->e2 : v->e1);
		get_definitions_of_constrained_vertices(e, t, I, D);

		// We consider the polygon located between I[0] and I[1], the one that includes v.
		// We generate one constrained vertex, by intersection of e and I[1].
		// When scheduling its list of events, we take none of the lines of I into account

		std::list<Intersection_Line*> I_discarded;
		std::copy(I.begin(), I.end(), std::back_inserter(I_discarded));

		Polygon* P_init = v->get_polygon();
		const int seed = P_init->seed;
		Signature S = P_init->get_cell()->get_signature();

		Polygon_Vertex_R* v_1;
		P_init->redirect_constrained_vertex(I[1], v, e, t, D, C, I_discarded, v_1);
		init_1_uplets_of_unidirectional_segments_of_null_length(I[1], V_t, v_1, C_0_ts, t);

		// Propagation in fractions of planes determined by I[k] and I[k + 1]

		for (int k = 1; k < n; k++) {

			// Before building each polygon we have to determine if the polygon can cross the line I[k]
			std::list<Constraint> C_limits_k;
			get_stop_constraints_at_multiple_intersections(C, k, true, C_limits_k);
			polygon_set->get_signature_of_adjacent_cell(S, I[k]);

			if (propagation_continues_at_intersection(I[k], v, false, V_t, W_t, C_limits_k, t, S, seed)) {
				++v->crossed_lines;

				// v crosses the line I[k] : we build a triangle
				// and initialize two segments led by v_k and v_l
				Polygon_Vertex_R *v_k, *v_l;
				Polygon* P = Polygon::build_as_constrained_vertex_propagates(v, t, seed, k, D, C, I_discarded, v_k, v_l);

				init_1_uplets_of_unidirectional_segments_of_null_length(I[k], V_t, v_k, C_0_ts, t);
				if (k + 1 != n) {
					// If k + 1 == n, v_l == v and the segment initially carried by it continues its way
					init_1_uplets_of_unidirectional_segments_of_null_length(I[k + 1], V_t, v_l, C_0_ts, t);
				}

				// Inserts P in the set
				polygon_set->insert(S, P);

			} else {
				v->stop_segments(I[1], seed, t);
				v->stop(t);
				delete v;
				return false;
			}
		}

		// This point is reached when the the vertex v was able to propagate other the n lines of I.
		// The returned boolean determines if the polygon should propagate laterally,
		// which is the case when there is no segment capturing the point V_t on line I[0].

		for (int k = 1; k < n; k++) v->extend_segments(I[k]);

		// If necessary, computes more events for v.
		v->reschedule_events();

		std::list<Constraint> C_limits_0;
		get_stop_constraints_at_multiple_intersections(C, 0, true, C_limits_0);

		polygon_set->get_signature_of_adjacent_cell(S, I[0]);
		return propagation_continues_at_intersection(I[0], nullptr, false, V_t, W_t, C_limits_0, t, S, seed);
	}



	void Support_Plane::constrained_vertex_intersects_multiple_lines_os(Polygon_Vertex_R* v, const FT t, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t,
		const std::vector<Intersection_Line*> & I, const std::vector<std::pair<Constraint, Constraint> > & C)
	{
		int n = int(I.size());
		const Constraint & C_0_ts = C[0].first, C_0_os = C[0].second;
		const Constraint & C_1_ts = C[1].first, C_1_os = C[1].second;

		// Initializes the signature with the address of the cell located on the other side of I[0].

		std::vector<bool> S = v->get_polygon()->get_cell()->get_signature();
		polygon_set->get_signature_of_adjacent_cell(S, I[0]);

		// Initializes a Polygon_Tree, centered in V_t.
		// This polygon has the same shape as the one whose propagation resulted in the creation of v.

		const int seed = v->get_polygon()->seed;
		CGAL_Vector_2 OV_t = V_t - polygon_directions[seed]->get_barycenter();
		Polygon_Tree* T_R = new Polygon_Tree(v->id_plane, seed, polygon_directions[seed], OV_t, t);

		T_R->split(I[0], t, 0, false);
		Polygon_Tree* T = (C_0_ts.second == PLUS ? T_R->subtree_minus : T_R->subtree_plus);
		std::list<Polygon*> L_P;

		for (int k = 0; k < n - 1; k++) {

			// Splits the subtree T with line I[k + 1]. We get two subtrees :
			// - one corresponds to the polygon that we insert in the cell that represents the fraction of plane (I[k - 1], I[k])
			// - another is for iterating.

			const Constraint & C_k_os = C[k + 1].second;
			T->split(I[k + 1], t, 0, false);

			Polygon_Tree* T_poly = (C_k_os.second == PLUS ? T->subtree_plus : T->subtree_minus);
			Polygon_Tree* T_iter = (C_k_os.second == PLUS ? T->subtree_minus : T->subtree_plus);

			// Gets the polygon that corresponds to the current iteration
			// Inserts it at address S in the set of polygons

			Polygon* P = T_poly->remove_reference_to_polygon();
			polygon_set->insert(S, P);
			L_P.push_back(P);

			// Initialiazes segments
			// Note : I[0] can be used as init condition for a segment, except when the segment is on I[0] itself

			Polygon_Vertex_R* v_k_curr, *v_k_next;
			P->get_active_constrained_vertices(I[k], I[k + 1], v_k_curr, v_k_next);
			P->shift_vertices(v_k_curr, v_k_next, V_t, t);

			init_1_uplets_of_unidirectional_segments_of_null_length(I[k], V_t, v_k_curr, (k == 0 ? C_1_os : C_0_os), t);
			init_1_uplets_of_unidirectional_segments_of_null_length(I[k + 1], V_t, v_k_next, C_0_os, t);

			// Iterates on S
			polygon_set->get_signature_of_adjacent_cell(S, I[k + 1]);

			// If there is no segment including V_t on halfline I[k + 1], iterates on T
			// Otherwise, exits loops

			std::list<Constraint> C_limits;
			get_stop_constraints_at_multiple_intersections(C, k + 1, false, C_limits);

			if (propagation_continues_at_intersection(I[k + 1], nullptr, false, V_t, W_t, C_limits, t, S, seed)) {
				T = T_iter;

			} else {
				// Shifts remaining vertices of the tree and deletes its root
				std::list<Intersection_Line*> L_I;
				std::copy(I.begin(), I.end(), std::back_inserter(L_I));
				Polygon::shift_remaining_vertices_and_schedule_events(L_P, L_I, V_t, t);

				delete T_R;
				return;
			}
		}

		// Reaching this points means that there was no segment on halflines I_1, I_2, ... I_n
		// There is one last polygon to build : the one located in the fraction of plane (I[n], I[0])

		Polygon* P = T->remove_reference_to_polygon();
		polygon_set->insert(S, P);
		L_P.push_back(P);

		Polygon_Vertex_R* v_n, *v_0;
		P->get_active_constrained_vertices(I[n - 1], I[0], v_n, v_0);
		P->shift_vertices(v_n, v_0, V_t, t);

		init_1_uplets_of_unidirectional_segments_of_null_length(I[n - 1], V_t, v_n, C_0_os, t);
		init_1_uplets_of_unidirectional_segments_of_null_length(I[0], V_t, v_0, C_1_ts, t);

		// Shifts remaining vertices of the tree and deletes its root
		std::list<Intersection_Line*> L_I;
		std::copy(I.begin(), I.end(), std::back_inserter(L_I));
		Polygon::shift_remaining_vertices_and_schedule_events(L_P, L_I, V_t, t);

		delete T_R;
	}



	void Support_Plane::get_sorted_lines(Polygon_Vertex_R* v, const CGAL_Point_2 & V_t, const std::list<Event_Vertex_Line*> & E, std::vector<Intersection_Line*> & I)
	{
		// At time t, a vertex v, constrained by I_0, intersects a set of lines I = {I_1, I_2, ... I_n}.
		// We first need to sort all the lines by angle made with the line I_0, to determine in what order polygons should be added to the support plane.

		const Constraint C_0 = v->get_constraint();
		Intersection_Line* I_0 = C_0.first;
		Sign epsilon = C_0.second;

		// We consider a parallel line to I_0, on the same side as v. Let's call it J.
		// We compute all intersections of J with lines of I, and we sort these intersections.

		const CGAL_Line_2 & L_0 = I_0->line;
		const CGAL_Vector_2 n_0 (L_0.a(), L_0.b());

		CGAL_Point_2 M_ref = V_t + n_0;
		if (I_0->sign(M_ref) != epsilon) {
			M_ref = V_t - n_0;
		}
		const CGAL_Vector_2 & dM_ref = v->dM;

		typedef std::pair<FT, Intersection_Line*> Abscissa;
		std::vector<Abscissa> A;

		struct _Abscissa_Comparator {
			bool operator() (const Abscissa & L, const Abscissa & R) {
				return L.first < R.first;
			}
		} Abscissa_Comparator;

		for (std::list<Event_Vertex_Line*>::const_iterator it_e = E.begin(); it_e != E.end(); it_e++) {
			Intersection_Line* I_k = get_line_by_identifier((*it_e)->intersected);
			const FT & a_k = I_k->a(), &b_k = I_k->b(), &c_k = I_k->c();

			FT t_k = -(a_k * M_ref.x() + b_k * M_ref.y() + c_k) / (a_k * dM_ref.x() + b_k * dM_ref.y());
			A.push_back(std::make_pair(t_k, I_k));
		}
		
		std::sort(A.begin(), A.end(), Abscissa_Comparator);

		// Returns the result.

		I.reserve(1 + E.size());
		I.push_back(I_0);
		for (int k = 0; k < int(A.size()); k++) {
			Intersection_Line* I_k = A[k].second;
			I.push_back(I_k);
		}
	}


#if 0
	void Support_Plane::get_sorted_lines(Polygon_Vertex_R* v, const std::list<Event_Vertex_Line*> & E, std::vector<Intersection_Line*> & I)
	{
		// WARNING / Here we need to switch to doubles

		// At time t, a vertex v, constrained by I_0, intersects a set of lines I = {I_1, I_2, ... I_n}.
		// We first need to sort all the lines by angle made with the line I_0, to determine in what order polygons should be added to the support plane.

		const Constraint C_0 = v->get_constraint();
		Intersection_Line* I_0 = C_0.first;
		Sign epsilon = C_0.second;

		double theta_inf, theta_sup;

		const FT a = I_0->line.a(), b = I_0->line.b();
		const CGAL_Vector_2 n_0 = (epsilon == PLUS ? CGAL_Vector_2(a, b) : CGAL_Vector_2(-a, -b));
		const double alpha = atan2(to_double(n_0.y()), to_double(n_0.x()));

		theta_inf = alpha - PI / 2.0;
		theta_sup = alpha + PI / 2.0;

		// Step 1.
		// We consider the interval [theta_inf, theta_sup]
		// We are going to consider values of theta_k that represent the orientations of the direction vectors of lines I_k.
		// This way we build couples of (theta_k, I_k) that we sort on theta_k

		typedef std::pair<double, Intersection_Line*> Orientation;
		std::vector<Orientation> J;
		J.reserve(E.size());

		for (std::list<Event_Vertex_Line*>::const_iterator it_e = E.begin(); it_e != E.end(); it_e++) {
			Intersection_Line* I_k = get_line_by_identifier((*it_e)->intersected);
			const CGAL_Vector_2 u_k = I_k->line.to_vector();
			double theta_k = atan2(to_double(u_k.y()), to_double(u_k.x()));

			// Sets theta_k in [theta_0, theta_0 + PI]
			if (theta_k < theta_inf) {
				int p = int(floor((theta_sup - theta_k) / PI));
				theta_k += p * PI;
			} else if (theta_k > theta_sup) {
				int p = int(floor((theta_k - theta_inf) / PI));
				theta_k -= p * PI;
			}
			assert(theta_inf <= theta_k && theta_k <= theta_sup);

			J.push_back(std::make_pair(theta_k, I_k));
		}

		// Step 2.
		// Now, sorts elements of J
		// To this end we design an adapted functor, that sorts elements by increasing or decreasing orientation
		// We sort elements by increasing (resp. decreasing) order of the vertex v goes in the direction of theta_sup (resp. theta_inf)

		struct _Orientation_Increasing_Order {
			bool operator() (const Orientation & OL, const Orientation & OR) {
				return OL.first < OR.first;
			}
		} Orientation_Increasing_Order;

		struct _Orientation_Decreasing_Order {
			bool operator() (const Orientation & OL, const Orientation & OR) {
				return OL.first > OR.first;
			}
		} Orientation_Decreasing_Order;

		// Computes the direction of v, expressed as an angle phi in [-1/2 pi, 3/2 pi]
		// Indeed that's the interval to which theta_sup belongs

		const CGAL_Vector_2 & dM = v->dM;
		double phi = atan2(to_double(dM.y()), to_double(dM.x()));
		if (phi < -PI / 2) phi += 2 * PI;

		if (fabs(theta_sup - phi) < 1e-6) {
			std::sort(J.begin(), J.end(), Orientation_Increasing_Order);
		} else {
			std::sort(J.begin(), J.end(), Orientation_Decreasing_Order);
		}

		// Step 3.
		// Finally gets elements of I

		I.reserve(1 + E.size());
		I.push_back(I_0);
		for (int k = 0; k < int(J.size()); k++) {
			Intersection_Line* I_k = J[k].second;
			I.push_back(I_k);
		}
	}
#endif


	void Support_Plane::get_sorted_constraints(Polygon_Vertex_R* v, const std::vector<Intersection_Line*> & I, std::vector<std::pair<Constraint, Constraint> > & C)
	{
		C.reserve(I.size());

		// The first entry of C represents the line I_0, by which v is constrained.
		// We insert a pair of constraints according to the usual order (first = this side, second = other side).

		Constraint C_0_ts = v->get_constraint();
		Constraint C_0_os = Constraint(C_0_ts.first, C_0_ts.second == PLUS ? MINUS : PLUS);
		C.push_back(std::make_pair(C_0_ts, C_0_os));

		for (int k = 1; k < I.size(); k++) {

			// For other entries of C, we compute the signs of constraints according to the location of v->M.
			Intersection_Line* I_k = I[k];
			Sign epsilon_ts = I_k->sign(v, 0);

			Constraint C_k_ts = Constraint(I_k, epsilon_ts);
			Constraint C_k_os = Constraint(I_k, epsilon_ts == PLUS ? MINUS : PLUS);

			C.push_back(std::make_pair(C_k_ts, C_k_os));
		}
	}



	void Support_Plane::get_definitions_of_constrained_vertices(Polygon_Edge* e, const FT t, const std::vector<Intersection_Line*> & I, std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & D)
	{
		D.reserve(I.size());

		// Inserts a dummy entry for line I[0]
		// After that, computes the initial point and the direction of each constrained vertex generated by intersection of I[k] and e

		D.push_back(std::make_pair(CGAL_Point_2(0, 0), CGAL_Vector_2(0, 0)));

		for (int k = 1; k < int(I.size()); k++) {
			CGAL_Point_2 M;
			CGAL_Vector_2 dM;
			e->intersection_pt_dir(I[k], t, M, dM);
			D.push_back(std::make_pair(M, dM));
		}
	}



	void Support_Plane::get_stop_constraints_at_multiple_intersections(const std::vector<std::pair<Constraint, Constraint> > & C, const int k, const bool H_ts, std::list<Constraint> & C_k)
	{
		// As a vertex v intersects the line I[k] in V_t, it may be stopped by a segment supported by line I[k]
		// If this segment starts or stops in V_t, then we define the conditions at limits that would make v part of this segment
		// TODO : illustration

		const Constraint & C_0_ts = C[0].first, &C_0_os = C[0].second;
		if (k != 0) {
			C_k.push_back(H_ts ? C_0_ts : C_0_os);
		}

		for (int l = 1; l < k; l++) {
			const Constraint & C_ts = C[l].first, &C_os = C[l].second;
			C_k.push_back(H_ts ? C_os : C_ts);
		}

		for (int l = k + 1; l < C.size(); l++) {
			const Constraint & C_ts = C[l].first, &C_os = C[l].second;
			C_k.push_back(H_ts ? C_ts : C_os);
		}
	}



	//
	// HANDLERS CASE D1
	//



	bool Support_Plane::edge_is_propagating_frontally(Intersection_Line* I, Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, const FT & t,
		bool lookup, const Constraint & C_1, const Constraint & C_2,
		const CGAL_Point_2 & V_1, const CGAL_Point_3 & W_1, const CGAL_Point_2 & V_2, const CGAL_Point_3 & W_2)
	{
		// An edge e = (v1 v2), intersecting I at time t,
		// crosses the line I if v1 or v2 can propagate

		Polygon* P = v1->get_polygon();
		Signature S = P->get_cell()->get_signature();
		const int seed = P->seed;

		//const Constraint C_1 = (v1->unconstrained() ? Constraint() : v1->get_constraint());
		//const Constraint C_2 = (v2->unconstrained() ? Constraint() : v2->get_constraint());

		Signature S_os = polygon_set->get_adjacent_polygons_signature(P, I);

		bool v1_keeps_propagating, v2_keeps_propagating;
		if (C_1.first == nullptr) {
			v1_keeps_propagating = propagation_continues_outside_intersections(I, v1, lookup, V_1, W_1, t, S_os, seed);
		} else {
			v1_keeps_propagating = propagation_continues_at_intersection(I, v1, lookup, V_1, W_1, C_1, t, S_os, seed);
		}

		if (v1_keeps_propagating) return true;

		if (C_2.first == nullptr) {
			v2_keeps_propagating = propagation_continues_outside_intersections(I, v2, lookup, V_2, W_2, t, S_os, seed);
		} else {
			v2_keeps_propagating = propagation_continues_at_intersection(I, v2, lookup, V_2, W_2, C_2, t, S_os, seed);
		}

		return v2_keeps_propagating;
		//return (v1_keeps_propagating || v2_keeps_propagating);
	}



	bool Support_Plane::edge_is_propagating_laterally(Intersection_Line* I, Intersection_Line* I_l, Polygon_Vertex_R* v, const FT & t,
		const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t)
	{
		// At time t, a constrained vertex v moving on line I_l,
		// intersects the line I and induces a collision between an edge and I.
		// We want to determine if the propagation should also be performed laterally from V_t = v->pt(t).

		Polygon* P = v->get_polygon();
		Signature S_lat = polygon_set->get_adjacent_polygons_signature(P, I, I_l);
		Constraint C_I_os = Constraint(I, I->sign(v, 0) == PLUS ? MINUS : PLUS);

		// We are going to determine if the polygon must or must not propagate on the other sides of I and I_l.
		// - If the answer is 'must not', then we're done, and we return false for saying that no, 
		// this edge doesn't have to propagate laterally.
		// - If the answer is 'must', then we should make sure that somehow a polygon will propagate in this area.

		bool must_continue = propagation_continues_at_intersection(I_l, nullptr, false, V_t, W_t, C_I_os, t, S_lat, P->seed);
		if (!must_continue) return false;

		// So from now on, we know we have to propagate laterally.
		// Two solutions :
		// - Either we do it now, as we process the current event,
		// - Or another event in an adjacent cell will cause the propagation to happen.

		// If there is no polygon, or if no vertex v_os is about to intersect I : propagate now.

		Polygon* P_os = polygon_set->get_adjacent_polygon(P, I_l);
		if (P_os == nullptr) return true;

		Polygon_Vertex_R* v_os = P_os->get_tied_vertex(I_l, v, t, V_t);
		if (v_os == nullptr) return true;

		// v_os may be inactive, this means that there was an event (v_os intersects I) which has been processed
		// and which obviously didn't result in the creation of a polygon beyond I.

		// v_os is constrained by C_l_os
		// Necessarily, one of its edges is also constrained by C_l_os : we take the other one
		// and check if it is parallel with line I. If not, we propagate

		Polygon_Edge* e = (v_os->e1->is_constrained_by(I_l) ? v_os->e2 : v_os->e1);

		const CGAL_Point_2 V_1 = e->v1->pt(t), V_2 = e->v2->pt(t);
		const CGAL_Point_3 W_1 = backproject(V_1), W_2 = backproject(V_2);

		const CGAL_Vector_2 V_12 = V_2 - V_1;
		if (CGAL::determinant(V_12, I->line.to_vector()) != 0) return true;

		// If the edge is parallel, we determine if it propagates frontally
		// If so, then no need to propagate laterally because this event will be later processed

		Polygon_Vertex_R* e_v1_r = e->v1->to_r();
		Polygon_Vertex_R* e_v2_r = e->v2->to_r();
		Constraint C_e_v1 = (e_v1_r->unconstrained() ? Constraint() : e_v1_r->get_constraint());
		Constraint C_e_v2 = (e_v2_r->unconstrained() ? Constraint() : e_v2_r->get_constraint());

		return (e_v1_r == nullptr || e_v2_r == nullptr 
			|| !edge_is_propagating_frontally(I, e_v1_r, e_v2_r, t, true, C_e_v1, C_e_v2, V_1, W_1, V_2, W_2));
	}



	void Support_Plane::edge_intersects_line(const std::list<Event_Vertex_Line*> & E_1, const std::list<Event_Vertex_Line*> & E_2, 
		Polygon_Vertex_R* v1, Polygon_Vertex_R* v2)
	{
#ifdef PRINT_EVENTS
		append(E_1.front(), "D");
#endif

		/*
		// Illustration :
		// 
		//           I_1
		//            |                      
		//   (Q31)    |        (Q2)            
		//            |                      
		//   --------+---------------------- I
		//            |                     
		//   (Q41)    |  ^                ^ 
		//            |  |               /
		//            |  v1 ---------- v2 
		//            |  |              \ 
		//            |  |       P       \ 
		//            |  |                \ 
		//            |  
		//
		*/

		int E_size_1 = int(E_1.size());
		int E_size_2 = int(E_2.size());

		Event_Vertex_Line* e_vl = E_1.front();
		Intersection_Line* I = get_line_by_identifier(e_vl->intersected);
		const FT t = e_vl->t_intersectant;

		Constraint C_1, C_2;
		CGAL_Point_2 V_1, V_2;
		
		get_constraint_and_location_at_time_t(v1, I, t, C_1, V_1);
		get_constraint_and_location_at_time_t(v2, I, t, C_2, V_2);
		
		// First of all, pops events related to v2
		// We discard the case when v1 or v2 intersect more than one line at a time.

		if (E_size_1 > 1 || E_size_2 > 1) {
			std::cout << v1->id_object << " " << v2->id_object << " " << v1->pt(t) << " " << v2->pt(t) << std::endl;
			throw std::logic_error("Error : an edge intersects more than one line at time t : unhandled case.");
		}

		/*Event_Vertex_Line* e_vl_2;
		std::list<Event_Vertex_Line*> E_2;

		get_simultaneous_events_as_edge_intersects_line(v2, t, V_2, I, E_2, e_vl_2);
		v2->decrement_queued_events(int(E_2.size()));

		if (e_vl_2 == nullptr) {
			std::cout << v1->id_object << " " << v2->id_object << " " << v1->pt(t) << " " << v2->pt(t) << std::endl;
			throw std::logic_error("Error : cannot find simultaneous event");
		} else if (E_1.size() > 1 || E_2.size() > 1) {
			std::cout << v1->id_object << " " << v2->id_object << " " << v1->pt(t) << " " << v2->pt(t) << std::endl;
			throw std::logic_error("Error : an edge intersects more than one line at time t : unhandled case.");
		}*/

		v1->delete_events(E_1);
		v2->delete_events(E_2);
		v1->decrement_queued_events(E_size_1);
		v2->decrement_queued_events(E_size_2);

		// Step 1.
		// At time t, the Polygon_Edge e = (v1 v2) intersects the line I.
		// It can propagate in three possible directions : frontally, but also laterally,
		// which is the case if v1 and v2 are constrained.
		// Over a first phase, we determine the different propagation directions.

		const CGAL_Point_3 W_1 = backproject(V_1);
		const CGAL_Point_3 W_2 = backproject(V_2);

		bool propagate_frontally = edge_is_propagating_frontally(I, v1, v2, t, false, C_1, C_2, V_1, W_1, V_2, W_2);
		bool propagate_laterally_1 = false, propagate_laterally_2 = false;

		if (propagate_frontally) {
			if (!v1->unconstrained()) {
				Intersection_Line* I_1 = v1->get_constraint().first;
				propagate_laterally_1 = edge_is_propagating_laterally(I, I_1, v1, t, V_1, W_1);
			}

			if (!v2->unconstrained()) {
				Intersection_Line* I_2 = v2->get_constraint().first;
				propagate_laterally_2 = edge_is_propagating_laterally(I, I_2, v2, t, V_2, W_2);
			}
		}

		// Step 2.
		// We perform all the necessary modifications to P,
		// and propagate, if needed, the edge (v1 v2) on the other side of I.
		// We get two booleans that indicate the necessity to propagate the polygon on the other sides 
		// of lines I_1 and I_2, that are the possible constraints of v1 and v2

		edge_propagates_frontally(I, v1, v2, t, V_1, W_1, V_2, W_2, propagate_frontally);

		// Step 3
		// If required, propagates polygons laterally

		if (propagate_laterally_1) {
			edge_propagates_laterally(I, v1, t, V_1);
		}

		if (propagate_laterally_2) {
			edge_propagates_laterally(I, v2, t, V_2);
		}
}



	void Support_Plane::edge_propagates_frontally(Intersection_Line* I, Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, const FT & t,
		const CGAL_Point_2 & V_1, const CGAL_Point_3 & W_1, const CGAL_Point_2 & V_2, const CGAL_Point_3 & W_2, const bool propagate_frontally)
	{
		// Replaces the vertices v1 and v2 by constrained or still vertices inside their polygon

		Polygon_Edge* e;
		Polygon_Vertex *v1_ts, *v1_os, *v2_ts, *v2_os;

		Constraint C_1, C_2, C_ts;

		Polygon* P_ts = v1->get_polygon();
		const int seed = P_ts->seed;
		Signature S_ts = v1->get_polygon()->get_cell()->get_signature();
		P_ts->remove_intersectant_edge(I, v1, v2, t, V_1, V_2, e, v1_ts, v2_ts, C_1, C_2, C_ts);

		if (propagate_frontally) {

			// Builds another polygon
			Polygon* P_os = Polygon::build_as_edge_propagates(v1, v2, v1_ts, v2_ts, e, t, seed, C_ts, v1_os, v2_os);
			Signature S_os = polygon_set->get_adjacent_polygons_signature(P_ts, I);
			polygon_set->insert(S_os, P_os);

			// Creates segments on each side of I
			// The number of segments depends on the active/inactive state of v1_ts, v1_os, v2_ts and v2_os
			build_n_uplets_of_segments_as_edge_collides_line(I, C_1, C_2, V_1, W_1, V_2, W_2, v1_ts, v1_os, v2_ts, v2_os, t);

			if (!v1->unconstrained()) v1->extend_segments(I);
			if (!v2->unconstrained()) v2->extend_segments(I);

			// If necessary, computes more events for v1 and v2.
			v1->reschedule_events(); ++v1->crossed_lines;
			v2->reschedule_events(); ++v2->crossed_lines;

		} else {
			// If v1 and v2 are not allowed to cross I, then the segments they may carry should be stopped.
			if (!v1->unconstrained()) v1->stop_segments(I, seed, t);
			if (!v2->unconstrained()) v2->stop_segments(I, seed, t);

			// Create segments that correspond to [v1_ts(t), v2_ts(t)]
			build_n_uplets_of_segments_as_edge_collides_line(I, C_1, C_2, V_1, W_1, V_2, W_2, v1_ts, nullptr, v2_ts, nullptr, t);

			// Deletes the edge and its vertices
			delete e;
			v1->stop(t);
			v2->stop(t);
			delete v1;
			delete v2;
		}

		if (P_ts->running_vertices == 0) {
			P_ts->forget();
		}
	}



	void Support_Plane::edge_propagates_laterally(Intersection_Line* I, Polygon_Vertex_R* v, const FT & t, const CGAL_Point_2 & V_t)
	{
		// We assume that the constrained vertex v is part of an edge, that has propagated frontally.
		// It is now in the quarter of plane (Q2).

		// Now, it's all about propagating the polygon in (Q3) delimited by constraints (C_v_os, C_os),
		// and maybe in (Q4) delimited by (C_v_os, C_ts).

		const Constraint C_v_ts = v->get_constraint();
		Intersection_Line* I_v = C_v_ts.first;
		const Constraint C_v_os = Constraint(I_v, C_v_ts.second == PLUS ? MINUS : PLUS);

		const Constraint C_ts = Constraint(I, I->sign(v, 0));
		const Constraint C_os = Constraint(I, C_ts.second == PLUS ? MINUS : PLUS);

		// Introduction
		// Initializes a Polygon_Tree, centered in V_t.

		const int seed = v->get_polygon()->seed;
		CGAL_Vector_2 OV_t = V_t - polygon_directions[seed]->get_barycenter();
		Polygon_Tree* T_R = new Polygon_Tree(v->id_plane, seed, polygon_directions[seed], OV_t, t);

		T_R->split(I_v, t, 0, false);
		Polygon_Tree* T = (C_v_ts.second == PLUS ? T_R->subtree_minus : T_R->subtree_plus);
		std::list<Polygon*> L_P;

		// Part 1
		// Splits T to obtain a polygon in (Q3)

		T->split(I, t, 0, false);
		Polygon_Tree* T_3 = (C_os.second == PLUS ? T->subtree_plus : T->subtree_minus);
		Polygon_Tree* T_4 = (C_os.second == PLUS ? T->subtree_minus : T->subtree_plus);

		// Builds P_3 and inserts it
		Signature S = polygon_set->get_adjacent_polygons_signature(v->get_polygon(), I_v);
		Polygon* P_3 = T_3->remove_reference_to_polygon();
		polygon_set->insert(S, P_3);
		L_P.push_back(P_3);

		// Initializes segments

		Polygon_Vertex_R *v3_I, *v3_I_v;
		P_3->get_active_constrained_vertices(I, I_v, v3_I, v3_I_v);
		P_3->shift_vertices(v3_I, v3_I_v, V_t, t);

		init_1_uplets_of_unidirectional_segments_of_null_length(I_v, V_t, v3_I_v, C_os, t);
		init_1_uplets_of_unidirectional_segments_of_null_length(I, V_t, v3_I, C_v_os, t);

		// Part 2
		// Determines if we should create a polygon in (Q4) and if so, inserts it

		// So, we already know that the edge propagates laterally because there is no edge from a polygon in (Q4)
		// that is going to hit I at time t. Now, the only condition for instanciating a polygon in (Q4)
		// is that, actually, no polygon of the same seed exists.

		polygon_set->get_signature_of_adjacent_cell(S, I);
		Polygon* P_4;

		if (polygon_set->exists(S, seed)) {
			P_4 = nullptr;
		} else {
			// Gets the polygon
			P_4 = T_4->remove_reference_to_polygon();
			polygon_set->insert(S, P_4);
			L_P.push_back(P_4);

			// Initializes segments
			Polygon_Vertex_R *v4_I, *v4_I_v;
			P_4->get_active_constrained_vertices(I, I_v, v4_I, v4_I_v);
			P_4->shift_vertices(v4_I, v4_I_v, V_t, t);

			if (v3_I->id_object < v4_I->id_object) {
				Polygon_Vertex_R::set_paired_vertices(v3_I, v4_I);
			} else {
				Polygon_Vertex_R::set_paired_vertices(v4_I, v3_I);
			}

			init_1_uplets_of_unidirectional_segments_of_null_length(I_v, V_t, v4_I_v, C_ts, t);
			init_1_uplets_of_unidirectional_segments_of_null_length(I, V_t, v4_I, C_v_os, t);
		}

		// Eventually
		// Shifts remaining vertices of the tree and deletes its root

		Polygon::shift_remaining_vertices_and_schedule_events(P_3, P_4, I, I_v, V_t, t);
		delete T_R;
	}



	//
	// HANDLERS CASE D2
	//



	void Support_Plane::two_edges_intersect_two_lines(const std::list<Event_Vertex_Line*> & E_1, const std::list<Event_Vertex_Line*> & E, const std::list<Event_Vertex_Line*> & E_2,
		Polygon_Vertex_R* v1, Polygon_Vertex_R* v, Polygon_Vertex_R* v2)
	{
		// As this function is called, we simultaneously process events of lists E_1, E and E_2.
		// Three consecutive vertices of the same subpolygon intersect two lines of the plane.
		// v is unconstrained, but v1 and v2 can be constrained.
		// 
		//                  I_1               I_v2
		//                  /                   |  
		//         (Q3)    /       (Q2)         |  (Q2a_os)
		//                /                     |
		// I_2 --------- + -------------------- + -----------
		//              / ^             ^       |
		//      (Q1)   /   \            |       |  (Q2a_ts)
		//            /     v ------------- v2  |
		//           /     /    e2          |   |
		//          /     / e1              |   |
		//         /     /                  |   |
		//        /  <- v1        (Q0)     ...  |
		//       /       \                      |
		//     ...       ...
		//

		// Here we strictly follow the represented geometric configuration,
		// by assuming that |E_1| = |E_2| = 1 and |E| = 2. 
		// If it is not the case, then we raise an exception.

		int E_size_1 = int(E_1.size()), E_size_2 = int(E_2.size()), E_size = int(E.size());

		if (E_size_1 != 1 || E_size_2 != 1 || E_size != 2) {
			throw std::logic_error("Error : three vertices intersect more than two lines : not implemented yet.");
		}

		// Part 1
		// Preparation of the update of the kinetic data structure.

		// Step 1.1
		// Exhibits I_1 and I_2, and the sides of these lines where v is before intersecting them.

		FT t = E.front()->t_intersectant;
		Intersection_Line* I_1 = get_line_by_identifier(E_1.front()->intersected);
		Intersection_Line* I_2 = get_line_by_identifier(E_2.front()->intersected);

		Sign eps_1 = I_1->sign(v, 0), eps_2 = I_2->sign(v, 0);
		Constraint C_1 = std::make_pair(I_1, eps_1);
		Constraint C_2 = std::make_pair(I_2, eps_2);

		v->delete_events(E);
		v1->delete_events(E_1);
		v2->delete_events(E_2);

		v->decrement_queued_events(E_size);
		v1->decrement_queued_events(E_size_1);
		v2->decrement_queued_events(E_size_2);

		// Step 1.2
		// Gets the 2D and 3D locations of the simultaneous events.

		CGAL_Point_2 V_t = (Universe::params->use_landmarks ? get_landmark(I_1, I_2) : get_intersection_point(I_1, I_2));

		Constraint C_v1, C_v2;
		CGAL_Point_2 V1_t, V2_t;
		get_constraint_and_location_at_time_t(v1, I_1, t, C_v1, V1_t);
		get_constraint_and_location_at_time_t(v2, I_2, t, C_v2, V2_t);

		CGAL_Point_3 W_t, W1_t, W2_t;
		if (Universe::params->stopping_condition != 0) {
			W_t = backproject(V_t);
			W1_t = backproject(V1_t);
			W2_t = backproject(V2_t);
		} else {
			W_t = W1_t = W2_t = CGAL::ORIGIN;
		}
		
		// Step 1.3
		// Identifies regions where the polygon should propagate (7 regions, potentially).

		bool propagates_frontally_1 = edge_is_propagating_frontally(I_1, v, v1, t, false, C_2, C_v1, V_t, W_t, V1_t, W1_t);
		bool propagates_frontally_2 = edge_is_propagating_frontally(I_2, v, v2, t, false, C_1, C_v2, V_t, W_t, V2_t, W2_t);

		bool propagates_laterally_1 = false, propagates_laterally_2 = false;
		bool propagates_diagonally = false;

		if (propagates_frontally_1) {
			if (C_v1.first != nullptr) {
				Intersection_Line* I_v1 = C_v1.first;
				propagates_laterally_1 = edge_is_propagating_laterally(I_1, I_v1, v1, t, V1_t, W1_t);
			}
		}

		if (propagates_frontally_2) {
			if (C_v2.first != nullptr) {
				Intersection_Line* I_v2 = C_v2.first;
				propagates_laterally_2 = edge_is_propagating_laterally(I_2, I_v2, v2, t, V2_t, W2_t);
			}
		}

		if (propagates_frontally_1 || propagates_frontally_2) {
			throw std::logic_error("Error in handler D2 : case not implemented yet.");
		}


		// Part 2.
		// Performs the update of the kinetic data structure.

		// Part 2.0
		// Updates the original polygon (region Q0).

		Polygon_Vertex_S* v_de;
		Polygon_Vertex* v1_ts, *v2_ts;

		Polygon* P_0 = v->get_polygon();
		P_0->remove_intersectant_edges(I_1, I_2, v1, v, v2, t, V1_t, V_t, V2_t, v1_ts, v_de, v2_ts);

		if (propagates_frontally_1 || propagates_frontally_2 || propagates_diagonally) {

			// Part 2.1
			// If necessary, propagates frontally the edge e1 in region Q1.
			// ...

			// Part 2.2
			// Same as before, with edge e2.
			// ...

			// Part 2.3
			// If necessary, propagates diagonally the vertex v.
			// ...

			throw std::logic_error("Error in handler D2 : case not implemented yet.");

		} else {
			build_n_uplets_of_segments_as_edge_collides_line(I_2, C_1, C_v2, V_t, W_t, V2_t, W2_t, v_de, nullptr, v2_ts, nullptr, t);
			build_n_uplets_of_segments_as_edge_collides_line(I_1, C_2, C_v1, V_t, W_t, V1_t, W1_t, v_de, nullptr, v1_ts, nullptr, t);

			v1->stop(t); v2->stop(t); v->stop(t); 
			delete v1;
			delete v2;
			delete v;
		}

		if (P_0->running_vertices == 0) {
			P_0->forget();
		}
	}   
}