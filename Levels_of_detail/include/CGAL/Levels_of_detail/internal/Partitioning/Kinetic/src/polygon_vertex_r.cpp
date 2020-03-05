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
	Polygon_Vertex_R::Polygon_Vertex_R(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const CGAL_Vector_2 & _dM, const int _K, Event_Flags flags)
		: Polygon_Vertex(_id_plane)
	{
		type = RUNNING_VERTEX;

		Support_Plane* SP = Universe::map_of_planes[id_plane];
		SP->vertices_r[id_object] = this;
		++SP->running_vertices;
		
		M = _M;
		dM = _dM;
		K = _K;

		assert(dM.x() != 0 || dM.y() != 0);

		t_init = _t_init;
		//t_stop = _t_init;

		k_th_interval = 0;
		queued_events = 0;

		// Sets the lower and upper bounds for the time tau,
		// which is the time required by this vertex to move by a distance D
		if (Universe::params->D_is_set) set_tau();

		// Adds this vertex to the table of moving objects.
		//Universe::map_of_objects[id_object] = this;
		++Universe::moving_objects;

		if (Universe::moving_objects > KP_Stats::peak_vertices) {
			KP_Stats::peak_vertices = Universe::moving_objects;
		}

		crossed_lines = 0;

		guided_segments = std::list<Polygon_Segment_R*>();

		paired_vertex = Companion();

		if (flags & SCHEDULE) schedule_events();
	}



	Polygon_Vertex_R::Polygon_Vertex_R(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const CGAL_Vector_2 & _dM, const Constraint & C, Intersection_Line* I_discarded, Polygon_Vertex_R* v_n, const int _K, Event_Flags flags)
		: Polygon_Vertex_R(_id_plane, _t_init, _M, _dM, _K, NO_SCHEDULE)
	{
		constraints.push_back(C);

		if (flags & SCHEDULE) schedule_events(I_discarded, v_n);
	}



	Polygon_Vertex_R::Polygon_Vertex_R(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const CGAL_Vector_2 & _dM, const Constraint & C, const std::list<Intersection_Line*> & I_discarded, Polygon_Vertex_R* v_n, const int _K, Event_Flags flags)
		: Polygon_Vertex_R(_id_plane, _t_init, _M, _dM, _K, NO_SCHEDULE)
	{
		constraints.push_back(C);

		if (flags & SCHEDULE) schedule_events(I_discarded, v_n);
	}



	Polygon_Vertex_R::Polygon_Vertex_R(Polygon_Vertex_R* v_ts, Event_Flags flags)
		: Polygon_Vertex_R(v_ts->id_plane, v_ts->t_init, v_ts->M, v_ts->dM, v_ts->K, NO_SCHEDULE)
	{
		// Constructs a vertex v_os which is the same vertex as v_ts,
		// except that it is on the other side of line I = C_ts.first

		Constraint C_ts = v_ts->get_constraint();
		Constraint C_os = Constraint(C_ts.first, (C_ts.second == PLUS ? MINUS : PLUS));
		constraints.push_back(C_os);

		set_paired_vertices(v_ts, this);

		// The value of the flags will not differ from the moment when v_ts has been itself created
		// So, if intersections are computed, we first copy t_border and I_border from v_ts,
		// and we duplicate the events involving v_ts.

		if (flags & SCHEDULE) schedule_events(v_ts);

		// Guided segments and owning polygon aren't known so far
	}


	Polygon_Vertex_R::~Polygon_Vertex_R()
	{
		++KP_Stats::life_expectancy_terms;
		KP_Stats::life_expectancy_lines += crossed_lines;

		++KP_Stats::schedule_events_vertices;
		KP_Stats::schedule_events_lines += int(lines_to_events.size());

		--Universe::moving_objects;

		reset_paired_vertex();

		deschedule_events();

		Support_Plane* SP = Universe::map_of_planes[id_plane];

		std::map<int, Polygon_Vertex_R*>::iterator it_v = SP->vertices_r.find(id_object);
		SP->vertices_r.erase(it_v);
		--SP->running_vertices;
	}


	void Polygon_Vertex_R::set_polygon(Polygon* _polygon)
	{
		if (_polygon != nullptr) {
			polygon = _polygon;
			++polygon->running_vertices;
		} else {
			--polygon->running_vertices;
			polygon = nullptr;
		}
	}



	void Polygon_Vertex_R::set_M(const CGAL_Point_2 & _M)
	{
		M = _M;
	}


	CGAL_Point_2 Polygon_Vertex_R::get_M() const
	{
		return M;
	}


	
	CGAL_Point_2 Polygon_Vertex_R::pt(const FT & t) const
	{
		return M + t * dM;
	}


	void Polygon_Vertex_R::stop(const FT & _t_stop)
	{
		CGAL_Vector_2 D = (_t_stop - t_init) * dM;
		KP_Stats::life_expectancy_distance += sqrt(to_double(D.squared_length()));
	}


	bool Polygon_Vertex_R::constrained_vertices_meet(const std::list<Event_Vertex_Line*> & E_1, const std::list<Event_Vertex_Line*> & E_2, 
		Polygon_Vertex_R* v1, Polygon_Vertex_R* v2)
	{
		if (v1->unconstrained() || v2->unconstrained()) return false;

		// v1 and v2 simultaneously intersect
		// iff v1, constrained by L_1, intersects (L_2, L_3 ... L_N)
		// and v2, constrained by L_2, intersects (L_1, L_3 ... L_N) 

		int L_1 = v1->get_constraint().first->id_object;
		int L_2 = v2->get_constraint().first->id_object;
		int id_1 = v1->id_object, id_2 = v2->id_object;

		std::list<int> I_1, I_2;

		for (std::list<Event_Vertex_Line*>::const_iterator it_e1 = E_1.begin() ; it_e1 != E_1.end() ; ++it_e1) {
			I_1.push_back((*it_e1)->intersected);
		}

		for (std::list<Event_Vertex_Line*>::const_iterator it_e2 = E_2.begin() ; it_e2 != E_2.end() ; ++it_e2) {
			I_2.push_back((*it_e2)->intersected);
		}

		// Evaluates condition

		return (std::find(I_1.begin(), I_1.end(), L_2) != I_1.end()) && (std::find(I_2.begin(), I_2.end(), L_1) != I_2.end());
	}


	Polygon_Vertex_R* Polygon_Vertex_R::get_constrained_neighbor(Intersection_Line* I_0, const std::list<Intersection_Line*> & I) const
	{
		// This function is called as we try to factorize events of the queue.
		// Typically, this vertex intersects a set of lines I and we try to determine
		// if it meets a neighbor, constrained by one of the lines of I, at the same time.

		// I_0 represents the possible constraint of v. 
		// If e1 or e2, the edges of v, lie on the line I_0 then we already know the vertices at their ends
		// are constrained by I_0, so it's not necessary to perform the test with elements of I.

		if (!e1->is_constrained_by(I_0)) {
			Polygon_Vertex_R* v1 = e1->other_vertex(this)->to_r();
			if (v1 != nullptr && v1->is_constrained_neighbor(I)) return v1;
		}

		if (!e2->is_constrained_by(I_0)) {
			Polygon_Vertex_R* v2 = e2->other_vertex(this)->to_r();
			if (v2 != nullptr && v2->is_constrained_neighbor(I)) return v2;
		}

		return nullptr;
	}


	Polygon_Vertex_R* Polygon_Vertex_R::get_neighbor_intersecting_identical_line(Intersection_Line* I_0, const std::list<Intersection_Line*> & I, const FT & t) const
	{
		// This function is called as we try to factorize events of the queue.
		// Typically this vertex intersects a set of lines I, and here we try to determine
		// if one of the neighbors v1 and v2 (ends of e1 and e2) intersect one of the lines in I.

		// Assumption : v1 and v2 are NOT constrained by any of the lines of I,
		// and therefore, edges e1 and e2 are not of length zero.

		if (!e1->is_constrained_by(I_0)) {
			Polygon_Vertex* v1_base = e1->other_vertex(this);
			if (Polygon_Vertex_R* v1 = dynamic_cast<Polygon_Vertex_R*>(v1_base)) {
				CGAL_Point_2 M = v1->pt(t);

				for (std::list<Intersection_Line*>::const_iterator it_l = I.begin(); it_l != I.end(); it_l++) {
					if ((*it_l)->includes(M)) return v1;
				}
			}
		}

		if (!e2->is_constrained_by(I_0)) {
			Polygon_Vertex* v2_base = e2->other_vertex(this);
			if (Polygon_Vertex_R* v2 = dynamic_cast<Polygon_Vertex_R*>(v2_base)) {
				CGAL_Point_2 M = v2->pt(t);

				for (std::list<Intersection_Line*>::const_iterator it_l = I.begin(); it_l != I.end(); it_l++) {
					if ((*it_l)->includes(M)) return v2;
				}
			}
		}

		return nullptr;
	}



	void Polygon_Vertex_R::transfer_segments(Polygon_Vertex_R* v_dest, const FT & t, const CGAL_Point_2 & V_t)
	{
		// Gets the 3D position (a) and the speed (b - a) of the vertex v_dest, which has just been created
		// const CGAL_Point_2 V_t = v_dest->pt(t);
		const CGAL_Point_2 V_u = V_t + v_dest->dM;

		const CGAL_Point_3 a = Universe::map_of_planes[v_dest->id_plane]->backproject(V_t);
		const CGAL_Point_3 b = Universe::map_of_planes[v_dest->id_plane]->backproject(V_u);

		for (std::list<Polygon_Segment_R*>::iterator it_s = guided_segments.begin(); it_s != guided_segments.end(); it_s++) {
			Polygon_Segment_R* s = (*it_s);

			Support_Plane* SP = Universe::map_of_planes[s->id_plane];

			// Updates the speed of the segment
			const CGAL_Point_2 A = SP->project(a);
			const CGAL_Point_2 B = SP->project(b);
			s->update_translation(t, A, B - A);

			// From now on, s is now going to move along with v_dest
			v_dest->set_as_guided_segment(s);
		}

		guided_segments.clear();
	}



	void Polygon_Vertex_R::stop_segments(Intersection_Line* L_ik, const int seed, const FT & t)
	{
		// Segments carried by this vertex no longer propagate after time t, when v(t) intersects a line I.
		// This vertex v belongs to plane i, and I represents the plane k : it is denoted as L_ik.
		// A segment lives on plane j, on line L_ji. We compute the its constraint C_jk.

		const int i = id_plane;
		const int k = L_ik->planes.front();

		for (std::list<Polygon_Segment_R*>::iterator it_s = guided_segments.begin(); it_s != guided_segments.end(); it_s++) {
			Polygon_Segment_R* s1 = (*it_s);

			const int j = s1->id_plane;
			Support_Plane* P_j = Universe::map_of_planes[j];

			Intersection_Line* L_jk = P_j->get_line_for_plane(k);
			Sign eps_jk = L_jk->sign(s1->origin());
			assert(eps_jk != ZERO);

			Constraint C_jk(L_jk, eps_jk);
			s1->stop(C_jk, t);

			// Once segments have stopped, the idea is to convert them into lighter Polygon_Segment_S objects.

			Polygon_Segment_R* s2 = s1->get_opposite();
			if (s2 == nullptr) {

				// Case 1.
				// Unidirectional segments.

				// We build a Polygon_Segment_S that collects all the information about s :
				// the initial, final constraints and the crossed intersection lines.

				Intersection_Line* I_supp = s1->support;
				Sign eps = s1->get_sign();

				Constraint C_supp (I_supp, eps);
				Polygon_Segment_S* s_const = new Polygon_Segment_S(s1->id_plane, seed, s1->C_init, C_supp, s1->C_crossed, s1->C_stop, true);

				// Erases and deletes s.

				std::list<Segment*> & I_supp_segments = (eps == PLUS ? I_supp->segments_plus : I_supp->segments_minus);
				I_supp_segments.erase(s1->iterator);
				delete s1;

			} else {
			
				// Case 2.
				// Bidirectional segments.

				// We may perform the same operation as before, but iff the two segments are stopped.

				if (s2->stopped()) {
					Intersection_Line* I_supp = s2->support;
					Sign eps = s1->get_sign();

					Constraint C_supp (I_supp, eps);
					Polygon_Segment_S* s_const = new Polygon_Segment_S(s1->id_plane, seed, s1->C_stop, C_supp, s1->C_crossed, s2->C_crossed, s2->C_stop);

					// Erases and deletes s and s_opposite.
					
					std::list<Segment*> & I_supp_segments = (eps == PLUS ? I_supp->segments_plus : I_supp->segments_minus);
					I_supp_segments.erase(s2->iterator);
					I_supp_segments.erase(s1->iterator);
					delete s2;
					delete s1;

				} else if (s2->C_init.first == nullptr) {
					// s_opposite still runs,
					// But the stop constraint for s may be used for delimiting it.
					s2->set_pseudo_init_constraint(C_jk);
				}
			}
		}

		guided_segments.clear();
	}



	void Polygon_Vertex_R::extend_segments(Intersection_Line* L_ik) const
	{
		// This is the opposite function to before.
		// A constraint vertex v intersects a line I but it can propagate beyond that line,
		// so all segments carried by this segment also keep propagating.
		// This vertex v belongs to plane i, and I represents the plane k : it is denoted as L_ik.

		const int i = id_plane;

		for (std::list<Polygon_Segment_R*>::const_iterator it_s = guided_segments.begin(); it_s != guided_segments.end(); it_s++) {
			Polygon_Segment_R* s = (*it_s);

			const int j = s->id_plane;
			Support_Plane* P_j = Universe::map_of_planes[j];

			// We always insert the line that corresponds to the plane k in the list of lines intersected by s

			std::set<Intersection_Line*> _CL;
			for (std::list<int>::iterator it_p = L_ik->planes.begin(); it_p != L_ik->planes.end(); it_p++) {
				int k = L_ik->planes.front();
				Intersection_Line* L_jk = P_j->get_line_for_plane(k);
				_CL.insert(L_jk);
			}

			std::list<Intersection_Line*> CL(_CL.begin(), _CL.end());
			for (std::list<Intersection_Line*>::iterator it_l = CL.begin(); it_l != CL.end(); it_l++) {
				s->insert_as_crossed_line(*it_l);
			}
			Intersection_Line* L_jk = CL.front();

			// If the segment is bidirectional,
			// L_jk may be used to delimit the opposite segment and set a pseudo-init constraint.
			if (s->exists_opposite_segment_without_initinal_constraint()) {
				Polygon_Segment_R* s_opposite = s->get_opposite();
				Sign eps_jk = L_jk->sign(s_opposite->origin());
				Constraint C_jk(L_jk, eps_jk);

				s_opposite->set_pseudo_init_constraint(C_jk);
			}
		}
	}



	void Polygon_Vertex_R::indicate_line_initially_crossed_by_segments(Intersection_Line* L_ik) const
	{
		// This function is called during the initialization process, and is very similar to before, 
		// except that we don't try to set an initial constraint to the unidirectional segments
		// because such a constraint already exists.

		const int i = id_plane;

		for (std::list<Polygon_Segment_R*>::const_iterator it_s = guided_segments.begin(); it_s != guided_segments.end(); it_s++) {
			Polygon_Segment_R* s = (*it_s);

			const int j = s->id_plane;
			Support_Plane* P_j = Universe::map_of_planes[j];

			// We always insert the line that corresponds to the plane k in the list of lines intersected by s

			std::set<Intersection_Line*> _CL;
			for (std::list<int>::iterator it_p = L_ik->planes.begin(); it_p != L_ik->planes.end(); it_p++) {
				int k = L_ik->planes.front();
				Intersection_Line* L_jk = P_j->get_line_for_plane(k);
				_CL.insert(L_jk);
			}

			std::list<Intersection_Line*> CL(_CL.begin(), _CL.end());
			for (std::list<Intersection_Line*>::iterator it_l = CL.begin(); it_l != CL.end(); it_l++) {
				s->insert_as_crossed_line(*it_l);
			}
			Intersection_Line* L_jk = CL.front();
		}
	}



	void Polygon_Vertex_R::copy_crossed_lines(Polygon_Vertex_R* v_os) const
	{
		// This function is called during the initialization process.
		// We suppose that we know the lines intersected by segments 
		// carried by a vertex that is on one side of a line.

		// Now, we want to copy these informations to segments carried by the same vertex,
		// but that is located on the other side of the same line.

		std::list<Polygon_Segment_R*>::const_iterator it_s1 = guided_segments.begin();
		std::list<Polygon_Segment_R*>::iterator it_s2 = v_os->guided_segments.begin();

		while (it_s1 != guided_segments.end() && it_s2 != v_os->guided_segments.end()) {
			Polygon_Segment_R* s1 = (*it_s1);
			Polygon_Segment_R* s2 = (*it_s2);

			for (std::list<Intersection_Line*>::const_iterator it_l = s1->crossed_lines_begin(); it_l != s1->crossed_lines_end(); ++it_l) {
				s2->insert_as_crossed_line(*it_l);
			}

			++it_s1;
			++it_s2;
		}
	}



	void Polygon_Vertex_R::set_as_guided_segment(Polygon_Segment_R* s)
	{
		guided_segments.push_back(s);
	}


	void Polygon_Vertex_R::set_tau()
	{
		// Sets the lower and upper bounds for the time tau,
		// which is the time required by this vertex to move by a distance D
		// defined in the parameters.

		const FT & D_inf_2 = Universe::params->D_inf_2;
		const FT & D_sup_2 = Universe::params->D_sup_2;
		const FT & D_2 = Universe::params->D_2;

		const FT dm_sq = dM.squared_length();
		const FT tau_2 = D_2 / dm_sq;

		double prec = 1.0 / (1 << 30) / (1 << 10);
		FT::set_relative_precision_of_to_double(prec);

		std::stringstream stream;
		double f_tau_inf = sqrt(CGAL::to_double(D_inf_2 / dm_sq));
		double f_tau_sup = sqrt(CGAL::to_double(D_sup_2 / dm_sq));
		stream << f_tau_inf << " " << f_tau_sup;
		stream >> tau_inf >> tau_sup;

		assert(tau_inf * tau_inf <= tau_2 && tau_2 <= tau_sup * tau_sup);
	}


	bool Polygon_Vertex_R::has_events() const
	{
		return (queued_events != 0);
	}


	void Polygon_Vertex_R::schedule_events()
	{
		std::list<Intersection_Line*> discarded;
		schedule_events(discarded, nullptr);
	}



	void Polygon_Vertex_R::schedule_events(Intersection_Line* I_redirect, Polygon_Vertex_R* v_n)
	{
		std::list<Intersection_Line*> discarded(1, I_redirect);
		schedule_events(discarded, v_n);
	}



	void Polygon_Vertex_R::schedule_events(const std::list<Intersection_Line*> & discarded, Polygon_Vertex_R* v_n)
	{
		// Each vertex has a map that associates lines to events.
		// A null event indicates that the events has already been processed, or cannot happen.
		// For this reason, we insert associate null events to the lines that must be discarded
		// from the algorithm.

		if (!discarded.empty()) {
			for (std::list<Intersection_Line*>::const_iterator it_l = discarded.begin(); it_l != discarded.end(); it_l++) {
				Intersection_Line* I = (*it_l);
				if (I != nullptr) lines_to_events[I->id_object] = nullptr;
			}
		}

		reschedule_events(v_n);
	}



	bool Polygon_Vertex_R::meets_constrained_neighbor_at_next_event(Polygon_Vertex_R* v_n, Intersection_Line* & I_n, std::list<Event_Vertex_Line*> & E) const
	{
		// In this function, we determine if this vertex is about to meet v at the next intersection,
		// in an event of type C1, which would spare the algorithm the full scheduling of events for this vertex :
		// this is the only case when we return true.

		if (unconstrained() || v_n->unconstrained()) return false;

		Intersection_Line* I = constraints.front().first;
		v_n->get_upcoming_event(E);

		for (std::list<Event_Vertex_Line*>::const_iterator it_e = E.begin() ; it_e != E.end() ; ++it_e) {
			if ((*it_e)->intersected == I->id_object) {
				I_n = v_n->get_constraint().first;
				return true;
			}
		}

		return false;
	}



	void Polygon_Vertex_R::get_upcoming_event(std::list<Event_Vertex_Line*> & E) const
	{
		E.clear();
		if (queued_events == 0) return;

		for (std::map<int, Event_Vertex_Line*>::const_iterator it_e = lines_to_events.begin() ; it_e != lines_to_events.end() ; ++it_e) {
			Event_Vertex_Line* e_vl = it_e->second;
			if (e_vl != nullptr) {

				// If E is empty, then e_vl is the first valid event we process. It can be used as reference.
				// Otherwise, we use an invariant property according to what E always contains the event with the lowest t.
				// E contains two elements or more iff there are simultaneous events.

				if (E.empty()) {
					E.push_back(e_vl);
				} else {
					Event_Vertex_Line* e = E.front();
					if (e_vl->t_intersectant < e->t_intersectant) {
						E.clear();
						E.push_back(e_vl);
					} else if (e_vl->t_intersectant == e->t_intersectant) {
						E.push_back(e_vl);
					}
				}
			}
		}
	}



	void Polygon_Vertex_R::reschedule_events()
	{
		reschedule_events(nullptr);
	}



	void Polygon_Vertex_R::reschedule_events(Polygon_Vertex_R* v_n)
	{
		// We assume that this vertex is active and there are no events associated to this vertex in the queue.
		// If it is not the case, returns.

		if (queued_events != 0) return;

		std::list<Event_Vertex_Line*> E_VL;

		// Part 1.
		// If the vertex is constrained, then it may meet one of its neighbors at the next event popped from the queue.
		// If so, we don't want to run the expensive algorithm that searches for intersectable lines in the plane :
		// we just build events after such neighbors.

		reschedule_events_after_neighbors(v_n, E_VL);

		// Part 2.
		// If we couldn't make a schedule after one of this vertex's neighbors,
		// then we run the normal procedure.

		if (E_VL.empty()) {
			reschedule_events_by_finding_lines(E_VL);
		}

		// Inserts events in the map.
		Universe::event_queue->push(E_VL);
	}



	void Polygon_Vertex_R::reschedule_events_after_neighbors(Polygon_Vertex_R* v_n, std::list<Event_Vertex_Line*> & E_VL)
	{
		// This is a particular way of scheduling events.
		// If this vertex is constrained by a line I, 
		// and has a neighbor v_n constrained by a line I_n whose upcoming event is (v_n intersects I),
		// then we schedule only one event (v intersects I_n) : both vertices are going to meet and get deleted.
		
		if (constraints.empty()) return;

		Intersection_Line* I_ref = constraints.front().first;

		Intersection_Line* I_neighbor = nullptr;
		std::list<Event_Vertex_Line*> E_neighbor;
		
		Polygon_Vertex_R* v1 = (e1 == nullptr ? nullptr : dynamic_cast<Polygon_Vertex_R*>(e1->other_vertex(this)));
		Polygon_Vertex_R* v2 = (e2 == nullptr ? nullptr : dynamic_cast<Polygon_Vertex_R*>(e2->other_vertex(this)));

		bool meets_constrained_neighbor = ((v_n != nullptr && meets_constrained_neighbor_at_next_event(v_n, I_neighbor, E_neighbor))
			|| (v1 != nullptr && meets_constrained_neighbor_at_next_event(v1, I_neighbor, E_neighbor))
			|| (v2 != nullptr && meets_constrained_neighbor_at_next_event(v2, I_neighbor, E_neighbor)));

		if (meets_constrained_neighbor) {

			// We copy all the events of E_n, in which the constrained neighbor is involved.

			for (std::list<Event_Vertex_Line*>::iterator it_e = E_neighbor.begin() ; it_e != E_neighbor.end() ; ++it_e) {
				Event_Vertex_Line* e_ref = (*it_e);
				Event_Vertex_Line* e_vl = nullptr;

				if (e_ref->intersected == I_ref->id_object) {
					e_vl = new Event_Vertex_Line(id_plane, id_object, I_neighbor->id_object, e_ref->t_intersectant);
					lines_to_events[I_neighbor->id_object] = e_vl;
				} else {
					e_vl = new Event_Vertex_Line(id_plane, id_object, e_ref->intersected, e_ref->t_intersectant);
					lines_to_events[e_ref->intersected] = e_vl;
				}

				++queued_events;
				E_VL.push_back(e_vl);
			}
		}
	}



	void Polygon_Vertex_R::reschedule_events_by_finding_lines(std::list<Event_Vertex_Line*> & E_VL)
	{
		// We assume there are still lines that it can intersect.
		// In this function, we are going to compute a set of events occuring in I_k = [tk_min, tk_max],
		// where tk_min = t_init + k * tau, tk_max = t_init + (k + 1) * tau.

		// The variable tau is the time when the vertex moves by a distance of D (a parameter).
		// As it is impossible to compute it exactly, we use two bounds tau_inf and tau_sup
		// to set t_min and t_max.

		// Given our initial assumption, we expect to compute at least one event that we insert in E_VL.
		// We iteratively loop on intervals I_k until E_VL gets no longer empty.

		Support_Plane* SP = Universe::map_of_planes[id_plane];

		while (E_VL.empty()) {

			// We define the current interval I_k.

			FT t_min = t_init + k_th_interval * tau_inf;
			FT t_max = t_init + (k_th_interval + 1) * tau_sup;

			// We compute v->pt(t_min) and v->pt(t_max).
			// We define J, a list of lines that might have been crossed in the meanwhile.

			std::list<std::tuple<Intersection_Line*, bool, FT> > J;
			SP->search_lines_in_neighborhood(this, t_min, t_max, J);

			for (std::list<std::tuple<Intersection_Line*, bool, FT> >::iterator it_l = J.begin(); it_l != J.end(); it_l++) {
				const std::tuple<Intersection_Line*, bool, FT> & J_tuple = (*it_l);

				Intersection_Line* I = std::get<0>(J_tuple);
				bool t_already_known = std::get<1>(J_tuple);
				Event_Vertex_Line* e_vl = nullptr;

				// For each line I, we query the look-up table lines_to_events.
				// The non-existence of an entry implies that the line is not discarded
				// or that no intersection with this line has been computed before.
				// When we compute an event, we increment by anticipation the number
				// of queued events involving this vertex.

				if (lines_to_events.find(I->id_object) == lines_to_events.end()) {
					if (t_already_known) {
						e_vl = make_event(I, std::get<2>(J_tuple));
					} else {
						e_vl = make_event(I);
					}
					lines_to_events[I->id_object] = e_vl;

					if (e_vl != nullptr) {
						// e_vl can be null is the intersection takes place behind the vertex
						E_VL.push_back(e_vl);
						++queued_events;
					}
				}
			}

			// Iterates on k.
			++k_th_interval;

			// If k takes high values (multiples of 32768),
			// we check if it is still inside the bounding polygon.
			// If not it means that we are stuck in an infinite loop and that the program should abort.

			if ((k_th_interval << 17) == 0) {
				std::cout << "Strange vertex : [" << id_plane << "] " << id_object << std::endl;
				if (!SP->assert_inside(this, t_min)) {
					throw std::logic_error("Error : infinite loop detected, please rerun with option --no-fast-schedule");
				}
			}
		}
	}


	
	void Polygon_Vertex_R::decrement_queued_events(const int n)
	{
		queued_events -= n;
	}



	Event_Vertex_Line* Polygon_Vertex_R::get_event_for_line(const int I_object) const
	{
		// This function is called to access an event that has been computed before

		std::map<int, Event_Vertex_Line*>::const_iterator it_e = lines_to_events.find(I_object);
		assert(it_e != lines_to_events.end());

		return it_e->second;
	}



	void Polygon_Vertex_R::delete_event(Event_Vertex_Line* e_vl)
	{
		// This function is called to delete a specified event, typically when it is processed,
		// and set the local reference to such event to nullptr.

		std::map<int, Event_Vertex_Line*>::iterator it_e = lines_to_events.find(e_vl->intersected);
		assert(it_e != lines_to_events.end() && it_e->second == e_vl);

		delete e_vl;
		it_e->second = nullptr;
	}



	void Polygon_Vertex_R::delete_events(const std::list<Event_Vertex_Line*> & E_VL)
	{
		for (std::list<Event_Vertex_Line*>::const_iterator it_e = E_VL.begin(); it_e != E_VL.end(); it_e++) {
			delete_event(*it_e);
		}
	}



	void Polygon_Vertex_R::copy_events(Polygon_Vertex_R* v_ts)
	{
		// Copies local events that were not existing before

		if (lines_to_events.size() != v_ts->lines_to_events.size()) {
			std::map<int, Event_Vertex_Line*>::iterator it_e_ts;

			for (it_e_ts = v_ts->lines_to_events.begin(); it_e_ts != v_ts->lines_to_events.end(); it_e_ts++) {
				int I_id_object = it_e_ts->first;
				Event_Vertex_Line* e_ts = it_e_ts->second;

				// Missing event : makes a local copy if e_ts is not nullptr
				if (lines_to_events.find(I_id_object) == lines_to_events.end()) {
					if (e_ts == nullptr) {
						lines_to_events[I_id_object] = nullptr;
					} else {
						lines_to_events[I_id_object] = new Event_Vertex_Line(id_plane, id_object, I_id_object, e_ts->t_intersectant);
					}
				}
			}
		}
	}



	void Polygon_Vertex_R::schedule_events(Polygon_Vertex_R* v_ts)
	{
		// We assume that some events have been computed for the master vertex v_ts.
		// We just read the look-up table of events and make a deep copy of all vertices.

		for (std::map<int, Event_Vertex_Line*>::iterator it_e = v_ts->lines_to_events.begin(); it_e != v_ts->lines_to_events.end(); ++it_e) {
			int index = it_e->first;
			Event_Vertex_Line* e_ts = it_e->second;

			if (lines_to_events.find(index) == lines_to_events.end()) {
				if (e_ts == nullptr) {
					lines_to_events[index] = nullptr;
				} else {
					Event_Vertex_Line* e_os = new Event_Vertex_Line(id_plane, id_object, e_ts->intersected, e_ts->t_intersectant);
					lines_to_events[index] = e_os;
					Universe::event_queue->push(e_os, e_ts->queue_iterator);
				}
			}
		}

		queued_events = v_ts->queued_events;
		k_th_interval = v_ts->k_th_interval;
	}



	void Polygon_Vertex_R::deschedule_events()
	{
		// Removes events from the queue before destroying them.

		for (std::map<int, Event_Vertex_Line*>::iterator it_e = lines_to_events.begin(); it_e != lines_to_events.end(); it_e++) {
			Event_Vertex_Line* e_vl = it_e->second;
			if (e_vl != nullptr) {
				Universe::event_queue->erase(e_vl);
				delete it_e->second;
			}
		}

		lines_to_events.clear();
		k_th_interval = 0;
		queued_events = 0;
	}



	std::pair<FT, bool> Polygon_Vertex_R::get_intersection_time(Intersection_Line* I) const
	{
		// We assume that v is active.
		// We are going to determine the time t when the vertex v(t) = M + t * dM intersects I.

		// When the direction of the line is collinear to the direction of the vertex, 
		// there is no intersection and we return FLT_MAX.

		const FT & a = I->a(), b = I->b(), c = I->c();

		FT t_den = a * dM.x() + b * dM.y();
		if (t_den == 0) {
			// The direction of the vector is collinear to the line
			return std::make_pair(0, false);
		}

		FT t_num = -(a * M.x() + b * M.y() + c);

		FT t = t_num / t_den;
		return std::make_pair(t, true);
	}



	Event_Vertex_Line* Polygon_Vertex_R::make_event(Intersection_Line* I) const
	{
		std::pair<FT, bool> R = get_intersection_time(I);

		if (!R.second || R.first < t_init) {
			return nullptr;
		} else {
			return new Event_Vertex_Line(id_plane, id_object, I->id_object, R.first);
		}
	}



	Event_Vertex_Line* Polygon_Vertex_R::make_event(Intersection_Line* I, const FT & t) const
	{
		if (t < t_init) {
			return nullptr;
		} else {
			return new Event_Vertex_Line(id_plane, id_object, I->id_object, t);
		}
	}



	void Polygon_Vertex_R::set_paired_vertices(Polygon_Vertex_R* v_ts, Polygon_Vertex_R* v_os)
	{
		// We create a link between two vertices v_ts and v_os.
		// v_os is considered as a copy of v_ts, located on the other side of a line I.

		// In this function, we provide the algorithm the information according to what
		// events of v_os don't have to be computed : they are copied from v_ts as well
		// since v_ts and v_os propagate at the same speed.

		// The duplication will take place at the insertion of events (v_ts, I) in
		// the queue and matrices of events.

		v_ts->set_paired_vertex(Companion(v_os, true));  // v_ts is the master
		v_os->set_paired_vertex(Companion(v_ts, false)); // v_os is the slave
	}



	bool Polygon_Vertex_R::are_paired_vertices(Polygon_Vertex_R* v_ts, Polygon_Vertex_R* v_os)
	{
		if (!v_ts->has_paired_vertex() || !v_os->has_paired_vertex()) return false;

		return (v_ts->paired_vertex.first == v_os && v_os->paired_vertex.first == v_ts);
	}



	void Polygon_Vertex_R::set_paired_vertex(const Companion & _paired_vertex)
	{
		paired_vertex = _paired_vertex;
	}



	void Polygon_Vertex_R::reset_paired_vertex()
	{
		Polygon_Vertex_R* v_os = paired_vertex.first;
		if (v_os != nullptr) {
			v_os->set_paired_vertex(Companion());
			paired_vertex = Companion();
		}
	}



	bool Polygon_Vertex_R::has_paired_vertex() const
	{
		return (paired_vertex.first != nullptr);
	}



	bool Polygon_Vertex_R::is_independent() const
	{
		return (paired_vertex.first == nullptr || paired_vertex.second);
	}



	Polygon_Vertex_R* Polygon_Vertex_R::get_paired_vertex() const
	{
		return paired_vertex.first;
	}


	Polygon_Vertex* Polygon_Vertex_R::to_base()
	{
		return dynamic_cast<Polygon_Vertex*>(this);
	}
}