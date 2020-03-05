#include "../include/support_plane.h"
#include "../include/parameters.h"
#include "../include/intersection_line.h"
#include "../include/polygon_vertex.h"
#include "../include/segment.h"
#include "../include/universe.h"


namespace Skippy {

	void Support_Plane::init_1_uplets_of_unidirectional_segments_of_null_length(Intersection_Line* I, const CGAL_Point_2 & V_t,
		Polygon_Vertex_R* v, const Constraint C, const FT & t)
	{
		// This function performs exactly the same process as the next one
		// We call it by replacing v_os by nullptr

		init_2_uplets_of_unidirectional_segments_of_null_length(I, V_t, v, nullptr, C, t);
	}



	void Support_Plane::init_2_uplets_of_unidirectional_segments_of_null_length(Intersection_Line* I, const CGAL_Point_2 & V_t,
		Polygon_Vertex_R* v_ts, Polygon_Vertex_R* v_os, const Constraint C, const FT & t)
	{
		// We would like to initialize a 2-uplet of segments on each of the planes represented by I
		// These segments are guided by v_ts and v_os, which are at time t, at the intersection of lines I and I_1
		// We assume that v_ts and v_os move at constant speed

		assert(v_ts->get_polygon() != nullptr);
		int seed = v_ts->get_polygon()->seed;

		Sign v_ts_eps = v_ts->sign_of_constraint(I);
		Sign v_os_eps = (v_ts_eps == PLUS ? MINUS : PLUS);

		const CGAL_Point_3 b = backproject(V_t + v_ts->dM);
		//const CGAL_Point_3 & a = W_t, b = backproject(V_t + v_ts->dM);

		Intersection_Line* I_1 = C.first;
		int p_0 = v_ts->id_plane, p_1 = I_1->planes.front();

		for (std::list<int>::iterator it_p = I->planes.begin(); it_p != I->planes.end(); it_p++) {
			Support_Plane* SP = Universe::map_of_planes[*it_p];

			// Identifies the lines of interest :
			// J_0 is the line on which we are going to initialize the segment associated to s_ts
			// J_1 is the line that represents its initial constraint
			Intersection_Line *J_0 = SP->get_line_for_plane(p_0), *J_1 = SP->get_line_for_plane(p_1);
			if (J_0->reject_segments) continue;

			// Projection of the 3D coordinates onto the 2D plane SP
			CGAL_Point_2 A;
			if (Universe::params->use_landmarks) {
				A = SP->get_landmark(J_0, J_1);
			} else {
				A = SP->get_intersection_point(J_0, J_1);
			}

			CGAL_Point_2 B = SP->project(b);
			// CGAL_Point_2 A = SP->project(a), B = SP->project(b);
			CGAL_Vector_2 dA = B - A;

			// Transposes C to get the initial constraint of the segment
			Sign eps_1 = J_1->sign(B);
			assert(eps_1 != ZERO);
			Constraint C_init(J_1, eps_1);

			// Gets support constraints
			Constraint C_ts(J_0, v_ts_eps);
			Polygon_Segment_R* s_ts = new Polygon_Segment_R(SP->id, seed, t, C_init, C_ts, v_ts, A, dA);

			if (v_os != nullptr) {
				Constraint C_os(J_0, v_os_eps);
				Polygon_Segment_R* s_os = new Polygon_Segment_R(SP->id, seed, t, C_init, C_os, v_os, A, dA);
			}
		}
	}



	void Support_Plane::init_1_uplets_of_unidirectional_segments(Intersection_Line* I, const CGAL_Point_2 & V_c, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t,
		Polygon_Vertex_R* v, const Constraint C_init, const FT & t)
	{
		// This function performs exactly the same process as the next one
		// We call it by replacing vc_os and v_os by nullptr

		init_2_uplets_of_unidirectional_segments(I, V_c, V_t, W_t, v, nullptr, C_init, t);
	}



	void Support_Plane::init_2_uplets_of_unidirectional_segments(Intersection_Line* I,
		const CGAL_Point_2 & V_c, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t,
		Polygon_Vertex_R* v_ts, Polygon_Vertex_R* v_os, const Constraint _C_init, const FT & t)
	{
		// We could like to initialize 2-uplets of segments in each of the planes represented by I.
		// These segments are guided by v_ts and v_os, are are delimited by vc_ts and vc_os which are constant.

		assert(v_ts->get_polygon() != nullptr);
		int seed = v_ts->get_polygon()->seed;

		const Sign v_ts_eps = v_ts->sign_of_constraint(I);
		const Sign v_os_eps = (v_ts_eps == PLUS ? MINUS : PLUS);

		//const CGAL_Point_3 & o = W_c, &a = W_t, b = backproject(V_t + v_ts->dM);
		const CGAL_Point_3 &a = W_t, b = backproject(V_t + v_ts->dM);

		Intersection_Line* I_1 = _C_init.first;
		int p_0 = v_ts->id_plane, p_1 = I_1->planes.front();

		for (std::list<int>::iterator it_p = I->planes.begin(); it_p != I->planes.end(); it_p++) {
			Support_Plane* SP = Universe::map_of_planes[*it_p];

			// For each instanciated segment, there are three elements we need to compute :
			// the initial points O and A, the propagation speed dA.

			Intersection_Line *J_0 = SP->get_line_for_plane(p_0), *J_1 = SP->get_line_for_plane(p_1);
			if (J_0->reject_segments) continue;

			//CGAL_Point_2 O = SP->project(o), A = SP->project(a), B = SP->project(b);
			CGAL_Point_2 O;
			if (Universe::params->use_landmarks) {
				O = SP->get_landmark(J_0, J_1);			
			} else {
				O = SP->get_intersection_point(J_0, J_1);
			}

			CGAL_Point_2 A = SP->project(a), B = SP->project(b);
			CGAL_Vector_2 dA = B - A;

			// Transposes C_init to get the initial constraint of the segment
			Sign eps_1 = J_1->sign(B);
			assert(eps_1 != ZERO);
			Constraint C_init(J_1, eps_1);

			// Gets support constraints
			Constraint C_ts(J_0, v_ts_eps);
			Polygon_Segment_R* s_ts = new Polygon_Segment_R(SP->id, seed, t, C_init, C_ts, v_ts, O, A, dA);

			if (v_os != nullptr) {
				Constraint C_os(J_0, v_os_eps);
				Polygon_Segment_R* s_os = new Polygon_Segment_R(SP->id, seed, t, C_init, C_os, v_os, O, A, dA);
			}
		}
	}



	void Support_Plane::init_2_uplets_of_bidirectional_segments_of_null_length(Intersection_Line* I, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t,
		Polygon_Vertex_R* v1_ts, Polygon_Vertex_R* v2_ts, const FT & t)
	{
		init_4_uplets_of_bidirectional_segments_of_null_length(I, V_t, W_t, v1_ts, nullptr, v2_ts, nullptr, t);
	}



	void Support_Plane::init_4_uplets_of_bidirectional_segments_of_null_length(Intersection_Line* I, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t,
		Polygon_Vertex_R* v1_ts, Polygon_Vertex_R* v1_os, Polygon_Vertex_R* v2_ts, Polygon_Vertex_R* v2_os, const FT & t)
	{
		// Just like the other segment initializers, we want to add segments on each of the plane represented by I
		// Here we consider the case when the bidirectional segments are points at the time of their creation

		assert(v1_ts->get_polygon() != nullptr);
		int seed = v1_ts->get_polygon()->seed;

		Sign v_ts_eps = v1_ts->sign_of_constraint(I);
		Sign v_os_eps = (v_ts_eps == PLUS ? MINUS : PLUS);
		int p_0 = v1_ts->id_plane;

		const CGAL_Point_3 & a = W_t;
		const CGAL_Point_3 b_1 = backproject(V_t + v1_ts->dM), b_2 = backproject(V_t + v2_ts->dM);

		for (std::list<int>::iterator it_p = I->planes.begin(); it_p != I->planes.end(); it_p++) {
			Support_Plane* SP = Universe::map_of_planes[*it_p];

			// Builds a support constraint
			Intersection_Line *J_0 = SP->get_line_for_plane(p_0);
			if (J_0->reject_segments) continue;
			Constraint C_ts(J_0, v_ts_eps), C_os(J_0, v_os_eps);

			// Projects coordinates of v1 and v2 relatively to the frame of SP
			const CGAL_Point_2 A = SP->project(a);
			const CGAL_Point_2 B_1 = SP->project(b_1), B_2 = SP->project(b_2);
			const CGAL_Vector_2 dA_1 = B_1 - A, dA_2 = B_2 - A;

			// We are now ready to build segments.
			Polygon_Segment_R *s1_ts, *s1_os, *s2_ts, *s2_os;

			// If a single Polygon_Vertex intersects an Intersection_Line, the obtained segments
			// are, initially, two points {A_1} and {A_2} (A_1 = A_2 = A)
			s1_ts = new Polygon_Segment_R(SP->id, seed, t, C_ts, v1_ts, A, dA_1);
			s2_ts = new Polygon_Segment_R(SP->id, seed, t, C_ts, v2_ts, A, dA_2);
			Polygon_Segment_R::set_as_opposite_bidirectional_segments(s1_ts, s2_ts);

			if (v1_os != nullptr && v2_os != nullptr) {
				s1_os = new Polygon_Segment_R(SP->id, seed, t, C_os, v1_os, A, dA_1);
				s2_os = new Polygon_Segment_R(SP->id, seed, t, C_os, v2_os, A, dA_2);
				Polygon_Segment_R::set_as_opposite_bidirectional_segments(s1_os, s2_os);
			}
		}
	}



	void Support_Plane::init_2_uplets_of_bidirectional_segments(Intersection_Line* I, const CGAL_Point_2 & V1_t, const CGAL_Point_3 & W1_t, const CGAL_Point_2 & V2_t, const CGAL_Point_3 & W2_t,
		Polygon_Vertex_R* v1, Polygon_Vertex_R* v2, const FT & t)
	{
		init_4_uplets_of_bidirectional_segments(I, V1_t, W1_t, V2_t, W2_t, v1, nullptr, v2, nullptr, t);
	}



	void Support_Plane::init_4_uplets_of_bidirectional_segments(Intersection_Line* I, const CGAL_Point_2 & V1_t, const CGAL_Point_3 & W1_t, const CGAL_Point_2 & V2_t, const CGAL_Point_3 & W2_t,
		Polygon_Vertex_R* v1_ts, Polygon_Vertex_R* v1_os, Polygon_Vertex_R* v2_ts, Polygon_Vertex_R* v2_os, const FT & t)
	{
		// Here we consider the case when the bidirectional segments are real segments at the time of their creation

		assert(v1_ts->get_polygon() != nullptr);
		int seed = v1_ts->get_polygon()->seed;

		Sign v_ts_eps = v1_ts->sign_of_constraint(I);
		Sign v_os_eps = (v_ts_eps == PLUS ? MINUS : PLUS);
		int p_0 = v1_ts->id_plane;

		const CGAL_Point_3 & a_1 = W1_t, &a_2 = W2_t;
		const CGAL_Point_3 b_1 = backproject(V1_t + v1_ts->dM), b_2 = backproject(V2_t + v2_ts->dM);

		for (std::list<int>::iterator it_p = I->planes.begin(); it_p != I->planes.end(); it_p++) {
			Support_Plane* SP = Universe::map_of_planes[*it_p];

			// Builds a support constraint
			Intersection_Line *J_0 = SP->get_line_for_plane(p_0);
			Constraint C_ts(J_0, v_ts_eps), C_os(J_0, v_os_eps);
			if (J_0->reject_segments) continue;

			// Projects coordinates of v1 and v2 relatively to the frame of SP
			const CGAL_Point_2 A_1 = SP->project(a_1), A_2 = SP->project(a_2);
			const CGAL_Point_2 B_1 = SP->project(b_1), B_2 = SP->project(b_2);
			const CGAL_Vector_2 dA_1 = B_1 - A_1, dA_2 = B_2 - A_2;

			Polygon_Segment_R *s1_ts, *s1_os, *s2_ts, *s2_os;

			// If the intersection of two polygons is not empty at the initialization, 
			// or when a Polygon_Edge intersects a parallel Intersection_Line,
			// the created segments are not points {A_1} or {A_2} but segments [OA_1] and [OA_2]
			// where O is the mid-point of [A_1 A_2]
			const CGAL_Point_2 O = CGAL::midpoint(A_1, A_2);
			s1_ts = new Polygon_Segment_R(SP->id, seed, t, C_ts, v1_ts, O, A_1, dA_1);
			s2_ts = new Polygon_Segment_R(SP->id, seed, t, C_ts, v2_ts, O, A_2, dA_2);
			Polygon_Segment_R::set_as_opposite_bidirectional_segments(s1_ts, s2_ts);

			if (v1_os != nullptr && v2_os != nullptr) {
				s1_os = new Polygon_Segment_R(SP->id, seed, t, C_os, v1_os, O, A_1, dA_1);
				s2_os = new Polygon_Segment_R(SP->id, seed, t, C_os, v2_os, O, A_2, dA_2);
				Polygon_Segment_R::set_as_opposite_bidirectional_segments(s1_os, s2_os);
			}
		}
	}



	void Support_Plane::init_1_uplets_of_constant_segments(Intersection_Line* I,
		Polygon_Vertex* v1_cts, Polygon_Vertex* v2_cts, const Constraint _C_1, const Constraint _C_2, const FT & t)
	{
		// Same function as the next one, that we can directly reuse 
		// by setting v1_cos and v2_cos to nullptr

		std::list<Intersection_Line*> C_crossed_empty;
		init_2_uplets_of_constant_segments(I, v1_cts, nullptr, v2_cts, nullptr, _C_1, _C_2, C_crossed_empty, t);
	}



	void Support_Plane::init_1_uplets_of_constant_segments(Intersection_Line* I,
		Polygon_Vertex* v1_cts, Polygon_Vertex* v2_cts, const Constraint _C_1, const Constraint _C_2, 
		const std::list<Intersection_Line*> & C_crossed, const FT & t)
	{
		// Same function as the next one, that we can directly reuse 
		// by setting v1_cos and v2_cos to nullptr

		init_2_uplets_of_constant_segments(I, v1_cts, nullptr, v2_cts, nullptr, _C_1, _C_2, C_crossed, t);
	}



	void Support_Plane::init_2_uplets_of_constant_segments(Intersection_Line* I,
		Polygon_Vertex* v1_cts, Polygon_Vertex* v1_cos, Polygon_Vertex* v2_cts, Polygon_Vertex* v2_cos,
		const Constraint _C_1, const Constraint _C_2, const FT & t)
	{
		std::list<Intersection_Line*> C_crossed_empty;
		init_2_uplets_of_constant_segments(I, v1_cts, nullptr, v2_cts, nullptr, _C_1, _C_2, C_crossed_empty, t);
	}



	void Support_Plane::init_2_uplets_of_constant_segments(Intersection_Line* I,
		Polygon_Vertex* v1_cts, Polygon_Vertex* v1_cos, Polygon_Vertex* v2_cts, Polygon_Vertex* v2_cos,
		const Constraint _C_1, const Constraint _C_2, const std::list<Intersection_Line*> & C_crossed, const FT & t)
	{
		// We suppose that the aforementioned vertices are constant,
		// we want to initialize a 2-uplet of constant segments located on both sides of the line

		assert(v1_cts->get_polygon() != nullptr);
		int seed = v1_cts->get_polygon()->seed;

		Sign v_ts_eps = v1_cts->sign_of_constraint(I);
		Sign v_os_eps = (v_ts_eps == PLUS ? MINUS : PLUS);
		//const CGAL_Point_3 & a = W1_t, &b = W2_t;

		// v1 and v2 are constrained by I_1 and I_2, 
		// that represent the intersection of v1 and v2's support plane (p_0) and other planes p_1 and p_2
		// More than one plane may be associated to I_1 and I_2 but this doesn't matter

		Intersection_Line *I_1 = _C_1.first, *I_2 = _C_2.first;
		int p_0 = v1_cts->id_plane, p_1 = I_1->planes.front(), p_2 = I_2->planes.front();

		// We are going to create a segment in the planes represented by I
		// Coordinates of [AB] = [v1 v2] are computed as intersections of lines representing p_0, p_1 and p_2

		for (std::list<int>::iterator it_p = I->planes.begin(); it_p != I->planes.end(); it_p++) {
			Support_Plane* SP = Universe::map_of_planes[*it_p];
			Intersection_Line *J_0 = SP->get_line_for_plane(p_0), *J_1 = SP->get_line_for_plane(p_1), *J_2 = SP->get_line_for_plane(p_2);
			if (J_0->reject_segments) continue;

			// A = v1 (resp. B = v2), intersection of J_0 and J_1 (resp. J_2) on plane SP
			CGAL_Point_2 A, B;
			if (Universe::params->use_landmarks) {
				A = SP->get_landmark(J_0, J_1);
				B = SP->get_landmark(J_0, J_2);
			} else {
				A = SP->get_intersection_point(J_0, J_1);
				B = SP->get_intersection_point(J_0, J_2);
			}
			//const CGAL_Point_2 A = SP->project(a), B = SP->project(b);

			// Initializes constraints at limits for A and B
			Sign eps_1 = J_1->sign(B), eps_2 = J_2->sign(A);
			assert(eps_1 != ZERO && eps_2 != ZERO);
			const Constraint C_1(J_1, eps_1), C_2(J_2, eps_2);

			// Initializes support constraints for these segments
			const Constraint C_ts(J_0, v_ts_eps);

			if (C_crossed.empty()) {
				Polygon_Segment_S* s_ts = new Polygon_Segment_S(SP->id, seed, C_1, C_ts, C_2, A, B);
				if (v1_cos != nullptr && v2_cos != nullptr) {
					const Constraint C_os(J_0, v_os_eps);
					Polygon_Segment_S* s_os = new Polygon_Segment_S(SP->id, seed, C_1, C_os, C_2, A, B);
				}

			} else {
				std::list<Intersection_Line*> J_crossed;
				for (Intersection_Line* I : C_crossed) {
					int p = I->planes.front();
					Intersection_Line* J = SP->get_line_for_plane(p);
					J_crossed.push_back(J);
				}

				Polygon_Segment_S* s_ts = new Polygon_Segment_S(SP->id, seed, C_1, C_ts, J_crossed, C_2, false);
				if (v1_cos != nullptr && v2_cos != nullptr) {
					const Constraint C_os(J_0, v_os_eps);
					Polygon_Segment_S* s_os = new Polygon_Segment_S(SP->id, seed, C_1, C_os, J_crossed, C_2, false);
				}
			}

			
		}
	}



	void Support_Plane::build_n_uplets_of_segments_as_edge_collides_line(Intersection_Line* I, const Constraint & C_1, const Constraint & C_2,
		const CGAL_Point_2 & V1_t, const CGAL_Point_3 & W1_t, const CGAL_Point_2 & V2_t, const CGAL_Point_3 & W2_t,
		Polygon_Vertex* v1_ts, Polygon_Vertex* v1_os, Polygon_Vertex* v2_ts, Polygon_Vertex* v2_os, const FT & t)
	{
		// Convenience function for procedure 'edge_intersects_parallel_line'
		// We initialize different n-uplets of segments depending on the constrained/unconstrained state of v1 and v2
		// which intersected I at time t, resulting in the initialization of vertices v1_tos, v2_tos

		Polygon_Vertex_R *v1_ts_r = v1_ts->to_r(), *v2_ts_r = v2_ts->to_r();
		
		if (v1_os != nullptr && v2_os != nullptr) {
			// If the edge (v1 v2) could cross I
			Polygon_Vertex_R *v1_os_r = v1_os->to_r(), *v2_os_r = v2_os->to_r();

			if (v1_ts_r != nullptr && v2_ts_r != nullptr) {
				// Similar situation to an Event_Vertex_Line of type A,
				// except that the Polygon_Segment is initially a segment instead of a singleton
				init_4_uplets_of_bidirectional_segments(I, V1_t, W1_t, V2_t, W2_t, v1_ts_r, v1_os_r, v2_ts_r, v2_os_r, t);

			} else if (v1_ts_r != nullptr && v2_ts_r == nullptr) {
				// We initialize two segment that start from v2_tos and go towards the directions of v1_tos
				init_2_uplets_of_unidirectional_segments(I, V2_t, V1_t, W1_t, v1_ts_r, v1_os_r, C_2, t);

			} else if (v1_ts_r == nullptr && v2_ts_r != nullptr) {
				// Symmetrical case
				init_2_uplets_of_unidirectional_segments(I, V1_t, V2_t, W2_t, v2_ts_r, v2_os_r, C_1, t);

			} else {
				// v1_tos and v2_tos don't move, neither do the generated segments
				init_2_uplets_of_constant_segments(I, v1_ts, v1_os, v2_ts, v2_os, C_1, C_2, t);
			}

		} else {
			if (v1_ts_r != nullptr && v2_ts_r != nullptr) {
				init_2_uplets_of_bidirectional_segments(I, V1_t, W1_t, V2_t, W2_t, v1_ts_r, v2_ts_r, t);
			} else if (v1_ts_r != nullptr && v2_ts_r == nullptr) {
				init_1_uplets_of_unidirectional_segments(I, V2_t, V1_t, W1_t, v1_ts_r, C_2, t);
			} else if (v1_ts_r == nullptr && v2_ts_r != nullptr) {
				init_1_uplets_of_unidirectional_segments(I, V1_t, V2_t, W2_t, v2_ts_r, C_1, t);
			} else {
				init_1_uplets_of_constant_segments(I, v1_ts, v2_ts, C_1, C_2, t);
			}
		}
	}

}