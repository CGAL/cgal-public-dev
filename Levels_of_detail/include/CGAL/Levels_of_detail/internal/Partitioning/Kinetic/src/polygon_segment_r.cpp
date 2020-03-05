#include "../include/segment.h"
#include "../include/segment_translation.h"
#include "../include/polygon_vertex.h"
#include "../include/support_plane.h"
#include "../include/universe.h"
#include "../include/parameters.h"
#include "../include/intersection_line.h"


namespace Skippy
{
	Polygon_Segment_R::Polygon_Segment_R(const int _id_plane, const int _seed, const Constraint & C_support)
		: Polygon_Segment(_id_plane, _seed, C_support)
	{
		sign = C_support.second;

		A = nullptr;

		Tr = nullptr;
		Tr_previous = std::list<Segment_Translation *>();

		opposite = nullptr;
	}


	Polygon_Segment_R::Polygon_Segment_R(const int _id_plane, const int _seed, const FT & _t_init, const Constraint & C_support, Polygon_Vertex_R* & v, const CGAL_Point_2 & _O, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA)
		: Polygon_Segment_R(_id_plane, _seed, C_support)
	{
		// An object of class Polygon_Segment can be divided into three different parts :
		//
		// a) An initial segment [O A_0] where A_0 = Tr_previous.front()->A. [O A_0] is null if the segment is created at t > 0.
		// b) A list of segments [A_i A_{i + 1}] where each segment corresponds to the path followed by the ray 
		//    between t->int_start and t->int_end starting from t->A, where t is an element of Tr_previous.
		// c) Finally, the half-line [Tr->A (Tr->A + t * Tr->dA)). It doesn't exist if the segment is not active anymore.

		A = new CGAL_Point_2(_O);
		Tr = new Segment_Translation(_t_init, _A, _dA);
		Tr_previous.push_back(new Segment_Translation(_t_init, _O, _A));

		// This construtor is typically called during the initialization process,
		// or when a Polygon_Edge with two unconstrained vertices intersects an Intersection_Line.

		v->set_as_guided_segment(this);
	}



	Polygon_Segment_R::Polygon_Segment_R(const int _id_plane, const int _seed, const FT & _t_init, const Constraint & C_support, Polygon_Vertex_R* & v, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA)
		: Polygon_Segment_R(_id_plane, _seed, C_support)
	{
		// This constructor is typically called as an unconstrained Polygon_Vertex intersects an Intersection_Line.

		A = new CGAL_Point_2(_A);
		Tr = new Segment_Translation(_t_init, _A, _dA);

		v->set_as_guided_segment(this);
	}



	Polygon_Segment_R::Polygon_Segment_R(const int _id_plane, const int _seed, const FT & _t_init, const Constraint & _C_init, const Constraint & C_support, Polygon_Vertex_R* & v,
		const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA)
		: Polygon_Segment_R(_id_plane, _seed, C_support)
	{
		// This constructor is typically called as an constrained Polygon_Vertex intersects an Intersection_Line.

		//Universe::map_of_planes[id_plane]->set_landmark(C_support.first, _C_init.first, _A);
		Tr = new Segment_Translation(_t_init, _A, _dA);

		C_init = _C_init;

		v->set_as_guided_segment(this);
	}



	Polygon_Segment_R::Polygon_Segment_R(const int _id_plane, const int _seed, const FT & _t_init, const Constraint & _C_init, const Constraint & C_support, Polygon_Vertex_R* & v,
		const CGAL_Point_2 & _O, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA)
		: Polygon_Segment_R(_id_plane, _seed, C_support)
	{
		// This constructor builds Polygon_Segments which are created at time t.
		// Contrary to before, at the time of its creation, the Polygon_Segment is a segment [OA] and not a single point {A}.
		// It is typically called when a Polygon_Edge with one constrained vertex intersects an Intersection_Line.

		//Universe::map_of_planes[id_plane]->set_landmark(_C_init.first, C_support.first, _O);

		Tr = new Segment_Translation(_t_init, _A, _dA);
		Tr_previous.push_back(new Segment_Translation(_t_init, _O, _A));

		C_init = _C_init;

		v->set_as_guided_segment(this);
	}


	Polygon_Segment_R::~Polygon_Segment_R()
	{
		if (Tr != nullptr) {
			delete Tr;
			Tr = nullptr;
		}

		for (std::list<Segment_Translation*>::iterator it = Tr_previous.begin(); it != Tr_previous.end(); it++) delete (*it);
		Tr_previous.clear();
	}


	Sign Polygon_Segment_R::get_sign() const
	{
		return sign;
	}


	CGAL_Point_2 Polygon_Segment_R::origin() const
	{
		if (A != nullptr) {
			return *A;
		} else {
			Support_Plane* SP = Universe::map_of_planes[id_plane];
			if (Universe::params->use_landmarks) {
				return SP->get_landmark(support, C_init.first);
			} else {
				return SP->get_intersection_point(support, C_init.first);
			}
		}
	}


	/*const CGAL_Point_2 & Polygon_Segment_R::origin() const
	{
		if (A != nullptr) {
			return *A;
		} else {
			//return Universe::map_of_planes[id_plane]->get_landmark(support, C_init.first);
			return Universe::map_of_planes[id_plane]->get_intersection_point(support, C_init.first);
		}
	}*/


	bool Polygon_Segment_R::stopped() const
	{
		return (Tr == nullptr);
	}


#pragma warning(push)
#pragma warning(disable:4715)
	CGAL_Point_2 Polygon_Segment_R::pt(const FT & t) const
	{
		if (stopped()) return end();

		// We first check if the requested t corresponds to the current interval
		if (Tr != nullptr && Tr->t_int_start <= t && t <= Tr->t_int_end) {
			return Tr->A + (t - Tr->t_int_start) * Tr->dA;
		}

		// If it not the case, then we loop on all the previous intervals 
		// until finding the one that contains t 
		else {
			assert(Tr_previous.size() > 0);
			assert(t >= Tr_previous.front()->t_int_start);

			// Loops on different intervals
			std::list<Segment_Translation*>::const_iterator it;
			for (it = Tr_previous.begin(); it != Tr_previous.end(); it++) {
				Segment_Translation* T = (*it);
				if (T->type == PROGRESSIVE) {
					// If the translation is progressive, 
					// i.e. if a polygon progressively collides with another polygon
					if (T->t_int_start <= t && t <= T->t_int_end) {
						return T->A + (t - T->t_int_start) * T->dA;
					}
				} else {
					// If the translation is instantaneous,
					// i.e. if this segment has been created at t = 0,
					// or if an edge of polygon has intersected orthogonally another polygon at t > 0
					// By convention it seems more logical to take the last point of the interval,
					// the one that is the closest to the trajectory of the segment in the next 
					// temporal interval
					if (T->t_int_end == t) return T->B;
				}
			}

			assert(it != Tr_previous.end());
		}
	}
#pragma warning(pop)


	void Polygon_Segment_R::update_translation(const FT & t, const CGAL_Point_2 & A_t, const CGAL_Vector_2 & dA_t)
	{
		// Indicates that the current translation is no longer valid after time t,
		// and inserts Tr into the list of previous translations
		Tr->set_end(t);
		Tr_previous.push_back(Tr);

		// Creates a new translation
		Tr = new Segment_Translation(t, A_t, dA_t);
	}


	void Polygon_Segment_R::stop(const Constraint C, const FT & t)
	{
		// Indicates that the current translation is no longer valid after time t,
		// and inserts Tr into the list of previous translations

		//CGAL_Point_2 B = Tr->A + (t - Tr->t_int_start) * Tr->dA;
		Support_Plane* SP = Universe::map_of_planes[id_plane];
		Intersection_Line* I = C.first;

		if (Universe::params->use_landmarks) {
			SP->get_landmark(I, support);
		}

		//const CGAL_Line_2 & L = I->line;
		//const CGAL_Line_2 & L_S = support->line;
		//CGAL_Point_2 B = SP->get_landmark(I, support);

		/*if (SP->exists_landmark(I, support)) {
			B = SP->get_landmark(I, support);
		} else {
			CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object = CGAL::intersection(L, L_S);
			if (object) {
				if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object)) {
					B = *ptr;
				}
			}
			SP->set_landmark(I, support, B);
		}*/

		delete Tr;
		Tr = nullptr;

		for (std::list<Segment_Translation*>::iterator it_t = Tr_previous.begin(); it_t != Tr_previous.end(); ++it_t) {
			delete (*it_t);
		}
		Tr_previous.clear();

		C_stop = C;
	}



	// The next function only concern opposite bidirectional segments,
	// which are initialized with a C_init = nullptr.

	// In order to improve the later performances of the algorithm,
	// we would like to find the closest lines to the origins of bidirectional segments,
	// and set C_init anyway to a consistent value.



	void Polygon_Segment_R::set_as_opposite_bidirectional_segments(Polygon_Segment_R* s_1, Polygon_Segment_R* s_2)
	{
		s_1->opposite = s_2;
		s_2->opposite = s_1;
	}



	bool Polygon_Segment_R::exists_opposite_segment_without_initinal_constraint() const
	{
		// This function is for preventing useless and expensive computations.

		return (opposite != nullptr && opposite->C_init.first == nullptr);
	}



	Polygon_Segment_R* Polygon_Segment_R::get_opposite() const
	{
		return opposite;
	}



	void Polygon_Segment_R::set_pseudo_init_constraint(const CGAL_Point_2 & A, std::list<Intersection_Line*> & L)
	{
		assert(C_init.first == nullptr);

		// We consider two opposite bidirectional segments :
		// - A is the extremum of the first one (this object),
		// - L the list of Intersection_Lines that cut the second.

		// We can set an initial constraint for the first segment,
		// using the line that is the closest to it.
		// To detect this line, we use the coordinates of A.

		FT dist_min = FLT_MAX;
		FT absolute_dist_min = FLT_MAX;
		Intersection_Line* argmin = nullptr;

		for (std::list<Intersection_Line*>::iterator it_l = L.begin(); it_l != L.end(); it_l++) {
			Intersection_Line* J = (*it_l);

			const CGAL_Line_2 & J_line = J->line;
			FT dist = J_line.a() * A.x() + J_line.b() * A.y() + J_line.c();
			FT absolute_dist = CGAL::abs(dist);
			if (absolute_dist < absolute_dist_min) {
				dist_min = dist;
				absolute_dist_min = absolute_dist;
				argmin = J;
			}
		}

		C_init = Constraint(argmin, dist_min < 0 ? MINUS : PLUS);
	}



	void Polygon_Segment_R::set_pseudo_init_constraint(const Constraint C_pseudo_init)
	{
		assert(C_init.first == nullptr);

		C_init = C_pseudo_init;

		delete A;
		A = nullptr;

		if (Universe::params->use_landmarks) {
			Universe::map_of_planes[id_plane]->get_landmark(C_init.first, support);
		}

		/*if (!Universe::map_of_planes[id_plane]->exists_landmark(C_init.first, support)) {
			const CGAL_Line_2 & L_init = C_init.first->line;
			const CGAL_Line_2 & L_supp = support->line;

			FT det = L_init.a() * L_supp.b() - L_supp.a() * L_init.b();
			FT x = (L_init.b() * L_supp.c() - L_init.c() * L_supp.b()) / det;
			FT y = (L_init.c() * L_supp.a() - L_init.a() * L_supp.c()) / det;
			CGAL_Point_2 M = CGAL_Point_2(x, y);

			Universe::map_of_planes[id_plane]->set_landmark(C_init.first, support, M);
		}*/
	}



	Constraint Polygon_Segment_R::get_pseudo_init_constraint() const
	{
		return C_init;
	}


	bool Polygon_Segment_R::includes_point_on_support_line(const CGAL_Point_2 & M, const FT & t) const
	{
		// We assume that M is a point of the support line of this segment.
		// Here, we determine if at time t, M is contained by this segment,
		// which is the case if A.x <= M.x <= B.x, or A.y <= M.y <= B.y,
		// where A and B represent the segments' ends.

		const CGAL_Point_2 & A = origin();
		CGAL_Point_2 B = pt(t);
		return closed_segment_includes(M, A, B);
	}



	bool Polygon_Segment_R::includes_point_at_intersection(const CGAL_Point_2 & M, const Constraint & C, const FT & t) const
	{
		// Like before, we assume that M is a point of the support line of this segment.
		// Here, we know that M is the intersection of the 'support' line and C.first,
		// and we also know on which side (or halfplane) M is located, via C.second

		Intersection_Line* J = C.first;
		Sign J_eps = C.second;

		// Test (1).
		// If J is crossed by the segment, then M belongs to it
		std::list<Intersection_Line*>::const_iterator it_l = std::find(C_crossed.cbegin(), C_crossed.cend(), J);
		if (it_l != C_crossed.cend()) return true;

		// Test (2).
		// Compares J to the lines that delimit the segment.
		// If they are not nullptr, then the segment includes M 
		// iff the initial or stop constraint have the same sign as C

		Intersection_Line *I_init = C_init.first, *I_stop = C_stop.first;
		if (I_init == J) {
			return (C_init.second == J_eps);
		} else if (I_stop == J) {
			return (C_stop.second == J_eps);
		}

		// Test (3).
		// If the initial and stop conditions of the segment are determined,
		// then M is necessarily located outside the point because J is not I_init, I_stop or in C_crossed
		if (I_init != nullptr && I_stop != nullptr) return false;

		// // Unfortunately we have to perform a numerical computation
		// // M belongs to the line, so if A and B represent the segments' ends,
		// // we should have A.x <= M.x <= B.x (or A.y <= M.y <= B.y)
		// return includes_point_on_support_line(M, t);

		// Test (4).
		// All the remaining cases imply numerical computations.

		// If s has stopped :
		// (4.1) We return the result of the test (M in s).
		// Else :
		// (4.2) We consider A = origin(), B = pt(t).
		//       If M is in [AB[, then the segment includes M : we return true.
		//       Else, if M = B, then determine on which side of the line J is pt(t).
		//       Else, we return false (M is outside [AB]).

		const CGAL_Point_2 & A = origin();

		// Test (4.2)
		CGAL_Point_2 B = pt(t);
		if (half_closed_segment_includes(M, A, B)) {
			return true;
		} else if (M == B) {
			return (J_eps == J->sign(A));
		} else {
			return false;
		}
	}


	void Polygon_Segment_R::who_is() const
	{
		who_is();
		std::cout << "opposite = " << opposite << std::endl;
	}


	const Segment_Translation* Polygon_Segment_R::get_current_translation() const
	{
		return Tr;
	}
}