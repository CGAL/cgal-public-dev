#include "../include/segment.h"
#include "../include/support_plane.h"
#include "../include/universe.h"
#include "../include/intersection_line.h"
#include "../include/parameters.h"



namespace Skippy
{
	Polygon_Segment_S::Polygon_Segment_S(const int _id_plane, const int _seed, const Constraint & _C_init, const Constraint & C_support, const Constraint & _C_stop,
		const CGAL_Point_2 & _A, const CGAL_Point_2 & _B)
		: Polygon_Segment(_id_plane, _seed, C_support)
	{
		// This constructor builds Polygon_Segments that are constant.
		// At the time of their creation, their initial and final points and constraints are already determined.
		// It is typically called when a Polygon_Edge with two constrained vertices intersects an Intersection_Line.

		/*if (Universe::params->use_landmarks) {
			Universe::map_of_planes[id_plane]->set_landmark(_C_init.first, C_support.first, _A);
			Universe::map_of_planes[id_plane]->set_landmark(_C_stop.first, C_support.first, _B);
		}*/

		std::list<Intersection_Line*> empty_list;
		init_by_merging(_C_init, C_support, empty_list, _C_stop);
	}


	Polygon_Segment_S::Polygon_Segment_S(const int _id_plane, const int _seed, const Constraint & _C_init, const Constraint & C_support, const std::list<Intersection_Line*> & _C_crossed, const Constraint & _C_stop, const bool merge_at_boundaries)
		: Polygon_Segment(_id_plane, _seed, C_support)
	{
		if (merge_at_boundaries) {
			// This constructor is called as a unidirectional Polygon_Segment_R stops,
			// giving us full knowledge of its initial and final constraints, as well as the list of edges it crossed.

			// No need to call the set_landmark function : we assume that it has been done in Polygon_Segment_R::stop().
			init_by_merging(_C_init, C_support, _C_crossed, _C_stop);

		} else {
			// Contrary to before we already know the lines between I_init and I_stop
			// and there is no need to call the function init_by_merging
			C_init = _C_init;
			C_stop = _C_stop;
			C_crossed = _C_crossed;
		}
	}


	Polygon_Segment_S::Polygon_Segment_S(const int _id_plane, const int _seed, const Constraint & _C_init, const Constraint & C_support, const std::list<Intersection_Line*> & C_crossed_1, const std::list<Intersection_Line*> & C_crossed_2, const Constraint & _C_stop)
		: Polygon_Segment(_id_plane, _seed, C_support)
	{
		// This constructor is called as two opposite bidirectional Polygon_Segment_R objects have stopped.
		// We now obtain a big segment (s1, s2) = (s1->C_stop, s1->crossed, s2->crossed, s2->C_stop).

		// Like before, no need to call the set_landmark function.

		std::list<Intersection_Line*> L;
		std::copy(C_crossed_1.begin(), C_crossed_1.end(), std::back_inserter(L));
		std::copy(C_crossed_2.begin(), C_crossed_2.end(), std::back_inserter(L));

		init_by_merging(_C_init, C_support, L, _C_stop);
	}



	Polygon_Segment_S::~Polygon_Segment_S()
	{
		
	}


	Constraint Polygon_Segment_S::get_other_constraint(const Constraint & C) const
	{
		return (C_init == C ? C_stop : C_init);
	}


	void Polygon_Segment_S::init_by_merging(const Constraint & _C_init, const Constraint & C_support, std::list<Intersection_Line*> _C_crossed, const Constraint & _C_stop)
	{
		C_init = _C_init;
		C_stop = _C_stop;
		C_crossed = _C_crossed;

		// This Polygon_Segment may be adjacent to one or two Polygon_Segment_S objects,
		// created by the same primitive during the propagation phase.

		// If it is the case, these adjacent objects will be merged with the current one,
		// whose initial and final constraints will be changed.

		Intersection_Line* I = C_support.first;
		Sign eps = C_support.second;
		std::list<Segment*> & L = (eps == PLUS ? I->segments_plus : I->segments_minus);

		Constraint C_1 = std::make_pair(C_init.first, C_init.second == PLUS ? MINUS : PLUS);
		Constraint C_2 = std::make_pair(C_stop.first, C_stop.second == PLUS ? MINUS : PLUS);

		Polygon_Segment_S* s_1 = get_adjacent_segment(L, C_1);
		Polygon_Segment_S* s_2 = get_adjacent_segment(L, C_2);

		if (s_1 != nullptr && s_2 == nullptr) {
			
			// The initial constraint gets equal to the other constraint of s_1
			// We add I_init and s_1->C_crossed to C_crossed
			// Finally s_1 is deleted

			C_init = s_1->get_other_constraint(C_1);
			C_crossed.push_back(C_1.first);
			std::copy(s_1->C_crossed.begin(), s_1->C_crossed.end(), std::back_inserter(C_crossed));
			L.erase(s_1->iterator);
			delete s_1;

		} else if (s_1 == nullptr && s_2 != nullptr) {

			// Similar reasoning with s_2 and C_stop
			
			C_stop = s_2->get_other_constraint(C_2);
			C_crossed.push_back(C_2.first);
			std::copy(s_2->C_crossed.begin(), s_2->C_crossed.end(), std::back_inserter(C_crossed));
			L.erase(s_2->iterator);
			delete s_2;

		} else if (s_1 != nullptr && s_2 != nullptr) {
			
			// The initial and final constraints gets equal to those of s_1 and s_2
			// I_init, I_stop, s_1->C_crossed, s_2->C_crossed are added to C_crossed
			// s_1 and s_2 are finally deleted.

			C_init = s_1->get_other_constraint(C_1);
			C_stop = s_2->get_other_constraint(C_2);
			C_crossed.push_back(C_1.first);
			C_crossed.push_back(C_2.first);
			std::copy(s_1->C_crossed.begin(), s_1->C_crossed.end(), std::back_inserter(C_crossed));
			std::copy(s_2->C_crossed.begin(), s_2->C_crossed.end(), std::back_inserter(C_crossed));
			L.erase(s_1->iterator);
			L.erase(s_2->iterator);
			delete s_1;
			delete s_2;
		}

		if (Universe::params->cl_flush != -1) {
			if (C_crossed.size() > Universe::params->cl_flush) C_crossed.clear();
		}
	}


	Polygon_Segment_S* Polygon_Segment_S::get_adjacent_segment(const std::list<Segment*> & L, const Constraint & C) const
	{
		for (std::list<Segment*>::const_iterator it_l = L.begin() ; it_l != L.end() ; ++it_l) {
			if (it_l == iterator) continue;

			if (Polygon_Segment_S* po_s = dynamic_cast<Polygon_Segment_S*>(*it_l)) {
				if (po_s->seed == seed && (po_s->C_init == C || po_s->C_stop == C)) return po_s;
			}
		}

		return nullptr;
	}


	CGAL_Point_2 Polygon_Segment_S::origin() const
	{
		if (Universe::params->use_landmarks) {
			return Universe::map_of_planes[id_plane]->get_landmark(support, C_init.first);
		} else {
			return Universe::map_of_planes[id_plane]->get_intersection_point(support, C_init.first);
		}
	}


	/*const CGAL_Point_2 & Polygon_Segment_S::origin() const
	{
		//return Universe::map_of_planes[id_plane]->get_landmark(support, C_init.first);
		return Universe::map_of_planes[id_plane]->get_intersection_point(support, C_init.first);
	}*/



	bool Polygon_Segment_S::includes_point_on_support_line(const CGAL_Point_2 & M, const FT & t) const
	{
		// We assume that M is a point of the support line of this segment.
		// Here, we determine if at time t, M is contained by this segment,
		// which is the case if A.x <= M.x <= B.x, or A.y <= M.y <= B.y,
		// where A and B represent the segments' ends.

		const CGAL_Point_2 & A = origin();
		const CGAL_Point_2 & B = end();
		return closed_segment_includes(M, A, B);
	}



	bool Polygon_Segment_S::includes_point_at_intersection(const CGAL_Point_2 & M, const Constraint & C, const FT & t) const
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
		if (Universe::params->cl_flush == -1) {
			if (I_init != nullptr && I_stop != nullptr) return false;
		}
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

		// Test (4.1)
		const CGAL_Point_2 & A = origin();
		const CGAL_Point_2 & B = end();
		return closed_segment_includes(M, A, B);
	}
}