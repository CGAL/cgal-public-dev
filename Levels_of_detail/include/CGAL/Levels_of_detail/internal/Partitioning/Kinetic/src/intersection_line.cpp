#include "../include/intersection_line.h"
#include "../include/polygon_vertex.h"
#include "../include/polygon_edge.h"
#include "../include/segment.h"
#include "../include/support_plane.h"
#include "../include/universe.h"
#include "../include/parameters.h"

#include <list>


namespace Skippy {

	Intersection_Line::Intersection_Line(const int _id_plane, const CGAL_Line_2 & _line, int intersected)
		: Support_Plane_Object(_id_plane),
		line(_line),
		_a(line.a()),
		_b(line.b()),
		_c(line.c())
	{
		double prec = 1.0 / (1 << 30) / (1 << 10);
		FT::set_relative_precision_of_to_double(prec);
		
		hint_a = to_double(_a);
		hint_b = to_double(_b);
		hint_c = to_double(_c);

		// If the id of the plane represented by the line is lower than 5,
		// then it corresponds to a facet of the bounding box.
		is_border = (intersected <= 5);

		// Considers by default that the line is inside the bounding polygon of the associated plane
		is_inside = true;

		// For space optimization purposes, it may become useless to build segments on some lines,
		// once a polygon has stopped propagating on a given support plane. However, when lines are
		// created, it is never the case.
		reject_segments = false;

		planes = std::list<int>(1, intersected);
		segments_minus = std::list<Segment *>();
		segments_plus = std::list<Segment *>();
	}


	Intersection_Line::~Intersection_Line()
	{
		if (!is_border) {
			clear_segments();
		}
	}


	void Intersection_Line::clear_segments()
	{
		std::list<Segment*> S;
		combine_signs(S);

		for (std::list<Segment*>::iterator it_s = S.begin(); it_s != S.end(); ++it_s) {
			if (Polygon_Segment_S* po_s = dynamic_cast<Polygon_Segment_S*>(*it_s)) {
				delete po_s;
			} else if (Polygon_Segment_R* po_r = dynamic_cast<Polygon_Segment_R*>(*it_s)) {
				assert(po_r != nullptr);
				delete po_r;
			}
		}

		segments_plus.clear();
		segments_minus.clear();
	}


	void Intersection_Line::combine_signs(std::list<Segment*> & S) const
	{
		S.clear();
		std::copy(segments_minus.begin(), segments_minus.end(), std::back_inserter(S));
		std::copy(segments_plus.begin(), segments_plus.end(), std::back_inserter(S));
	}


	void Intersection_Line::mark_as_intersected(int intersected)
	{
		planes.push_back(intersected);
	}



	bool Intersection_Line::intersects(const int id_plane) const
	{
		for (std::list<int>::const_iterator it_p = planes.begin(); it_p != planes.end(); it_p++) {
			if ((*it_p) == id_plane) {
				return true;
			}
		}
		return false;
	}



	void Intersection_Line::set_inside(const bool _is_inside)
	{
		is_inside = _is_inside;
	}



	bool Intersection_Line::includes(const CGAL_Point_2 & M) const
	{
		return (_a * M.x() + _b * M.y() + _c == 0);
	}



	Sign Intersection_Line::sign(const CGAL_Point_2 & pt) const
	{
		FT t = _a * pt.x() + _b * pt.y() + _c;
		if (t > 0) {
			return PLUS;
		} else if (t < 0) {
			return MINUS;
		} else {
			return ZERO;
		}
	}



	Sign Intersection_Line::sign(Polygon_Vertex* v, const FT & t) const
	{
		int K = 0;

		CGAL_Point_2 M = v->pt(t);
		FT dist = _a * M.x() + _b * M.y() + _c;

		while (dist == 0) {
			// While it is not clear whether the vertex v, at time t, 
			// is clearly on one or another side of I, we exponentially increment
			// the argument of function v->pt()

			K += 1;
			M = v->pt(t - int(pow(10, K)));
			dist = _a * M.x() + _b * M.y() + _c;
		}

		if (dist > 0) {
			return PLUS;
		} else {
			return MINUS;
		}
	}



	bool Intersection_Line::is_parallel(Intersection_Line* I) const
	{
		assert(this != I);

		const CGAL_Line_2 & I_line = I->line;
		return (CGAL::determinant(line.to_vector(), I_line.to_vector()) == 0);
	}



	bool Intersection_Line::exist_segments_including_point_outside_intersections(const CGAL_Point_2 & V_t, const FT & t) const
	{
		// We suppose that this line is not a border of the bounding box.
		// Here, we determine if there exists a couple of segments, of positive and negative signs, that include V_t at time t.

		return (exists_segment_including_point_outside_intersections(V_t, t, segments_plus)
			&& exists_segment_including_point_outside_intersections(V_t, t, segments_minus));
	}



	bool Intersection_Line::exists_segment_including_point_outside_intersections(const CGAL_Point_2 & V_t, const FT & t, const std::list<Segment*> & segments) const
	{
		// We loop on all the segments of a given sign
		// If one includes V_t at time t, we return early from the function

		for (std::list<Segment*>::const_iterator it_s = segments.begin(); it_s != segments.end(); it_s++) {
			if (Polygon_Segment* po_s = dynamic_cast<Polygon_Segment*>(*it_s)) {
				if (po_s->includes_point_on_support_line(V_t, t)) {
					return true;
				}
			}
		}

		return false;
	}



	bool Intersection_Line::exist_segments_including_point_at_intersection(const CGAL_Point_2 & V_t, const Constraint & C, const FT & t) const
	{
		// This is the same code as above,
		// except that we do not call the same method of the Polygon_Segment :
		// this time we make use of the fact that V_t is at the intersection this line and C.first

		return (exists_segment_including_point_at_intersection(V_t, C, t, segments_plus)
			&& exists_segment_including_point_at_intersection(V_t, C, t, segments_minus));
	}



	bool Intersection_Line::exists_segment_including_point_at_intersection(const CGAL_Point_2 & V_t, const Constraint & C, const FT & t, const std::list<Segment*> & segments) const
	{
		for (std::list<Segment*>::const_iterator it_s = segments.begin(); it_s != segments.end(); it_s++) {
			if (Polygon_Segment* po_s = dynamic_cast<Polygon_Segment*>(*it_s)) {
				if (po_s->includes_point_at_intersection(V_t, C, t)) {
					return true;
				}
			}
		}

		return false;
	}



	bool Intersection_Line::exist_segments_including_point_at_intersection(const CGAL_Point_2 & V_t, const std::list<Constraint> & C_limits, const FT & t) const
	{
		// This is the same code as above
		// It takes into account the case of a multi-line intersection

		return exists_segment_including_point_at_intersection(V_t, C_limits, t, segments_plus)
			&& exists_segment_including_point_at_intersection(V_t, C_limits, t, segments_minus);
	}



	bool Intersection_Line::exists_segment_including_point_at_intersection(const CGAL_Point_2 & V_t, const std::list<Constraint> & C_limits, const FT & t, const std::list<Segment*> & segments) const
	{
		for (std::list<Segment*>::const_iterator it_s = segments.begin(); it_s != segments.end(); it_s++) {
			if (Polygon_Segment* po_s = dynamic_cast<Polygon_Segment*>(*it_s)) {
				if (po_s->includes_point_at_intersection(V_t, C_limits, t)) {
					return true;
				}
			}
		}

		return false;
	}


	bool Intersection_Line::exists_segment_adjacent_to_edge(const Constraint & C_1, const Constraint & C, const Constraint & C_2) const
	{
		if (is_border) return true;

		// e = (v1 v2), where v1 = (C_1, C) and v2 = (C, C_2)

		Support_Plane* SP = Universe::map_of_planes[this->id_plane];

		CGAL_Point_2 V_1, V_2;
		if (Universe::params->use_landmarks) {
			V_1 = SP->get_landmark(C_1.first, C.first);
			V_2 = SP->get_landmark(C_2.first, C.first);
		} else {
			V_1 = SP->get_intersection_point(C_1.first, C.first);
			V_2 = SP->get_intersection_point(C_2.first, C.first);
		}
		

		//const CGAL_Point_2 & V_1 = SP->get_landmark(C_1.first, C.first);
		//const CGAL_Point_2 & V_2 = SP->get_landmark(C_2.first, C.first);

		// Gets lines that are concurrent with this object, C_1.first and C_2.first
		std::list<Intersection_Line*> CL_1, CL_2, CL;

		SP->get_concurrent_lines(this, C_1.first, CL_1);
		SP->get_concurrent_lines(this, C_2.first, CL_2);

		std::copy(CL_1.begin(), CL_1.end(), std::back_inserter(CL));
		std::copy(CL_2.begin(), CL_2.end(), std::back_inserter(CL));

		// Loops on all segments until find one that contains the edge

		for (std::list<Segment*>::const_iterator it_s = segments_plus.begin(); it_s != segments_plus.end(); it_s++) {
			if (Polygon_Segment* po_s = dynamic_cast<Polygon_Segment*>(*it_s)) {
				if (po_s->includes_edge(V_1, V_2, C_1, C_2, CL)) {
					// Early exit
					return true;
				}
			}
		}

		for (std::list<Segment*>::const_iterator it_s = segments_minus.begin(); it_s != segments_minus.end(); it_s++) {
			if (Polygon_Segment* po_s = dynamic_cast<Polygon_Segment*>(*it_s)) {
				if (po_s->includes_edge(V_1, V_2, C_1, C_2, CL)) {
					return true;
				}
			}
		}

		return false;
	}



	/*bool Intersection_Line::exists_segment_adjacent_to_edge(Polygon_Edge* e) const
	{
		if (is_border) return true;

		assert(e->v1->type == STILL_VERTEX && e->v2->type == STILL_VERTEX);

		// e = (v1 v2), where v1 and v2 are double-constrained vertices
		// We get the constraints C_1 and C_2 that correspond to such vertices,
		// in which this object is not involved

		Polygon_Vertex_S* e_v1 = e->v1->to_s();
		Constraint C_1 = e_v1->get_constraint();
		if (C_1.first == this) C_1 = e_v1->get_second_constraint();

		Polygon_Vertex_S* e_v2 = e->v2->to_s();
		Constraint C_2 = e_v2->get_constraint();
		if (C_2.first == this) C_2 = e_v2->get_second_constraint();

		// Gets vertices

		const CGAL_Point_2 & V_1 = e_v1->get_M();
		const CGAL_Point_2 & V_2 = e_v2->get_M();

		// Gets lines that are concurrent with this object, C_1.first and C_2.first
		std::list<Intersection_Line*> CL_1, CL_2, CL;

		Support_Plane* SP = Universe::map_of_planes[this->id_plane];
		SP->get_concurrent_lines(this, C_1.first, CL_1);
		SP->get_concurrent_lines(this, C_2.first, CL_2);

		std::copy(CL_1.begin(), CL_1.end(), std::back_inserter(CL));
		std::copy(CL_2.begin(), CL_2.end(), std::back_inserter(CL));

		// Loops on all segments until find one that contains the edge

		for (std::list<Segment*>::const_iterator it_s = segments.begin(); it_s != segments.end(); it_s++) {
			if (Polygon_Segment* po_s = dynamic_cast<Polygon_Segment*>(*it_s)) {
				if (po_s->includes_edge(V_1, V_2, C_1, C_2, CL)) {
					// Early exit
					return true;
				}
			}
		}

		return false;
	}*/


	bool Intersection_Line::represents_same_intersection(const std::pair<Constraint, Constraint> & v1, const std::pair<Constraint, Constraint> & v2)
	{
		const Constraint & C_11 = v1.first, & C_12 = v1.second;
		const Constraint & C_21 = v2.first, & C_22 = v2.second;
		
		// We get the indices of the lines that correspond to such constraints.
		const Intersection_Line *I_11 = C_11.first, *I_12 = C_12.first;
		const Intersection_Line *I_21 = C_21.first, *I_22 = C_22.first;
		int p = I_11->id_plane;

		// If (C_1, C_2) and (C_v1, C_v2) represent the same lines, then this vertex and v represent the same intersection
		if ((I_11 == I_21 && I_12 == I_22) || (I_11 == I_22 && I_12 == I_21)) return true;

		// If we couldn't simultaneously match (C_1, C_2) and (C_v1, C_v2), 
		// then we get triplets of concurrent lines in which we find (C_1 and C_2)
		// If C_v1 or C_v2 is part of such a triplet, then we return true.
		std::list<Intersection_Line *> I_L;
		Universe::map_of_planes[p]->get_concurrent_lines(I_11, I_12, I_L);

		for (Intersection_Line* L : I_L) {
			if (I_21 == L || I_22 == L) return true;
		}

		// The search failed
		return false;
	}


	const FT & Intersection_Line::a() const
	{
		return _a;
	}


	const FT & Intersection_Line::b() const
	{
		return _b;
	}


	const FT & Intersection_Line::c() const
	{
		return _c;
	}


	
	void Intersection_Line::reject_any_further_segment()
	{
		reject_segments = true;

		if (!is_border) {
			reject_any_further_segment(segments_plus);
			reject_any_further_segment(segments_minus);
		}
	}


	void Intersection_Line::reject_any_further_segment(std::list<Segment*> & S)
	{
		for (std::list<Segment*>::iterator it_s = S.begin(); it_s != S.end(); ++it_s) {
			if (Polygon_Segment_S* po_s = dynamic_cast<Polygon_Segment_S*>(*it_s)) {
				delete po_s;
				it_s = S.erase(it_s);
			} else if (Polygon_Segment_R* po_r = dynamic_cast<Polygon_Segment_R*>(*it_s)) {
				std::cout << "TODO" << std::endl;
			}
		}
	}
}