#include "../include/segment.h"
#include "../include/intersection_line.h"



namespace Skippy {

	Segment::Segment(const int _id_plane, const Constraint & C_support)
		: Support_Plane_Object(_id_plane)
	{
		support = C_support.first;
		if (C_support.second == PLUS) {
			support->segments_plus.push_back(this);
		} else {
			support->segments_minus.push_back(this);
		}
	}


	Segment::~Segment()
	{

	}


	bool Segment::closed_segment_includes(const CGAL_Point_2 & M, const CGAL_Point_2 & A, const CGAL_Point_2 & B)
	{
		// Returns the result of (M in [AB]).

		FT x_a = A.x(), x_b = B.x(), x = M.x();
		if (x_a != x_b) {
			return ((x_a <= x && x <= x_b) || (x_b <= x && x <= x_a));
		} else {
			FT y_a = A.y(), y_b = B.y(), y = M.y();
			return ((y_a <= y && y <= y_b) || (y_b <= y && y <= y_a));
		}
	}


	bool Segment::half_closed_segment_includes(const CGAL_Point_2 & M, const CGAL_Point_2 & A, const CGAL_Point_2 & B)
	{
		// Returns the result of (M in [AB[).

		FT x_a = A.x(), x_b = B.x(), x = M.x();
		if (x_a != x_b) {
			return ((x_a <= x && x < x_b) || (x_b < x && x <= x_a));
		} else {
			FT y_a = A.y(), y_b = B.y(), y = M.y();
			return ((y_a <= y && y < y_b) || (y_b < y && y <= y_a));
		}
	}
}