#include "../include/segment.h"
#include "../include/support_plane.h"
#include "../include/universe.h"


namespace Skippy {
	Planar_Segment::Planar_Segment(const int _id_plane, const Constraint & C_support, CGAL_Point_2 & _A, CGAL_Point_2 & _B)
		: Segment(_id_plane, C_support)
	{
		Universe::map_of_planes[_id_plane]->borders[id_object] = this;

		A = _A;
		B = _B;
	}


	Planar_Segment::~Planar_Segment()
	{
	}


	bool Planar_Segment::checks_if_belongs(const CGAL_Point_2 & V) const
	{
		// Returns true if V is on [AB], false otherwise

		CGAL_Segment_2 AB(A, B);
		return AB.has_on(V);
	}
}