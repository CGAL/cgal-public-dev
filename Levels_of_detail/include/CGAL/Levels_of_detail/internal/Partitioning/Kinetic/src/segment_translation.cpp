#include "../include/segment_translation.h"


namespace Skippy {

	Segment_Translation::Segment_Translation(const FT & t, const CGAL_Point_2 & _A, const CGAL_Point_2 & _B)
		: type(INSTANTANEOUS),
		t_int_start(t),
		t_int_end(t),
		A(_A),
		dA(CGAL_Vector_2(0, 0)),
		B(_B)
	{

	}


	Segment_Translation::Segment_Translation(const FT & t, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA)
		: type(PROGRESSIVE),
		t_int_start(t),
		t_int_end(FLT_MAX),
		A(_A),
		dA(_dA),
		B(CGAL_Point_2(0, 0))
	{

	}


	Segment_Translation::~Segment_Translation()
	{

	}


	void Segment_Translation::set_end(const FT & t)
	{
		assert(type == PROGRESSIVE);

		t_int_end = t;
		B = A + (t_int_end - t_int_start) * dA;
	}
}