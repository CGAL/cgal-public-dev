#pragma once
#include "octree_base_vertex.h"
#include "defs_cgal.h"


namespace Skippy
{
	class Point_Cloud_Vertex : public Octree_Base_Vertex
	{
	public:
		Point_Cloud_Vertex(const CGAL_Point_3 & _M, const CGAL_Inexact_Point_3 & _hint_M, const CGAL_Inexact_Vector_3 & _N);

		~Point_Cloud_Vertex();

	public:
		CGAL_Inexact_Vector_3 hint_normal;
	};
}