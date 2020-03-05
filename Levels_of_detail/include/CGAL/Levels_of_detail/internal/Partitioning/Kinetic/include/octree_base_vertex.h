#pragma once
#include "defs_cgal.h"


namespace Skippy 
{
	class Octree_Base_Vertex
	{
	public:
		Octree_Base_Vertex(const CGAL_Point_3 & _M, const CGAL_Inexact_Point_3 & _hint_M);

		virtual ~Octree_Base_Vertex();

	public:
		CGAL_Point_3 M;
		CGAL_Inexact_Point_3 hint_M;
	};
}