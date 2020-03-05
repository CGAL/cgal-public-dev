#include "../include/octree_base_vertex.h"


namespace Skippy 
{
	Octree_Base_Vertex::Octree_Base_Vertex(const CGAL_Point_3 & _M, const CGAL_Inexact_Point_3 & _hint_M)
	{
		M = _M;
		hint_M = _hint_M;
	}


	Octree_Base_Vertex::~Octree_Base_Vertex()
	{
	}
}