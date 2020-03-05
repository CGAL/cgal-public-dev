#pragma once
#include "defs_cgal.h"
#include "octree_base_vertex.h"


namespace Skippy 
{
	class Polygon_Vertex_Octree_Data : public Octree_Base_Vertex
	{
	public:
		Polygon_Vertex_Octree_Data(const CGAL_Point_3 & _M, const CGAL_Inexact_Point_3 & _hint_M, const int & _id);

		~Polygon_Vertex_Octree_Data();

	public:
		int id;
	};
}