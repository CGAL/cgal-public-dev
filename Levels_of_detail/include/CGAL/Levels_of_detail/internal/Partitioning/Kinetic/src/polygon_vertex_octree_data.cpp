#include "../include/polygon_vertex_octree_data.h"



namespace Skippy 
{
	Polygon_Vertex_Octree_Data::Polygon_Vertex_Octree_Data(const CGAL_Point_3 & _M, const CGAL_Inexact_Point_3 & _hint_M, const int & _id)
		: Octree_Base_Vertex(_M, _hint_M)
	{
		id = _id;
	}


	Polygon_Vertex_Octree_Data::~Polygon_Vertex_Octree_Data()
	{
	
	}
}