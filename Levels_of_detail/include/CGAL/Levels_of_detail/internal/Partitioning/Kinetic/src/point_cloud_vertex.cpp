#include "../include/point_cloud_vertex.h"


namespace Skippy 
{
	Point_Cloud_Vertex::Point_Cloud_Vertex(const CGAL_Point_3 & _M, const CGAL_Inexact_Point_3 & _hint_M, const CGAL_Inexact_Vector_3 & _N)
		: Octree_Base_Vertex (_M, _hint_M)
	{
		hint_normal = _N;
	}


	Point_Cloud_Vertex::~Point_Cloud_Vertex()
	{
	}
}