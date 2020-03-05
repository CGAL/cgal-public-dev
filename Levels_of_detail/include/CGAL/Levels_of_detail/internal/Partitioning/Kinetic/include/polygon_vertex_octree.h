#pragma once
#include <list>
#include <vector>
#include <cassert>
#include <iterator>
#include "defs.h"
#include "octree_base_vertex.h"
#include "octree_base.h"
#include "polygon_vertex_octree_data.h"


namespace Skippy 
{
	class Polygon_Vertex_Octree : public Octree_Base
	{
	public:
		Polygon_Vertex_Octree(const std::vector<double> & dims);

		~Polygon_Vertex_Octree();

		int get_identifier(const CGAL_Point_3 & M, const CGAL_Inexact_Point_3 & hint_M);

		void get_all_sorted_vertices(std::vector<Polygon_Vertex_Octree_Data*> & R) const;

	public:
		static int id_vertex;
	};
}