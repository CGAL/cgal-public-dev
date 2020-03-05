#pragma once
#include "polygon.h"


namespace Skippy 
{
	// G_Edge[1] = edge constraint in the Polygon_Node
	// G_Edge[1] x G_Edge[0] = first vertex
	// G_Edge[1] x G_Edge[2] = second vertex
	typedef std::tuple<Constraint, Constraint, Constraint> G_Edge;
	typedef std::pair<Constraint, Constraint> G_Vertex;


	class Polygon_Group 
	{
	public:
		Polygon_Group(const std::list<Polygon_Node*> & _nodes);

		~Polygon_Group();

		void make_borders(bool get_edges = true);

		size_t nodes_size() const;

		std::list<G_Edge>::const_iterator borders_begin();

		std::list<G_Edge>::const_iterator borders_end();

		std::list<G_Vertex>::const_iterator borders_v_begin();

		std::list<G_Vertex>::const_iterator borders_v_end();

	protected:
		std::list<Polygon_Node*> nodes;
		std::list<G_Edge> borders;
		std::list<G_Vertex> borders_v;
	};
}