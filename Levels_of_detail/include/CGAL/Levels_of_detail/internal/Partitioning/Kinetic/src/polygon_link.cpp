#include "../include/polygon_node.h"

namespace Skippy 
{
	Polygon_Link::Polygon_Link(Polygon_Node* _v1, Polygon_Node* _v2)
	{
		v1 = _v1;
		v2 = _v2;

		v1->insert_link(this);
		v2->insert_link(this);
	}


	Polygon_Link::~Polygon_Link()
	{
	
	}


	Polygon_Node* Polygon_Link::other_node(Polygon_Node* v) const
	{
		return (v1 == v ? v2 : v1);
	}
}