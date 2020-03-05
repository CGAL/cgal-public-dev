#include "../include/polygon_group.h"
#include "../include/support_plane_objects.h"
#include "../include/polygon_vertex.h"
#include "../include/polygon_edge.h"
#include "../include/intersection_line.h"
#include "../include/polygon_node.h"



namespace Skippy 
{
	Polygon_Group::Polygon_Group(const std::list<Polygon_Node*> & _nodes)
	{
		nodes = _nodes;
		make_borders();
	}


	Polygon_Group::~Polygon_Group()
	{

	}


	void Polygon_Group::make_borders(bool get_edges)
	{
		// Each Polygon_Node has :
		// - a list of delimiter constraints, 
		// - a list of booleans that indicate if they stand on the borders of the group of connected nodes.
		// Here, we loop on these borders to get a complete sequence.

		// Step 1.
		// Gets an unordered sequence of edges.

		std::list<G_Edge> unordered;
		for (std::list<Polygon_Node*>::const_iterator it_n = nodes.begin() ; it_n != nodes.end() ; ++it_n) {
			Polygon_Node* n = (*it_n);
			n->append_contours(unordered);
		}
		
		// Step 2.
		// Ranks points.

		G_Edge e_prev = unordered.front();
		G_Vertex v_anchor = std::make_pair(std::get<1>(e_prev), std::get<2>(e_prev));

		if (get_edges) {
			borders.push_back(e_prev);
		} else {
			borders_v.push_back(v_anchor);
		}
		unordered.pop_front();

		while (!unordered.empty()) {
			std::list<G_Edge>::iterator it_e = unordered.begin();

			// Loops on the remaining edges until finding the one that is connected to v_anchor
			// Exits the loop when such an element is found
			while (true) {
				G_Edge e = (*it_e);
				G_Vertex e_v1 = std::make_pair(std::get<0>(e), std::get<1>(e));
				G_Vertex e_v2 = std::make_pair(std::get<1>(e), std::get<2>(e));
				if (Intersection_Line::represents_same_intersection(v_anchor, e_v1)) {
					v_anchor = e_v2;
					break;
				} else if (Intersection_Line::represents_same_intersection(v_anchor, e_v2)) {
					v_anchor = e_v1;
					break;
				} else {
					++it_e;
				}
			}

			// Adds the element to the list of ordered borders
			// Removes it from the list of unordered elements
			if (get_edges) {
				borders.push_back(*it_e);
			} else {
				borders_v.push_back(v_anchor);
			}
			unordered.erase(it_e);
		}
	}


	size_t Polygon_Group::nodes_size() const
	{
		return nodes.size();
	}


	std::list<G_Edge>::const_iterator Polygon_Group::borders_begin()
	{
		return borders.cbegin();
	}


	std::list<G_Edge>::const_iterator Polygon_Group::borders_end()
	{
		return borders.cend();
	}


	std::list<G_Vertex>::const_iterator Polygon_Group::borders_v_begin()
	{
		return borders_v.cbegin();
	}


	std::list<G_Vertex>::const_iterator Polygon_Group::borders_v_end()
	{
		return borders_v.cend();
	}
}