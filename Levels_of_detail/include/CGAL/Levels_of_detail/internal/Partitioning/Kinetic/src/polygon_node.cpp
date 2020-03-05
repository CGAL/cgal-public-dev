#include "../include/polygon_node.h"
#include "../include/intersection_line.h"


namespace Skippy 
{
	Polygon_Node::Polygon_Node(const Signature & S, Polygon* P)
	{
		signature = S;
		push(P);

		borders = nullptr;

		visited_node = false;
		queued_node = false;
	}


	Polygon_Node::~Polygon_Node()
	{
		for (std::list<Polygon*>::iterator it_p = polygons.begin(); it_p != polygons.end(); it_p++) {
			delete (*it_p);
		}
		polygons.clear();

		if (borders != nullptr) delete[] borders;
	}


	Signature Polygon_Node::make_signature(const CGAL_Point_2 & O, const std::vector<Intersection_Line*> & L)
	{
		// Initializes a vector
		Signature S(L.size(), false);

		int entry = -1;
		for (std::vector<Intersection_Line*>::const_iterator it_l = L.begin(); it_l != L.end(); it_l++) {
			const CGAL_Line_2 & l = (*it_l)->line;
			S[++entry] = (l.a() * O.x() + l.b() * O.y() + l.c() > 0);
		}

		return S;
	}


	Signature Polygon_Node::get_signature() const
	{
		return signature;
	}


	Polygon* Polygon_Node::get_one_polygon() const
	{
		return polygons.front();
	}


	void Polygon_Node::push(Polygon* P)
	{
		polygons.push_back(P);
		P->set_cell(this);
	}


	size_t Polygon_Node::size() const
	{
		return polygons.size();
	}


	std::list<Polygon*>::const_iterator Polygon_Node::polygons_begin() const
	{
		return polygons.cbegin();
	}


	std::list<Polygon*>::const_iterator Polygon_Node::polygons_end() const
	{
		return polygons.cend();
	}


	bool Polygon_Node::contours_are_empty() const
	{
		return contours.empty();
	}


	void Polygon_Node::set_contours(const std::list<Constraint> & _contours)
	{
		// The node is delimited by the contours passed as argument
		contours = std::vector<Constraint>(_contours.begin(), _contours.end());

		// For now, we consider that these contours also delimit the borders
		// of the group of adjacent nodes in which it is contained
		size_t n = contours.size();
		borders = new bool[n];
		for (size_t i = 0; i < n; ++i) borders[i] = true;
	}


	size_t Polygon_Node::contours_size() const
	{
		return contours.size();
	}


	Constraint Polygon_Node::get_contour(size_t i) const
	{
		assert(i < contours.size());
		return contours[i];
	}


	void Polygon_Node::set_contour_as_not_on_border(size_t i)
	{
		borders[i] = false;
	}


	void Polygon_Node::set_contour_as_not_on_border(Intersection_Line* I)
	{
		size_t n = contours.size();

		for (size_t i = 0 ; i < n ; ++i) {
			if (contours[i].first == I) {
				borders[i] = false;
				return;
			}
		}
	}


	void Polygon_Node::append_contours(std::list<std::tuple<Constraint, Constraint, Constraint> > & T) const
	{
		// This function is called as Nodes have been inserted in Groups.
		// By constructing groups, we know that some edges delimiting the nodes
		// cannot also define the group. Here, we just collect the remaining elements.

		size_t n = contours.size();

		for (size_t i = 0 ; i < n ; ++i) {
			if (borders[i]) {
				size_t i_prev = (i == 0 ? n - 1 : i - 1);
				size_t i_next = (i == n - 1 ? 0 : i + 1);
				T.push_back(std::make_tuple(contours[i_prev], contours[i], contours[i_next]));
			}
		}
	}


	void Polygon_Node::insert_link(Polygon_Link* l)
	{
		links.push_back(l);
	}


	void Polygon_Node::clear_links()
	{
		links.clear();
	}



	bool Polygon_Node::is_visited() const
	{
		return visited_node;
	}


	void Polygon_Node::set_visited()
	{
		visited_node = true;
	}


	bool Polygon_Node::is_queued() const
	{
		return queued_node;
	}


	void Polygon_Node::set_queued()
	{
		queued_node = true;
	}


	std::list<Polygon_Link*>::const_iterator Polygon_Node::links_begin() const
	{
		return links.cbegin();
	}


	std::list<Polygon_Link*>::const_iterator Polygon_Node::links_end() const
	{
		return links.cend();
	}
}