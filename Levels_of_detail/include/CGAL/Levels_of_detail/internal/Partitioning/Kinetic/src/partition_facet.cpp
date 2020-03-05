#include "../include/partition_objects.h"
#include "../include/vars.h"
#include "../include/universe.h"
#include "../include/support_plane.h"


namespace Skippy 
{
	Partition_Facet::Partition_Facet(const int _P, const std::list<Partition_Edge*> & E)
		: id(++Counters::id_partition_facet),
		p(_P),
		polyhedron_1(nullptr),
		polyhedron_2(nullptr)
	{
		edges.clear();
		std::copy(E.begin(), E.end(), std::back_inserter(edges));

		for (std::list<Partition_Edge*>::iterator it_e = edges.begin(); it_e != edges.end(); it_e++) {
			(*it_e)->push(this);
		}
	}


	Partition_Facet::~Partition_Facet()
	{
		for (std::list<Partition_Edge*>::iterator it_e = edges.begin(); it_e != edges.end(); it_e++) {
			(*it_e)->pop(this);
		}
	}


	Partition_Vertex* Partition_Facet::get_projectible_vertex(Partition_Edge* e) const
	{
		// We assume that e is contained in the definition of the current facet
		// We search for a vertex that is not in e = (v1 v2) and that is not aligned with e

		Partition_Vertex* v_1 = e->source(true);
		Partition_Vertex* v_2 = e->target(true);

		const CGAL_Point_3 & M_1 = v_1->M, &M_2 = v_2->M;

		for (std::list<Partition_Edge*>::const_iterator it_e = edges.begin(); it_e != edges.end(); it_e++) {
			Partition_Edge* e_curr = (*it_e);
			if (e_curr != e) {
				Partition_Vertex* v_3 = e_curr->source(true);
				Partition_Vertex* v_4 = e_curr->target(true);

				if (v_3 != v_1 && v_3 != v_2) {
					if (!CGAL::collinear(M_1, M_2, v_3->M)) return v_3;
				}

				if (v_4 != v_1 && v_4 != v_2) {
					if (!CGAL::collinear(M_1, M_2, v_4->M)) return v_4;
				}
			}
		}

		return nullptr;
	}


	void Partition_Facet::push(Partition_Edge* e)
	{
		edges.push_back(e);
		e->push(this);
	}


	void Partition_Facet::pop(Partition_Edge* e)
	{
		std::list<Partition_Edge*>::iterator it_e = edges.begin();
		while (it_e != edges.end()) {
			if ((*it_e) == e) {
				break;
			} else {
				++it_e;
			}
		}
		assert(it_e != edges.end());
		it_e = edges.erase(it_e);
	}


	void Partition_Facet::push(Partition_Polyhedron* P)
	{
		if (polyhedron_1 == nullptr) {
			polyhedron_1 = P;
		} else {
			assert(polyhedron_2 == nullptr);
			polyhedron_2 = P;
		}
	}


	void Partition_Facet::pop(Partition_Polyhedron* P)
	{
		if (polyhedron_1 == P) {
			polyhedron_1 = nullptr;
		} else {
			assert(polyhedron_2 == P);
			polyhedron_2 = nullptr;
		}
	}


	void Partition_Facet::get_circular_sequence_of_vertices(const std::vector<CGAL_Plane> & P, std::list<Partition_Vertex*> & V, bool orientation) const
	{
		V.clear();
		std::map<Partition_Vertex*, std::pair<Partition_Edge*, Partition_Edge*> > adjacent_edges;

		// Part 1.
		// Loops on the sequence of edges and assigns two adjacent edges to each vertex

		for (std::list<Partition_Edge*>::const_iterator it_e = edges.begin(); it_e != edges.end(); it_e++) {
			Partition_Edge* e = (*it_e);

			Partition_Vertex *v_s = e->source(true), *v_t = e->target(true);
			std::map<Partition_Vertex*, std::pair<Partition_Edge*, Partition_Edge*> >::iterator it_v;

			it_v = adjacent_edges.find(v_s);
			if (it_v == adjacent_edges.end()) {
				// If the vertex hasn't been encountered yet
				adjacent_edges[v_s] = std::make_pair(e, nullptr);
			} else {
				// If the vertex has already been met once
				it_v->second.second = e;
			}

			it_v = adjacent_edges.find(v_t);
			if (it_v == adjacent_edges.end()) {
				// If the vertex hasn't been encountered yet
				adjacent_edges[v_t] = std::make_pair(e, nullptr);
			} else {
				// If the vertex has already been met once
				it_v->second.second = e;
			}
		}

		// Part 2.
		// Loops on vertices

		Partition_Vertex* v_init = adjacent_edges.begin()->first, *v_curr = v_init;
		Partition_Edge* e_prev = nullptr;

		do {
			// Gets the first of the two edges associated to v_curr
			// Compares it to the last processed edge, to loop on a new vertex
			Partition_Edge* e = adjacent_edges[v_curr].first;
			if (e == e_prev) {
				e = adjacent_edges[v_curr].second;
			}

			// Iterates
			v_curr = e->second_vertex(v_curr);
			V.push_back(v_curr);
			e_prev = e;

			if (V.size() > (1 << 20)) {
				throw std::logic_error("Alert : infinite loop in get_circular_sequence_of_vertices");
			}

		} while (v_curr != v_init);

		// Part 3.
		// Determines if the sequence of vertices should been reverted, to match the desired orientation

		// Indeed, if v1, v2 and v3 represent the three first vertices of V,
		// The variable orientation is set to true v_12 ^ v_13 has the same direction as the support plane's normal vector

		std::list<Partition_Vertex*>::iterator it_v = V.begin();
		Partition_Vertex *v1 = (*it_v); ++it_v;
		Partition_Vertex *v2 = (*it_v);
		Partition_Vertex *v3 = nullptr;

		CGAL_Vector_3 v_12 = v2->M - v1->M;

		while (++it_v != V.end()) {
			// Chooses v3 as the first non colinear vertex found
			v3 = (*it_v);
			if (!CGAL::collinear(v1->M, v2->M, v3->M)) break;
		}
		assert(v3 != nullptr);

		CGAL_Vector_3 v_13 = v3->M - v1->M;
		CGAL_Vector_3 v_n = CGAL::cross_product(v_12, v_13);

		const CGAL_Plane & H_p = P[this->p];
		CGAL_Vector_3 N(H_p.a(), H_p.b(), H_p.c());
		FT lambda = (CGAL::abs(N.x()) > 0 ? v_n.x() / N.x() : (CGAL::abs(N.y()) > 0 ? v_n.y() / N.y() : v_n.z() / N.z()));
		assert(lambda != 0);

		if (lambda > 0) {
			// v_n has the same orientation as the normal
			if (!orientation) V.reverse();
		} else {
			// v_n has the opposite orientation to the normal
			if (orientation) V.reverse();
		}

		for (std::list<Partition_Vertex*>::iterator it_v = V.begin(); it_v != V.end(); it_v++) {
			assert((*it_v)->id >= 0);
		}
	}



	void Partition_Facet::get_circular_sequence_of_halfedges(const std::vector<CGAL_Plane> & P, std::list<Partition_HalfEdge> & H, bool orientation) const
	{
		H.clear();

		std::list<Partition_Vertex*> V;
		get_circular_sequence_of_vertices(P, V, orientation);

		std::list<Partition_Vertex*>::iterator it_v1 = V.begin(), it_v2 = ++V.begin();
		while (it_v1 != V.end()) {
			Partition_Vertex *v1 = (*it_v1), *v2 = (*it_v2);

			// Finds the edge between v1 and v2...
			for (std::list<Partition_Edge*>::iterator it_e = v1->edges_begin(); it_e != v1->edges_end(); ++it_e) {
				Partition_Edge* e = (*it_e);
				if (e->reaches(v2)) {
					// Then the halfedge
					bool sign = (e->source(true) == v1);
					H.push_back(std::make_pair(e, sign));
					break;
				}
			}

			// Iterates
			++it_v1;
			++it_v2;
			if (it_v2 == V.end()) it_v2 = V.begin();
		}
	}



	std::list<Partition_Edge*>::const_iterator Partition_Facet::edges_begin() const
	{
		return edges.cbegin();
	}


	std::list<Partition_Edge*>::const_iterator Partition_Facet::edges_end() const
	{
		return edges.cend();
	}


	std::list<Partition_Edge*>::iterator Partition_Facet::edges_begin()
	{
		return edges.begin();
	}


	std::list<Partition_Edge*>::iterator Partition_Facet::edges_end()
	{
		return edges.end();
	}


	Partition_Polyhedron* Partition_Facet::get_polyhedron_1() const
	{
		return polyhedron_1;
	}


	Partition_Polyhedron* Partition_Facet::get_polyhedron_2() const
	{
		return polyhedron_2;
	}
}