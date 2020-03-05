#include "../include/partition_objects.h"
#include "../include/vars.h"

namespace Skippy {
	using CGAL::to_double;


	Partition_Vertex::Partition_Vertex(const CGAL_Point_3 & _M, std::list<int> & P)
		: M(_M),
		hint_M(CGAL_Inexact_Point_3(to_double(M.x()), to_double(M.y()), to_double(M.z()))),
		id(++Counters::id_partition_vertex)
	{
		for (std::list<int>::iterator it_p = P.begin(); it_p != P.end(); it_p++) {
			local_ids[*it_p] = ++Counters::par_v_local_ids[*it_p];
		}
	}



	Partition_Vertex::~Partition_Vertex()
	{

	}



	bool Partition_Vertex::definitions_match(const std::list<int> & Q) const
	{
		// We assume that elements of Q are sorted
		// That's why this test consists in looping on all elements of local_ids and Q
		// and making sure that local_ids->first == q at each iteration
		if (local_ids.size() != Q.size()) return false;

		std::map<int, int>::const_iterator it_p = local_ids.begin();
		std::list<int>::const_iterator it_q = Q.begin();

		while (it_p != local_ids.end() && it_q != Q.end()) {
			if (it_p->first != (*it_q)) return false;
			++it_p, ++it_q;
		}

		return true;
	}



	bool Partition_Vertex::definitions_partially_match(const std::list<int> & Q) const
	{
		// We assume that elements of Q are sorted
		// We loop on elements of local_ids and Q and check that planes listed in Q are represented in local_ids


		std::list<int>::const_iterator it_q = Q.begin();
		for (std::map<int, int>::const_iterator it_p = local_ids.begin(); it_p != local_ids.end(); it_p++) {
			if (it_p->first == *it_q) {
				++it_q;
			}
		}

		return (it_q == Q.end());
	}



	void Partition_Vertex::get_definition(std::list<int> & P) const
	{
		for (std::map<int, int>::const_iterator it_p = local_ids.begin(); it_p != local_ids.end(); it_p++) {
			P.push_back(it_p->first);
		}
	}


	void Partition_Vertex::switch_to_plane(const int i, const int j)
	{
		int index = -1;

		std::map<int, int>::iterator it_p = local_ids.begin();
		while (it_p != local_ids.end()) {
			if (it_p->first == i) {
				index = it_p->second;
				break;
			} else {
				++it_p;
			}
		}

		local_ids.erase(it_p);
		local_ids[j] = index;
	}


	void Partition_Vertex::push(Partition_Edge* e)
	{
		edges.push_back(e);
	}



	void Partition_Vertex::pop(Partition_Edge* e)
	{
		std::list<Partition_Edge*>::iterator it_e;
		for (it_e = edges.begin(); it_e != edges.end(); it_e++) {
			if ((*it_e) == e) break;
		}
		assert(it_e != edges.end());
		it_e = edges.erase(it_e);
	}



	bool Partition_Vertex::belongs_to_plane(const int id) const
	{
		return local_ids.find(id) != local_ids.end();
	}



	bool Partition_Vertex::belongs_to_planes(const int i, const int j) const
	{
		return (local_ids.find(i) != local_ids.end() && local_ids.find(j) != local_ids.end());
	}



	int Partition_Vertex::get_local_id(const int id) const
	{
		std::map<int, int>::const_iterator it = local_ids.find(id);
		return (it != local_ids.end() ? it->second : -1);
	}



	std::map<int, int>::const_iterator Partition_Vertex::local_ids_begin() const
	{
		return local_ids.cbegin();
	}


	std::map<int, int>::const_iterator Partition_Vertex::local_ids_end() const
	{
		return local_ids.cend();
	}


	std::map<int, int>::iterator Partition_Vertex::local_ids_begin()
	{
		return local_ids.begin();
	}


	std::map<int, int>::iterator Partition_Vertex::local_ids_end()
	{
		return local_ids.end();
	}


	size_t Partition_Vertex::local_ids_size() const
	{
		return local_ids.size();
	}


	std::list<Partition_Edge*>::const_iterator Partition_Vertex::edges_begin() const
	{
		return edges.cbegin();
	}


	std::list<Partition_Edge*>::const_iterator Partition_Vertex::edges_end() const
	{
		return edges.cend();
	}


	std::list<Partition_Edge*>::iterator Partition_Vertex::edges_begin()
	{
		return edges.begin();
	}


	std::list<Partition_Edge*>::iterator Partition_Vertex::edges_end()
	{
		return edges.end();
	}


	Partition_Edge* Partition_Vertex::edges_front() const
	{
		return edges.front();
	}


	Partition_Edge* Partition_Vertex::edges_back() const
	{
		return edges.back();
	}


	int Partition_Vertex::connectivity() const
	{
		return int(edges.size());
	}


	bool Partition_Vertex::delimits_duplicate_paths(const int slice_id, std::list< std::pair<Partition_Edge*, Partition_Edge*> > & D)
	{
		D.clear();

		// This function is called as we exhibit duplicate paths of edges on the slicing plane slice_id.

		std::map<int, bool> processed_edges;
		for (std::list<Partition_Edge*>::const_iterator it_e = edges.begin() ; it_e != edges.end() ; ++it_e) {
			processed_edges[(*it_e)->id] = false;
		}

		std::list<Partition_Edge*>::const_iterator it_e1, it_e2;

		for (it_e1 = edges.begin() ; it_e1 != edges.end() ; ++it_e1) {
			Partition_Edge* e1 = (*it_e1);
			processed_edges[e1->id] = true;

			// If e1 is not on the slicing plane, passes.
			if (e1->get_local_id(slice_id) == -1) continue;

			it_e2 = it_e1;
			++it_e2;

			// We loop on the other edges of this object.
			// We are looking for one that results from the intersection of the same planes as e1.

			while (it_e2 != edges.end()) {
				Partition_Edge* e2 = (*it_e2);
				if (e2->get_local_id(slice_id) == -1 || !e2->belongs_to_same_planes(e1)) {
					++it_e2;
					continue;
				}

				// Are e1 and e2 duplicate ?
				// We compute the dot product of the vectors OM_1 and OM_2.
				// If it is positive, then e1 and e2 are duplicate.

				CGAL_Vector_3 OM_1 = e1->second_vertex(this)->M - M;
				CGAL_Vector_3 OM_2 = e2->second_vertex(this)->M - M;
				if (OM_1 * OM_2 > 0) {
					D.push_back(std::make_pair(e1, e2));
					processed_edges[e2->id] = true;
					break;
				}

				++it_e2;
			}
		}

		return (!D.empty());
	}



	Partition_Edge* Partition_Vertex::prolongate(Partition_Edge* e_ref)
	{
		for (std::list<Partition_Edge*>::const_iterator it_e = edges.begin() ; it_e != edges.end() ; ++it_e) {
			Partition_Edge* e = (*it_e);
			if (e_ref != e && e_ref->belongs_to_same_planes(e)) {

				// We must test if e_ref and e are duplicate.
				// We set M as reference point, and measure the angle made by e and e_ref.
				// If the inner product of the edges is negative, then e and e_ref go to opposite directions
				// which is what we want...

				CGAL_Vector_3 u_ref = e_ref->second_vertex(this)->M - M;
				CGAL_Vector_3 u = e->second_vertex(this)->M - M;
				if (u * u_ref < 0) {
					return e;
				}
			}
		}

		// Default result
		return nullptr;
	}
}