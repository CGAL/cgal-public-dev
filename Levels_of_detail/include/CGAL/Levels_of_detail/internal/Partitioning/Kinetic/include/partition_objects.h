#pragma once
#include "defs.h"
#include "support_plane_objects.h"

namespace Skippy {
	class Partition_Vertex;
	class Partition_Edge;
	class Partition_Facet;
	class Partition_Polyhedron;

	typedef std::pair<Partition_Edge*, bool> Partition_HalfEdge;
	typedef std::pair<Partition_Facet*, bool> Partition_Side;


	class Partition_Vertex
	{
	public:
		Partition_Vertex(const CGAL_Point_3 & _M, std::list<int> & P);

		~Partition_Vertex();

		bool definitions_match(const std::list<int> & Q) const;

		bool definitions_partially_match(const std::list<int> & Q) const;

		void get_definition(std::list<int> & P) const;

		void switch_to_plane(const int i, const int j);

		bool belongs_to_plane(const int id) const;

		bool belongs_to_planes(const int i, const int j) const;

		int get_local_id(const int id) const;

		void push(Partition_Edge* e);

		void pop(Partition_Edge* e);



		KINETIC_PARTITION_API std::map<int, int>::const_iterator local_ids_begin() const;

		KINETIC_PARTITION_API std::map<int, int>::const_iterator local_ids_end() const;

		KINETIC_PARTITION_API std::map<int, int>::iterator local_ids_begin();

		KINETIC_PARTITION_API std::map<int, int>::iterator local_ids_end();

		KINETIC_PARTITION_API size_t local_ids_size() const;

		KINETIC_PARTITION_API std::list<Partition_Edge*>::const_iterator edges_begin() const;

		KINETIC_PARTITION_API std::list<Partition_Edge*>::const_iterator edges_end() const;

		KINETIC_PARTITION_API std::list<Partition_Edge*>::iterator edges_begin();

		KINETIC_PARTITION_API std::list<Partition_Edge*>::iterator edges_end();

		Partition_Edge* edges_front() const;

		Partition_Edge* edges_back() const;

		int connectivity() const;

	public:
		bool delimits_duplicate_paths(const int slice_id, std::list< std::pair<Partition_Edge*, Partition_Edge*> > & D);

		Partition_Edge* prolongate(Partition_Edge* e);

	public:
		int id;
		const CGAL_Point_3 M;
		const CGAL_Inexact_Point_3 hint_M;

	private:
		std::map<int, int> local_ids;
		std::list<Partition_Edge*> edges;
	};


	inline std::ostream & operator<< (std::ostream & os, Partition_Vertex & v)
	{
		os << "[" << v.id << "] M : (" << v.M.x() << ", " << v.M.y() << ", " << v.M.z() << "), P : {";
		std::map<int, int>::iterator it_l_id = v.local_ids_begin();
		while (it_l_id != v.local_ids_end()) {
			os << "(" << it_l_id->first << ", " << it_l_id->second << ")";
			if (++it_l_id != v.local_ids_end()) os << ", ";
		}
		os << " }";
		return os;
	}



	class Partition_Edge
	{
	public:
		Partition_Edge(Partition_Vertex* _v1, Partition_Vertex* _v2);

		Partition_Edge(Partition_Vertex* _v1, Partition_Vertex* _v2, int F_i, int F_j);

		~Partition_Edge();

		int get_index_of_another_plane(const int excluded_id) const;

		int get_local_id(const int id) const;

		int is_edge_of_another_bounding_facet(int F_i) const;

		KINETIC_PARTITION_API bool reaches(Partition_Vertex* v) const;

		KINETIC_PARTITION_API Partition_Vertex* source(bool v1_v2) const;

		KINETIC_PARTITION_API Partition_Vertex* target(bool v1_v2) const;

		KINETIC_PARTITION_API Partition_Vertex* second_vertex(Partition_Vertex* v) const;

		void push(Partition_Facet* f);

		void pop(Partition_Facet* f);

		KINETIC_PARTITION_API std::map<int, int>::const_iterator local_ids_begin() const;

		KINETIC_PARTITION_API std::map<int, int>::const_iterator local_ids_end() const;

		KINETIC_PARTITION_API std::map<int, int>::iterator local_ids_begin();

		KINETIC_PARTITION_API std::map<int, int>::iterator local_ids_end();

		KINETIC_PARTITION_API size_t local_ids_size() const;

		KINETIC_PARTITION_API std::list<Partition_Facet*>::const_iterator facets_begin() const;

		KINETIC_PARTITION_API std::list<Partition_Facet*>::const_iterator facets_end() const;

		KINETIC_PARTITION_API std::list<Partition_Facet*>::iterator facets_begin();

		KINETIC_PARTITION_API std::list<Partition_Facet*>::iterator facets_end();

		void built_while_processing_negative_side_of_slicing_plane();

		bool get_negative_side_of_slicing_plane() const;

		bool belongs_to_same_planes(Partition_Edge* e) const;

		void substitute_in_facets(std::list<Partition_Edge*> S);

		static void definition_at_intersection(Partition_Edge* e1, Partition_Edge* e2, std::list<int> & D);

	public:
		int id;

	protected:
		Partition_Vertex* v1;
		Partition_Vertex* v2;
		std::map<int, int> local_ids;

		std::list<Partition_Facet*> facets;
		bool negative_side_of_slicing_plane;
	};


	class Partition_Facet
	{
	public:
		Partition_Facet(const int P, const std::list<Partition_Edge*> & E);

		~Partition_Facet();

		Partition_Vertex* get_projectible_vertex(Partition_Edge* e) const;

		void push(Partition_Edge* e);

		void pop(Partition_Edge* e);

		void push(Partition_Polyhedron* P);

		void pop(Partition_Polyhedron* P);

		KINETIC_PARTITION_API void get_circular_sequence_of_vertices(const std::vector<CGAL_Plane> & P, std::list<Partition_Vertex*> & V, bool orientation) const;

		KINETIC_PARTITION_API void get_circular_sequence_of_halfedges(const std::vector<CGAL_Plane> & P, std::list<Partition_HalfEdge> & H, bool orientation) const;

		KINETIC_PARTITION_API std::list<Partition_Edge*>::const_iterator edges_begin() const;

		KINETIC_PARTITION_API std::list<Partition_Edge*>::const_iterator edges_end() const;

		KINETIC_PARTITION_API std::list<Partition_Edge*>::iterator edges_begin();

		KINETIC_PARTITION_API std::list<Partition_Edge*>::iterator edges_end();

		KINETIC_PARTITION_API Partition_Polyhedron* get_polyhedron_1() const;

		KINETIC_PARTITION_API Partition_Polyhedron* get_polyhedron_2() const;

	public:
		int id;
		int p;

	protected:
		std::list<Partition_Edge*> edges;
		Partition_Polyhedron* polyhedron_1;
		Partition_Polyhedron* polyhedron_2;
	};


	inline std::ostream & operator<< (std::ostream & os, Partition_Facet & f)
	{
		os << "** Face " << f.id << std::endl;

		std::list<Partition_Vertex*> V;
		for (std::list<Partition_Edge*>::const_iterator it_e = f.edges_begin(); it_e != f.edges_end(); it_e++) {

			std::list<Partition_Vertex*>::iterator it_v;
			for (it_v = V.begin(); it_v != V.end(); it_v++) {
				if ((*it_v) == (*it_e)->source(true)) {
					break;
				}
			}
			if (it_v == V.end()) V.push_back((*it_e)->source(true));

			for (it_v = V.begin(); it_v != V.end(); it_v++) {
				if ((*it_v) == (*it_e)->target(true)) {
					break;
				}
			}
			if (it_v == V.end()) V.push_back((*it_e)->target(true));
		}

		for (std::list<Partition_Vertex*>::iterator it_v = V.begin(); it_v != V.end(); it_v++) {
			os << *(*it_v) << std::endl;
		}
		return os;
	}


	class Partition_Polyhedron
	{
	public:
		Partition_Polyhedron(const std::set<Partition_Side> & P);

		~Partition_Polyhedron();

		KINETIC_PARTITION_API std::list<Partition_Side>::const_iterator facets_begin() const;

		KINETIC_PARTITION_API std::list<Partition_Side>::const_iterator facets_end() const;

		KINETIC_PARTITION_API std::list<Partition_Side>::iterator facets_begin();

		KINETIC_PARTITION_API std::list<Partition_Side>::iterator facets_end();

		KINETIC_PARTITION_API size_t facets_size() const;

		void reverse_sides(const int id);

	public:
		int id;

	protected:
		std::list<Partition_Side> facets;
	};


	inline bool sort_by_vertex_id(Partition_Vertex* v1, Partition_Vertex* v2) { return (v1->id < v2->id); }
	inline bool sort_by_x_coordinate(Partition_Vertex* v1, Partition_Vertex* v2) { return (v1->M.x() < v2->M.x()); }
	inline bool sort_by_y_coordinate(Partition_Vertex* v1, Partition_Vertex* v2) { return (v1->M.y() < v2->M.y()); }
	inline bool sort_by_z_coordinate(Partition_Vertex* v1, Partition_Vertex* v2) { return (v1->M.z() < v2->M.z()); }
}