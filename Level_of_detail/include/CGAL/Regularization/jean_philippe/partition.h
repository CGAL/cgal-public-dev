#pragma once
#pragma once
#include "parameters.h"
#include "partition_elements.h"
#include "indexed_event.h"
#include "quadtree.h"
#include "matrix.h"
#include <set>

using std::list;
using std::pair;
using std::vector;
using std::set;

using cv::Vec3b;
using cv::Vec3d;


class Partition
{
public:
	Partition();

    Partition(uint _rows, uint _cols, Parameters* _params);

	Partition(const Partition & P);

	Partition & operator= (const Partition & P);

	~Partition();

private:
	void clear();
	
	void copy(const Partition & P);

public:
	int & get_id_vertices();

	void init_edges(vector<Ray *> & rays, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
		vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges, vector<Vertex *> & initial_vertices_of_rays);

	void build_edge(vector<Ray *> & rays, IndexedEvent* current_event, Vertex* intersection_point, bool vertex_is_corner, bool vertex_is_new, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
		vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges, vector<Vertex *> & initial_vertices_of_rays);

	void build_faces(Size2i & size);

	void merge_containers(list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices, vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges);

	void merge_thin_facets(vector<Ray *> & rays, int verbose_level);

	void find_pixels_inside_facets();

	void seal_and_remove_bivalent_vertices();

	void reset_indices();

	void split_facet(Face* f, double a, double b, double c);

private:
	Face* third_facet(Vertex* v_i, HalfEdge* h_ij, Face* f_i, Face* f_j);
	
	void determine_best_edge_for_joining(Face* f_i, list<pair<uint, uint> > & possibilities, pair<uint, uint> & selected);

	void join_thin_facet(map<int, list<Face *>::iterator> & map_facets, Face* & f_i, list<pair<uint, uint> > & edges_joinable_facets);

	void merge_two_facets(map<int, list<Face *>::iterator> & map_facets, Face* & f_i, Face* & f_j, pair<uint, uint> & h_ij,
		bool exclude_f_i_from_list, bool exclude_f_j_from_list);

	void determine_best_edge_for_slicing(vector<HalfEdge *> & f_i_edges, vector<pair<Vertex *, Vertex_Type> > & f_i_vertices, 
		list<pair<uint, uint> > & possibilites, pair<uint, uint> & selected);

	void determine_intersection_point_when_slicing(Vertex* v_curr, Vec2d & sliced_edge_direction, list<HalfEdge *> & intersectable_halfedges, 
	IndexedEvent* & intersection, HalfEdge* & extended_halfedge, HalfEdge* & intersected_halfedge);

	void update_graph_when_slicing(vector<Ray *> & rays, list<HalfEdge *> & halfedges_to_trivalent_vertex, Vertex* trivalent_vertex, IndexedEvent* event,
		HalfEdge* extended_halfedge, HalfEdge* intersected_halfedge, list<HalfEdge *> & remaining_halfedges_of_f_i,
		map<HalfEdge*, list<HalfEdge *>::iterator> & access_remaining_halfedges_of_f_i, list<HalfEdge *> & sub_facet_f_i);
		
	void slice_thin_facet(vector<Ray *> & rays, map<int, list<Face *>::iterator> & map_facets, Face* & f_i, list<pair<uint, uint> > & edges_to_slice);

	void loop_inside_sub_facet(HalfEdge* h_e, list<HalfEdge*> & sub_facet);

	void classify_along_sliced_edge(vector<HalfEdge *> & big_sliced_edge, vector<list<HalfEdge *> > & sequences_halfedges_to_trivalent, vector<list<Vertex *> > & sequences_bivalent_vertices);

	void split_edge_and_update_adjacent_facet(HalfEdge* h, Vertex* v);
public:
    //void fill_matrix_of_iterators();

    int count_invalid_faces();

    void save_edges(Matrix<uchar> & I, std::string & directory, std::string & basename, double f = 1.0);

	void save_faces(std::string & directory, std::string & basename, bool color_after_thinness);

    void save_graph_definition(std::string & directory, std::string & basename);

    void save_boundaries(std::string & directory, std::string & basename);

	void save_parameters(std::string &directory, std::string &basename, int N);

    void save_labels(std::string & directory, std::string & basename);

    void save_liuyuns_input(std::string & directory, std::string & basename, vector<Segment *> & segments);

	//void read_graph_definition(std::string & filename);

    //void stats();


private:
	Vertex* erase(Vertex* v);

	Edge* erase(Edge* e, bool destroy);

	void erase(list<Edge *> & l_e, bool destroy);

	void push_back(Vertex* v);

	void push_back(Edge* e);

public:
    void draw_edges(Matrix<uchar> & I, Matrix<uchar> & J, double f = 1.0);

	void draw_faces(Matrix<uchar> & J, bool color_after_thinness);

	void draw_intermediate_graph(Matrix<uchar> & I, Matrix<uchar> & J, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
		vector<list<Edge *> > & outer_edges, vector<list<Edge *> > & inner_edges);

public:
	uint rows;
	uint cols;

    Parameters* params;

	Quadtree* quadtree;

	bool is_valid;
	int id_vertices;
	int id_edges;
	int id_faces;

	Vertex* vertices_head;
	Vertex* vertices_tail;
	uint v_size;

	Edge* edges_head;
	Edge* edges_tail;
	uint e_size;

	list<Face *> faces;
    //Matrix<list<Face *>::iterator> F;
};


inline bool sort_by_thinness(pair<int, double> & first_pair, pair<int, double> & second_pair)
{
	return (first_pair.second < second_pair.second);
}