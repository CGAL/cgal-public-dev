#include "partition.h"
#include <iostream>
#include <fstream>
#include <queue>
#include <random>
#include <ctime>
#include <set>
#include <map>
#include "defs.h"
#include "trace.h"
#include "colormap.h"
#include "geometry.h"
#include "svg.h"
#include <cassert>

using std::set;
using std::map;


Partition::Partition() : 
	rows (0), 
	cols (0), 
	params (nullptr), 
	quadtree (nullptr), 
	is_valid(false), 
	id_vertices(0), 
	id_edges (0), 
	id_faces(0), 
	vertices_head (nullptr), 
	vertices_tail(nullptr), 
	v_size(0), 
	edges_head(nullptr), 
	edges_tail(nullptr), 
	e_size(0)
{
	faces = list<Face *>();
}


Partition::Partition(uint _rows, uint _cols, Parameters *_params) :
	rows (_rows),
	cols (_cols),
    params(_params)
{
	quadtree = new Quadtree(0, cols, 0, rows);
	
	is_valid = false;
	id_vertices = 0;
	id_edges = 0;
	id_faces = 0;

	vertices_head = NULL;
	vertices_tail = NULL;
	v_size = 0;

	edges_head = NULL;
	edges_tail = NULL;
	e_size = 0;

    faces = list<Face *>();
}


Partition::~Partition()
{
	clear();
}


void Partition::clear()
{
	// Erases all facets
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		delete (*it_f);
	}
	faces.clear();

	// Erases all edges
	Edge *e_curr = edges_head, *e_next = nullptr;
	while (e_curr != nullptr) {
		e_next = e_curr->e_next;
		delete e_curr;
		e_curr = e_next;
	}
	edges_head = edges_tail = nullptr;
	e_size = 0;
	
	// Erases all vertices and deletes quadtree
	Vertex *v_curr = vertices_head, *v_next = nullptr;
	while (v_curr != nullptr) {
		v_next = v_curr->v_next;
		delete v_curr;
		v_curr = v_next;
	}
	vertices_head = vertices_tail = nullptr;
	v_size = 0;

	if (quadtree != nullptr) delete quadtree;
	quadtree = nullptr;

	// params is set to nullptr (not destroyed here)
	params = nullptr;

	// Finally, resets the number of rows and columns
	rows = cols = 0;
}


void Partition::copy(const Partition & P)
{
	rows = P.rows;
	cols = P.cols;
	params = P.params;
	is_valid = P.is_valid;

	// Copies vertices
	// First initializes a quadtree
	quadtree = new Quadtree(0, cols, 0, rows);

	// Initializes a vector of vertices that we will later use to initialize our edges
	vector<Vertex *> v_vertices;
	v_vertices.reserve(P.v_size);

	id_vertices = 0;
	vertices_head = vertices_tail = nullptr;
	v_size = 0;

	Vertex* p_v = P.vertices_head;
	while (p_v != nullptr) {
		// Deep copy of V
		Vertex* v = new Vertex(id_vertices, p_v->pt.x, p_v->pt.y);
		for (list<IndexedEvent *>::iterator it_e = p_v->events.begin() ; it_e != p_v->events.end() ; it_e++) {
			v->events.push_back(new IndexedEvent(**it_e));
		}
		quadtree->add(v);
		push_back(v);
		v_vertices.push_back(v);

		// Iterates
		p_v = p_v->v_next;
	}

	// Initializes a vector of edges that we will later use to build copies of P's facets
	vector<Edge *> v_edges;
	v_edges.reserve(P.e_size);

	id_edges = 0;
	edges_head = edges_tail = nullptr;
	e_size = 0;

	Edge* p_e = P.edges_head;
	while (p_e != nullptr) {
		double alpha_v1 = p_e->v1->incidence_angle(p_e);
		double alpha_v2 = p_e->v2->incidence_angle(p_e);
		Vertex* v1 = v_vertices[p_e->v1->id_vertex];
		Vertex* v2 = v_vertices[p_e->v2->id_vertex];
		Edge* e = nullptr;
		if (p_e->type == INNER_EDGE) {
			Inner_Edge* casted_p_e = static_cast<Inner_Edge *>(p_e);
			set<Ray *> support;
			casted_p_e->get_supporting_rays(support);
			e = new Inner_Edge(id_edges, support, 0, v1, v2, alpha_v1, alpha_v2);
		} else if (p_e->type == OUTER_EDGE) {
			Outer_Edge* casted_p_e = static_cast<Outer_Edge *>(p_e);
			Image_Boundary boundary = casted_p_e->boundary;
			e = new Outer_Edge(id_edges, boundary, v1, v2, alpha_v1, alpha_v2);
		}
		push_back(e);
		v_edges.push_back(e);
		p_e = p_e->e_next;
	}

	// And now here come the facets
	id_faces = 0;
	for (list<Face *>::const_iterator it_f = P.faces.begin() ; it_f != P.faces.end() ; it_f++) {
		Face* p_f = (*it_f);
		list<HalfEdge *> & p_f_halfedges = p_f->edges;

		list<HalfEdge *> f_halfedges;
		for (list<HalfEdge *>::iterator it_h = p_f_halfedges.begin() ; it_h != p_f_halfedges.end() ; ++it_h) {
			HalfEdge* p_h = (*it_h);
			HalfEdge* h = (p_h->v1_v2 ? v_edges[p_h->e->id_edge]->v1_v2 : v_edges[p_h->e->id_edge]->v2_v1);
			f_halfedges.push_back(h);
		}
		faces.push_back(new Face(id_faces, f_halfedges));
	}

	v_edges.clear();
	v_vertices.clear();
}


Partition::Partition(const Partition & P)
{
	copy(P);
}


Partition & Partition::operator= (const Partition & P)
{
	if (&P != this) 
	{
		clear();
		copy(P);
	}
	return *this;
}


int & Partition::get_id_vertices()
{
	return id_vertices;
}



Vertex* Partition::erase(Vertex* v)
{	
	Vertex *v_prev = v->v_prev, *v_next = v->v_next;
	if (v_prev != NULL) {
		v_prev->v_next = v_next;
	} else {
		vertices_head = v_next;
	}
	if (v_next != NULL) {
		v_next->v_prev = v_prev;
	} else {
		vertices_tail = v_prev;
	}
	--v_size;
	return v_next;
}



Edge* Partition::erase(Edge* e, bool destroy)
{
	Edge *e_prev = e->e_prev, *e_next = e->e_next;
	if (e_prev != NULL) {
		e_prev->e_next = e_next;
	} else {
		edges_head = e_next;
	}
	if (e_next != NULL) {
		e_next->e_prev = e_prev;
	} else {
		edges_tail = e_prev;
	}
	if (destroy) delete e;
	--e_size;
	return e_next;
}



void Partition::erase(list<Edge *> & l_e, bool destroy)
{
	for (list<Edge *>::iterator it_e = l_e.begin(); it_e != l_e.end(); it_e++) erase(*it_e, destroy);
	l_e.clear();
}



void Partition::push_back(Vertex* v)
{
	if (vertices_head == NULL) {
		vertices_head = v;
		vertices_tail = v;
		v->v_prev = v->v_next = NULL;
	} else {
		Vertex* v_last = vertices_tail;
		v_last->v_next = v;
		v->v_prev = v_last;
		vertices_tail = v;
	}
	++v_size;
}



void Partition::push_back(Edge* e)
{
	if (edges_head == NULL) {
		edges_head = e;
		edges_tail = e;
		e->e_prev = e->e_next = NULL;
	} else {
		Edge* e_last = edges_tail;
		e_last->e_next = e;
		e->e_prev = e_last;
		edges_tail = e;
	}
	++e_size;
}



void Partition::init_edges(vector<Ray *> & rays, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
	vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges, vector<Vertex *> & initial_vertices_of_rays)
{
	// To init the graph, we first build a set of vertices that represent the corners of the image.
	Vertex* top_left = new Vertex(id_vertices, 0, rows);
	Vertex* top_right = new Vertex(id_vertices, cols, rows);
	Vertex* bottom_right = new Vertex(id_vertices, cols, 0);
	Vertex* bottom_left = new Vertex(id_vertices, 0, 0);

	outer_vertices.push_back(top_left);
	outer_vertices.push_back(top_right);
	outer_vertices.push_back(bottom_right);
	outer_vertices.push_back(bottom_left);
	for (list<Vertex *>::iterator it_v = outer_vertices.begin(); it_v != outer_vertices.end(); it_v++) quadtree->add(*it_v);

	// We also build the edges that connect such vertices.
	Outer_Edge* top = new Outer_Edge(id_edges, TOP_IMAGE, top_left, top_right, 0, PI, 0, cols);
	Outer_Edge* bottom = new Outer_Edge(id_edges, BOTTOM_IMAGE, bottom_left, bottom_right, 0, PI, 0, cols);
	Outer_Edge* left = new Outer_Edge(id_edges, LEFT_IMAGE, bottom_left, top_left, PI / 2, -PI / 2, 0, rows);
	Outer_Edge* right = new Outer_Edge(id_edges, RIGHT_IMAGE, bottom_right, top_right, PI / 2, -PI / 2, 0, rows);

	outer_edges.reserve(4);
	outer_edges.push_back(list<Outer_Edge *>(1, top));
	outer_edges.push_back(list<Outer_Edge *>(1, bottom));
	outer_edges.push_back(list<Outer_Edge *>(1, left));
	outer_edges.push_back(list<Outer_Edge *>(1, right));

	// Finally, we build initial vertices that correspond to the centers of the segments.
	inner_edges.reserve(rays.size());
	initial_vertices_of_rays.reserve(rays.size());
	for (int i = 0; i < rays.size() / 2; i++) {
		Point2d & O = rays[2 * i]->O;
		Vertex* v = new Vertex(id_vertices, O.x, O.y);
		//v->events.push_back(new IndexedEvent(2 * i, 2 * i + 1, -rays[2 * i]->initial_length, -rays[2 * i + 1]->initial_length, true));
		quadtree->add(v);
		inner_vertices.push_back(v);
		// Then we insert two empty lists in the vector of edges.
		// Indeed, each entry of this vector is associated to a ray.
		inner_edges.push_back(list<Inner_Edge *>());
		inner_edges.push_back(list<Inner_Edge *>());
		// Remember the vertices we created in a vector
		initial_vertices_of_rays.push_back(v);
		initial_vertices_of_rays.push_back(v);
	}
}



void Partition::build_edge(vector<Ray *> & rays, IndexedEvent* current_event, Vertex* vertex, bool vertex_is_corner, bool vertex_is_new, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
	vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges, vector<Vertex *> & initial_vertices_of_rays)
{
	Ray* r_i = rays[current_event->intersectant];
	Ray* r_j = (current_event->intersected >= 0 ? rays[current_event->intersected] : NULL);

	double abscissa;

	if (!current_event->is_colliding_ray && !vertex_is_corner && vertex_is_new) {

		// On one hand, if the interesected ray is an image boundary
		outer_vertices.push_back(vertex);

		// First, we determine which edge must be split
		Image_Boundary intersected_boundary = Image_Boundary(current_event->intersected);
		if (intersected_boundary == TOP_IMAGE || intersected_boundary == BOTTOM_IMAGE) {
			abscissa = vertex->pt.x;
		} else if (intersected_boundary == LEFT_IMAGE || intersected_boundary == RIGHT_IMAGE) {
			abscissa = vertex->pt.y;
		}

		// We loop on the list of edges associated to this boundary
		list<Outer_Edge *> & edges_for_this_corner = outer_edges[-1 - current_event->intersected];
		list<Outer_Edge *>::iterator it1 = edges_for_this_corner.begin();
		while (it1 != edges_for_this_corner.end()) {
			double t1, t2;
			(*it1)->time_range(t1, t2);
			if (t1 <= abscissa && abscissa < t2) break;
			// if ((*it1)->t1 <= abscissa && abscissa < (*it1)->t2) break;
			it1++;
		}

		// The loop has stopped on the edge to split
		Outer_Edge* split_edge = *it1;
		Vertex* v1 = split_edge->v1, *v2 = split_edge->v2;
		double t1, t2;
		split_edge->time_range(t1, t2);
		//double t1 = split_edge->t1, t2 = split_edge->t2;
		double alpha_1 = v1->incidence_angle(split_edge), alpha_2 = v2->incidence_angle(split_edge);
		Outer_Edge* v1_v = new Outer_Edge(id_edges, intersected_boundary, v1, vertex, alpha_1, alpha_2, t1, abscissa);
		Outer_Edge* v_v2 = new Outer_Edge(id_edges, intersected_boundary, vertex, v2, alpha_1, alpha_2, abscissa, t2);

		it1 = edges_for_this_corner.erase(it1);
		edges_for_this_corner.insert(it1, v1_v);
		edges_for_this_corner.insert(it1, v_v2);
		delete split_edge;

	} else if (current_event->is_colliding_ray) {

		// On the other hand, if the intersected ray is a real ray
		if (vertex_is_new && !vertex_is_corner) inner_vertices.push_back(vertex);

		abscissa = current_event->t_intersected;
		list<Inner_Edge *> & edges_for_this_intersected_ray = inner_edges[current_event->intersected];
		if (edges_for_this_intersected_ray.size() == 0) {

			// In case no edge modelizes the path followed by the intersected ray, we create one
			if (initial_vertices_of_rays[current_event->intersected] != vertex) {
				//int tag = (current_event->t_intersected > r_j->t_swap ? 1 : 0);
				double alpha_1 = r_j->opposite()->incidence_angle, alpha_2 = r_j->incidence_angle;

				Inner_Edge* v1_v = new Inner_Edge(id_edges, r_j, 0, initial_vertices_of_rays[current_event->intersected], vertex, alpha_1, alpha_2, -r_j->initial_length, abscissa);
				edges_for_this_intersected_ray.push_back(v1_v);
			}

		} else {

			// Loops on the list of edges associated to the intersected ray
			// Exhibits the one that contains the intersection point, if it does exist
			Vertex* v1 = NULL;
			Vertex* v2 = NULL;
			double last_abscissa = -FLT_MAX;
			list<Inner_Edge *>::iterator it2 = edges_for_this_intersected_ray.begin();
			while (it2 != edges_for_this_intersected_ray.end()) {
				double t1, t2;
				(*it2)->time_range(r_j, t1, t2);
				if (t1 <= abscissa && abscissa < t2) break;
				// if ((*it2)->t1 <= abscissa && abscissa < (*it2)->t2) break;
				v2 = (*it2)->v2;
				last_abscissa = t2; // (*it2)->t2;
				it2++;
			}

			if (it2 == edges_for_this_intersected_ray.end()) {
				// Adds a new edge at the end of the list
				if (v2 != vertex) {
					//int tag = (current_event->t_intersected > rays[current_event->intersected]->t_swap ? 1 : 0);
					double alpha_1 = r_j->opposite()->incidence_angle, alpha_2 = r_j->incidence_angle;

					Inner_Edge* v2_v = new Inner_Edge(id_edges, r_j, 0, v2, vertex, alpha_1, alpha_2, last_abscissa, abscissa);
					edges_for_this_intersected_ray.push_back(v2_v);
				}
			} else {
				// Splits this edge
				Inner_Edge* split_edge = *it2;
				Vertex* v1 = split_edge->v1, *v2 = split_edge->v2;
				double t1, t2;
				split_edge->time_range(r_j, t1, t2);
				//double t1 = split_edge->t1, t2 = split_edge->t2;
				double alpha_1 = v1->incidence_angle(split_edge), alpha_2 = v2->incidence_angle(split_edge);

				if (v1 != vertex && v2 != vertex) {
					Inner_Edge* v1_v = new Inner_Edge(id_edges, r_j, 0, v1, vertex, alpha_1, alpha_2, t1, abscissa);
					Inner_Edge* v_v2 = new Inner_Edge(id_edges, r_j, 0, vertex, v2, alpha_1, alpha_2, abscissa, t2);

					it2 = edges_for_this_intersected_ray.erase(it2);
					edges_for_this_intersected_ray.insert(it2, v1_v);
					edges_for_this_intersected_ray.insert(it2, v_v2);
					delete split_edge;
				}
			}
		}
	}

	// Second, the intersectant ray : we create an edge that comes from the last intersection point
	// in which it was actively involved, to the interesection point we just created
	
	Vertex* v3 = NULL;
	Inner_Edge* v3_v = NULL;
	abscissa = current_event->t_intersectant;
	double alpha_v3 = r_i->opposite()->incidence_angle;
	double alpha_v = r_i->incidence_angle;

	list<Inner_Edge *> & edges_for_this_intersectant_ray = inner_edges[current_event->intersectant];
	list<Inner_Edge *>::iterator it2 = edges_for_this_intersectant_ray.begin();
	if (edges_for_this_intersectant_ray.size() == 0) {
		v3 = initial_vertices_of_rays[current_event->intersectant];
		if (v3 != vertex) {
			//int tag = int(rays[current_event->intersectant]->primary_condition ? 0 : rays[current_event->intersectant]->secondary_condition);
			v3_v = new Inner_Edge(id_edges, r_i, 0, v3, vertex, alpha_v3, alpha_v, -r_i->initial_length, abscissa);
		}
	} else {
		double last_abscissa = -FLT_MAX;
		while (it2 != edges_for_this_intersectant_ray.end()) {
			v3 = (*it2)->v2;
			double t1, t2;
			(*it2)->time_range(r_i, t1, t2);
			//last_abscissa = (*it2)->t2;
			it2++;
		}
		if (v3 != vertex) {
			//int tag = int(rays[current_event->intersectant]->primary_condition ? 0 : rays[current_event->intersectant]->secondary_condition);
			v3_v = new Inner_Edge(id_edges, r_i, 0, v3, vertex, alpha_v3, alpha_v, last_abscissa, abscissa);
		}
	}
	if (v3_v != NULL) edges_for_this_intersectant_ray.push_back(v3_v);
}



void Partition::find_pixels_inside_facets()
{
	for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
		(*it_f)->find_pixels_inside_facet();
	}
}



void Partition::merge_containers(list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices, vector<list<Outer_Edge *> > & outer_edges, vector<list<Inner_Edge *> > & inner_edges)
{
	// We are going to insert elements contained in outer_vertices, inner_vertices in a double-linked list
	for (list<Vertex *>::iterator it_v = outer_vertices.begin(); it_v != outer_vertices.end(); it_v++) push_back(*it_v);
	outer_vertices.clear();
	for (list<Vertex *>::iterator it_v = inner_vertices.begin(); it_v != inner_vertices.end(); it_v++) push_back(*it_v);
	inner_vertices.clear();

	// We also merge the different lists of edges
	// Step 1 : inner_edges
	for (int i = 0 ; i < inner_edges.size() ; i++) {
		for (list<Inner_Edge *>::iterator it_e = inner_edges[i].begin(); it_e != inner_edges[i].end(); it_e++) push_back(*it_e);
		inner_edges[i].clear();
	}
	inner_edges.clear();
	// Step 2 : outer_edges
	for (int i = 0; i < outer_edges.size(); i++) {
		for (list<Outer_Edge *>::iterator it_e = outer_edges[i].begin(); it_e != outer_edges[i].end(); it_e++) push_back(*it_e);
		outer_edges[i].clear();
	}
	outer_edges.clear();
}



void Partition::build_faces(Size2i & size)
{
	id_edges = 0;
	int nb_edges = int(e_size);
	bool** queued = new bool*[nb_edges];

	Edge* e = edges_head;
	while (e != NULL) {
		queued[id_edges] = new bool[2];
		// The algorithm of face retrieval consists in crossing over the edges of the
		// graph in clockwise order to find cycles. This leads us to consider virtual
		// half-edges : for every e = (v1 v2), there "exists" the half-edge (v1 v2)
		// and the half-edge (v2 v1), except for the outer edges for which only one
		// half-edge exists. We eliminate such half-edges from the algorithm by
		// considering that they have been already processed.

		bool top_image = false, bottom_image = false, right_image = false, left_image = false;
		if (e->type == OUTER_EDGE) {
			Outer_Edge* oe = static_cast<Outer_Edge *>(e);
			top_image = (oe->boundary == TOP_IMAGE);
			bottom_image = (oe->boundary == BOTTOM_IMAGE);
			right_image = (oe->boundary == RIGHT_IMAGE);
			left_image = (oe->boundary == LEFT_IMAGE);
		}

		if (top_image || left_image) {
			// Entry 0 correspond to the direction (v1 v2), entry 1 to (v2 v1)
			// Values are determined according to the way outer edges are created
			queued[id_edges][0] = true;
			queued[id_edges][1] = false;
		} else if (bottom_image || right_image) {
			queued[id_edges][0] = false;
			queued[id_edges][1] = true;
		} else {
			// General case : inner edge
			queued[id_edges][0] = false;
			queued[id_edges][1] = false;
		}
		e->id_edge = id_edges++;
		e = e->e_next;
	}

	// Defines a queue and initializes with the top left half-edge of the figure
	std::queue<HalfEdge *> queue;
	Vertex* top_left = vertices_head;
	HalfEdge* h = top_left->directions[top_left->directions.size() - 1].second;
	queue.push(h);

	// While there are elements left in the queue
	while (queue.size() != 0) {

		// We define a reference to the first element of the queue,
		// which is the first element of a new face
		HalfEdge* h_e = queue.front();
		HalfEdge* h_f = h_e;
		list<HalfEdge *> current_face;

		// If the element hasn't been already processed in the meanwhile
		if (!queued[h_f->e->id_edge][h_f->v1_v2]) {
			queued[h_f->e->id_edge][h_f->v1_v2] = true;

			do {
				// Finds the next half-edge in clockwise order
				// Adds it to the face we are currently defining
				Vertex* v = (h_f->v1_v2 ? h_f->e->v2 : h_f->e->v1);
				int path = 0;
				for (int i = 0; i < v->directions.size(); i++) {
					if (v->directions[i].second->e == h_f->e) {
						path = (i + 1) % v->directions.size();
						break;
					}
				}
				h_f = v->directions[path].second;
				assert((h_f->v1_v2 && h_f->e->v1 == v) || (!h_f->v1_v2 && h_f->e->v2 == v));
				current_face.push_back(h_f);

				// The current half-edge shouldn't be queued
				queued[h_f->e->id_edge][h_f->v1_v2] = true;

				// Inserts the opposite half-edge into the queue
				if (!queued[h_f->e->id_edge][!h_f->v1_v2]) {
					if (h_f->e->v1_v2 == h_f) {
						queue.push(h_f->e->v2_v1);
					} else {
						queue.push(h_f->e->v1_v2);
					}
				}

			} while (h_f != h_e);
			faces.push_back(new Face(id_faces, current_face));
		}

		// Removes the first element of the queue, since we're done with it
		queue.pop();
	}

	// Deletes the structure that has previously been allowed
	for (int i = 0; i < int(e_size); i++) {
		delete[] queued[i];
	}
	delete[] queued;

	// The graph can be used to compute statistics
    is_valid = (count_invalid_faces() == 0);
}



void Partition::seal_and_remove_bivalent_vertices()
{
	for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
		Face* f = (*it_f);
		for (list<HalfEdge *>::iterator it_h = f->edges.begin() ; it_h != f->edges.end() ; it_h++) {
			HalfEdge* h = (*it_h);
			int i = h->e->id_edge;
		}
	}

	Vertex* v = vertices_head;
	while (v != NULL) {
		assert(v->connectivity() > 1);

		// The 4 first vertices are image corners
		if (v->connectivity() == 2 && v->id_vertex > 3) {

			list<Edge *> edges_to_delete;
			list<Vertex *> vertices_to_delete;
			set<Ray *> rays;
			Image_Boundary boundary = INVALID_BORDER;

			// The vertex v is bivalent : that's why we find the vertices to which v is connected
			// However, maybe the connected vertices are themselves bivalent : we need to iterate on v1 and v2 until
			// finding trivalent vertices. Until then, we define vertices and edges to delete
			Edge* e1 = v->directions[0].second->e;
			Edge* e2 = v->directions[1].second->e;

			Edge_Type e_type = e1->type;
			Inner_Edge *ie = nullptr;
			Outer_Edge* oe = nullptr;

			Vertex* v1 = (e1->v1 == v ? e1->v2 : e1->v1);
			Vertex* v2 = (e2->v1 == v ? e2->v2 : e2->v1);
			double alpha_v1 = v1->incidence_angle(e1);
			double alpha_v2 = v2->incidence_angle(e2);
			
			// We define a sequence of vertices (v1 .. v .. v2)
			// Obtains the two facets adjacent to the big edge (v1 v2) : to this end, we can use e2
			Face *f_v1_v2 = nullptr, *f_v2_v1 = nullptr;
			HalfEdge* h = e2->v1_v2;
			if (h->e->v1 == v && h->e->v2 == v2) {
				f_v1_v2 = h->f;
				f_v2_v1 = h->opposite()->f;
			} else {
				f_v1_v2 = h->opposite()->f;
				f_v2_v1 = h->f;
			}

			if (e_type == INNER_EDGE) {
				ie = static_cast<Inner_Edge*>(e1);
				ie->get_supporting_rays(rays);
				ie = static_cast<Inner_Edge*>(e2);
				ie->get_supporting_rays(rays);
			} else {
				oe = static_cast<Outer_Edge*>(e1);
				boundary = oe->boundary;
			}

			edges_to_delete.push_back(e1);
			edges_to_delete.push_back(e2);

			// Finds the edge we have not explored so far on v1's side
			while (v1->connectivity() == 2 && v1->id_vertex > 3) {
				vertices_to_delete.push_back(v1);
				if (v1->directions[0].second->e == e1) {
					e1 = v1->directions[1].second->e;
				} else {
					e1 = v1->directions[0].second->e;
				}
				edges_to_delete.push_back(e1);
				if (e_type == INNER_EDGE) {
					ie = static_cast<Inner_Edge*>(e1);
					ie->get_supporting_rays(rays);
				}
				v1 = (e1->v1 == v1 ? e1->v2 : e1->v1);
			}

			// Same operation for v2
			while (v2->connectivity() == 2 && v2->id_vertex > 3) {
				vertices_to_delete.push_back(v2);
				if (v2->directions[0].second->e == e2) {
					e2 = v2->directions[1].second->e;
				} else {
					e2 = v2->directions[0].second->e;
				}
				edges_to_delete.push_back(e2);

				if (e_type == INNER_EDGE) {
					ie = static_cast<Inner_Edge*>(e2);
					ie->get_supporting_rays(rays);
				}
				v2 = (e2->v1 == v2 ? e2->v2 : e2->v1);
			}

			// Creates the merged edge
			Edge* e_12 = nullptr;
			if (e_type == INNER_EDGE) {
				e_12 = new Inner_Edge(id_edges, rays, 0, v1, v2, alpha_v1, alpha_v2);
			} else {
				e_12 = new Outer_Edge(id_edges, boundary, v1, v2, alpha_v1, alpha_v2);
			}

			if (f_v1_v2 != nullptr) {
				f_v1_v2->remove_non_corners(v1, v2, e_12->v1_v2);

				for (list<HalfEdge *>::iterator it_h = f_v1_v2->edges.begin() ; it_h != f_v1_v2->edges.end() ; it_h++) {
					int i = (*it_h)->e->id_edge;
				}
			}

			if (f_v2_v1 != nullptr) {
				f_v2_v1->remove_non_corners(v2, v1, e_12->v2_v1);

				for (list<HalfEdge *>::iterator it_h = f_v2_v1->edges.begin() ; it_h != f_v2_v1->edges.end() ; it_h++) {
					int i = (*it_h)->e->id_edge;
				}
			}

			// Deletes the references to the splitted edges and the bivalent edges
			erase(edges_to_delete, true);
			for (list<Vertex *>::iterator it_v = vertices_to_delete.begin(); it_v != vertices_to_delete.end(); it_v++) {
				erase(*it_v);
				quadtree->remove((*it_v), true);
			}

			// Deletes the last bivalent edge
			Vertex* v_next = erase(v);
			quadtree->remove(v, true);
			v = v_next;

			// Adds the merged edge to the list of vertices
			push_back(e_12);

		} else {
			v = v->v_next;
		}
	}

	for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
		Face* f = (*it_f);
		for (list<HalfEdge *>::iterator it_h = f->edges.begin() ; it_h != f->edges.end() ; it_h++) {
			HalfEdge* h = (*it_h);
			int i = h->e->id_edge;
		}
	}
}



void Partition::merge_thin_facets(vector<Ray *> & rays, int verbose_level)
{
	//clock_t t_begin = clock();

	// Computes pre-statistics
#if NOT_MEASURING_PERFORMANCES
	uint convex_union_exists = 0;
	uint convex_cut_exists = 0;
	uint convex_union_and_cut_exist = 0;
	uint other_cells = 0;
	uint thin_cells = 0;
	for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; ++it_f) {

		Face* f_i = (*it_f);
		if (f_i->thinness <= params->merge_min_thinness) {
			list<pair<uint, uint> > edges_to_remove;
			list<pair<uint, uint> > edges_to_slice;
			f_i->convex_union(edges_to_remove);
			f_i->convex_cut(edges_to_slice);
			if (!edges_to_remove.empty() && !edges_to_slice.empty()) {
				convex_union_and_cut_exist++;
			} else if (!edges_to_remove.empty() && edges_to_slice.empty()) {
				convex_union_exists++;
			} else if (edges_to_remove.empty() && !edges_to_slice.empty()) {
				convex_cut_exists++;
			} else {
				other_cells++;
			}
			thin_cells++;
		}
	}

	trace(verbose_level, 15, "** Pre-statistics for the merge of thin cells :");
	trace(verbose_level, 15, "** Convex union    : " + std::to_string(convex_union_exists) + "/" + std::to_string(thin_cells) + " (" + std::to_string(double(100 * convex_union_exists) / thin_cells) + " %)");
	trace(verbose_level, 15, "** Convex cut      : " + std::to_string(convex_cut_exists) + "/" + std::to_string(thin_cells) + " (" + std::to_string(double(100 * convex_cut_exists) / thin_cells) + " %)");
	trace(verbose_level, 15, "** Both applicable : " + std::to_string(convex_union_and_cut_exist) + "/" + std::to_string(thin_cells) + " (" + std::to_string(double(100 * convex_union_and_cut_exist) / thin_cells) + " %)");
	trace(verbose_level, 15, "** Other cells     : " + std::to_string(other_cells) + "/" + std::to_string(thin_cells) + " (" + std::to_string(double(100 * other_cells) / thin_cells) + " %)");
#endif
	map<int, list<Face *>::iterator> map_facets;
	for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
		map_facets[(*it_f)->id_face] = it_f;
	}

	list<Face *>::iterator it_f = faces.begin();
	while (it_f != faces.end()) {

		Face* f_i = (*it_f);
		if (f_i->thinness > params->merge_min_thinness) {
			++it_f;
			continue;
		}

		// We try to merge tiny facets while keeping the convexity of the cells. Two operations can be performed :
		// - merging a facet with adjacent facets, whose edges are prolongating the current facet's edges.
		// - diving a facet into several subfacets along the adjacent facet's edges and merging them with the adjacent cells.

		list<pair<uint, uint> > edges_to_remove;
		list<pair<uint, uint> > edges_to_slice;

		if (f_i->convex_union(edges_to_remove)) {
			join_thin_facet(map_facets, f_i, edges_to_remove);
			it_f = faces.erase(it_f);
		} else if (f_i->convex_cut(edges_to_slice)) {
			slice_thin_facet(rays, map_facets, f_i, edges_to_slice);
			it_f = faces.erase(it_f);
		} else {
			++it_f;
		}
	}

	clock_t t_end = clock();

	// Resets indices
	id_faces = -1;
	for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
		(*it_f)->id_face = ++id_faces;
	}

	//double elapsed_time = double(t_end - t_begin) / CLOCKS_PER_SEC;
	//trace(verbose_level, 5, "** Merged thin facets (remains : " + std::to_string(faces.size()) + ") in " + std::to_string(elapsed_time) + " s.");
}



void Partition::reset_indices()
{
	id_vertices = id_edges = id_faces = 0;

	Vertex* v = vertices_head;
	while (v != nullptr) {
		v->id_vertex = id_vertices++;
		v = v->v_next;
	}

	Edge* e = edges_head;
	while (e != nullptr) {
		e->id_edge = id_edges++;
		e = e->e_next;
	}

	for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
		(*it_f)->id_face = id_faces++;
	}
}



Face* Partition::third_facet(Vertex* v_i, HalfEdge* h_ij, Face* f_i, Face* f_j)
{
	for (uint i = 0; i < v_i->directions.size(); i++) {
		HalfEdge* h_i = v_i->directions[i].second;
		if (h_i == h_ij) {
			continue;
		} else {
			if (h_i->f == f_i || h_i->f == f_j) {
				return h_i->opposite()->f;
			} else {
				return h_i->f;
			}
		}
	}
	return nullptr;
}



void Partition::determine_best_edge_for_joining(Face* f_i, list<pair<uint, uint> > & possibilities, pair<uint, uint> & selected)
{
	vector<HalfEdge *> f_i_edges(f_i->edges.begin(), f_i->edges.end());

	// If there are several facets adjacent to f_i with which this f_i can be merged, we are going to select the thinnest of them
	double th_min = FLT_MAX;
	for (list<pair<uint, uint> >::iterator it_u = possibilities.begin() ; it_u != possibilities.end() ; it_u++) {
		int seq_start = it_u->first, seq_length = it_u->second;
		Edge* e = f_i_edges[seq_start]->e;
		Face* f_adj = (e->v1_v2->f == f_i ? e->v2_v1->f : e->v1_v2->f);
		if (f_adj->thinness < th_min) {
			th_min = f_adj->thinness;
			selected = *it_u;
		}
	}
}



void Partition::join_thin_facet(map<int, list<Face *>::iterator> & map_facets, Face* & f_i, list<pair<uint, uint> > & edges_joinable_facets)
{
	uint fi_id = f_i->id_face;

	pair<uint, uint> selected;
	if (edges_joinable_facets.size() == 1) {
		selected = edges_joinable_facets.front();
	} else {
		determine_best_edge_for_joining(f_i, edges_joinable_facets, selected);
	}

	list<HalfEdge *>::iterator it_h = f_i->edges.begin();
	for (uint i = 0; i < selected.first; ++i) ++it_h;
	Face* f_j = (*it_h)->opposite()->f;
	uint fj_id = f_j->id_face;

	merge_two_facets(map_facets, f_i, f_j, selected, false, true);
}



void Partition::merge_two_facets(map<int, list<Face *>::iterator> & map_facets, Face* & f_i, Face* & f_j, pair<uint, uint> & h_ij, 
	bool exclude_f_i_from_list, bool exclude_f_j_from_list)
{
	typedef pair<Vertex*, Vertex_Type> Vertex_Element;

	uint id_i = f_i->id_face;
	uint id_j = f_j->id_face;

	uint n_i = f_i->edges.size();
	vector<HalfEdge *> edges_i(f_i->edges.begin(), f_i->edges.end());
	vector<Vertex_Element> vertices_i(f_i->vertices.begin(), f_i->vertices.end());

	// Anticipates the merge of f_i and f_j, by listing the edges of the new facet (f_i U f_j)
	uint seq_start = h_ij.first;
	uint seq_length = h_ij.second;
	Vertex* v_1 = vertices_i[seq_start].first;
	Vertex* v_2 = vertices_i[(seq_start + seq_length) % n_i].first;

	list<HalfEdge *> halfedges_ij;
	if (seq_length == 1) {
		HalfEdge* h = edges_i[seq_start];
		Face::merged_halfedges(f_i, f_j, h->e, halfedges_ij);
		
	} else {
		HalfEdge* h_1 = edges_i[seq_start];
		HalfEdge* h_2 = edges_i[(seq_start + seq_length - 1) % n_i];
		Face::merged_halfedges(f_i, f_j, h_1->e, h_2->e, halfedges_ij);
	}

	if (v_1->connectivity() == 3) {
		Face* f_1 = third_facet(v_1, edges_i[seq_start], f_i, f_j);
		if (f_1 != nullptr) f_1->set_as_bivalent(v_1);
	}
	if (v_2->connectivity() == 3) {
		Face* f_2 = third_facet(v_2, edges_i[(seq_start + seq_length - 1) % n_i], f_i, f_j);
		if (f_2 != nullptr) f_2->set_as_bivalent(v_2);
	}

	// Erases the facets f_i and f_j
	if (exclude_f_i_from_list) faces.erase(map_facets[id_i]);	
	if (exclude_f_j_from_list) faces.erase(map_facets[id_j]);
	map_facets.erase(id_i);
	map_facets.erase(id_j);
	delete f_j;
	delete f_i;

	// Erases the edges between f_i and f_j
	for (uint i = 0 ; i < seq_length ; i++) {
		HalfEdge* h_i = edges_i[(seq_start + i) % n_i];
		erase(h_i->e, true);
	}

	// Erases the bivalent vertices found along the edge (h_1, h_2)
	if (seq_length > 1) {
		for (uint i = 0 ; i < seq_length - 1 ; i++) {
			Vertex* v_i = vertices_i[(seq_start + i + 1) % n_i].first;
			erase(v_i);
			quadtree->remove(v_i, true);
		}
	}

	// Builds the new facet
	Face* f_ij = new Face(id_faces, halfedges_ij);
	faces.push_back(f_ij);
	map_facets[f_ij->id_face] = --faces.end();
}



void Partition::determine_best_edge_for_slicing(vector<HalfEdge *> & f_i_edges, vector<pair<Vertex *, Vertex_Type> > & f_i_vertices, 
	list<pair<uint, uint> > & possibilites, pair<uint, uint> & selected)
{
	// In order to determine the best line of cut for slicing f_i, we use two criteria
	// The first one is the length of the edge : it is useless to merge f_i with other thin faces
	uint n = uint(f_i_edges.size());

	double e_longest = -FLT_MAX;
	list<pair<uint, uint> > e_longest_argmax;
	for (list<pair<uint, uint> >::iterator it_p = possibilites.begin() ; it_p != possibilites.end() ; it_p++) {

		pair<uint, uint> sequence = (*it_p);
		uint seq_start = sequence.first, seq_length = sequence.second;

		HalfEdge* h_first = f_i_edges[seq_start];
		Vertex* v1 = (h_first->v1_v2 ? h_first->e->v1 : h_first->e->v2);

		HalfEdge* h_last = f_i_edges[(seq_start + seq_length) % n];
		Vertex* v2 = (h_last->v1_v2 ? h_last->e->v2 : h_last->e->v1);
		double e_length = (v2->pt.x - v1->pt.x) * (v2->pt.x - v1->pt.x) + (v2->pt.y - v1->pt.y) * (v2->pt.y - v1->pt.y);

		if (fabs(e_longest - e_length) < 1e-3) {
			e_longest_argmax.push_back(sequence);
		} else if (e_length > e_longest) {
			e_longest_argmax.clear();
			e_longest_argmax.push_back(sequence);
		}
	}

	if (e_longest_argmax.size() == 1) {
		selected = e_longest_argmax.front();
		return;
	}

	// If we have to choose between several lines of cut of equal length
	// then we select the line of cut which offers, in average, the most orthogonal direction of cut
	double alpha_mean_best = -FLT_MAX;
	for (list<pair<uint, uint> >::iterator it_p = e_longest_argmax.begin() ; it_p != e_longest_argmax.end() ; it_p++) {
	
		double alpha_mean = 0;
		pair<uint, uint> sequence = (*it_p);
		uint seq_start = sequence.first, seq_length = sequence.second;

		// Reference angle, determines the orientation of the sequence of vertices
		Vertex* corner_1 = f_i_vertices[seq_start].first;
		Vertex* corner_2 = f_i_vertices[(seq_start + seq_length) % n].first;
		Vec2d corner_12 = Vec2d(corner_2->pt.x - corner_1->pt.x, corner_2->pt.y - corner_1->pt.y);
		corner_12 = cv::normalize(corner_12);

		// For each vertex which is not a corner, we evaluate the angle made by its most orthogonal halfedge (to corner_12)
		int trivalent_vertices = 0;
		for (uint k = 0 ; k < seq_length - 1 ; k++) {

			if (f_i_vertices[(seq_start + k + 1) % n].second == NON_CORNER_TRIVALENT) {
				Vertex* v = f_i_vertices[(seq_start + k + 1) % n].first;
				double min_abs_dot_product = FLT_MAX;

				// Gets vertices surrouding v on line (corner_1, corner_2)
				Vertex* v_prev = f_i_vertices[(seq_start + k) % n].first;
				Vertex* v_next = f_i_vertices[(seq_start + k + 2) % n].first;

				// For each other direction that (corner_1, corner_2) evaluates the minimal dot product
				// The dot product corresponds to a angle that we want as close to PI/2 as possible
				double alpha = 0;
				HalfEdge* h;
				v->most_orthogonal_halfedge(v_prev, v_next, corner_12, alpha, h);

				alpha_mean += alpha;
				++trivalent_vertices;
			}
		}

		// If in average, lines of cut are more orthogonal for this sequence, then we choose it
		alpha_mean /= trivalent_vertices;
		if (alpha_mean < alpha_mean_best) {
			alpha_mean_best = alpha_mean;
			selected = (*it_p);
		}
	}
}


void Partition::determine_intersection_point_when_slicing(Vertex* v_curr, Vec2d & sliced_edge_direction, list<HalfEdge *> & intersectable_halfedges, 
	IndexedEvent* & intersection, HalfEdge* & extended_halfedge, HalfEdge* & intersected_halfedge)
{
	typedef pair<double, double> Time_Range;

	// Gets the extended halfedge and ray
	HalfEdge* h_cut = v_curr->most_orthogonal_halfedge(sliced_edge_direction);
	extended_halfedge = h_cut->opposite();

	Inner_Edge* extended_edge = static_cast<Inner_Edge *>(extended_halfedge->e);
	Ray* r_extended = extended_edge->get_front_supporting_ray();

	intersected_halfedge = nullptr;
	list<Ray *> intersectable_rays;
	list<Image_Boundary> intersectable_boundaries;

	Ray* r_intersected = nullptr;
	Image_Boundary b_intersected = INVALID_BORDER;

	// Determines valid times of execution for each of the intersectable rays
	map<int, Time_Range> intersection_times;
	map<HalfEdge *, Time_Range> intersection_times_per_halfedge;
	for (list<HalfEdge *>::iterator it_h = intersectable_halfedges.begin() ; it_h != intersectable_halfedges.end() ; ++it_h) {
		HalfEdge* h = (*it_h);
		int ray_index;

		double t_1, t_2;
		if (h->e->type == INNER_EDGE) {
			// If the intersected halfedge is inside the image, it has a support ray
			Inner_Edge* e = static_cast<Inner_Edge *>(h->e);
			r_intersected = e->get_front_supporting_ray();
			e->time_range(r_intersected, t_1, t_2);
			ray_index = int(r_intersected->index);
		} else {
			// If the intersected halfedge is a border, gets the index of the boundary
			Outer_Edge* e = static_cast<Outer_Edge *>(h->e);
			b_intersected = e->boundary;
			e->time_range(t_1, t_2);
			ray_index = int(b_intersected);
		}

		double t_inf, t_sup;
		if (t_1 < t_2) { 
			t_inf = t_1; t_sup = t_2;
		} else {
			t_inf = t_2; t_sup = t_1;
		}
		intersection_times_per_halfedge[h] = std::make_pair(t_inf, t_sup);

		// Updates the map of valid intersection times
		if (intersection_times.find(ray_index) == intersection_times.end()) {
			intersection_times[ray_index] = std::make_pair(t_inf, t_sup);
			if (ray_index >= 0) {
				intersectable_rays.push_back(r_intersected);
			} else {
				intersectable_boundaries.push_back(b_intersected);
			}
		} else {
			Time_Range & range = intersection_times[ray_index];
			if (t_inf < range.first) range.first = t_inf;
			if (t_sup > range.second) range.second = t_sup;
		}
	}

	int index_intersected_ray = -INT_MAX;
	double t_i = FLT_MAX, t_j = FLT_MAX;
	
	// See notes below
	double dist_required_interval = FLT_MAX;
	double hidden_t_i = FLT_MAX, hidden_t_j = FLT_MAX;
	int ray_argmin_dist_interval = -INT_MAX;

	// Now, for each of the rays that delimit the intersectable halfedges
	// We compute the times of intersection with r_extended, and check if the result is valid
	for (list<Ray *>::iterator it_r = intersectable_rays.begin() ; it_r != intersectable_rays.end() ; it_r++) {
		r_intersected = (*it_r);
		Geometry::direct_intersection(r_extended, r_intersected, t_i, t_j);
		
		// There are two ways to determine when we hit the correct ray
		// The first method is when t_j is inside the segment [range.first, range.second]

		Time_Range & range = intersection_times[r_intersected->index];
		if (range.first <= t_j && t_j <= range.second) {
			index_intersected_ray = r_intersected->index;
			break;
		} else {

			// However, it happens that the algorithm misses the interval by an epsilon
			// But on the other hand, we don't want to return the wrong segment
			// So we compute a distance that must be as small as possible
			// and that we use if we are enable to find a value that clearly falls within a interval
			double d = (t_j < range.first ? range.first - t_j : t_j - range.second);
			if (d < dist_required_interval) {
				dist_required_interval = d;
				hidden_t_i = t_i;
				hidden_t_j = t_j;
				ray_argmin_dist_interval = r_intersected->index;
			}
		}
	}

	if (index_intersected_ray == -INT_MAX) {
		for (list<Image_Boundary>::iterator it_b = intersectable_boundaries.begin() ; it_b != intersectable_boundaries.end() ; it_b++) {
			b_intersected = (*it_b);
			Geometry::direct_intersection(r_extended, b_intersected, rows, cols, t_i, t_j);

			Time_Range & range = intersection_times[int(b_intersected)];
			if (range.first <= t_j && t_j <= range.second) {
				index_intersected_ray = int(b_intersected);
				break;
			} else {
				double d = (t_j < range.first ? range.first - t_j : t_j - range.second);
				if (d < dist_required_interval) {
					dist_required_interval = d;
					hidden_t_i = t_i;
					hidden_t_j = t_j;
					ray_argmin_dist_interval = int(b_intersected);
				}
			}
		}
	}

	// Case when we have to handle the round precision problem,
	// and use the distances to the measured intervals
	if (index_intersected_ray == -INT_MAX) {
		if (dist_required_interval >= 1e-6) {
			std::cout << intersectable_boundaries.size() << " " << intersectable_halfedges.size() << std::endl;
			std::cout << dist_required_interval << std::endl;
		}
		assert(dist_required_interval < 1e-6);
		t_i = hidden_t_i;
		t_j = hidden_t_j;
		index_intersected_ray = ray_argmin_dist_interval;
		if (index_intersected_ray >= 0) {
			for (list<Ray *>::iterator it_r = intersectable_rays.begin() ; it_r != intersectable_rays.end() ; it_r++) {
				if (index_intersected_ray == (*it_r)->index) {
					r_intersected = (*it_r);
					break;
				}
			}
		} else {
			for (list<Image_Boundary>::iterator it_b = intersectable_boundaries.begin() ; it_b != intersectable_boundaries.end() ; it_b++) {
				if (index_intersected_ray == int(*it_b)) {
					b_intersected = (*it_b);
					break;
				}
			}
		}
	}

	// Returns an event and the halfedge intersected by the ray
	dist_required_interval = FLT_MAX;
	
	HalfEdge* halfedge_argmin_dist_interval = nullptr;

	intersected_halfedge = nullptr;
	for (map<HalfEdge*, Time_Range>::iterator it_m = intersection_times_per_halfedge.begin() ; it_m != intersection_times_per_halfedge.end() ; it_m++) {
		HalfEdge* h = it_m->first;
		Time_Range & range = it_m->second;
		if ((index_intersected_ray >= 0) && (h->e->type == INNER_EDGE)) {
			Inner_Edge* e = static_cast<Inner_Edge *>(h->e);
			if (e->get_front_supporting_ray() == r_intersected) { 
				if (range.first <= t_j && t_j <= range.second) {
					intersected_halfedge = h;
					intersection = new IndexedEvent(r_extended->index, r_intersected->index, t_i, t_j, false);
					return;

				} else {
					double d = (t_j < range.first ? range.first - t_j : t_j - range.second);
					if (d < dist_required_interval) {
						dist_required_interval = d;
						halfedge_argmin_dist_interval = h;
					}
				}
			}
		} else if ((index_intersected_ray < 0) && (h->e->type == OUTER_EDGE)) {
			Outer_Edge* e = static_cast<Outer_Edge *>(h->e);
			if (e->boundary == b_intersected) {
				if (range.first <= t_j && t_j <= range.second) {
					intersected_halfedge = h;
					intersection = new IndexedEvent(r_extended->index, b_intersected, t_i, t_j);
					return;

				} else {
					double d = (t_j < range.first ? range.first - t_j : t_j - range.second);
					if (d < dist_required_interval) {
						dist_required_interval = d;
						halfedge_argmin_dist_interval = h;
					}
				}
			}
		}
	}

	// If we reach this point we are facing the problem of round precision
	intersected_halfedge = halfedge_argmin_dist_interval;
	if (index_intersected_ray >= 0) {
		intersection = new IndexedEvent(r_extended->index, r_intersected->index, t_i, t_j, false);
	} else {
		intersection = new IndexedEvent(r_extended->index, b_intersected, t_i, t_j);
	}
}


void Partition::classify_along_sliced_edge(vector<HalfEdge *> & big_sliced_edge, vector<list<HalfEdge *> > & sequences_halfedges_to_trivalent, vector<list<Vertex *> > & sequences_bivalent_vertices)
{
	list<HalfEdge *> curr_sequence_halfedges_to_trivalent;
	list<Vertex *> curr_sequence_bivalent_vertices;

	for (uint i = 0 ; i < big_sliced_edge.size() - 1 ; i++) {

		HalfEdge* h_i = big_sliced_edge[i];
		Vertex* v_i = (h_i->v1_v2 ? h_i->e->v2 : h_i->e->v1);

		curr_sequence_halfedges_to_trivalent.push_back(h_i);
		if (v_i->connectivity() == 2) {
			// If we reach a bivalent vertex, the current sequence is not over
			// and we can add another vertex we should later delete
			curr_sequence_bivalent_vertices.push_back(v_i);

		} else {
			// If we reach a non-bivalent vertex the sequence ends
			sequences_halfedges_to_trivalent.push_back(curr_sequence_halfedges_to_trivalent);
			sequences_bivalent_vertices.push_back(curr_sequence_bivalent_vertices);

			// We clear the containers to start we can start a new sequence
			curr_sequence_halfedges_to_trivalent.clear();
			curr_sequence_bivalent_vertices.clear();
		}
	}

	// The last halfedge found concludes the current sequence, and we know it is a corner
	curr_sequence_halfedges_to_trivalent.push_back(big_sliced_edge.back());
	
	sequences_halfedges_to_trivalent.push_back(curr_sequence_halfedges_to_trivalent);
	sequences_bivalent_vertices.push_back(curr_sequence_bivalent_vertices);
}



void Partition::slice_thin_facet(vector<Ray *> & rays, map<int, list<Face *>::iterator> & map_facets, Face* & f_i, list<pair<uint, uint> > & edges_to_slice)
{
	uint f_id = f_i->id_face;

	vector<HalfEdge *> f_i_edges (f_i->edges.begin(), f_i->edges.end());
	vector<pair<Vertex *, Vertex_Type> > f_i_vertices (f_i->vertices.begin(), f_i->vertices.end());
	uint n = uint(f_i->vertices.size());

	// There may exist several possibilities for slicing the image (several pairs of corners separated by at least one non-corner and trivalent vertex)
	// Therefore and as a preliminar operation, we must select the best edge to use for slicing the image
	pair<uint, uint> big_sliced_edge;
	if (edges_to_slice.size() == 1) {
		big_sliced_edge = edges_to_slice.front();
	} else {
		determine_best_edge_for_slicing(f_i_edges, f_i_vertices, edges_to_slice, big_sliced_edge);
	}

	// We have obtained a edge which is a sequence of halfedges. 
	// The idea is to destroy all them, and associate the right parts of f_i to the adjacent facets by extending the rays.
	uint seq_start = big_sliced_edge.first, seq_length = big_sliced_edge.second;

	// Keeps in memory a list of edges and vertices along the sliced edge
	// Gets its direction
	vector<HalfEdge *> sliced_edges;
	vector<Vertex *> sliced_vertices;
	sliced_edges.reserve(seq_length);
	sliced_vertices.reserve(seq_length + 1);

	for (uint i = 0 ; i < seq_length ; i++) {
		sliced_edges.push_back(f_i_edges[(seq_start + i) % n]);
		sliced_vertices.push_back(f_i_vertices[(seq_start + i) % n].first);
	}
	sliced_vertices.push_back(f_i_vertices[(seq_start + seq_length) % n].first);

	Inner_Edge* sliced_edge_as_inner_edge = static_cast<Inner_Edge*>(sliced_edges.front()->e);
	Segment* sliced_edge_seed = sliced_edge_as_inner_edge->get_front_supporting_segment();
	Vec2d sliced_edge_direction = Vec2d(-sliced_edge_seed->b, sliced_edge_seed->a);

	// For each non-corner vertex, from corner_1 to corner_2, we get the most orthogonal direction to corner_12 and cut the cell into pieces.
	// To this end we need to maintain a list of edges that the extended rays can intersect :
	// At first it consists in edges of the cell that are not included in the big, sliced edge.
	
	list<HalfEdge *> remaining_halfedges_of_f_i;
	map<HalfEdge*, list<HalfEdge *>::iterator> access_remaining_halfedges_of_f_i;

	for (uint i = 0 ; i < n - seq_length ; i++) {
		HalfEdge* h_i = f_i_edges[(seq_start + seq_length + i) % n];
		remaining_halfedges_of_f_i.push_back(h_i);
		access_remaining_halfedges_of_f_i[h_i] = --remaining_halfedges_of_f_i.end();
	}
	
	map_facets.erase(f_i->id_face);
	delete f_i;

	// Classifies the halfedges between the two corners
	// We want sequences of halfedges between non corners, trivalent vertices
	vector<list<HalfEdge *> > sequences_halfedges_to_trivalent;
	vector<list<Vertex *> > sequences_bivalent_vertices;
	classify_along_sliced_edge(sliced_edges, sequences_halfedges_to_trivalent, sequences_bivalent_vertices);
	uint subfacets = uint(sequences_halfedges_to_trivalent.size());
	
	// For each sequence of halfedges between corner_1 and corner_2, that ends with a trivalent vertex
	for (uint i = 0; i < subfacets - 1; i++) {

		// Gets the current trivalent vertex from which a new ray is coming
		list<HalfEdge *> current_sequence = sequences_halfedges_to_trivalent[i];
		HalfEdge* h = current_sequence.back();
		//int h_support = (*(static_cast<Inner_Edge *>(h->e)->rays.begin()))->index;
		Vertex* v = (h->v1_v2 ? h->e->v2 : h->e->v1);
		assert(v->directions.size() > 2);

		// The first part of the process consists in determining the path which is our line of cut
		IndexedEvent* event = nullptr;
		HalfEdge *extended_edge = nullptr, *intersected_edge = nullptr;
		determine_intersection_point_when_slicing(v, sliced_edge_direction, remaining_halfedges_of_f_i, event, extended_edge, intersected_edge);
		//int intersected_edge_support = (*(static_cast<Inner_Edge *>(intersected_edge->e)->rays.begin()))->index;
		//assert(intersected_edge_support != h_support);

		// The second part of the process consists in modifying the graph.
		// A vertex may be added : the intersection point
		// At least one edge is created, the path of the extended ray.
		list<HalfEdge *> sub_facet_f_i;
		update_graph_when_slicing(rays, current_sequence, v, event, extended_edge, intersected_edge, remaining_halfedges_of_f_i, access_remaining_halfedges_of_f_i, sub_facet_f_i);

		// Creates a subfacet
		Face* sub_f_i = new Face(id_faces, sub_facet_f_i);
		faces.push_back(sub_f_i);
		map_facets[sub_f_i->id_face] = --faces.end();
	}

	// Last sequence of halfedges : no need to extend a ray
	HalfEdge* last_sliced_edge = sliced_edges[seq_length - 1];
	list<HalfEdge *> sub_f_i_last_halfedges;
	loop_inside_sub_facet(last_sliced_edge, sub_f_i_last_halfedges);

	Face* sub_f_i_last = new Face(id_faces, sub_f_i_last_halfedges);
	faces.push_back(sub_f_i_last);
	map_facets[sub_f_i_last->id_face] = --faces.end();

	// Last step of the algorithm : we merge all the subfacets created before, by associating each of them
	// to the corresponding facet on the other side of the big sliced edge
	for (uint i = 0 ; i < subfacets ; i++) {

		list<HalfEdge *> halfedges = sequences_halfedges_to_trivalent[i];
		HalfEdge* h_i = halfedges.front();

		uint ind_i = 0;
		list<HalfEdge *>::iterator it_h = h_i->f->edges.begin();
		while ((*it_h) != h_i) {
			++it_h;
			++ind_i;
		}

		pair<uint, uint> ref_halfedges = std::make_pair(ind_i, uint(halfedges.size()));
		merge_two_facets(map_facets, h_i->f, h_i->opposite()->f, ref_halfedges, true, true);
	}
}




void Partition::update_graph_when_slicing(vector<Ray *> & rays, list<HalfEdge *> & halfedges_to_trivalent_vertex, Vertex* trivalent_vertex, IndexedEvent* event,
		HalfEdge* extended_halfedge, HalfEdge* intersected_halfedge, list<HalfEdge *> & remaining_halfedges_of_f_i,
		map<HalfEdge*, list<HalfEdge *>::iterator> & access_remaining_halfedges_of_f_i, list<HalfEdge *> & sub_facet_f_i)
{
	// The first step consists in initializing a Vertex that represents the intersection point between the extended and the intersected edges.
	Vertex* v = nullptr;
	bool v_just_created;
	Ray* r_extended = rays[event->intersectant];
	Ray* r_intersected = (event->intersected >= 0 ? rays[event->intersected] : nullptr);

	Point2d pt;
	Vertex::approximate_coordinates(r_extended, nullptr, event->t_intersectant, rows, cols, pt);

	list<Vertex *> neighbors;
	quadtree->search(pt, 1, neighbors);
	for (list<Vertex *>::iterator it_vn = neighbors.begin() ; it_vn != neighbors.end() ; it_vn++) {
		Vertex* v_n = (*it_vn);
		
		bool verbose = false;
		if ((r_intersected == nullptr && v_n->outer_vertex_is_very_close(pt))
			|| (r_intersected != nullptr && v_n->inner_vertex_created_by_colinear_ray(r_extended, r_intersected, rays, verbose))
			|| (r_intersected != nullptr && v_n->same_intersection_point(r_extended, r_intersected))) {
			v = v_n;
			v->events.push_back(event);
			/*if (!v->has_no_duplicated_event()) {
				v->events.pop_back();
			}*/
			break;
		}
	}

	v_just_created = (v == nullptr);
	if (v == nullptr) {
		v = new Vertex(id_vertices, event, pt);
		quadtree->add(v);
		push_back(v);
	}

	// The second step consists in creating new edges
	// At this point, we create the edge (v3 = trivalent_vertex, v)
	double alpha_v3 = r_extended->opposite()->incidence_angle;
	double alpha_v = r_extended->incidence_angle;
	Inner_Edge* v3_v = nullptr;

	// Normally, we should add v3_v bind v3_v to v3 (resp. v) with angle alpha_v3 (resp. alpha_v)
	// However, in some cases it seems necessary to invert the angle, because the extended rays may run in its opposite direction
	// which is for example the case if we recreate an edge which has previously been deleted
	Vec2d v3_v_dir = v->pt - trivalent_vertex->pt;
	if (v3_v_dir.ddot(r_extended->OA) > 0) {
		v3_v = new Inner_Edge(id_edges, r_extended, 1, trivalent_vertex, v, alpha_v3, alpha_v);
	} else {
		v3_v = new Inner_Edge(id_edges, r_extended, 1, trivalent_vertex, v, alpha_v, alpha_v3);
	}

	push_back(v3_v);
	remaining_halfedges_of_f_i.push_back(v3_v->v2_v1);
	access_remaining_halfedges_of_f_i[v3_v->v2_v1] = --remaining_halfedges_of_f_i.end();

	if (v_just_created) {
		// Then the intersected edge needs to be split
		bool intersected_halfedge_v1_v2 = intersected_halfedge->v1_v2;
		Vertex* v1 = intersected_halfedge->e->v1;
		Vertex* v2 = intersected_halfedge->e->v2;
		double alpha_1 = v1->incidence_angle(intersected_halfedge->e);
		double alpha_2 = v2->incidence_angle(intersected_halfedge->e);

		Edge *v1_v = nullptr, *v_v2 = nullptr;
		if (intersected_halfedge->e->type == INNER_EDGE) {
			v1_v = new Inner_Edge(id_edges, r_intersected, 2, v1, v, alpha_1, alpha_2);
			v_v2 = new Inner_Edge(id_edges, r_intersected, 2, v, v2, alpha_1, alpha_2);
			/*if (v1_v->length() < 1e-9 || v_v2->length() < 1e-9) {
				std::cout << "BREAK" << std::endl;
			}*/
			Face* f_adj = intersected_halfedge->opposite()->f;
			if (intersected_halfedge_v1_v2) {
				// The halfedge to replace in the adjacent facet is : (v2 v1) -> (v2 vi + vi v1)
				f_adj->add_non_corner(intersected_halfedge->opposite(), v, v_v2->v2_v1, v1_v->v2_v1);
			} else {
				// The halfedge to replace in the adjacent facet is : (v1 v2) -> (v1 vi + vi v2)
				f_adj->add_non_corner(intersected_halfedge->opposite(), v, v1_v->v1_v2, v_v2->v1_v2);
			}
		} else {
			Image_Boundary boundary = Image_Boundary(event->intersected);
			v1_v = new Outer_Edge(id_edges, boundary, v1, v, alpha_1, alpha_2);
			v_v2 = new Outer_Edge(id_edges, boundary, v, v2, alpha_1, alpha_2);
		}

		// Replaces the intersected edge in the graph by the newly created edges
		remaining_halfedges_of_f_i.erase(access_remaining_halfedges_of_f_i[intersected_halfedge]);
		access_remaining_halfedges_of_f_i.erase(intersected_halfedge);
		erase(intersected_halfedge->e, true);

		push_back(v1_v);
		HalfEdge* intersectable_v1_v = (intersected_halfedge_v1_v2 ? v1_v->v1_v2 : v1_v->v2_v1);
		remaining_halfedges_of_f_i.push_back(intersectable_v1_v);
		access_remaining_halfedges_of_f_i[intersectable_v1_v] = --remaining_halfedges_of_f_i.end();

		push_back(v_v2);
		HalfEdge* intersectable_v_v2 = (intersected_halfedge_v1_v2 ? v_v2->v1_v2 : v_v2->v2_v1);
		remaining_halfedges_of_f_i.push_back(intersectable_v_v2);
		access_remaining_halfedges_of_f_i[intersectable_v_v2] = --remaining_halfedges_of_f_i.end();

	} else if (v->connectivity() == 3) {
		// If we happen to have turned a bivalent vertex into a trivalent vertex by creating v3_v,
		// then we need to update vertex info in the adjacent facet
		if (!v->faces.empty()) {
			// v used to be adjacent to two facets, f_i and the one we want, but we have just deleted f_i
			if (v->faces.size() != 1) {
				std::cout << *v << std::endl;
				std::cout << "faces.size() = " << v->faces.size() << std::endl;
				for (int v_d = 0; v_d < v->directions.size() ; v_d++) {
					Edge* e_d = v->directions[v_d].second->e;
					if (e_d->type == INNER_EDGE) {
						Inner_Edge* ie_d = static_cast<Inner_Edge*>(e_d);
						std::cout << "Edge " << v_d << " : tag = " << ie_d->tag << std::endl;
					}
					Vertex* v_r = (e_d->v1 == v ? e_d->v2 : e_d->v1);
					std::cout << *v_r << std::endl;
				}
			}
			assert(v->faces.size() == 1);
			Face* f_adj = (*v->faces.begin());
			f_adj->set_as_trivalent(v);
		}
	}

	// The third step consists in determining the subfacet that will be later joined to the facet that is
	// on the other side of 'halfedges_to_trivalent_vertex'
	sub_facet_f_i.clear();
	loop_inside_sub_facet(v3_v->v1_v2, sub_facet_f_i);

	// This sub-facet can be used to update the list of remaining edges of f_i
	// The idea is to remove all elements we find in sub_facet_f_i from remaining_edges_of_f_i, 
	// except 'halfedges_to_trivalent_vertex' and v3_v
	for (list<HalfEdge *>::iterator it_h1 = sub_facet_f_i.begin() ; it_h1 != sub_facet_f_i.end() ; it_h1++) {
		HalfEdge* h = (*it_h1);

		// First case to discard
		if (h->e == v3_v) continue;

		// Second case to discard
		list<HalfEdge *>::iterator it_h2 = halfedges_to_trivalent_vertex.begin();
		while (it_h2 != halfedges_to_trivalent_vertex.end()) {
			if ((*it_h2) == h) break;
			++it_h2;
		}
		if (it_h2 != halfedges_to_trivalent_vertex.end()) continue;

		// It's now OK to remove h
		assert(access_remaining_halfedges_of_f_i.find(h) != access_remaining_halfedges_of_f_i.end());
		remaining_halfedges_of_f_i.erase(access_remaining_halfedges_of_f_i[h]);
		access_remaining_halfedges_of_f_i.erase(h);
	}
}



void Partition::loop_inside_sub_facet(HalfEdge* h_e, list<HalfEdge*> & sub_facet)
{
	HalfEdge* h_f = h_e;
	do {
		// Finds the next half-edge in clockwise order
		// Adds it to the facet we are currently defining
		Vertex* v = (h_f->v1_v2 ? h_f->e->v2 : h_f->e->v1);
		uint v_n = uint(v->directions.size());
		int path = 0;
		for (uint i = 0 ; i < v_n ; i++) {
			if (v->directions[i].second->e == h_f->e) {
				path = (i + 1) % v_n;
				break;
			}
		}
		h_f = v->directions[path].second;
		sub_facet.push_back(h_f);
	} while (h_f != h_e);
}


void Partition::split_facet(Face* f, double a, double b, double c)
{
	typedef pair<double, double> Time_Range;

	uint f_id = f->id_face;

	// We would like to split f into two subfacets, separated by a line of equation (D) ax + by + c.
	vector<HalfEdge *> f_edges (f->edges.begin(), f->edges.end());
	vector<pair<Vertex *, Vertex_Type> > f_vertices (f->vertices.begin(), f->vertices.end());
	uint n = uint(f->vertices.size());

	delete f;

	// First we need to find two intersection points
	Vertex *v1 = nullptr, *v2 = nullptr;
	HalfEdge* h1 = nullptr, *h2 = nullptr;
	for (vector<HalfEdge *>::iterator it_h = f_edges.begin() ; it_h != f_edges.end() ; it_h++) {
		HalfEdge* h = (*it_h);
		
		// The current edge is delimited by A and B.
		// We would like to check if there exists a point M(t) = A + t * (B - A) that belongs to (D) with 0 <= t <= 1		
		Vertex* A = h->e->v1;
		Vertex* B = h->e->v2;
		double x_a = A->pt.x, y_a = A->pt.y, x_b = B->pt.x, y_b = B->pt.y;
		double d = a * (x_b - x_a) + b * (y_b - y_a);
		if (fabs(d) < 1e-6) continue;

		double t = -(a * x_a + b * y_a + c) / d;
		Point2d M = Point2d(x_a + t * (x_b - x_a), y_a + t * (y_b - y_a));
		Vec2d tmppp;
		Vec2d delta = (v1 != nullptr ? static_cast<Vec2d>(v1->pt - M) : tmppp);
		if (0 <= t && t <= 1) {
			if (v1 == nullptr) {
				v1 = new Vertex(id_vertices, nullptr, M);
				quadtree->add(v1);
				push_back(v1);
				h1 = h;
			} else {
				if (cv::norm(delta) > 1e-2) {
					v2 = new Vertex(id_vertices, nullptr, M);
					quadtree->add(v2);
					push_back(v2);
					h2 = h;
				}
			}
		} else {
			Vec2d AM = M - A->pt, BM = M - B->pt;
			if (cv::norm(AM) < 1e-2) {
				if (v1 == nullptr) {
					v1 = A;
				} else {
					if (cv::norm(delta) > 1e-2) {
						v2 = A;
					}
				}
			} else if (cv::norm(BM) < 1e-2) {
				if (v1 == nullptr) {
					v1 = B;
				} else {
					if (cv::norm(delta) > 1e-2) {
						v2 = B;
					}
				}
			}
		}

		if (v1 != nullptr && v2 != nullptr) break;
	}

	assert(v1 != nullptr && v2 != nullptr);

	// At this stage of the algorithm, we identified two vertices that represent the intersections of f and D.
	// In case v1 and v2 don't correspond to preexisting points of the facet, it will be necessary to split the halfedges h1 and h2.
	// Anyway, first of all we create the artificial edge (v1 v2).
	double alpha_v1 = atan2(v2->pt.y - v1->pt.y, v2->pt.x - v1->pt.x);
	double alpha_v2 = (alpha_v1 <= 0 ? alpha_v1 + PI : alpha_v1 - PI);
	double e_a = -sin(alpha_v1), e_b = cos(alpha_v1), e_c = -e_a * v1->pt.x - e_b * v1->pt.y;
	Edge* e = new Artificial_Edge(id_edges, e_a, e_b, e_c, v1, v2, alpha_v1, alpha_v2);
	push_back(e);

	// Splits the intersected halfedges
	if (h1 != nullptr) {
		split_edge_and_update_adjacent_facet(h1, v1);
	}

	if (h2 != nullptr) {
		split_edge_and_update_adjacent_facet(h2, v2);
	}

	// Creates the two subfacets.
	list<HalfEdge *> halfedges_1, halfedges_2;
	loop_inside_sub_facet(e->v1_v2, halfedges_1);
	loop_inside_sub_facet(e->v2_v1, halfedges_2);

	Face* f1 = new Face(id_faces, halfedges_1);
	Face* f2 = new Face(id_faces, halfedges_2);
	faces.push_back(f1);
	faces.push_back(f2);
}


void Partition::split_edge_and_update_adjacent_facet(HalfEdge* h, Vertex* v)
{
	// We want to access the opposite facet of h which is going to take one more vertex in its definition : v.
	Face* f = h->opposite()->f;

	bool v1_v2 = h->v1_v2;
	Vertex *v1 = h->e->v1, *v2 = h->e->v2;
	double alpha_v1 = v1->incidence_angle(h->e), alpha_v2 = v2->incidence_angle(h->e);
	Edge *v1_v = nullptr, *v_v2 = nullptr;

	if (h->e->type == INNER_EDGE) {
		Inner_Edge* e = static_cast<Inner_Edge *>(h->e);
		set<Ray *> support;
		e->get_supporting_rays(support);
		v1_v = new Inner_Edge(id_edges, support, 3, v1, v, alpha_v1, alpha_v2);
		v_v2 = new Inner_Edge(id_edges, support, 3, v, v2, alpha_v1, alpha_v2);
	} else if (h->e->type == OUTER_EDGE) {
		Outer_Edge* e = static_cast<Outer_Edge *>(h->e);
		Image_Boundary boundary = e->boundary;
		v1_v = new Outer_Edge(id_edges, boundary, v1, v, alpha_v1, alpha_v2);
		v_v2 = new Outer_Edge(id_edges, boundary, v, v2, alpha_v1, alpha_v2);
	} else if (h->e->type == ARTIFICIAL_EDGE) {
		Artificial_Edge* e = static_cast<Artificial_Edge *>(h->e);
		double ae_a, ae_b, ae_c;
		e->line_equation(ae_a, ae_b, ae_c);
		v1_v = new Artificial_Edge(id_edges, ae_a, ae_b, ae_c, v1, v, alpha_v1, alpha_v2);
		v_v2 = new Artificial_Edge(id_edges, ae_a, ae_b, ae_c, v, v2, alpha_v1, alpha_v2);
	}

	if (f != nullptr) {
		if (v1_v2) {
			// Then the order or the vertices in the adjacent facet is (w2 v w1)
			f->add_non_corner(h->opposite(), v, v_v2->v2_v1, v1_v->v2_v1);
		} else {
			// Reverted order : (w1 v w2)
			f->add_non_corner(h->opposite(), v, v1_v->v1_v2, v_v2->v1_v2);
		}
	}
	push_back(v1_v);
	push_back(v_v2);
	erase(h->e, true);
}



#if 0
void Partition::fill_matrix_of_iterators()
{
    for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end() ; it_f++) {
        Face* f = (*it_f);
        for (list<Point2i>::iterator it_p = f->pixels.begin() ; it_p != f->pixels.end() ; it_p++) {
            F(I.rows - it_p->y, it_p->x) = it_f;
        }
    }
}
#endif


int Partition::count_invalid_faces()
{
    int count = 0;
    for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
        if (!(*it_f)->is_simple(false)) {
            ++count;
        }
    }
    return count;
}


void Partition::draw_intermediate_graph(Matrix<uchar> & I, Matrix<uchar> & J, list<Vertex *> & outer_vertices, list<Vertex *> & inner_vertices,
    vector<list<Edge *> > & outer_edges, vector<list<Edge *> > & inner_edges)
{
    // Copies an enlarged version of the image
    uint f = 4;
    J = Matrix<uchar>(f * I.rows, f * I.cols, 3);
    for (uint i = 0 ; i < I.rows ; i++) {
        for (uint j = 0 ; j < I.cols ; j++) {
            for (uint k = 0 ; k < f ; k++) {
                for (uint l = 0 ; l < f ; l++) {
                    for (uint c = 0 ; c < 3 ; c++) J(f * i + k, f * j + l, c) = I(i, j, c);
                }
            }
        }
    }

    uchar blue[3] = {0, 0, 255};
    uchar green[3] = {0, 255, 0};
    uchar red[3] = {255, 0, 0};
    uchar yellow[3] = {255, 255, 0};

    // Outputs edges
    for (uint i = 0; i < outer_edges.size(); i++) {
        for (list<Edge *>::iterator it = outer_edges[i].begin(); it != outer_edges[i].end(); it++) {
            Edge* e = *it;
            Point2d & pt1 = e->v1->pt;
            Point2d & pt2 = e->v2->pt;
            Point2i pi1 = Point2i(int(round(jclamp(0, f * pt1.x, J.cols - 1))), int(round(jclamp(0, J.rows - f * pt1.y, J.rows - 1))));
            Point2i pi2 = Point2i(int(round(jclamp(0, f * pt2.x, J.cols - 1))), int(round(jclamp(0, J.rows - f * pt2.y, J.rows - 1))));
            J.line(pi1.y, pi1.x, pi2.y, pi2.x, red);
        }
    }
    for (uint i = 0; i < inner_edges.size(); i++) {
        for (list<Edge *>::iterator it = inner_edges[i].begin(); it != inner_edges[i].end(); it++) {
            Edge* e = *it;
            Point2d & pt1 = e->v1->pt;
            Point2d & pt2 = e->v2->pt;
            Point2i pi1 = Point2i(int(round(jclamp(0, f * pt1.x, J.cols - 1))), int(round(jclamp(0, J.rows - f * pt1.y, J.rows - 1))));
            Point2i pi2 = Point2i(int(round(jclamp(0, f * pt2.x, J.cols - 1))), int(round(jclamp(0, J.rows - f * pt2.y, J.rows - 1))));
            J.line(pi1.y, pi1.x, pi2.y, pi2.x, red);
        }
    }

    // Draws vertices
    for (int l = 0 ; l < 2 ; l++) {
        list<Vertex *> & reference_to_list = (l == 0 ? inner_vertices : outer_vertices);
        for (list<Vertex *>::iterator it = reference_to_list.begin(); it != reference_to_list.end(); it++) {
            Vertex* v = *it;
            Point2d & pt = v->pt;
            Point2i pi = Point2i(int(round(jclamp(0, f * pt.x, J.cols - 1))), int(round(jclamp(0, J.rows - f * pt.y, J.rows - 1))));
            if (v->connectivity() == 2) {
                for (int c = 0 ; c < 3 ; c++) J(pi.y, pi.x, c) = green[c];
            } else if (v->connectivity() == 3) {
                for (int c = 0 ; c < 3 ; c++) J(pi.y, pi.x, c) = blue[c];
            } else {
                for (int c = 0 ; c < 3 ; c++) J(pi.y, pi.x, c) = yellow[c];
            }
        }
    }
}


void Partition::draw_edges(Matrix<uchar> & I, Matrix<uchar> & J, double f)
{
	if (f == 1.0) {
		J = I;
	} else {
		J = Matrix<uchar>(I, f);
	}

    uchar blue[3] = {0, 0, 255};
    uchar dark_green[3] = {0, 128, 0};
    uchar red[3] = {255, 0, 0};
    uchar yellow[3] = {255, 255, 0};

    Edge* e = edges_head;
    while (e != NULL) {
        Point2d & pt1 = e->v1->pt;
        Point2d & pt2 = e->v2->pt;
        Point2i pi1 = Point2i(int(jclamp(0, f * pt1.x, J.cols - 1)), int(jclamp(0, J.rows - f * pt1.y, J.rows - 1)));
        Point2i pi2 = Point2i(int(jclamp(0, f * pt2.x, J.cols - 1)), int(jclamp(0, J.rows - f * pt2.y, J.rows - 1)));
        J.line(pi1.y, pi1.x, pi2.y, pi2.x, red);
        e = e->e_next;
    }

    Vertex* v = vertices_head;
    while (v != NULL) {
        Point2d & pt = v->pt;
        Point2i pi = Point2i(int(jclamp(0, f * pt.x, J.cols - 1)), int(jclamp(0, J.rows - f * pt.y, J.rows - 1)));
        if (v->connectivity() == 2) {
            for (int c = 0 ; c < 3 ; c++) J(pi.y, pi.x, c) = yellow[c];
        } else if (v->connectivity() == 3) {
            for (int c = 0 ; c < 3 ; c++) J(pi.y, pi.x, c) = blue[c];
        } else {
            for (int c = 0 ; c < 3 ; c++) J(pi.y, pi.x, c) = dark_green[c];
        }
        v = v->v_next;
    }
}


void Partition::draw_faces(Matrix<uchar> & J, bool color_after_thinness)
{
    Matrix<int> B(rows, cols, 1);
    for (uint i = 0 ; i < rows ; i++) {
        for (uint j = 0 ; j < cols ; j++) {
            B(i, j) = -1;
        }
    }
    // Assigns a facet index to all pixels
    for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
        Face* f = (*it_f);
		f->find_pixels_inside_facet();
        int k = f->id_face;
        for (list<Point2i>::iterator it_p = f->pixels.begin(); it_p != f->pixels.end(); it_p++) {
            uint i = uint(int(B.rows) - it_p->y);
            uint j = uint(it_p->x);
            assert(B(i, j) == -1);
            B(i, j) = k;
        }
    }

    // Prepares palette
    vector<Vec3b> colors;
    colors.reserve (faces.size());
	if (color_after_thinness) {
		for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
			double thinness = (*it_f)->thinness;
			colors.push_back(ColorMap::blue_white_red(thinness / 10));
		}
	} else {
	    std::default_random_engine generator;
	    std::uniform_int_distribution<int> uniform_dist (128, 255);
		for (int i = 0; i < faces.size(); i++) {
			int r = uniform_dist(generator);
			int g = uniform_dist(generator);
			int b = uniform_dist(generator);
			colors.push_back(Vec3b(r, g, b));
		}
	}
    // Colors facets
    J = Matrix<uchar>(rows, cols, 3);
    for (uint i = 0; i < B.rows; i++) {
        for (uint j = 0; j < B.cols; j++) {
            int k = B(i, j);
            assert(k != -1);
            for (int c = 0 ; c < 3 ; c++) J(i, j, c) = uchar(colors[k][c]);
        }
    }
}


void Partition::save_edges(Matrix<uchar> & I, std::string &directory, std::string &basename, double f)
{
    Matrix<uchar> J;
    draw_edges(I, J, f);
	std::string filename;

    if (params->lsd_create_additional) {
        filename += directory + '\\' + basename + "_edges_cross_propagation" + (params->suffix == "" ? "" : "_" + params->suffix) + ".tiff";
    } else {
        filename += directory + '\\' + basename + "_edges" + (params->suffix == "" ? "" : "_" + params->suffix) + ".tiff";
    }

    J.write_uchar(filename);
}


void Partition::save_faces(std::string &directory, std::string &basename, bool color_after_thinness)
{
    Matrix<uchar> J;
	draw_faces(J, color_after_thinness);
	std::string filename = directory + '\\' + basename + "_faces" + (params->suffix == "" ? "" : "_" + params->suffix) + ".tiff";
    J.write_uchar(filename);
}


void Partition::save_graph_definition(std::string &directory, std::string &basename)
{
	int id_vertex = 0;
	int id_edge = 0;
	int id_face = 0;
    std::string filename = directory + '\\' + basename + "_graph.txt";

	FILE* file = fopen(filename.c_str(), "w");
	if (file != NULL) {
		// If we have been able to successfully open a file, then we can start writing in it
		// First of all we set the size of the associated image
		// fprintf(file, "size\n");
        fprintf(file, "%i %i\n", rows, cols);
		fprintf(file, "%i %i %i\n", int(v_size), int(e_size), int(faces.size()));

		// We list the vertices
		fprintf(file, "vertices\n");
		Vertex* v = vertices_head;
		while (v != NULL) {
			v->id_vertex = id_vertex++;
			fprintf(file, "%i %lf %lf\n", v->id_vertex, v->pt.x, v->pt.y);
			v = v->v_next;
		}

		// We list the edges
		fprintf(file, "edges\n");
		Edge* e = edges_head;
		while (e != NULL) {
			e->id_edge = id_edge++;
			fprintf(file, "%i %i %i\n", e->id_edge, e->v1->id_vertex, e->v2->id_vertex);
			e = e->e_next;
		}

		// We list the faces
		fprintf(file, "faces\n");
		for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
			Face* f = (*it_f);
			f->id_face = id_face++;
			fprintf(file, "%i %i", f->id_face, int(f->edges.size()));
			for (list<HalfEdge *>::iterator it_h = f->edges.begin(); it_h != f->edges.end(); it_h++) {
				HalfEdge* h = (*it_h);
				fprintf(file, " %i %i", h->e->id_edge, int(h->v1_v2));
			}
			fprintf(file, "\n");
		}
		fclose(file);
	}
}


void Partition::save_liuyuns_input(std::string &directory, std::string &basename, vector<Segment *> & segments)
{
    std::string common_basename = directory + '\\' + basename;

    // First output : the segments
    FILE* file_segv = fopen((common_basename + ".segV").c_str(), "w");
    if (file_segv != NULL) {
        fprintf(file_segv, "%i\n", int(segments.size()));
        for (int i = 0 ; i < segments.size() ; i++) {
            Segment* s_i = segments[i];
            int index = i;
            int junction_type = -1;
            Point2d end1 = Point2d(s_i->finalEnd1.x, rows - s_i->finalEnd1.y);
            Point2d end2 = Point2d(s_i->finalEnd2.x, rows - s_i->finalEnd2.y);
            Point2d barycenter = (end1 + end2) / 2.0;
            Vec2d dir = normalize(Vec2d(end2 - end1));
            double angle = atan2(dir[1], dir[0]) * 180 / PI;
            double length = s_i->length;
            double fitting_quality = 0.8;
            double a = -dir[1], b = dir[0];
            double c = -a * barycenter.x - b * barycenter.y;
            int type = 1;
            fprintf(file_segv, "%i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %i\n",
                index, junction_type,
                end1.x, end1.y, end2.x, end2.y, barycenter.x, barycenter.y,
                dir[0], dir[1], length, fitting_quality, angle, a, b, c, type);
        }
        fclose(file_segv);
    }

    // Second output : the cells
    vector<int> num;
    vector<Point2d> vcell;
    vector<Point2d> off;
    vector<int> seedtype;
    double global_average_dist = 0;
    for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
        Face* f = (*it_f);
        if (f->pixels.size() == cols * rows) continue;

        // 'num' receives the number of edges of the cell
        num.push_back(int(f->edges.size()));

        // 'seedtype' receives 1
        seedtype.push_back(1);

        // We build a vector that lists the vertices of the facet in counter-clockwise order
        vector<Point2d> vertices;
        for (list<pair<Vertex *, Vertex_Type> >::iterator it_v = f->vertices.begin() ; it_v != f->vertices.end() ; it_v++) {
            vertices.push_back(Point2d(it_v->first->pt.x, rows - it_v->first->pt.y));
        }

        // 'v_off' receives the center of the cells, 'v_cell' the list of cells that defines the current cell
        Point2d center = Point2d(0, 0);
        for (int i = 0; i < vertices.size(); i++) center += vertices[i];
        center /= int(vertices.size());

        double average_dist = 0;
        for (int i = 0; i < vertices.size(); i++) average_dist += norm(vertices[i] - center);
        average_dist /= int(vertices.size());
        global_average_dist += average_dist;

        for (int i = 0 ; i < vertices.size() ; i++) vcell.push_back(vertices[i]);
        off.push_back(center);
    }
    global_average_dist /= faces.size();


    FILE* file_num = fopen((common_basename + ".num").c_str(), "w");
    if (file_num != NULL) {
        for (int i = 0 ; i < faces.size() ; i++) {
            fprintf(file_num, "%i\n", num[i]);
        }
        fclose(file_num);
    }

    FILE* file_vcell = fopen((common_basename + ".vcell").c_str(), "w");
    if (file_vcell != NULL) {
        int offset = 0;
        for (int i = 0 ; i < faces.size() ; i++) {
            for (int j = 0 ; j < num[i] ; j++) {
                fprintf(file_vcell, "%lf %lf\n", vcell[offset + j].x, vcell[offset + j].y);
            }
            offset += num[i];
        }
        fclose(file_vcell);
    }

    FILE* file_off = fopen((common_basename + ".off").c_str(), "w");
    if (file_off != NULL) {
        fprintf(file_off, "%i\n%i\n", cols, rows);
        fprintf(file_off, "%i\n", int(global_average_dist));
        for (int i = 0 ; i < faces.size() ; i++) {
            fprintf(file_off, "%lf %lf\n", off[i].x, off[i].y);
        }
        fclose(file_off);
    }

    FILE* file_seedtype = fopen((common_basename + ".seedtype").c_str(), "w");
    if (file_seedtype != NULL) {
        for (int i = 0; i < faces.size(); i++) {
            fprintf(file_seedtype, "%i\n", seedtype[i]);
        }
        fclose(file_seedtype);
    }

    // Third output : indices of segments that delimit the cells
    FILE* file_segments_indices = fopen((common_basename + ".segments_indices").c_str(), "w");
    if (file_segments_indices != NULL) {
        fprintf(file_segments_indices, "%i\n", int(faces.size()));
        for (list<Face *>::iterator it_f = faces.begin() ; it_f != faces.end() ; it_f++) {
            set<Segment *> s_segments;
            (*it_f)->get_supporting_segments(s_segments);
            fprintf(file_segments_indices, "%i", int(s_segments.size()));
            for (set<Segment *>::iterator it_s = s_segments.begin() ; it_s != s_segments.end() ; it_s++) {
                fprintf(file_segments_indices, " %i", (*it_s)->index);
            }
            fprintf(file_segments_indices, "\n");
        }
        fclose(file_segments_indices);
    }

    // Fourth output : neighboring faces
    FILE* file_neighboring_faces = fopen((common_basename + ".neighbourInd").c_str(), "w");
    if (file_neighboring_faces != NULL) {
        for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
            set<Face *> neighbors;
            (*it_f)->get_neighboring_faces(neighbors);
            fprintf(file_neighboring_faces, "%i\n", int(neighbors.size()));
            for (set<Face *>::iterator it_nf = neighbors.begin(); it_nf != neighbors.end(); it_nf++) {
                fprintf(file_neighboring_faces, "%i ", (*it_nf)->id_face);
            }
            fprintf(file_neighboring_faces, "\n");
        }
        fclose(file_neighboring_faces);
    }
}



void Partition::save_boundaries(std::string & directory, std::string & basename)
{
    const uchar white[3] = {255, 255, 255};

    Matrix<uchar> B(rows, cols, 3);
    B.set(0);

    Edge* e = edges_head;
    while (e != NULL) {
        Point2d & pt1 = e->v1->pt;
        Point2d & pt2 = e->v2->pt;
        Point2i pi1 = Point2i(int(jclamp(0, pt1.x, int(cols) - 1)), int(jclamp(0, int(rows) - pt1.y, int(rows) - 1)));
        Point2i pi2 = Point2i(int(jclamp(0, pt2.x, int(cols) - 1)), int(jclamp(0, int(rows) - pt2.y, int(rows) - 1)));
        bool is_boundary = (pi1.x == 0 && pi2.x == 0) || (pi1.x == int(cols) - 1 && pi2.x == int(cols) - 1)
                || (pi1.y == 0 && pi2.y == 0) || (pi1.y == int(rows) - 1 && pi2.y == int(rows) - 1);
        if (!is_boundary) {
            B.line(pi1.y, pi1.x, pi2.y, pi2.x, white);
        }
        e = e->e_next;
    }

    std::string filename = directory + '\\' + basename + "_boundaries" + (params->suffix == "" ? "" : "_" + params->suffix) + ".tiff";
    B.write_uchar(filename);
}



void Partition::save_labels(std::string &directory, std::string &basename)
{
    Matrix<int> B(rows, cols, 1);
    for (uint i = 0 ; i < rows ; i++) {
        for (uint j = 0 ; j < cols ; j++) {
            B(i, j) = -1;
        }
    }
    // Assigns a facet index to all pixels
    for (list<Face *>::iterator it_f = faces.begin(); it_f != faces.end(); it_f++) {
        Face* f = (*it_f);
        int k = f->id_face;
		f->find_pixels_inside_facet();
        for (list<Point2i>::iterator it_p = f->pixels.begin(); it_p != f->pixels.end(); it_p++) {
            uint i = uint(int(B.rows) - it_p->y);
            uint j = uint(it_p->x);
            //assert(B(i, j) == -1);
            B(i, j) = k;
        }
    }

	std::string filename = directory + '\\' + basename + "_labels" + (params->suffix == "" ? "" : "_" + params->suffix) + ".txt";
    FILE* file = fopen(filename.c_str(), "w");
    if (file != NULL) {
        for (uint i = 0 ; i < rows ; i++) {
            for (uint j = 0 ; j < cols ; j++) {
				assert(B(i, j) != -1);
                fprintf(file, "%i\n", B(i, j));
            }
        }
        fclose(file);
    }
}


void Partition::save_parameters(std::string &directory, std::string &basename, int N)
{
	std::string filename = directory + '\\' + basename + "_params" + (params->suffix == "" ? "" : "_" + params->suffix) + ".txt";
	std::filebuf fb;
	fb.open(filename, std::ios::out);
	std::ostream os(&fb);
	os << "** LSD" << std::endl;
	os << "Scale : " << params->lsd_scale << std::endl;
	os << "Sigma scale : " << params->lsd_sigma_scale << std::endl;
	os << "Angle tol : " << params->lsd_angle << std::endl;
	os << "Quant : " << params->lsd_quant << std::endl;
	os << "Log-eps : " << params->lsd_log_eps << std::endl;
	os << "Density : " << params->lsd_density << std::endl;
	os << std::endl;
	os << "** Regularization" << std::endl;
	os << "Enabled : " << params->rega_regp_enabled << std::endl;
	os << "Algorithm : " << params->rega_method << std::endl;
	os << "dtheta max : " << params->rega_angle_const << std::endl;
	os << "dt max :" << params->regp_trans_const << std::endl;
	os << std::endl;
	os << "** Propagation" << std::endl;
	os << "K : " << params->prop_ttl << std::endl;
	os << "Extra propagation : " << params->prop_extra_enabled << std::endl;
	os << std::endl;
	os << "** Merge" << std::endl;
	os << "Thinness : " << params->merge_min_thinness << std::endl;
	os << std::endl;
	os << "** Cells : " << N << std::endl;
	fb.close();
}


/*
void Partition::read_graph_definition(std::string & filename)
{
	std::vector<Vertex *> read_vertices;
	std::vector<Edge *> read_edges;
	std::vector<Face *> read_faces;

	try {
		// Opens the provided file : if it exists, we start parsing
		FILE* file = fopen(filename.c_str(), "r");
		if (file == NULL) throw std::ios_base::failure("Couldn't read the provided file graph.txt");

		// First of all we should find the number of vertices, edges and faces
		int index;
		char str[9];
		int nb_vertices, nb_edges, nb_faces;
		if (fscanf(file, "%i %i\n", &rows, &cols) != 2) {
			throw std::logic_error("Badly formatted file graph.txt : number of rows and columns of the image not found");
		}
		Quadtree* quadtree = new Quadtree(0, cols, 0, rows);

		if (fscanf(file, "%i %i %i\n", &nb_vertices, &nb_edges, &nb_faces) != 3) {
			throw std::logic_error("Badly formatted file graph.txt : number of elements not found");
		} else {
			read_vertices.reserve(nb_vertices);
			read_edges.reserve(nb_edges);
			read_faces.reserve(nb_faces);
		}

		// Reads the list of vertices
		if (fscanf(file, "%s\n", str) == 1) {
			std::cout << strcmp(str, "vertices") << std::endl;
			if (strcmp(str, "vertices")) {
				throw std::logic_error("Badly formatted file graph.txt : expected keyword 'vertices' not found");
			}
		}
		for (int i = 0; i < nb_vertices; i++) {
			double x, y;
			if (fscanf(file, "%i %lf %lf\n", &index, &x, &y) != 3) {
				throw std::logic_error("Badly formatted file graph.txt : expected a definition of a vertex");
			} else {
				read_vertices.push_back(new Vertex(index, NULL, x, y));
			}
		}

		// Reads the list of edges
		if (fscanf(file, "%s\n", str) == 1) {
			std::cout << strcmp(str, "edges") << std::endl;
			if (strcmp(str, "edges")) {
				throw std::logic_error("Badly formatted file graph.txt : expected keyword 'edges' not found");
			}
		}
		for (int i = 0; i < nb_edges; i++) {
			int v1, v2;
			if (fscanf(file, "%i %i %i\n", &index, &v1, &v2) != 3) {
				throw std::logic_error("Badly formatted file graph.txt : expected a definition of an edge");
			} else {
				read_edges.push_back(new Edge(index, read_vertices[v1], read_vertices[v2], rows, cols));
			}
		}

		// Reads the list of faces
		if (fscanf(file, "%s\n", str) == 1) {
			std::cout << strcmp(str, "faces") << std::endl;
			if (strcmp(str, "faces")) {
				throw std::logic_error("Badly formatted file graph.txt : expected keyword 'faces' not found");
			}
		}
		for (int i = 0; i < nb_faces; i++) {
			int lines;
			std::list<HalfEdge *> halfedges;
			if (fscanf(file, "%i %i", &index, &lines) != 2) {
				throw std::logic_error("Badly formatted file graph.txt : expected a definition of a face");
			}
			for (int j = 0; j < lines; j++) {
				int e, direction;
				if (((j < lines - 1) && (fscanf(file, "%i %i", &e, &direction) != 2)) || ((j == lines - 1) && (fscanf(file, "%i %i\n", &e, &direction) != 2))) {
					throw std::logic_error("Badly formatted file graph.txt : expected a list of halfedges");
				} else {
					halfedges.push_back((direction == 1 ? read_edges[e]->v1_v2 : read_edges[e]->v2_v1));
				}
			}
			read_faces.push_back(new Face(index, halfedges));
		}

		// Copies the elements we read to their final containers
		for (int i = 0; i < read_vertices.size(); i++) push_back(read_vertices[i]);
		v_size = uint(read_vertices.size());
		read_vertices.clear();
		for (int i = 0; i < read_edges.size(); i++) push_back(read_edges[i]);
		e_size = uint(read_edges.size());
		read_edges.clear();
		faces = std::list<Face *>(read_faces.begin(), read_faces.end());
		read_faces.clear();

		id_vertices = nb_vertices;
		id_edges = nb_edges;
		id_faces = nb_faces;
		is_valid = true;

		fclose(file);

		// We're done
		std::cout << "Successfully read the provided file graph.txt" << std::endl;

	} catch (const std::exception & e) {
		std::cout << "Error : " + std::string(e.what()) << std::endl;
		for (int i = 0; i < read_faces.size(); i++) delete read_faces[i];
		read_faces.clear();
		for (int i = 0; i < read_edges.size(); i++) delete read_edges[i];
		read_edges.clear();
		for (int i = 0; i < read_vertices.size(); i++) delete read_vertices[i];
		read_vertices.clear();
	}
}
*/

/*
void Partition::stats()
{
	std::cout << "** Statistics : computed over " << faces.size() << " faces" << std::endl;
	faces_sizes(params, faces);
	if (params.is_disparity_map) {
		disparity_mask_with_satellite_image(params, faces);
		disparity_mask_with_height_map(params, faces);
		disparity_least_squares(params, faces);
	}
}
*/
