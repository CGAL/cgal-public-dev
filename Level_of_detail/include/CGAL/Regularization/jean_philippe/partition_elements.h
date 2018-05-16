#pragma once
#include <list>
#include <vector>
#include <set>
#include <opencv2/core.hpp>
#include "indexed_event.h"
#include "quadtree_point.h"
#include "segment_ray.h"


using cv::Point2d;
using cv::Point2i;
using std::list;
using std::vector;
using std::pair;
using std::set;


class Vertex;

class HalfEdge;

class Face;


typedef enum { INNER_EDGE, OUTER_EDGE, ARTIFICIAL_EDGE } Edge_Type;


class Edge
{
protected:
	Edge(int & _id_edge, Vertex* _v1, Vertex* _v2, double & alpha_1, double & alpha_2);

	//Edge(const Edge & E);

public:
	virtual ~Edge();

public:
	void init_vertices_and_halfedges(Vertex* _v1, double & _alpha_1, Vertex* _v2, double & _alpha_2);

	double length();

    virtual void line_equation(double & a, double & b, double & c) = 0;

public:
	int id_edge;

	Edge_Type type;
	
	Vertex* v1;
	Vertex* v2;

private:
	double alpha_1;
	double alpha_2;

public:
	HalfEdge* v1_v2;
	HalfEdge* v2_v1;

	Edge* e_prev;
	Edge* e_next;

	//int tag;
};


class Artificial_Edge : public Edge
{
public:
	Artificial_Edge(int & _id_edge, double _a, double _b, double _c, Vertex* _v1, Vertex* _v2, double _alpha_1, double _alpha_2);

	//Artificial_Edge(const Artificial_Edge & E);

	~Artificial_Edge();

	void line_equation(double &a, double &b, double &c);
private:
	double a;
	double b;
	double c;
};


class Inner_Edge : public Edge
{
public:
	Inner_Edge(int & _id_edge, Ray* _r, int _tag, Vertex* _v1, Vertex* _v2, double alpha_1, double alpha_2, double _t1 = 0, double _t2 = 0);

	Inner_Edge(int & _id_edge, set<Ray*> & _support, int _tag, Vertex* _v1, Vertex* _v2, double alpha_1, double alpha_2);

	//Inner_Edge(const Inner_Edge & E);

	~Inner_Edge();

public:
	void add_supporting_ray(Ray* s);
	
	void add_supporting_rays(set<Ray *> & s);

	void get_supporting_rays(set<Ray *> & s);

	void get_supporting_segments(set<Segment *> & s);

	Ray* get_front_supporting_ray();

	Segment* get_front_supporting_segment();

	void line_equation(double &a, double &b, double &c);

	void time_range(Ray* r, double &t_1, double &t_2);
public:
	set<Ray *> rays;
	int tag;
};



class Outer_Edge : public Edge
{
public:
	Outer_Edge(int & _id_edge, Image_Boundary _boundary, Vertex* _v1, Vertex* _v2, double alpha_1, double alpha_2, double _t1 = 0, double _t2 = 0);

	//Outer_Edge(const Outer_Edge & E);

	~Outer_Edge();

	void line_equation(double &a, double &b, double &c);

	void time_range(double &t_1, double &t_2);
public:
	Image_Boundary boundary;
};



class HalfEdge
{
public:
	HalfEdge(Edge* _e, bool _v1_v2);

	~HalfEdge();

	void set(Face* _f = nullptr);

	HalfEdge* opposite();

	static void intersects_if_extended(HalfEdge* h, list<HalfEdge *> & intersectable_halfedges, Point2d & intersection_point, HalfEdge* & intersected);

public:
	Edge* e;
	bool v1_v2;
	Face* f;
};


class Vertex : public Quadtree_Point
{
public:
	Vertex(int & _id_vertex, double x, double y);

	Vertex(int & _id_vertex, IndexedEvent* _event, double x, double y);

	Vertex(int & _id_vertex, IndexedEvent* _event, Point2d & _pt);

	//Vertex(const Vertex & V);

	~Vertex();

	Point2d keypoint() { return pt; }

	double incidence_angle(const Edge* e);
	
	void add(double alpha, HalfEdge *h);

	void add(Face *f);

	void remove(HalfEdge *h);

	void remove(Face *f);

	int connectivity();

	bool outer_vertex_is_very_close(Point2d & _pt);

	bool inner_vertex_created_by_colinear_ray(Ray* r_intersectant, Ray* r_intersected, vector<Ray *> & rays, bool verbose);

	bool inner_vertex_barycenter_of_brother_segment(Ray* r_intersectant, Ray* r_intersected, vector<Ray *> & rays);

	bool same_intersection_point(Ray* r_intersectant, Ray* r_intersected);

	bool has_no_duplicated_event();

	HalfEdge* most_orthogonal_halfedge(Vec2d & dir);

	void most_orthogonal_halfedge(Vertex* v_prev, Vertex* v_next, Vec2d & v_prev_next, double & abs_angle, HalfEdge* & h);

	static bool right_turn(HalfEdge* h_prev, Vertex* v_curr, HalfEdge* h_next);

	static bool right_turn(Vertex* v_prev, Vertex* v_curr, Vertex* v_next);

	static bool approximate_coordinates(Ray* r_i, Ray* r_j, double t_i, int rows, int cols, Point2d & pt);

private:
	static bool round_coordinates(Point2d & pt, int rows, int cols, double eps);

public:
	int id_vertex;
	list<IndexedEvent*> events;
	Point2d pt;
	vector<pair<double, HalfEdge*> > directions;
	set<Face *> faces;

	Vertex* v_prev;
	Vertex* v_next;
};


typedef enum {
	UNDETERMINED = -1,
	CORNER = 0,
	NON_CORNER_TRIVALENT = 1,
	NON_CORNER_BIVALENT = 2
} Vertex_Type;


class Face
{
public:
	Face(int & _id_face, list<HalfEdge *> & _edges);

	~Face();

	bool is_simple(bool verbose);

	void find_pixels_inside_facet();
	
	void print_details();

private:
	void list_vertices();

	void find_leftmost_and_rightmost_edges(double & x_min, double & x_max, double & y_min, double & y_max, list<HalfEdge *>::iterator & it_l, list<HalfEdge *>::iterator & it_r);

	void find_pixels(double x_min, double x_max, double y_min, double y_max, list<HalfEdge *>::iterator it_l, list<HalfEdge *>::iterator it_r, bool verbose = false);

public:
	void get_supporting_segments(set<Segment *> & supporting_segments);

	void get_neighboring_faces(set<Face *> & neighboring_faces);

	void classify();

	void set_thinness();

	void add_non_corner(HalfEdge* h, Vertex* v, HalfEdge* h_1, HalfEdge* h_2);

	void remove_non_corners(Vertex* v1, Vertex* v2, HalfEdge* h);

	void set_as_bivalent(Vertex* v);

	void set_as_trivalent(Vertex* v);

	bool convex_union(list<pair<uint, uint> > & adjacent_edges);

	bool convex_cut(list<pair<uint, uint> > & lines_of_cut);

	static void merged_halfedges(Face* f_1, Face* f_2, Edge* e, list<HalfEdge *> & merged);

	static void merged_halfedges(list<HalfEdge *> & v_1, list<HalfEdge *> & v_2, Edge* e, list<HalfEdge *> & merged);

	static void merged_halfedges(Face* f_1, Face* f_2, Edge* e_1, Edge* e_2, list<HalfEdge *> & merged);

	static void merged_halfedges(list<HalfEdge *> & l_1, list<HalfEdge *> & l_2, Edge* e_1, Edge* e_2, list<HalfEdge *> & merged);
public:
	int id_face;
	list<HalfEdge *> edges;
	list<pair<Vertex *, Vertex_Type> > vertices;

	double thinness;
	list<Point2i> pixels;
};


inline std::ostream & operator<< (std::ostream & stream, Vertex & v)
{
	stream << "Vertex " << v.id_vertex << " " << v.pt << std::endl;
	for (int d = 0 ; d < v.directions.size() ; d++) {
		double angle = v.directions[d].first * 180 / CV_PI;
		Edge* e = v.directions[d].second->e;
		bool v1_v2 = v.directions[d].second->v1_v2;
		if (v1_v2) {
			stream << "   " << angle << " : Edge " << e->id_edge << " between vertices " << e->v1->id_vertex << " " << e->v1->pt << " and " << e->v2->id_vertex << " " << e->v2->pt << std::endl;
		} else {
			stream << "   " << angle << " : Edge " << e->id_edge << " between vertices " << e->v2->id_vertex << " " << e->v2->pt << " and " << e->v1->id_vertex << " " << e->v1->pt << std::endl;
		}
	}
	stream << "Events : " << std::endl;
	for (list<IndexedEvent *>::iterator it = v.events.begin() ; it != v.events.end() ; it++) {
		if ((*it) != NULL) {
			stream << "c_ray = " << (*it)->is_colliding_ray << ", c_col = " << (*it)->is_colliding_colinear << ", i = " << (*it)->intersectant << ", j = " << (*it)->intersected <<
				", t_i = " << (*it)->t_intersectant << ", t_j = " << (*it)->t_intersected << std::endl;
		} else {
			stream << "Null event" << std::endl;
		}
	}
	return stream;
}


inline std::ostream & operator << (std::ostream & stream, Face & f)
{
	stream << std::endl << "** Face " << f.id_face << " : " << std::endl;
	stream << "Half-edges : " << std::endl;
	for (list<HalfEdge *>::iterator it_h = f.edges.begin(); it_h != f.edges.end(); it_h++) {
		if ((*it_h)->v1_v2) {
			stream << "Vertex " << (*it_h)->e->v1->id_vertex << " " << (*it_h)->e->v1->pt << " "
				<< "to " << (*it_h)->e->v2->id_vertex << " " << (*it_h)->e->v2->pt << std::endl;
		} else {
			stream << "Vertex " << (*it_h)->e->v2->id_vertex << " " << (*it_h)->e->v2->pt << " "
				<< "to " << (*it_h)->e->v1->id_vertex << " " << (*it_h)->e->v1->pt << std::endl;
		}
	}
	stream << std::endl;
	stream << "Vertices, in clock-wise order : " << std::endl;
	for (list<pair<Vertex *, Vertex_Type> >::iterator it_v = f.vertices.begin(); it_v != f.vertices.end(); it_v++) {
		Vertex* v = it_v->first;
		stream << v->id_vertex << " " << v->pt << std::endl;
	}
	stream << std::endl;
	stream << "Vertices, details : " << std::endl;
	for (list<pair<Vertex *, Vertex_Type> >::iterator it_v = f.vertices.begin(); it_v != f.vertices.end(); it_v++) {
		stream << *(it_v->first) << std::endl;
	}
	return stream;
}
