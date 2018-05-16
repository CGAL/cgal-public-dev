#include "partition_elements.h"
#include <iostream>
#include <map>
#include <opencv2/core/core.hpp>

using cv::Mat;


Vertex::Vertex(int & _id_vertex, double x, double y)
{
	id_vertex = _id_vertex++;

	events = list<IndexedEvent *>();

	pt = Point2d(x, y);
	directions = vector<pair<double, HalfEdge *> >();

	v_prev = v_next = nullptr;
}


Vertex::Vertex(int & _id_vertex, IndexedEvent* event, double x, double y)
{
	id_vertex = _id_vertex++;

	events = list<IndexedEvent *>();
	if (event != nullptr) events.push_back(event);

	pt = Point2d(x, y);
	directions = vector<pair<double, HalfEdge *> >();

	v_prev = v_next = nullptr;
}


Vertex::Vertex(int & _id_vertex, IndexedEvent* event, Point2d & _pt)
{
	id_vertex = _id_vertex++;

	events = list<IndexedEvent *>();
	if (event != nullptr) events.push_back(event);
	
	pt = _pt;
	directions = vector<pair<double, HalfEdge*> >();

	v_prev = v_next = nullptr;
}


/*
Vertex::Vertex(const Vertex & V)
{
	// The copy constructor is only called when a partition is copied
	// And when a partition is copied, we first create the vertices, then the edges, then the facets
	id_vertex = V.id_vertex;

	// Deep copy of events
	events = list<IndexedEvent *>();
	for (list<IndexedEvent *>::const_iterator it_e = V.events.begin() ; it_e != V.events.end() ; it_e++) {
		events.push_back(new IndexedEvent(*(*it_e)));
	}

	pt = V.pt;

	// We do not copy the vector of directions, since we are going to call the copy constructor for directions later
	directions = vector<pair<double, HalfEdge *> >();

	v_prev = v_next = nullptr;
}
*/

Vertex::~Vertex()
{
	for (list<IndexedEvent*>::iterator it_e = events.begin() ; it_e != events.end() ; it_e++) {
		delete (*it_e);
	}
	events.clear();

	assert(directions.size() == 0);
	assert(faces.empty());
	faces.clear();
	directions.clear();
}


double Vertex::incidence_angle(const Edge* e)
{
	for (vector<pair<double, HalfEdge *> >::const_iterator it = directions.begin(); it != directions.end() ; ++it) {
		if (it->second->e == e) {
			return it->first;
		}
	}
	assert(false);
	return 0;
}


void Vertex::add(double alpha, HalfEdge *h)
{
	vector<pair<double, HalfEdge *> >::iterator it;
	for (it = directions.begin(); it != directions.end(); it++) {
		if (it->first > alpha) break;
	}
	directions.insert(it, std::make_pair(alpha, h));
}


void Vertex::remove(HalfEdge *h)
{
	vector<pair<double, HalfEdge *> >::iterator it;
	for (it = directions.begin() ; it != directions.end() ; it++) {
		if (it->second == h) break;
	}
	directions.erase(it);
}


int Vertex::connectivity()
{
	return int(directions.size());
}


void Vertex::add(Face *f)
{
	faces.insert(f);
}


void Vertex::remove(Face *f)
{
	set<Face *>::iterator it;
	for (it = faces.begin(); it != faces.end(); it++) {
		if (*it == f) break;
	}
	faces.erase(it);
}


bool Vertex::outer_vertex_is_very_close(Point2d & _pt)
{
	return (pt.x == _pt.x && pt.y == _pt.y);
}


bool Vertex::inner_vertex_created_by_colinear_ray(Ray* r_intersectant, Ray* r_intersected, vector<Ray *> & rays, bool verbose)
{
	// We are going to check if this vertex has been created by another ray that belongs
	// to the group of colinear segments identified by 'node'
	Node_Colinear_Segments* node = r_intersectant->parent->node_colinear;

	for (list<IndexedEvent *>::iterator it_e = events.begin(); it_e != events.end(); it_e++) {
		IndexedEvent* e = *it_e;
		int intersectant = e->intersectant, intersected = e->intersected;
		if (intersected >= 0 && rays[intersected] == r_intersected && rays[intersectant]->parent->node_colinear == node) {
			return true;
		}
	}

	return false;
}


bool Vertex::inner_vertex_barycenter_of_brother_segment(Ray* r_intersectant, Ray* r_intersected, vector<Ray *> & rays)
{
	int s_intersected = r_intersected->parent->index;
	int s_intersectant = r_intersectant->parent->index;

	for (list<IndexedEvent *>::iterator it_e = events.begin() ; it_e != events.end() ; it_e++) {
		IndexedEvent* e = (*it_e);
		int intersectant = e->intersectant, intersected = e->intersected;
		if (intersected == intersectant + 1) {
			if (rays[intersectant]->parent->index == s_intersectant || rays[intersectant]->parent->index == s_intersected) {
				return true;
			}
		}
	}

	return false;
}


bool Vertex::same_intersection_point(Ray* r_intersectant, Ray* r_intersected)
{
	int id_intersectant = r_intersectant->index, id_intersected = r_intersected->index;

	for (list<IndexedEvent *>::iterator it_e = events.begin() ; it_e != events.end() ; it_e++) {
		IndexedEvent* e = (*it_e);
		int intersectant = e->intersectant, intersected = e->intersected;
		if ((intersectant == id_intersectant && intersected == id_intersected)
			|| (intersectant == id_intersected && intersected == id_intersectant))
			return true;
	}

	return false;
}


bool Vertex::has_no_duplicated_event()
{
	if (events.size() <= 1) return true;

	list<IndexedEvent*>::iterator it_e1, it_e2;
	for (it_e1 = events.begin() ; it_e1 != events.end() ; it_e1++) {
		int e1_intersectant = (*it_e1)->intersectant;
		int e1_intersected = (*it_e1)->intersected;
		it_e2 = it_e1;
		it_e2++;
		while (it_e2 != events.end()) {
			int e2_intersectant = (*it_e2)->intersectant;
			int e2_intersected = (*it_e2)->intersected;
			if (e1_intersectant == e2_intersectant && e1_intersected == e2_intersected) return false;
			it_e2++;
		}
	}

	return true;
}


HalfEdge* Vertex::most_orthogonal_halfedge(Vec2d & direction)
{
	double min_abs_dot_product = FLT_MAX;
	HalfEdge* argmin_abs_dot_product = nullptr;

	for (uint d = 0 ; d < directions.size() ; d++) {
		double alpha_d = directions[d].first;
		Vec2d u = Vec2d(cos(alpha_d), sin(alpha_d));
		double abs_dot_product = fabs(direction.ddot(u));
		if (abs_dot_product < min_abs_dot_product) {
			min_abs_dot_product = abs_dot_product;
			argmin_abs_dot_product = directions[d].second;
		}
	}

	return argmin_abs_dot_product;
}


void Vertex::most_orthogonal_halfedge(Vertex* v_prev, Vertex* v_next, Vec2d & v_prev_next, double & abs_angle, HalfEdge* & h)
{
	double min_abs_dot_product = FLT_MAX;
	HalfEdge* argmin_abs_dot_product = nullptr;

	for (uint d = 0; d < directions.size(); d++) {
		HalfEdge* h_d = directions[d].second;
		if (h_d->v1_v2 && (h_d->e->v2 == v_prev || h_d->e->v2 == v_next)) {
			continue;
		} else if (!h_d->v1_v2 && (h_d->e->v1 == v_prev || h_d->e->v1 == v_next)) {
			continue;
		} else {
			Vec2d u = Vec2d(cos(directions[d].first), sin(directions[d].first));
			double abs_dot_product = fabs(v_prev_next.ddot(u));
			if (abs_dot_product < min_abs_dot_product) {
				min_abs_dot_product = abs_dot_product;
				argmin_abs_dot_product = h_d;
			}
		}
	}

	abs_angle = acos(min_abs_dot_product);
	h = argmin_abs_dot_product;
}


bool Vertex::right_turn(HalfEdge* h_prev, Vertex* v_curr, HalfEdge* h_next)
{
	// Handle colinearity cases
	if (h_prev->e->type == INNER_EDGE && h_next->e->type == INNER_EDGE) {
		Inner_Edge* e_prev = static_cast<Inner_Edge *>(h_prev->e);
		Inner_Edge* e_next = static_cast<Inner_Edge *>(h_next->e);
		Segment* s_prev = e_prev->get_front_supporting_segment();
		Segment* s_next = e_next->get_front_supporting_segment();
		if (s_prev->node_colinear != NULL && s_prev->node_colinear == s_next->node_colinear) {
			return true;
		}
	} else if (h_prev->e->type == OUTER_EDGE && h_next->e->type == OUTER_EDGE) {
		Outer_Edge* e_prev = static_cast<Outer_Edge *>(h_prev->e);
		Outer_Edge* e_next = static_cast<Outer_Edge *>(h_next->e);
		Image_Boundary b_prev = e_prev->boundary;
		Image_Boundary b_next = e_next->boundary;
		if (b_prev == b_next) {
			return true;
		}
	}

	double alpha_p = 0, alpha_n = 0;
	for (uint d = 0 ; d < uint(v_curr->directions.size()) ; ++d) {
		double & alpha_d = v_curr->directions[d].first;
		HalfEdge* h_d = v_curr->directions[d].second;
		if (h_d == h_prev) {
			alpha_p = alpha_d;
		} else if (h_d == h_next) {
			alpha_n = alpha_d;
		}
	}

	double eps = 1e-6;
	if (alpha_p <= 0) {
		// We should have alpha_n in [alpha_p, alpha_p + PI]
		return (fabs(alpha_p - alpha_n) < eps) || (fabs(alpha_p + PI - alpha_n) < eps) || (alpha_p < alpha_n && alpha_n < alpha_p + PI);
	} else {
		// We should have alpha_n in [-PI, alpha_p - PI] U [alpha_p, PI]
		bool first_interval = (fabs(alpha_n - (alpha_p - PI)) < eps) || (alpha_n < alpha_p - PI);
		if (first_interval) {
			return true;
		} else {
			bool second_interval = (fabs(alpha_n - alpha_p) < eps) || (alpha_p < alpha_n);
			return second_interval;
		}
	}
}


bool Vertex::right_turn(Vertex* v_prev, Vertex* v_curr, Vertex* v_next)
{
	HalfEdge *h_prev = nullptr, *h_next = nullptr;
	for (uint d = 0; d < v_curr->directions.size(); d++) {
		HalfEdge* h_d = v_curr->directions[d].second;
		if ((h_d->v1_v2 && h_d->e->v2 == v_prev) || (!h_d->v1_v2 && h_d->e->v1 == v_prev)) {
			h_prev = h_d;
		} else if ((h_d->v1_v2 && h_d->e->v2 == v_next) || (!h_d->v1_v2 && h_d->e->v1 == v_next)) {
			h_next = h_d;
		}
	}

	assert(h_prev != nullptr && h_next != nullptr);
	return right_turn(h_prev, v_curr, h_next);
}


bool Vertex::approximate_coordinates(Ray* r_i, Ray* r_j, double t_i, int rows, int cols, Point2d & pt)
{
	// The return value is true iff pt is a corner of the image
	double a = r_i->parent->a, b = r_i->parent->b, c = r_i->parent->c;
	double x = 0, y = 0;
	if (fabs(b) > fabs(a)) {
		x = r_i->A.x + t_i * r_i->OA[0];
		y = (-c - a * x) / b;
	} else {
		y = r_i->A.y + t_i * r_i->OA[1];
		x = (-c - b * y) / a;
	}
	pt = Point2d(x, y);
	
	if (r_j != nullptr) {
		return false;
	} else {
		// Case when a boundary gets intersected
		return Vertex::round_coordinates(pt, rows, cols, 0.001);
	}
}


bool Vertex::round_coordinates(Point2d & pt, int rows, int cols, double eps)
{
	// The return value is true iff pt is a corner of the image
	if (fabs(pt.x) < eps) {
		if (fabs(pt.y) < eps) {
			// Bottom left corner
			pt.x = 0;
			pt.y = 0;
			return true;
		} else if (fabs(double(rows) - pt.y) < eps) {
			// Top left corner
			pt.x = 0;
			pt.y = double(rows);
			return true;
		} else {
			// Left border ?
			if (fabs(pt.x) < eps) pt.x = 0;
			return false;
		}
	} else if (fabs(double(cols) - pt.x) < eps) {
		if (fabs(pt.y) < eps) {
			// Bottom right corner
			pt.x = double(cols);
			pt.y = 0;
			return true;
		} else if (fabs(double(rows) - pt.y) < eps) {
			// Top right corner
			pt.x = double(cols);
			pt.y = double(rows);
			return true;
		} else {
			// Right border ?
			if (fabs(double(cols) - pt.x) < eps) pt.x = double(cols);
			return false;
		}
	} else {
		// Bottom and top borders
		if (fabs(pt.y) < eps) {
			pt.y = 0;
		} else if (fabs(double(rows) - pt.y) < eps) {
			pt.y = double(rows);
		}
		return false;
	}
}


Edge::Edge(int & _id_edge, Vertex* _v1, Vertex* _v2, double & alpha_1, double & alpha_2)
{
	id_edge = _id_edge++;

	init_vertices_and_halfedges(_v1, alpha_1, _v2, alpha_2);

	e_prev = e_next = nullptr;
}


/*
Edge::Edge(const Edge & E)
{
	id_edge = E.id_edge;
	
	double alpha_1 = E.v1->incidence_angle(&E);
	double alpha_2 = E.v2->incidence_angle(&E);
	init_vertices_and_halfedges(E.v1, alpha_1, E.v2, alpha_2);

	e_prev = e_next = nullptr;
}
*/


Artificial_Edge::Artificial_Edge(int & _id_edge, double _a, double _b, double _c, Vertex* _v1, Vertex* _v2, double _alpha_1, double _alpha_2)
	: Edge (_id_edge, _v1, _v2, _alpha_1, _alpha_2)
{
	type = ARTIFICIAL_EDGE;
	a = _a;
	b = _b;
	c = _c;
}


/*
Artificial_Edge::Artificial_Edge(const Artificial_Edge & E)
	: Edge (E)
{
	type = ARTIFICIAL_EDGE;
	a = E.a;
	b = E.b;
	c = E.c;
}
*/


Outer_Edge::Outer_Edge(int & _id_edge, Image_Boundary _boundary, Vertex* _v1, Vertex* _v2, double alpha_1, double alpha_2, double _t1, double _t2)
	: Edge (_id_edge, _v1, _v2, alpha_1, alpha_2)
{
	type = OUTER_EDGE;
	boundary = _boundary;
}


/*
Outer_Edge::Outer_Edge(const Outer_Edge & E)
	: Edge(E)
{
	type = OUTER_EDGE;
	boundary = E.boundary;
}
*/


Inner_Edge::Inner_Edge(int & _id_edge, Ray* _r, int _tag, Vertex* _v1, Vertex* _v2, double alpha_1, double alpha_2, double _t1, double _t2)
	: Edge (_id_edge, _v1, _v2, alpha_1, alpha_2)
{
	type = INNER_EDGE;
	tag = _tag;
	rays.insert(_r);
}


Inner_Edge::Inner_Edge(int & _id_edge, set<Ray*> & _support, int _tag, Vertex* _v1, Vertex* _v2, double alpha_1, double alpha_2)
	: Edge(_id_edge, _v1, _v2, alpha_1, alpha_2)
{
	type = INNER_EDGE;
	tag = _tag;
	for (set<Ray *>::iterator it = _support.begin() ; it != _support.end() ; ++it) {
		rays.insert(*it);
	}
}


/*
Inner_Edge::Inner_Edge(const Inner_Edge & E)
	: Edge(E)
{
	type = INNER_EDGE;
	for (set<Ray *>::iterator it = E.rays.begin() ; it != E.rays.end() ; ++it) {
		rays.insert(*it);
	}
}
*/


Edge::~Edge()
{
	v1->remove(v1_v2);
	v2->remove(v2_v1);
	delete v1_v2;
	delete v2_v1;
}


Artificial_Edge::~Artificial_Edge()
{

}

Inner_Edge::~Inner_Edge()
{

}


Outer_Edge::~Outer_Edge()
{

}


void Edge::init_vertices_and_halfedges(Vertex* _v1, double & _alpha_1, Vertex* _v2, double & _alpha_2)
{
	v1 = _v1;
	v2 = _v2;

	v1_v2 = new HalfEdge(this, true);
	v2_v1 = new HalfEdge(this, false);
	
	alpha_1 = _alpha_1;
	alpha_2 = _alpha_2;

	v1->add(alpha_1, v1_v2);
	v2->add(alpha_2, v2_v1);
}


void Inner_Edge::add_supporting_ray(Ray* s)
{
	if (type == INNER_EDGE) {
		rays.insert(s);
	}
}


void Inner_Edge::add_supporting_rays(set<Ray *> & s)
{
	if (type == INNER_EDGE) {
		for (set<Ray *>::iterator it_s = s.begin(); it_s != s.end(); it_s++) {
			rays.insert(*it_s);
		}
	}
}


void Inner_Edge::get_supporting_rays(set<Ray *> & s)
{
	if (type == INNER_EDGE) {
		for (set<Ray *>::iterator it_s = rays.begin() ; it_s != rays.end() ; it_s++) {
			s.insert(*it_s);
		}
	}
}


void Inner_Edge::get_supporting_segments(set<Segment *> & s)
{
	if (type == INNER_EDGE) {
		for (set<Ray *>::iterator it_s = rays.begin() ; it_s != rays.end() ; it_s++) {
			s.insert((*it_s)->parent);
		}
	}
}


Ray* Inner_Edge::get_front_supporting_ray()
{
	if (type == INNER_EDGE) {
		Ray* r = *rays.begin();
		return r;
	} else {
		return nullptr;
	}
}


Segment* Inner_Edge::get_front_supporting_segment()
{
	if (type == INNER_EDGE) {
		Ray* r = *rays.begin();
		return r->parent;
	} else {
		return nullptr;
	}
}


void Inner_Edge::time_range(Ray* r, double & t_1, double & t_2)
{
	if (v1->events.empty()) {
		t_1 = -r->initial_length;
	} else {
		for (list<IndexedEvent *>::iterator it_e = v1->events.begin() ; it_e != v1->events.end() ; it_e++) {
			IndexedEvent* e = (*it_e);
			if (e->intersectant == r->index) {
				t_1 = e->t_intersectant;
				break;
			} else if (e->intersected == r->index) {
				t_1 = e->t_intersected;
				break;
			}
		}
	}

	if (v2->events.empty()) {
		t_2 = -r->initial_length;
	} else {
		for (list<IndexedEvent *>::iterator it_e = v2->events.begin() ; it_e != v2->events.end() ; it_e++) {
			IndexedEvent* e = (*it_e);
			if (e->intersectant == r->index) {
				t_2 = e->t_intersectant;
				break;
			} else if (e->intersected == r->index) {
				t_2 = e->t_intersected;
				break;
			}
		}
	}
}


double Edge::length()
{
	return sqrt((v2->pt.x - v1->pt.x) * (v2->pt.x - v1->pt.x) + (v2->pt.y - v1->pt.y) * (v2->pt.y - v1->pt.y));
}


void Inner_Edge::line_equation(double &a, double &b, double &c)
{
	if (type == INNER_EDGE) {
		Segment* s_ref = (*rays.begin())->parent;
		a = s_ref->a, b = s_ref->b, c = s_ref->c;
	}
}


void Outer_Edge::line_equation(double &a, double &b, double &c)
{
	if (type == OUTER_EDGE) {
		if (boundary == TOP_IMAGE || boundary == BOTTOM_IMAGE) {
			a = 0, b = 1, c = -v1->pt.y;
		} else if (boundary == LEFT_IMAGE || boundary == RIGHT_IMAGE) {
			a = 1, b = 0, c = -v1->pt.x;
		}
	}
}


void Outer_Edge::time_range(double & t_1, double & t_2)
{
	if (v1->events.empty()) {
		switch (boundary) {
		case TOP_IMAGE: t_1 = v1->pt.x; break;
		case BOTTOM_IMAGE: t_1 = v1->pt.x; break;
		case LEFT_IMAGE: t_1 = v1->pt.y; break;
		case RIGHT_IMAGE: t_1 = v1->pt.y; break;
		}
	} else {
		t_1 = v1->events.front()->t_intersected;
	}

	if (v2->events.empty()) {
		switch (boundary) {
		case TOP_IMAGE: t_2 = v2->pt.x; break;
		case BOTTOM_IMAGE: t_2 = v2->pt.x; break;
		case LEFT_IMAGE: t_2 = v2->pt.y; break;
		case RIGHT_IMAGE: t_2 = v2->pt.y; break;
		}
	} else {
		t_2 = v2->events.front()->t_intersected;
	}
}


void Artificial_Edge::line_equation(double & _a, double & _b, double & _c)
{
	_a = a; b = _b; c = _c;
}


HalfEdge::HalfEdge(Edge* _e, bool _v1_v2)
{
	e = _e;
	v1_v2 = _v1_v2;
	f = nullptr;
}


HalfEdge::~HalfEdge()
{

}


void HalfEdge::set(Face *_f)
{
	f = _f;
}


HalfEdge* HalfEdge::opposite()
{
	return (e->v1_v2 == this ? e->v2_v1 : e->v1_v2);
}



void HalfEdge::intersects_if_extended(HalfEdge* h, list<HalfEdge *> & intersectable_halfedges, Point2d & intersection_point, HalfEdge* & intersected)
{
	std::set<Segment *> support_extended;
	Inner_Edge* e = static_cast<Inner_Edge*>(h->e);
	e->get_supporting_segments(support_extended);
	Segment* s_extended = (*support_extended.begin());

	// Finds the point from which we extend the ray, finds its direction as well
	Point2d h_v = s_extended->finalBarycenter;
	Vec2d h_dir = Vec2d();

	Vec2d u = (h->v1_v2 ? h->e->v2->pt - h->e->v1->pt : h->e->v1->pt - h->e->v2->pt);
	Vec2d s_extended_dir = Vec2d(s_extended->finalEnd2 - s_extended->finalEnd1);
	if (u.ddot(s_extended->finalDirection) > 0) {
		h_dir = s_extended->finalDirection;
	} else {
		h_dir = -s_extended->finalDirection;
	}

	double eps = 1e-4;
	double t_min = FLT_MAX;
	for (list<HalfEdge *>::iterator it_h = intersectable_halfedges.begin() ; it_h != intersectable_halfedges.end(); ++it_h) {
	
		// Gets an edge
		HalfEdge* h_arg = (*it_h);
		Edge* arg = h_arg->e;

		// Gets its equation
		double a, b, c;
		arg->line_equation(a, b, c);

		// Extends rays and sees when its intersects the supporting line of the edge
		// Given that this function will only be used in the merging process...
		double den = (a * h_dir[0] + b * h_dir[1]);
		if (fabs(den) < 1e-6) {
			continue;
		} else {
			double t = -(a * h_v.x + b * h_v.y + c) / den;
			if (t < t_min) {
				intersection_point = Point2d(h_v.x + t * h_dir[0], h_v.y + t * h_dir[1]);
				
				Vertex *arg_v1 = arg->v1, *arg_v2 = arg->v2;
				double dx = arg_v2->pt.x - arg_v1->pt.x;
				double dy = arg_v2->pt.y - arg_v1->pt.y;
				bool check_x = true, check_y = true;
				
				if (fabs(dx) > 1e-6) {
					double u_x = (intersection_point.x - arg_v1->pt.x) / dx;
					check_x = (-eps <= u_x && u_x <= 1 + eps);
				}

				if (fabs(dy) > 1e-6) {
					double u_y = (intersection_point.y - arg_v1->pt.y) / dy;
					check_y = (-eps <= u_y && u_y <= 1 + eps);
				}

				if (check_x && check_y) {
					t_min = t;
					intersected = h_arg;
					return;
				}
			}
		}
	}
}



Face::Face(int & _id_face, list<HalfEdge *> & _edges)
{
	id_face = _id_face++;

	// Copies the list of edges
	edges = list<HalfEdge *>(_edges.begin(), _edges.end());

	// Reciprocally, for every edge we insert a reference to this face
	for (list<HalfEdge *>::iterator it = edges.begin(); it != edges.end() ; it++) {
		(*it)->set(this);
	}

	// We are now going to identify the list of vertices with the help of the half-edges
	list_vertices();

	classify();
	set_thinness();

    is_simple(true);
}



Face::~Face()
{
	// Removes references to this facet in every container
	for (list<HalfEdge *>::iterator it_h = edges.begin(); it_h != edges.end(); it_h++) {
		(*it_h)->set();
	}
	edges.clear();

	for (list<pair<Vertex *, Vertex_Type> >::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
		(*it_v).first->remove(this);
	}
	vertices.clear();
	
	pixels.clear();
}



bool Face::is_simple(bool verbose)
{
	assert(vertices.size() == edges.size());

	std::map<int, int> vertices_occurences;
	for (list<pair<Vertex *, Vertex_Type> >::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
		vertices_occurences.insert(std::make_pair((*it_v).first->id_vertex, 0));
	}

	for (list<HalfEdge *>::iterator it_h = edges.begin(); it_h != edges.end(); it_h++) {
		int v1 = (*it_h)->e->v1->id_vertex;
		int v2 = (*it_h)->e->v2->id_vertex;
		vertices_occurences[v1]++;
		vertices_occurences[v2]++;
	}

	for (std::map<int, int>::iterator it_m = vertices_occurences.begin(); it_m != vertices_occurences.end(); it_m++) {
		if (it_m->second != 2) {
			vector<double> x;
			vector<double> y;
			vector<int> l;
			
			for (list<pair<Vertex *, Vertex_Type> >::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
				Vertex* v = it_v->first;
				x.push_back(v->pt.x);
				y.push_back(v->pt.y);
				l.push_back(v->id_vertex);
			}
			FILE* file = fopen("error.R", "w");
			fprintf(file, "x <- c("); for (int i = 0; i < x.size() - 1; i++) fprintf(file, "%lf, ", x[i]); fprintf(file, "%lf)\n", x[x.size() - 1]);
			fprintf(file, "y <- c("); for (int i = 0; i < y.size() - 1; i++) fprintf(file, "%lf, ", y[i]); fprintf(file, "%lf)\n", y[y.size() - 1]);
			fprintf(file, "l <- c("); for (int i = 0; i < l.size() - 1; i++) fprintf(file, "%i, ", l[i]); fprintf(file, "%i)\n", l[l.size() - 1]);
			fclose(file);
            if (verbose) {
                std::cout << "****** ERROR : Face " << id_face << " is not simple" << std::endl;
                std::cout << *this << std::endl;
            }
            return false;
		}
	}

    return true;
}



void Face::list_vertices()
{
	list<HalfEdge *>::iterator it_h = edges.begin();
	HalfEdge* h_curr = (*it_h);
	HalfEdge* h_next = (*(++it_h));
	Vertex* v_init = nullptr;
	Vertex* v = nullptr;
	if (h_curr->e->v1 == h_next->e->v1 || h_curr->e->v1 == h_next->e->v2) {
		v_init = h_curr->e->v2;
		v = h_curr->e->v1;
	} else {
		v_init = h_curr->e->v1;
		v = h_curr->e->v2;
	}
	v_init->add(this);
	v->add(this);
	vertices.push_back(std::make_pair(v_init, Vertex_Type::UNDETERMINED));
	vertices.push_back(std::make_pair(v, Vertex_Type::UNDETERMINED));
	while (true) {
		h_curr = h_next;
		it_h++;
		if (it_h != edges.end()) {
			h_next = (*it_h);
			if (h_curr->e->v1 == h_next->e->v1 || h_curr->e->v1 == h_next->e->v2) {
				v = h_curr->e->v1;
			} else {
				v = h_curr->e->v2;
			}
			v->add(this);
			vertices.push_back(std::make_pair(v, Vertex_Type::UNDETERMINED));
		} else {
			break;
		}
	}
}



void Face::find_pixels_inside_facet()
{
	double x_min, x_max, y_min, y_max;
	list<HalfEdge *>::iterator it_l, it_r;

	find_leftmost_and_rightmost_edges(x_min, x_max, y_min, y_max, it_l, it_r);
	find_pixels(x_min, x_max, y_min, y_max, it_l, it_r);	
}



void Face::find_leftmost_and_rightmost_edges(double & x_min, double & x_max, double & y_min, double & y_max, list<HalfEdge *>::iterator & it_l, list<HalfEdge *>::iterator & it_r)
{
	// We are searching for the leftmost and rightmost half-edges of the face
	list<HalfEdge *>::iterator it_h = edges.begin();
	
	// Invariant : v1 and v2 represent the vertices of the previously crossed half-edge
	Vertex* v1 = nullptr;
	Vertex* v2 = nullptr;
	if ((*it_h)->v1_v2) {
		v1 = (*it_h)->e->v1;
		v2 = (*it_h)->e->v2;
	} else {
		v2 = (*it_h)->e->v1;
		v1 = (*it_h)->e->v2;
	}
	
	x_min = (v1->pt.x < v2->pt.x ? v1->pt.x : v2->pt.x);
	x_max = (v1->pt.x > v2->pt.x ? v1->pt.x : v2->pt.x);
	y_min = (v1->pt.y < v2->pt.y ? v1->pt.y : v2->pt.y);
	y_max = (v1->pt.y > v2->pt.y ? v1->pt.y : v2->pt.y);
	it_l = it_h, it_r = it_h;

	bool shift_l, shift_r;
	if (x_min == v2->pt.x) {
		shift_l = true;
		shift_r = false;
	} else if (x_max == v2->pt.x) {
		shift_l = false;
		shift_r = true;
	}
	it_h++;
	Vertex* v = nullptr;
	while (it_h != edges.end()) {
		HalfEdge* h = *it_h;
		if (h->e->v1 == v1 || h->e->v1 == v2) {
			// This vertex is common to the current half-edge and the previous one
			// This means that h->e->v2 is the targeted vertex
			// We are going to save it in the variable that doesn't contain h->e->v1
			if (v1 != h->e->v1) {
				v1 = h->e->v2;
				v = v1;
			} else {
				v2 = h->e->v2;
				v = v2;
			}
		} else {
			// h->e->v1 is the next new vertex, which means that h->e->v2 is the one
			// that is shared by the current and the previous half-edges. So we are
			// going to save h->e->v1 in the variable that doesn't contain h->e->v2
			if (h->e->v2 != v1) {
				v1 = h->e->v1;
				v = v1;
			} else {
				v2 = h->e->v1;
				v = v2;
			}
		}

		// If we have found the leftmost or the rightmost vertex so far,
		// we remember its coordinates and an iterator on the current half-edge
		if (v->pt.x < x_min) {
			x_min = v->pt.x;
			it_l = it_h;
			shift_l = true;
		}
		if (v->pt.x > x_max) {
			x_max = v->pt.x;
			it_r = it_h;
			shift_r = true;
		}
		y_min = MIN(y_min, v->pt.y);
		y_max = MAX(y_max, v->pt.y);
		it_h++;
	}

	if (shift_l) {
		if (++it_l == edges.end()) it_l = edges.begin();
	}
	if (shift_r) {
		if (++it_r == edges.end()) it_r = edges.begin();
	}
}



void Face::find_pixels(double x_min, double x_max, double y_min, double y_max, list<HalfEdge *>::iterator it_l, list<HalfEdge *>::iterator it_r, bool verbose)
{
	// Defines a grid covering the integers included in K = [x_min x_max[ x [y_min y_max[
	int x_0 = int(ceil(x_min));
	int y_0 = (y_min == ceil(y_min) ? int(y_min) + 1 : int(ceil(y_min)));
	int x_1 = (x_max == floor(x_max) ? int(x_max) - 1 : int(floor(x_max)));
	int y_1 = int(floor(y_max));
	int g_x = x_1 - x_0 + 1;
	int g_y = y_1 - y_0 + 1;
	Mat grid(g_y, g_x, CV_8U, cv::Scalar::all(1));

	// First half-loop : from left to right
	list<HalfEdge *>::iterator it_h = it_l;
	while (it_h != it_r) {
		Vertex* v1 = (*it_h)->e->v1;
		Vertex* v2 = (*it_h)->e->v2;
		// Finds the equation of the line (v1 v2)
		if (fabs(v2->pt.x - v1->pt.x) > 1e-6) {
			double a = (v2->pt.y - v1->pt.y) / (v2->pt.x - v1->pt.x);
			double b = v1->pt.y - a * v1->pt.x;
			// Discards every point above the line
			int xl_min = 0, xl_max = 0;
			if (v1->pt.x < v2->pt.x) {
				xl_min = int(ceil(v1->pt.x));
				xl_max = (v2->pt.x == floor(v2->pt.x) ? int(v2->pt.x) - 1 : int(floor(v2->pt.x)));
			} else {
				xl_min = int(ceil(v2->pt.x));
				xl_max = (v1->pt.x == floor(v1->pt.x) ? int(v1->pt.x) - 1 : int(floor(v1->pt.x)));
			}
			for (int x = xl_min ; x <= xl_max; x++) {
				int yl_min = int(a * x + b == ceil(a * x + b) ? int(a * x + b) + 1 : int(ceil(a * x + b)));
				int yl_max = y_1;
				for (int y = yl_min; y <= yl_max; y++) {
					grid.at<uchar>(y_1 - y, x - x_0) = 1 - grid.at<uchar>(y_1 - y, x - x_0);
				}
			}
		}
		// Iterates
		if (++it_h == edges.end()) it_h = edges.begin();
	}

	// Second half-loop : from right to left
	// Same algorithm as before, except that this time, we discard pixels below (v1 v2)
	it_h = it_r;
	while (it_h != it_l) {
		Vertex* v1 = (*it_h)->e->v1;
		Vertex* v2 = (*it_h)->e->v2;
		// Finds the equation of the line (v1 v2)
		if (fabs(v2->pt.x - v1->pt.x) > 1e-6) {
			double a = (v2->pt.y - v1->pt.y) / (v2->pt.x - v1->pt.x);
			double b = v1->pt.y - a * v1->pt.x;
			// Discards every point below the line
			int xl_min = 0, xl_max = 0;
			if (v1->pt.x < v2->pt.x) {
				xl_min = int(ceil(v1->pt.x));
				xl_max = (v2->pt.x == floor(v2->pt.x) ? int(v2->pt.x) - 1 : int(floor(v2->pt.x)));
			} else {
				xl_min = int(ceil(v2->pt.x));
				xl_max = (v1->pt.x == floor(v1->pt.x) ? int(v1->pt.x) - 1 : int(floor(v1->pt.x)));
			}
			for (int x = xl_min; x <= xl_max; x++) {
				int yl_min = y_0;
				int yl_max = int(floor(a * x + b));
				for (int y = yl_min; y <= yl_max ; y++) {
					grid.at<uchar>(y_1 - y, x - x_0) = 1 - grid.at<uchar>(y_1 - y, x - x_0);
				}
			}
		}
		// Iterates
		if (++it_h == edges.end()) it_h = edges.begin();
	}

	// Finally : uses the grid to find the pixels included in the face
	for (int i = 0; i < g_y; i++) {
		for (int j = 0; j < g_x; j++) {
			if (grid.at<uchar>(i, j) == 1) {
				pixels.push_back(Point2i(x_0 + j, y_1 - i));
			}
		}
	}
}



void Face::print_details()
{
	vector<double> x;
	vector<double> y;
	vector<int> l;
	for (list<pair<Vertex *, Vertex_Type> >::iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
		Vertex* v = it_v->first;
		x.push_back(v->pt.x);
		y.push_back(v->pt.y);
		l.push_back(v->id_vertex);
	}
	std::string name = "face_" + std::to_string(id_face) + ".R";
	FILE* file = fopen(name.c_str(), "w");
	fprintf(file, "x <- c("); for (int i = 0; i < x.size() - 1; i++) fprintf(file, "%lf, ", x[i]); fprintf(file, "%lf)\n", x[x.size() - 1]);
	fprintf(file, "y <- c("); for (int i = 0; i < y.size() - 1; i++) fprintf(file, "%lf, ", y[i]); fprintf(file, "%lf)\n", y[y.size() - 1]);
	fprintf(file, "l <- c("); for (int i = 0; i < l.size() - 1; i++) fprintf(file, "%i, ", l[i]); fprintf(file, "%i)\n", l[l.size() - 1]);
	fclose(file);
}



void Face::get_supporting_segments(set<Segment *> & supporting_segments)
{
	supporting_segments.clear();
	for (list<HalfEdge *>::iterator it_h = edges.begin(); it_h != edges.end(); it_h++) {
		Edge* e = (*it_h)->e;
		if (e->type == INNER_EDGE) {
			Inner_Edge* i_e = static_cast<Inner_Edge*>(e);
			i_e->get_supporting_segments(supporting_segments);
		}
	}
}


void Face::get_neighboring_faces(set<Face *> & neighboring_faces)
{
	neighboring_faces.clear();
	for (list<HalfEdge *>::iterator it_h = edges.begin() ; it_h != edges.end() ; it_h++) {
		Face* neighbor = (*it_h)->opposite()->f;
		if (neighbor != nullptr) {
			neighboring_faces.insert(neighbor);
		}
	}
}


void Face::classify()
{
	vector<HalfEdge *> v_edges(edges.begin(), edges.end());
	vector<pair<Vertex *, Vertex_Type> > v_vertices(vertices.begin(), vertices.end());
	uint n = uint(v_vertices.size());

	for (uint i = 0 ; i < vertices.size() ; i++) {
		Vertex* v = v_vertices[i].first;
		int v_c = v->connectivity();

		// Accesses the previous and next edges
		Edge* e_1 = v_edges[(i + n - 1) % n]->e;
		Edge* e_2 = v_edges[i]->e;

		// Exhibits their respective support segments
		// The idea is that an inner edge is created by extension of one segment or more
		// which is not the case of an outer edge which relies on the image borders and has no support

		if (e_1->type == Edge_Type::INNER_EDGE && e_2->type == Edge_Type::INNER_EDGE) {
			
			set<Segment *> support_1, support_2;
			Inner_Edge* ie_1 = static_cast<Inner_Edge*>(e_1);
			Inner_Edge* ie_2 = static_cast<Inner_Edge*>(e_2);
			ie_1->get_supporting_segments(support_1);
			ie_2->get_supporting_segments(support_2);

			// v is a corner if and only if s_1 and s_2 contain colinear segments
			bool detected_colinear_segments = false;
			for (set<Segment *>::iterator it_s1 = support_1.begin() ; it_s1 != support_1.end() ; it_s1++) {
				Segment* s_1 = (*it_s1);
				for (set<Segment *>::iterator it_s2 = support_2.begin() ; it_s2 != support_2.end() ; it_s2++) {
					Segment* s_2 = (*it_s2);
					if (s_1 == s_2 || (s_1->node_colinear != nullptr && s_1->node_colinear == s_2->node_colinear)) {
						detected_colinear_segments = true;
						break;
					}
				}
				if (detected_colinear_segments) break;
			}
			
			if (detected_colinear_segments) {
				v_vertices[i].second = (v_c == 2 ? NON_CORNER_BIVALENT : NON_CORNER_TRIVALENT);
			} else {
				v_vertices[i].second = CORNER;
			}

		} else if (e_1->type == Edge_Type::OUTER_EDGE && e_2->type == Edge_Type::OUTER_EDGE) {
			
			// Case when edges e_1 and e_2 are borders
			Outer_Edge* oe_1 = static_cast<Outer_Edge*>(e_1);
			Outer_Edge* oe_2 = static_cast<Outer_Edge*>(e_2);
			
			if (oe_1->boundary == oe_2->boundary) {
				v_vertices[i].second = (v_c == 2 ? NON_CORNER_BIVALENT : NON_CORNER_TRIVALENT);
			} else {
				v_vertices[i].second = CORNER;
			}

		} else {
			// v is the intersection of an inner and an other segment
			// It must be a corner
			v_vertices[i].second = CORNER;
		}
	}

	vertices = list<pair<Vertex*, Vertex_Type> > (v_vertices.begin(), v_vertices.end());
}


void Face::set_thinness()
{
	// Finds the main points of the facet, i.e. vertices which define a corner of the facet and not a flat angle
	list<Vertex *> corners;
	for (list<pair<Vertex *, Vertex_Type> >::iterator it_v = vertices.begin() ; it_v != vertices.end() ; it_v++) {
		if (it_v->second == CORNER) {
			corners.push_back(it_v->first);
		}
	}

	// Finds the center of the corners
	double x = 0, y = 0;
	for (list<Vertex *>::iterator it_v = corners.begin() ; it_v != corners.end() ; it_v++) {
		x += (*it_v)->pt.x;
		y += (*it_v)->pt.y;
	}
	x /= corners.size(); y /= corners.size();

	// Finds the principal inertia axis of the corners
	double I_xx = 0, I_yy = 0, I_xy = 0;
	for (list<Vertex *>::iterator it_v = corners.begin() ; it_v != corners.end() ; it_v++) {
		Point2d & pt = (*it_v)->pt;
		I_xx += (pt.y - y) * (pt.y - y);
		I_yy += (pt.x - x) * (pt.x - x);
		I_xy -= (pt.x - x) * (pt.y - y);
	}

	// Compute smallest eigen value and gets the angle
	double lambda = 0.5 * (I_xx + I_yy - sqrt((I_xx - I_yy) * (I_xx - I_yy) + 4 * I_xy * I_xy));
	double theta = fabs(I_xx) > fabs(I_yy) ? atan2(lambda - I_xx, I_xy) : atan2(I_xy, lambda - I_yy);

	double dx = cos(theta), dy = sin(theta);
	double l_min = 0, l_max = 0, w_min = 0, w_max = 0;

	// We compute l and w as the distance from the center of the corners points to the corners themselves
	// projected along the axes (dx, dy) and (-dy, dx).
	for (list<Vertex *>::iterator it_v = corners.begin(); it_v != corners.end() ; it_v++) {
		Point2d & pt = (*it_v)->pt;
		double l =  (pt.x - x) * dx + (pt.y - y) * dy;
		double w = -(pt.x - x) * dy + (pt.y - y) * dx;
		if (l > l_max) l_max = l;
		if (l < l_min) l_min = l;
		if (w > w_max) w_max = w;
		if (w < w_min) w_min = w;
	}

	double dl = l_max - l_min;
	double dw = w_max - w_min;

	// The thinness of a cell can also be given by dl / dw
	// Indeed if the facet is very thin, then dw should be very small compared to dl, and vice-versa
	thinness = dw;
}



void Face::add_non_corner(HalfEdge* h, Vertex* v, HalfEdge* h_1, HalfEdge* h_2)
{
	list<pair<Vertex*, Vertex_Type> >::iterator it_v = vertices.begin();
	list<HalfEdge *>::iterator it_h = edges.begin();
	
	// Localizes the edge to replace
	while (it_h != edges.end() && (*it_h) != h) {
		++it_v;
		++it_h;
	}
	++it_v;

	// Replaces h = (a b) by h_1 = (a v) and h_2 = (v b)
	it_h = edges.erase(it_h);
	
	edges.insert(it_h, h_1);
	edges.insert(it_h, h_2);
	h_1->f = this;
	h_2->f = this;

	// Inserts v
	vertices.insert(it_v, std::make_pair(v, NON_CORNER_TRIVALENT));
	v->add(this);
}



void Face::remove_non_corners(Vertex* v1, Vertex* v2, HalfEdge* h)
{
	// We assume that v1 and v2 are included in the list of vertices of the facet.
	// We assume that all vertices between v1 and v2 are non corner, bivalent vertices that we are going to remove.
	// Finally, the sequence of edges and vertices between v1 and v2 is replaced with a single halfedge h.
	list<pair<Vertex*, Vertex_Type> >::iterator it_v = vertices.begin();
	list<HalfEdge *>::iterator it_h = edges.begin();
	bool insert_last_position = false;

	// First of all we have to find v1
	while (it_v->first != v1) {
		++it_v; ++it_h;
	}

	// it_v is now pointing to v1, and it_h to the halfedge (v1 w)
	it_h = edges.erase(it_h);
	++it_v;
	if (it_v == vertices.end() && it_h == edges.end()) {
		it_h = edges.begin();
		it_v = vertices.begin();
		insert_last_position = true;
	}

	// We start removing all the vertices pointing by it_v, until finding v2
	while (it_v->first != v2) {
		it_v->first->remove(this);
		it_v = vertices.erase(it_v);
		it_h = edges.erase(it_h);
		if (it_v == vertices.end() && it_h == edges.end()) {
			it_v = vertices.begin();
			it_h = edges.begin();
			insert_last_position = true;
		}
	}
	
	// In the end, we insert h
	if (insert_last_position) {
		edges.insert(edges.end(), h);
	} else {
		edges.insert(it_h, h);
	}
	h->f = this;

	{
		list<pair<Vertex*, Vertex_Type> >::iterator it_vt = vertices.begin();
		list<HalfEdge *>::iterator it_ht = edges.begin();
		while (it_vt != vertices.end() && it_ht != edges.end()) {
			Vertex* v = it_vt->first;
			HalfEdge* h = (*it_ht);
			if (h->v1_v2) {
				assert(h->e->v1 == v);
			} else {
				assert(h->e->v2 == v);
			}
			++it_vt;
			++it_ht;
		}
	}
}


void Face::set_as_bivalent(Vertex* v)
{
	bool found_vertex = false;
	typedef pair<Vertex *, Vertex_Type> Vertex_Element;
	for (list<Vertex_Element>::iterator it_v = vertices.begin() ; it_v != vertices.end() ; it_v++) {
		if (it_v->first == v) {
			it_v->second = NON_CORNER_BIVALENT;
			found_vertex = true;
			break;
		}
	}
	assert(found_vertex);
}



void Face::set_as_trivalent(Vertex* v)
{
	bool found_vertex = false;
	typedef pair<Vertex *, Vertex_Type> Vertex_Element;
	for (list<Vertex_Element>::iterator it_v = vertices.begin() ; it_v != vertices.end() ; it_v++) {
		if (it_v->first == v) {
			it_v->second = NON_CORNER_TRIVALENT;
			found_vertex = true;
			break;
		}
	}
	assert(found_vertex);
}


bool Face::convex_union(list<pair<uint, uint> > & adjacent_edges)
{
	typedef pair<Vertex *, Vertex_Type> Vertex_Element;
	vector<Vertex_Element> v_vertices (vertices.begin(), vertices.end());
	vector<HalfEdge *> v_edges(edges.begin(), edges.end());
	uint n = uint(vertices.size());
	adjacent_edges.clear();

	for (uint i = 0 ; i < n ; i++) {
	
		// The current iteration is skipped if we can not link two corners of the facet,
		// between which there exists no other vertex, or bivalent vertices only
		Vertex_Element curr_1 = v_vertices[i];
		if (curr_1.second != CORNER) continue;

		uint j = 1;
		while (v_vertices[(i + j) % n].second == NON_CORNER_BIVALENT) j++;
		Vertex_Element curr_2 = v_vertices[(i + j) % n];
		if (curr_2.second != CORNER) {
			i += j;
			continue;
		}

		// We obtain the halfedge that preceeds v_curr_1 and the corner that follows v_curr_2
		Vertex *v_1 = curr_1.first, *v_2 = curr_2.first;
		Vertex_Element prev = v_vertices[(i + n - 1) % n];
		Vertex_Element next = v_vertices[(i + j + 1) % n];

		// We obtain the facet that is located on the other side of the big edge that binds v_curr_1 to v_curr_2
		Face* f_adj = v_edges[i]->opposite()->f;
		if (f_adj == nullptr) {
			i += j - 1;
			continue;
		}

		vector<Vertex_Element> adj_vertices (f_adj->vertices.begin(), f_adj->vertices.end());
		uint m = uint(adj_vertices.size());

		// We need to find the vertex that preceeds v_curr_2 and the one that follows v_curr_1 in f_adj
		uint k = 0, k_1 = -1, k_2 = -1;
		while (k < adj_vertices.size()) {
			if (adj_vertices[k].first == curr_1.first) {
				k_1 = k;
			} else if (adj_vertices[k].first == curr_2.first) {
				k_2 = k;
			}
			if (k_1 != -1 && k_2 != -1) {
				break;
			} else {
				k++;
			}
		}
		Vertex_Element a_prev_curr_2 = adj_vertices[(k_2 + m - 1) % m];
		Vertex_Element a_next_curr_1 = adj_vertices[(k_1 + 1) % m];

		// We are now going to check if (prev curr_1 a_next) and (a_prev curr_2 next) make a right turn
		// If so it means that the union of the two facets will not affect the convexity of the resulting cell
		// As output, we keep a trace of the only edge that delimits the two cells
		Vertex *v_prev = prev.first, *v_next = next.first;
		Vertex *v_a_prev = a_prev_curr_2.first, *v_a_next = a_next_curr_1.first;
		if (Vertex::right_turn(v_prev, v_1, v_a_next) && Vertex::right_turn(v_a_prev, v_2, v_next)) {
			adjacent_edges.push_back(std::make_pair(i, j));
		}
		
		i += j - 1;
	}

	return (!adjacent_edges.empty());
}


bool Face::convex_cut(list<pair<uint, uint> > & edges_to_slice)
{
	typedef pair<Vertex*, Vertex_Type> Vertex_Element;
	vector<Vertex_Element> v_vertices (vertices.begin(), vertices.end());
	vector<HalfEdge *> v_edges(edges.begin(), edges.end());

	uint n = uint(vertices.size());
	edges_to_slice.clear();

	for (uint i = 0; i < n; i++) {

		// First, we need to find a sequence of corners
		Vertex_Element corner_1 = v_vertices[i];
		if (corner_1.second != CORNER) continue;

		uint j = 1;
		while (v_vertices[(i + j) % n].second != CORNER) j++;
		Vertex_Element corner_2 = v_vertices[(i + j) % n];
		
		// There must exist at least one non corner, trivalent vertex between the two corners
		if (j == 1) continue;

		bool trivalent_vertex_found = false;
		for (uint k = 1; k < j; k++) {
			if (v_vertices[(i + k) % n].second == NON_CORNER_TRIVALENT) {
				trivalent_vertex_found = true;
				break;
			}
		}
		if (!trivalent_vertex_found) {
			i += j - 1;
			continue;
		}

		// We are going to consider the half-edge h_1 = (corner_1, vertices[i + 1])
		// We access to the facet that is located on the other side of h1, and check if the convexity properties are respected
		HalfEdge* h_1 = v_edges[i];
		Face* f_adj_1 = h_1->opposite()->f;
		if (f_adj_1 == nullptr) {
			i += j - 1;
			continue;
		}

		list<Vertex_Element> & adj_vertices_1 = f_adj_1->vertices;
		vector<Vertex_Element> v_adj_vertices_1(adj_vertices_1.begin(), adj_vertices_1.end());

		uint m_1 = uint(v_adj_vertices_1.size()), k_1 = 0;
		while (k_1 != v_adj_vertices_1.size()) {
			if (v_adj_vertices_1[k_1].first == corner_1.first) break;
			k_1++;
		}
		Vertex_Element prev = v_vertices[(i + n - 1) % n];
		Vertex_Element a_next_1 = v_adj_vertices_1[(k_1 + 1) % m_1];
		if (!Vertex::right_turn(prev.first, corner_1.first, a_next_1.first)) {
			i += j - 1;
			continue;
		}

		// Same operation for the half-edge h_2 = (vertices[i + j - 1], corner_2)
		HalfEdge* h_2 = v_edges[(i + j - 1) % n];
		Face* f_adj_2 = h_2->opposite()->f;

		list<Vertex_Element> & adj_vertices_2 = f_adj_2->vertices;
		vector<Vertex_Element> v_adj_vertices_2(adj_vertices_2.begin(), adj_vertices_2.end());

		uint m_2 = uint(v_adj_vertices_2.size()), k_2 = 0;
		while (k_2 != v_adj_vertices_2.size()) {
			if (v_adj_vertices_2[k_2].first == corner_2.first) break;
			k_2++;
		}
		Vertex_Element next = v_vertices[(i + j + 1) % n];
		Vertex_Element a_prev_2 = v_adj_vertices_2[(k_2 + m_2 - 1) % m_2];
		if (!Vertex::right_turn(a_prev_2.first, corner_2.first, next.first)) {
			i += j - 1;
			continue;
		}

		// If the function reachs this point, it means that all the conditions for slicing the cell between corner_1 and corner_2 are respected
		edges_to_slice.push_back(std::make_pair(i, j));

		// Iterates to save time
		i += j - 1;
	}

	return (!edges_to_slice.empty());
}



void Face::merged_halfedges(Face* f_1, Face* f_2, Edge* e, list<HalfEdge *> & merged)
{
	merged_halfedges(f_1->edges, f_2->edges, e, merged);
}


void Face::merged_halfedges(list<HalfEdge *> & v_1, list<HalfEdge *> & v_2, Edge* e, list<HalfEdge *> & merged)
{
	// We want to obtain a list of edges of a merged facet (f_1 U f_2) in such a way 
	// that the order of the edges is still clockwise, given e adjacent to f_1 and f_2.
	// To this end, a natural idea is to insert all the edges of f_2 in f_1's list
	// at the location of e.
	//
	// Here, f_1 and f_2 are symbolized by their respective vectors of edges.
	merged.clear();

	for (list<HalfEdge *>::iterator it_h1 = v_1.begin() ; it_h1 != v_1.end() ; it_h1++) {
		HalfEdge* h1 = (*it_h1);
		if (h1->e != e) {
			merged.push_back(h1);
		} else {
			// We insert all the halfedges of f_2, except one that we need to locate
			list<HalfEdge *>::iterator it_h2 = v_2.begin();
			while ((*it_h2)->e != e) it_h2++;

			// Now we loop on all other edges, until we are back to e
			list<HalfEdge *>::iterator it_h3 = it_h2;
			if (++it_h3 == v_2.end()) it_h3 = v_2.begin();
			while ((*it_h3)->e != e) {
				merged.push_back(*it_h3);
				if (++it_h3 == v_2.end()) {
					it_h3 = v_2.begin();
				}
			}
		}
	}
}


void Face::merged_halfedges(Face* f_1, Face* f_2, Edge* e_1, Edge* e_2, list<HalfEdge *> & merged)
{
	merged_halfedges(f_1->edges, f_2->edges, e_1, e_2, merged);
}


void Face::merged_halfedges(list<HalfEdge *> & l_1, list<HalfEdge *> & l_2, Edge* e_1, Edge* e_2, list<HalfEdge *> & merged)
{
	// We want to obtain a list of edges of a merged facet (f_1 U f_2) in such a way that the order of the edges is still clockwise, 
	// given a sequence of edges e = (e_1 e_2) common to f_1 and f_2.
	// Here, the first and the last element of e are given with respect to the sequence of edges in f_1.
	merged.clear();
	
	// We start by locating e_2 in l_1.
	list<HalfEdge*>::iterator it_h = l_1.begin();
	while ((*it_h)->e != e_2) {
		if (++it_h == l_1.end()) it_h = l_1.begin();
	}

	// We get an iterator that points to e_2.
	// We need to make it point to the next element.
	if (++it_h == l_1.end()) it_h = l_1.begin();
	
	// We can now iteratively copy all the edges of l_1 until find e_1.
	while ((*it_h)->e != e_1) {
		merged.push_back(*it_h);
		if (++it_h == l_1.end()) it_h = l_1.begin();
	}

	// Same process for the edges of l_2
	// Except that the roles played by e_1 and e_2 are inversed
	it_h = l_2.begin();
	while ((*it_h)->e != e_1) {
		if (++it_h == l_2.end()) it_h = l_2.begin();
	}

	// it_h is currently pointing to e_1, lets make it point to the next element
	if (++it_h == l_2.end()) it_h = l_2.begin();

	// We can now iteratively copy all the edges of l_2 until finding e_2.
	while ((*it_h)->e != e_2) {
		merged.push_back(*it_h);
		if (++it_h == l_2.end()) it_h = l_2.begin();
	}
}
