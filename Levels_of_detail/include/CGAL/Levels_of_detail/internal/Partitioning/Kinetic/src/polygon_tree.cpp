#include "../include/support_plane_objects.h"
#include "../include/intersection_line.h"
#include "../include/polygon_vertex.h"
#include "../include/polygon_edge.h"
#include "../include/support_plane.h"
#include "../include/polygon_group.h"
#include "../include/universe.h"
#include "../include/polygon.h"

namespace Skippy 
{
	Polygon_Tree::Polygon_Tree()
	{
		is_node = false;
		line = nullptr;
		parent = subtree_plus = subtree_minus = nullptr;
	}



	Polygon_Tree::Polygon_Tree(const int id_plane, const int seed, Polygon_Directions* D)
		: Polygon_Tree()
	{
		polygon = new Polygon(id_plane, seed, D);
	}


	Polygon_Tree::Polygon_Tree(const int id_plane, const int seed, Polygon_Directions* D, const CGAL_Vector_2 & OM, const FT & t)
		: Polygon_Tree()
	{
		polygon = new Polygon(id_plane, seed, D, OM, t);
	}



	Polygon_Tree::Polygon_Tree(Polygon_Tree* _parent, Polygon* _polygon)
		: Polygon_Tree()
	{
		parent = _parent;
		polygon = _polygon;
	}



	Polygon_Tree::~Polygon_Tree()
	{
		if (is_node) {
			delete subtree_plus;
			delete subtree_minus;
		} else {
			if (polygon != nullptr) delete polygon;
		}
	}



	bool Polygon_Tree::split(Intersection_Line* I, const FT & t, const int K, bool pair_constrained_vertices)
	{
		std::list<Polygon_Vertex*> V;
		return split(I, V, t, K, pair_constrained_vertices);
	}



	bool Polygon_Tree::split(Intersection_Line* I, std::list<Polygon_Vertex *> & V, const FT & t, const int K, bool pair_constrained_vertices)
	{
		if (is_node) {
			// If the tree is a node, then we apply the function to its two subtrees
			bool res_plus = subtree_plus->split(I, V, t, K, pair_constrained_vertices);
			bool res_minus = subtree_minus->split(I, V, t, K, pair_constrained_vertices);
			return (res_plus || res_minus);

		} else {
			// Otherwise, we check if one of the edges of the polygon is not cut by l
			std::map<int, Sign> signs;
			std::list<Polygon_Vertex *> vp, vn, vz;
			std::list<Polygon_Edge *> ep, en;
			Polygon_Edge *ep_front, *ep_back, *en_front, *en_back;

			if (polygon->intersection_exists(I, t, signs, vp, ep, ep_front, ep_back, vn, en, en_front, en_back, vz)) {

				// Given two sequences of vertices and edges strictly located on the positive (resp. negative) side of I,
				// our aim is to create the 4 vertices and 6 edges that would allow to split the current polygon into two subpolygons.
				// Let us consider the subpolygon on the positive side of I : we find the intersections of I with ep_front and ep_back.

				Polygon_Vertex *vp_front, *vp_back, *vn_front, *vn_back;

				vp_front = ep_front->intersection(I, PLUS, t, K);
				vp_back = ep_back->intersection(I, PLUS, t, K);
				const Constraint C_ep_front = ep_front->get_constraint(), C_ep_back = ep_back->get_constraint();
				V.push_back(vp_front);
				V.push_back(vp_back);

				// Same process for the other subpolygon
				// Avoids to recompute vn->M and vn->dM, when they are the same as before

				bool matches_vn_front_vp_front = false, matches_vn_front_vp_back = false;
				bool matches_vn_back_vp_front = false, matches_vn_back_vp_back = false;

				CGAL_Point_2 vp_front_M = vp_front->get_M();
				CGAL_Vector_2 vp_front_dM = CGAL::NULL_VECTOR;
				if (vp_front->type == RUNNING_VERTEX) vp_front_dM = vp_front->to_r()->dM;

				CGAL_Point_2 vp_back_M = vp_back->get_M();
				CGAL_Vector_2 vp_back_dM = CGAL::NULL_VECTOR;
				if (vp_back->type == RUNNING_VERTEX) vp_back_dM = vp_back->to_r()->dM;

				if (en_front == ep_front) {	
					vn_front = en_front->intersection(I, MINUS, t, vp_front_M, vp_front_dM, K);
					matches_vn_front_vp_front = true;
				} else if (en_front == ep_back) {
					vn_front = en_front->intersection(I, MINUS, t, vp_back_M, vp_back_dM, K);
					matches_vn_front_vp_back = true;
				} else {
					vn_front = en_front->intersection(I, MINUS, t, K);
				}

				if (en_back == ep_front) {
					vn_back = en_back->intersection(I, MINUS, t, vp_front_M, vp_front_dM, K);
					matches_vn_back_vp_front = true;
				} else if (en_back == ep_back) {
					vn_back = en_back->intersection(I, MINUS, t, vp_back_M, vp_back_dM, K);
					matches_vn_back_vp_back = true;
				} else {
					vn_back = en_back->intersection(I, MINUS, t, K);
				}

				if (pair_constrained_vertices) {
					if (Polygon_Vertex_R* vn_front_r = vn_front->to_r()) {
						if (matches_vn_front_vp_front) {
							assert(vp_front->to_r() != nullptr);
							Polygon_Vertex_R::set_paired_vertices(vp_front->to_r(), vn_front_r);
						} else if (matches_vn_front_vp_back) {
							assert(vp_back->to_r() != nullptr);
							Polygon_Vertex_R::set_paired_vertices(vp_back->to_r(), vn_front_r);
						}
					}

					if (Polygon_Vertex_R* vn_back_r = vn_back->to_r()) {
						if (matches_vn_back_vp_front) {
							assert(vp_front->to_r() != nullptr);
							Polygon_Vertex_R::set_paired_vertices(vp_front->to_r(), vn_back_r);
						} else if (matches_vn_back_vp_back) {
							assert(vp_back->to_r() != nullptr);
							Polygon_Vertex_R::set_paired_vertices(vp_back->to_r(), vn_back_r);
						}
					}
				}

				const Constraint C_en_front = en_front->get_constraint(), C_en_back = en_back->get_constraint();
				V.push_back(vn_front);
				V.push_back(vn_back);

				// Let us now create three edges for side of I
				// Before that we must find the vertices that link these vertices to the rest of the polygons
				Polygon_Vertex* vp_first = (signs[ep_front->v1->id_object] == PLUS ? ep_front->v1 : ep_front->v2);
				Polygon_Vertex* vp_last = (signs[ep_back->v1->id_object] == PLUS ? ep_back->v1 : ep_back->v2);
				Polygon_Vertex* vn_first = (signs[en_front->v1->id_object] == MINUS ? en_front->v1 : en_front->v2);
				Polygon_Vertex* vn_last = (signs[en_back->v1->id_object] == MINUS ? en_back->v1 : en_back->v2);

				// Deletes old edges
				std::set<Polygon_Edge *> to_delete;
				to_delete.insert(ep_front); to_delete.insert(ep_back); to_delete.insert(en_front); to_delete.insert(en_back);

				for (std::set<Polygon_Edge *>::iterator it_e = to_delete.begin(); it_e != to_delete.end(); it_e++) {
					std::list<Polygon_Edge *>::iterator it_me = std::find(polygon->edges.begin(), polygon->edges.end(), (*it_e));
					polygon->edges.erase(it_me);
					delete (*it_e);
				}

				// Deletes vertices that are actually located on I,
				// and which are not considered in the list of vertices vp and vn
				for (std::list<Polygon_Vertex*>::iterator it_v = vz.begin(); it_v != vz.end(); it_v++) {
					assert(std::find(vp.begin(), vp.end(), (*it_v)) == vp.end());
					assert(std::find(vn.begin(), vn.end(), (*it_v)) == vn.end());
					delete (*it_v);
				}

				// Creates new edges
				Polygon_Edge* shortened_ep_front;
				if (C_ep_front.first == nullptr) {
					shortened_ep_front = new Polygon_Edge(I->id_plane, vp_first, vp_front);
				} else {
					shortened_ep_front = new Polygon_Edge(I->id_plane, vp_first, vp_front, C_ep_front);
				}

				Polygon_Edge* shortened_ep_back;
				if (C_ep_back.first == nullptr) {
					shortened_ep_back = new Polygon_Edge(I->id_plane, vp_last, vp_back);
				} else {
					shortened_ep_back = new Polygon_Edge(I->id_plane, vp_last, vp_back, C_ep_back);
				}

				Polygon_Edge* shortened_en_front;
				if (C_en_front.first == nullptr) {
					shortened_en_front = new Polygon_Edge(I->id_plane, vn_first, vn_front);
				} else {
					shortened_en_front = new Polygon_Edge(I->id_plane, vn_first, vn_front, C_en_front);
				}

				Polygon_Edge* shortened_en_back;
				if (C_en_back.first == nullptr) {
					shortened_en_back = new Polygon_Edge(I->id_plane, vn_last, vn_back);
				} else {
					shortened_en_back = new Polygon_Edge(I->id_plane, vn_last, vn_back, C_en_back);
				}

				Polygon_Edge* vp_front_to_back = new Polygon_Edge(I->id_plane, vp_front, vp_back, I, PLUS);
				Polygon_Edge* vn_front_to_back = new Polygon_Edge(I->id_plane, vn_front, vn_back, I, MINUS);

				// We add the created objects to the lists of vertices and edges computed before
				// These lists will then be used to create the two subpolygons
				vp.push_back(vp_front); vp.push_back(vp_back);
				vn.push_back(vn_front); vn.push_back(vn_back);
				ep.push_back(shortened_ep_front); ep.push_back(shortened_ep_back); ep.push_back(vp_front_to_back);
				en.push_back(shortened_en_front); en.push_back(shortened_en_back); en.push_back(vn_front_to_back);

				int seed = polygon->seed;
				Polygon* subpolygon_plus = new Polygon(seed, vp, ep);
				Polygon* subpolygon_minus = new Polygon(seed, vn, en);

				// Last step of the algorithm : division of the node
				turn_into_node(I, subpolygon_plus, subpolygon_minus, true);

				return true;
			}

			return false;
		}
	}



	void Polygon_Tree::turn_into_node(Intersection_Line* I, Polygon* subpolygon_plus, Polygon* subpolygon_minus, bool destroy_polygon)
	{
		is_node = true;
		line = I;

		subtree_plus = new Polygon_Tree(this, subpolygon_plus);
		subtree_minus = new Polygon_Tree(this, subpolygon_minus);

		if (destroy_polygon) {
			polygon->vertices.clear();
			polygon->edges.clear();
			delete polygon;
		}

		polygon = nullptr;
	}



	void Polygon_Tree::get_polygons(std::list<Polygon*> & P)
	{
		// Gets all the leaves of a Polygon_Tree
		if (is_node) {
			subtree_plus->get_polygons(P);
			subtree_minus->get_polygons(P);
		} else {
			P.push_back(polygon);
		}
	}


	/*
	void Polygon_Tree::get_final_polygon_description(int id_plane, std::list<std::list<CGAL_Point_3> > & polygons, std::list<CGAL_Color> & colors)
	{
		// Step 1
		// Gets a list containing all the polygons listed in the tree
		std::list<Polygon*> P;
		get_polygons(P);
		int p_i = -1;

		// Step 2
		// Groups adjacent polygons, i.e. polygons separated by a single edge without segment on it
		std::list<Polygon_Group*> G;
		for (std::list<Polygon*>::iterator it_p = P.begin(); it_p != P.end(); it_p++) {
			std::list<std::list<Polygon_Group*>::iterator> G_adj;
			for (std::list<Polygon_Group*>::iterator it_g = G.begin() ; it_g != G.end() ; it_g++) {
				if ((*it_g)->is_adjacent_to(*it_p)) {
					G_adj.push_back(it_g);
				}
			}

			// If no adjacent polygon has been found, then we create a new group of adjacent polygons with a single element
			if (G_adj.empty()) {
				G.push_back(new Polygon_Group(*it_p));
			}

			// If a group of adjacent polygons is found, then we simply add P[i] to that group
			// In case there are more than one compatible groups, then we add P[i] to the first one
			// and merge all other groups to the first one
			else {
				std::list<std::list<Polygon_Group*>::iterator>::iterator it_1st_group = G_adj.begin();
				(*(*it_1st_group))->append(*it_p);

				std::list<std::list<Polygon_Group*>::iterator>::iterator it_kth_group = it_1st_group;
				while (++it_kth_group != G_adj.end()) {
					(*(*it_1st_group))->merge(*(*it_kth_group));
					G.erase(*it_kth_group);
				}
			}
		}

		// Step 3
		// Exports the result to the in-out arguments of this function

		for (std::list<Polygon_Group*>::iterator it_g = G.begin() ; it_g != G.end() ; it_g++) {
			std::vector<CGAL_Point_2> G_ch;
			(*it_g)->get_convex_hull(G_ch);

			// Gets the plane equation
			Support_Plane* SP = Universe::map_of_planes[id_plane];
			const CGAL_Plane & H = SP->plane;
			const CGAL_Vector_3 & h = H.orthogonal_vector();

			// Backprojects all points listed in the convex hull of the group of polygons
			std::vector<CGAL_Point_3> V;
			V.reserve(G_ch.size());
			for (size_t i = 0 ; i < G_ch.size() ; i++) {
				V.push_back(SP->backproject(G_ch[i]));
			}

			// Prints the result inside Q
			std::list<CGAL_Point_3> polygon(1, V[0]);
			bool forward = (CGAL::cross_product(V[1] - V[0], V[2] - V[0]) * h > 0);
			if (forward) {
				for (size_t i = 1; i < V.size(); i++) polygon.push_back(V[i]);
			} else {
				for (size_t i = V.size() - 1; i > 0; i--) polygon.push_back(V[i]);
			}
			polygons.push_back(polygon);
			colors.push_back(SP->color);
		}

		for (std::list<Polygon_Group*>::iterator it_g = G.begin(); it_g != G.end(); it_g++) delete (*it_g);
	}
	*/


	Polygon* Polygon_Tree::remove_reference_to_polygon()
	{
		assert(!is_node);

		Polygon* P = polygon;
		polygon = nullptr;

		return P;
	}


	void Polygon_Tree::remove_reference_to_polygons()
	{
		if (is_node) {
			subtree_plus->remove_reference_to_polygons();
			subtree_minus->remove_reference_to_polygons();
		} else {
			polygon = nullptr;
		}
	}



	void Polygon_Tree::sort_intersection_points(Intersection_Line* I, const std::list<Polygon_Vertex*> & V, const FT & t, std::vector<Polygon_Vertex*> & V_res)
	{
		// As we enter this function, we assume that points included in V
		// result from the split of a polygon into two subpolygons.
		// Here we are looking for a 4-uplet with the furthest points on I.

		const CGAL_Line_2 & L = I->line;
		const CGAL_Point_2 L_m = L.point();
		const CGAL_Vector_2 L_u = L.to_vector();
		bool L_horizontal = (L_u.x() == 0);

		// I = (L_m + u * L_u).
		// The idea is to find the abscissa u that correspond to each point of I,
		// and sort all elements by u.

		Support_Plane* SP = Universe::map_of_planes[I->id_plane];

		std::vector<std::pair<Polygon_Vertex*, FT> > R;
		for (std::list<Polygon_Vertex*>::const_iterator it_v = V.begin() ; it_v != V.end() ; ++it_v) {
			Polygon_Vertex* v = (*it_v);
			CGAL_Point_2 V_t;
			FT u;

			if (Polygon_Vertex_S* v_s = v->to_s()) {
				V_t = v_s->get_M();
			} else {
				V_t = v->to_r()->pt(t);
			}

			if (L_horizontal) {
				u = (V_t.y() - L_m.y()) / L_u.y();
			} else {
				u = (V_t.x() - L_m.x()) / L_u.x();
			}

			R.push_back(std::make_pair(v, u));
		}

		// Sorts by ascending order on u.
		// In the end, the desired points are V_res[0], V_res[1], V_res[n - 2] and V_res[n - 1]
		// where n is the size of the vectors R and V_res.

		struct _Comparator_Intersection_Points {
			bool operator() (const std::pair<Polygon_Vertex*, FT> & L, const std::pair<Polygon_Vertex*, FT> & R) {
				return L.second < R.second;
			}
		} Comparator_Intersection_Points;
		std::sort(R.begin(), R.end(), Comparator_Intersection_Points);

		V_res.clear();
		V_res.reserve(R.size());
		for (size_t i = 0 ; i < R.size() ; ++i) {
			V_res.push_back(R[i].first);
		}
	}
}