#include "../include/polygon_set.h"
#include "../include/polygon_vertex.h"
#include "../include/polygon_edge.h"
#include "../include/support_plane.h"
#include "../include/universe.h"
#include "../include/intersection_line.h"
#include <queue>

#include <list>
#include <iterator>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>


namespace Skippy 
{
	Polygon_Set::Polygon_Set(int _id_plane, const std::vector<Intersection_Line*> & L)
		: id_plane (_id_plane)
	{
		// Builds a map of correspondances between lines and bits of a signature
		int entry = -1;
		for (std::vector<Intersection_Line*>::const_iterator it_l = L.begin(); it_l != L.end(); it_l++) {
			dictionary[*it_l] = ++entry;
		}
	}


	Polygon_Set::~Polygon_Set()
	{
		dictionary.clear();

		for (std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::iterator it_s = nodes.begin(); it_s != nodes.end(); it_s++) {
			delete it_s->second;
		}
	}


	void Polygon_Set::insert(const Signature & S, Polygon* P)
	{
		// Searches for an entry with the same signature
		// If it does exist, we only add P to the cell
		// (This branch will only be used for initialization of coplanar polygons)
		// Otherwise we create a new cell containing the polygon

		std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::iterator it_s = nodes.find(S);

		if (it_s != nodes.end()) {
			it_s->second->push(P);
		} else {
			nodes[S] = new Polygon_Node(S, P);
		}
	}


	bool Polygon_Set::exists(const Signature & S)
	{
		std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::iterator it_s = nodes.find(S);
		return (it_s != nodes.end());
	}


	bool Polygon_Set::exists(const Signature & S, const int seed)
	{
		std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::iterator it_s = nodes.find(S);
		if (it_s == nodes.end()) return false;

		Polygon_Node* C = it_s->second;
		for (std::list<Polygon*>::const_iterator it_p = C->polygons_begin(); it_p != C->polygons_end(); it_p++) {
			Polygon* P = (*it_p);
			if (P->seed == seed) return true;
		}

		return false;
	}

	
	std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::const_iterator Polygon_Set::cells_begin() const
	{
		return nodes.cbegin();
	}
		

	std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::const_iterator Polygon_Set::cells_end() const
	{
		return nodes.cend();
	}


	Polygon* Polygon_Set::get_adjacent_polygon(Polygon* P_ts, Intersection_Line* I)
	{
		// Constructs the signature that corresponds to the polygon adjacent to the one containing v,
		// but located on the other side of the intersection line I, reached by v	
		std::vector<bool> S_os = get_adjacent_polygons_signature(P_ts, I);

		// Searches for such a polygon, returns it
		std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::iterator it_s = nodes.find(S_os);
		if (it_s == nodes.end()) {
			return nullptr;
		} else {
			for (std::list<Polygon*>::const_iterator it_p = it_s->second->polygons_begin(); it_p != it_s->second->polygons_end(); it_p++) {
				Polygon* P = (*it_p);
				if (P->seed == P_ts->seed) {
					return P;
				}
			}
			return nullptr;
			//return it_s->second->get_unique_polygon();
		}
	}


	std::vector<bool> Polygon_Set::get_adjacent_polygons_signature(Polygon* P, Intersection_Line* I)
	{
		// Gets original signature of P
		std::vector<bool> S_ts = P->get_cell()->get_signature();

		// Sets the bit that corresponds to line I
		std::vector<bool> S_os = S_ts;

		int b = dictionary[I];
		S_os[b] = !S_os[b];

		// Returns result
		return S_os;
	}


	std::vector<bool> Polygon_Set::get_adjacent_polygons_signature(Polygon* P, Intersection_Line* I_1, Intersection_Line* I_2)
	{
		// Same as before, except that we set two bits
		std::vector<bool> S_ts = P->get_cell()->get_signature();

		std::vector<bool> S_os = S_ts;

		int b_1 = dictionary[I_1];
		int b_2 = dictionary[I_2];
		S_os[b_1] = !S_os[b_1];
		S_os[b_2] = !S_os[b_2];

		return S_os;
	}


	void Polygon_Set::get_signature_of_adjacent_cell(std::vector<bool> & S, Intersection_Line* I)
	{
		int b = dictionary[I];
		S[b] = !S[b];
	}



	void Polygon_Set::get_polygons(std::list<Polygon*> & P)
	{
		P.clear();

		// There may be more than one polygon per cell
		// However, assuming that this function is called once all points have stopped propagating,
		// we only take one of those polygons into account (they all have the same definition in the end)

		for (std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::iterator it_s = nodes.begin(); it_s != nodes.end(); it_s++) {
			Polygon_Node* C = it_s->second;
			P.push_back(C->get_one_polygon());
		}
	}


#if 0
	void Polygon_Set::get_polygon_description(std::list<std::list<CGAL_Point_3> > & polygons, std::list<CGAL_Color> & colors, const FT & t)
	{
		typedef CGAL::Partition_traits_2<K>   Traits;
		typedef Traits::Point_2               Point_2;
		typedef Traits::Polygon_2             Polygon_2;
		typedef Polygon_2::Vertex_iterator    Vertex_iterator;

		Support_Plane* SP = Universe::map_of_planes[id_plane];

		for (size_t d = 0; d < SP->polygon_directions.size() ; ++d) {
		
			// Part 1.
			// At time t, the primitive d is split into a set of subpolygons.
			// We want to list all the edges that are part of the contours of the object.

			std::list<CGAL_Segment_2> unordered_ring;

			for (auto it_n1 = nodes.begin() ; it_n1 != nodes.end() ; ++it_n1) {
				Signature S_1 = it_n1->first;
				Polygon_Node* node_1 = it_n1->second;

				for (std::list<Polygon*>::const_iterator it_p1 = node_1->polygons_begin() ; it_p1 != node_1->polygons_end() ; ++it_p1) {
					Polygon* P = (*it_p1);
					if (P->seed == d) {
						if (P->edges.empty()) {

							// Case 1.1
							// The subpolygon is inactive.
							// We need to check if it is on the outer ring of polygons of the primitive at time t.
							// To this end, we loop on the edges of the cell, and check if there is a subpolygon
							// with a similar edge in the cell adjacent to the edge.

							size_t n = node_1->contours_size();
							for (size_t i = 0 ; i < n ; ++i) {
								size_t i_prev = (i == 0 ? n - 1 : i - 1);
								size_t i_next = (i == n - 1 ? 0 : i + 1);
								Constraint C = node_1->get_contour(i);
								Constraint C_prev = node_1->get_contour(i_prev), C_next = node_1->get_contour(i_next);
								
								CGAL_Point_2 A = SP->get_intersection_point(C_prev.first, C.first);
								CGAL_Point_2 B = SP->get_intersection_point(C.first, C_next.first);

								if (!exists_similar_edge_in_adjacent_cell(d, S_1, C, C_prev, C_next, A, B, t)) {
									unordered_ring.push_back(CGAL_Segment_2(A, B));
								}
							}
							
						} else {

							// Case 1.2
							// The subpolygon is propagating.

							for (auto it_e = P->edges.begin() ; it_e != P->edges.end() ; ++it_e) {
								Polygon_Edge* e = (*it_e);
								Polygon_Vertex *v1 = e->v1, *v2 = e->v2;
								Polygon_Vertex_R *v1_r = v1->to_r(), *v2_r = v2->to_r();
								Polygon_Vertex_S *v1_s = v1->to_s(), *v2_s = v2->to_s();

								if (v1_s != nullptr && v2_s != nullptr) {

									// Case 1.2.1
									// e = (v1 v2) is constant.
									// This is the same case as before.

									Constraint C = e->get_constraint();
									Constraint C_prev = v1_s->get_other_constraint(C);
									Constraint C_next = v2_s->get_other_constraint(C);

									CGAL_Point_2 A = SP->get_intersection_point(C_prev.first, C.first);
									CGAL_Point_2 B = SP->get_intersection_point(C.first, C_next.first);

									if (!exists_similar_edge_in_adjacent_cell(d, S_1, C, C_prev, C_next, A, B, t)) {										
										unordered_ring.push_back(CGAL_Segment_2(A, B));
									}

								} else if ((v1_s != nullptr && v2_s == nullptr) || (v1_s == nullptr && v2_s != nullptr)) {

									// Case 1.2.2
									// e = (v_s v_r), where v_s is still and v_r is propagating.
									// v_r and e are necessarily constrained.
									// Like before, we check if there exists a similar edge on the other side of e->constraint.

									Polygon_Vertex_S* v_s = (v1_s != nullptr ? v1_s : v2_s);
									Polygon_Vertex_R* v = (v1_r != nullptr ? v1_r : v2_r);

									Constraint C = e->get_constraint();
									Constraint C_s = v_s->get_other_constraint(C);

									if (!exists_similar_edge_in_adjacent_cell(d, S_1, C, v, v_s, t)) {
										CGAL_Point_2 A = v->pt(t);
										CGAL_Point_2 B = SP->get_intersection_point(C.first, C_s.first);
										unordered_ring.push_back(CGAL_Segment_2(A, B));
									}

								} else {

									// Case 1.2.3
									// e = (v1_r v2_r) is added to the list of contours,
									// except if e interesects a line at time t : in this case, we go for the usual check.

									std::list<Event_Vertex_Line*> E_1, E_2;
									v1_r->get_upcoming_event(E_1);
									v2_r->get_upcoming_event(E_2);
									CGAL_Point_2 A = v1_r->pt(t), B = v2_r->pt(t);
									if (E_1.front()->t_intersectant != t || E_1.front()->t_intersectant != t) {
										unordered_ring.push_back(CGAL_Segment_2(A, B));
									} else {
										assert(E_1.size() != 1);
										Intersection_Line* I = SP->get_line_by_identifier(E_1.front()->intersected);
										if (!exists_similar_edge_in_adjacent_cell(d, S_1, I, A, B, t)) {
											unordered_ring.push_back(CGAL_Segment_2(A, B));
										}
									}
								}
							}
						}
					}
				}
			}

			// Part 2.
			// Connects elements of the ring
			
			std::list<CGAL_Point_2> ordered_ring;

			CGAL_Point_2 anchor = unordered_ring.front().source();
			ordered_ring.push_back(unordered_ring.front().source());
			unordered_ring.pop_front();

			while (!unordered_ring.empty()) {
				std::list<CGAL_Segment_2>::iterator it_s = unordered_ring.begin();

				while (true) {
					if (it_s->source() == anchor) {
						anchor = it_s->target();
						break;
					} else if (it_s->target() == anchor) {
						anchor = it_s->source();
						break;
					} else {
						++it_s;
					}
				}
				
				ordered_ring.push_back(anchor);
				unordered_ring.erase(it_s);
			}

			// Part 3.
			// Removes collinear vertices

			std::list<CGAL_Point_2>::iterator it_prev, it_curr, it_next;
			it_prev = --ordered_ring.end(), it_curr = ordered_ring.begin(), it_next = ++ordered_ring.begin();

			std::list<CGAL_Point_2> simplified_ring;

			while (it_curr != ordered_ring.end()) {
				CGAL_Point_2 A = (*it_prev), B = (*it_curr), C = (*it_next);
				if (!CGAL::collinear(A, B, C)) {
					simplified_ring.push_back(B);
				}
				++it_prev, ++it_curr, ++it_next;
				if (it_prev == ordered_ring.end()) it_prev = ordered_ring.begin();
				if (it_next == ordered_ring.end()) it_next = ordered_ring.begin();
			}

			// Part 4.
			// Partitions the obtained polygon into a set of convex polygons.

			Polygon_2 inner_polygon (simplified_ring.begin(), simplified_ring.end());
			if (inner_polygon.orientation() == CGAL::CLOCKWISE) inner_polygon.reverse_orientation();

			std::list<Polygon_2> inner_subpolygons;
			CGAL::approx_convex_partition_2(inner_polygon.vertices_begin(), inner_polygon.vertices_end(), std::back_inserter(inner_subpolygons));

			for (std::list<Polygon_2>::const_iterator it_sp = inner_subpolygons.begin() ; it_sp != inner_subpolygons.end() ; ++it_sp) {
				Polygon_2 polygon = (*it_sp);
				std::list<CGAL_Point_3> P;
				get_sequence_of_3d_vertices(polygon, P);
				polygons.push_back(P);
				colors.push_back(SP->color);
			}
		}
	}


	bool Polygon_Set::exists_similar_edge_in_adjacent_cell(int d, const Signature & S, const Constraint & C, const Constraint & C_prev, const Constraint & C_next, const CGAL_Point_2 & A, const CGAL_Point_2 & B, const FT & t)
	{
		Signature S_adj = S;
		get_signature_of_adjacent_cell(S_adj, C.first);

		auto it_n = nodes.find(S_adj);
		if (it_n == nodes.end()) return false;

		Polygon_Node* node = it_n->second;

		for (std::list<Polygon*>::const_iterator it_p = node->polygons_begin() ; it_p != node->polygons_end() ; ++it_p) {
			Polygon* P = (*it_p);
			if (P->seed == d) {
				if (P->edges.empty()) {
					size_t n_c = node->contours_size.size();
					for (size_t i = 0 ; i < n_c ; ++i) {
						Constraint D_curr = node->get_contour(i);
						Constraint D_prev = node->get_contour(i == 0 ? n_c - 1 : i - 1);
						Constraint D_next = node->get_contour(i == n_c - 1 ? 0 : i + 1);
						if (C == D_curr) {
							return ((D_prev == C_prev && D_next == C_next) || (D_prev == C_next && D_next == C_prev));
						}
					}
				} else {
					for (std::list<Polygon_Edge*>::iterator it_e = P->edges.begin() ; it_e != P->edges.end() ; ++it_e) {
						Polygon_Edge* e = (*it_e);
						Polygon_Vertex *v1 = e->v1, *v2 = e->v2;
						Polygon_Vertex_R *v1_r = v1->to_r(), *v2_r = v2->to_r();
						Polygon_Vertex_S *v1_s = v2->to_s(), *v2_s = v2->to_s();
						CGAL_Point_2 S, T;
						// ...
					}
				}
			}
		}

		return false;
	}
#endif


	/*void Polygon_Set::get_sequence_of_3d_vertices(const CGAL::Partition_traits_2<K>::Polygon_2 & polygon, std::list<CGAL_Point_3> & P)
	{
		Support_Plane* SP = Universe::map_of_planes[id_plane];

		const CGAL_Plane & H = SP->plane;
		const CGAL_Vector_3 & h = H.orthogonal_vector();

		std::vector<CGAL_Point_3> V;
		size_t n = polygon.size();
		V.reserve(n);

		for (CGAL::Partition_traits_2<K>::Polygon_2::Vertex_iterator it_v = polygon.vertices_begin() ; it_v != polygon.vertices_end() ; ++it_v) {
			V.push_back(SP->backproject(*it_v));
		}

		P.push_back(V[0]);
		bool forward = (CGAL::cross_product(V[1] - V[0], V[2] - V[0]) * h > 0);
		if (forward) {
			for (size_t i = 1; i < n; i++) P.push_back(V[i]);
		} else {
			for (size_t i = n - 1; i > 0; i--) P.push_back(V[i]);
		}
	}*/


#if 0
	void Polygon_Set::get_polygon_description(std::list<std::list<CGAL_Point_3> > & polygons, std::list<CGAL_Color> & colors, const double t)
	{
		typedef CGAL::Partition_traits_2<K>                        Traits;
		typedef Traits::Point_2                                     Point_2;
		typedef Traits::Polygon_2                                   Polygon_2;
		typedef Polygon_2::Vertex_iterator                          Vertex_iterator;

		double prec = 1.0 / (1 << 30) / (1 << 10);
		FT::set_relative_precision_of_to_double(prec);

		Support_Plane* SP = Universe::map_of_planes[id_plane];

		for (size_t d = 0 ; d < SP->polygon_directions.size() ; ++d) {
			
			// For every primitive which has propagated,
			// we loop on its subpolygons which are either active or inactive.

			std::list<Polygon*> active_polygons;
			std::list<Polygon_Node*> inactive_nodes;

			for (std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::iterator it_c = nodes.begin(); it_c != nodes.end(); it_c++) {
				Polygon_Node* N = it_c->second;
				for (std::list<Polygon*>::const_iterator it_p = N->polygons_begin() ; it_p != N->polygons_end() ; ++it_p) {
					Polygon* P = (*it_p);
					if (P->seed == d) {
						if (P->vertices.empty()) {
							inactive_nodes.push_back(N);
						} else {
							active_polygons.push_back(P);
						}
					}
				}
			}

			// Active polygons are printed as they are.

			for (std::list<Polygon*>::const_iterator it_p = active_polygons.begin() ; it_p != active_polygons.end() ; ++it_p) {
				std::list<CGAL_Point_3> P;
				get_sequence_of_3d_vertices((*it_p), t, P);
				polygons.push_back(P);
				colors.push_back(SP->color);
			}

			// Inactive polygons are grouped together.

			std::list<CGAL_Point_2> V;
			Polygon_Group* G = new Polygon_Group(inactive_nodes);
			G->make_borders(false);
			for (std::list<G_Vertex>::const_iterator it_v = G->borders_v_begin() ; it_v != G->borders_v_end() ; ++it_v) {
				Constraint C_1 = it_v->first, C_2 = it_v->second;
				V.push_back(SP->get_intersection_point(C_1.first, C_2.first));
			}

			Polygon_2 inner_polygon (V.begin(), V.end());
			if (inner_polygon.orientation() == CGAL::CLOCKWISE) inner_polygon.reverse_orientation();

			std::list<Polygon_2> inner_subpolygons;
			CGAL::approx_convex_partition_2(inner_polygon.vertices_begin(), inner_polygon.vertices_end(), std::back_inserter(inner_subpolygons));

			for (std::list<Polygon_2>::const_iterator it_sp = inner_subpolygons.begin() ; it_sp != inner_subpolygons.end() ; ++it_sp) {
				Polygon_2 polygon = (*it_sp);
				std::list<CGAL_Point_3> P;
				get_sequence_of_3d_vertices(polygon, P);
				polygons.push_back(P);
				colors.push_back(SP->color);
			}
		}
	}
#endif


	void Polygon_Set::get_sequence_of_3d_vertices(Polygon* polygon, const double t, std::list<CGAL_Point_3> & P)
	{
		Support_Plane* SP = Universe::map_of_planes[id_plane];

		const CGAL_Plane & H = SP->plane;
		const CGAL_Vector_3 & h = H.orthogonal_vector();

		std::list<Polygon_Vertex*> & vertices = polygon->vertices;
		Polygon_Vertex *v_init = *(vertices.begin()), *v_curr = v_init, *v_next = nullptr;
		Polygon_Edge* e_prev = nullptr;

		std::vector<CGAL_Point_3> V;
		size_t n = vertices.size();
		V.reserve(n);

		// Gets a list of backprojected vertices
		while (true) {
			CGAL_Point_2 M = v_curr->pt(t);
			V.push_back(SP->backproject(M));
			Polygon_Edge* e = (v_curr->e1 == e_prev ? v_curr->e2 : v_curr->e1);
			Polygon_Vertex* v_next = (e->v1 == v_curr ? e->v2 : e->v1);
			if (v_next == v_init) {
				break;
			} else {
				e_prev = e;
				v_curr = v_next;
			}
		}

		// Builds facets

		P.clear();
		P.push_back(V[0]);
		bool forward = (CGAL::cross_product(V[1] - V[0], V[2] - V[0]) * h > 0);
		if (forward) {
			for (size_t i = 1; i < n; i++) P.push_back(V[i]);
		} else {
			for (size_t i = n - 1; i > 0; i--) P.push_back(V[i]);
		}
	}


	
	void Polygon_Set::get_polygon_description(std::list<std::list<CGAL_Point_3> > & polygons, std::list<CGAL_Color> & colors, const double t)
	{
		Support_Plane* SP = Universe::map_of_planes[id_plane];

		// Gets the description of all polygons at time t.
		// Their assigned colors are those of their support planes.

		for (std::map<Signature, Polygon_Node*, Vector_Bool_Comparator>::iterator it_c = nodes.begin(); it_c != nodes.end(); it_c++) {
			Polygon_Node* C = it_c->second;

			if (C->contours_are_empty()) {

				// If the contours of the cell aren't known so far,
				// then we get a description of the polygons that propagate.
				
				for (std::list<Polygon*>::const_iterator it_p = C->polygons_begin(); it_p != C->polygons_end(); it_p++) {
					Polygon* polygon = (*it_p);
					if (polygon->vertices.empty()) continue;

					std::list<CGAL_Point_3> P;
					get_sequence_of_3d_vertices(polygon, t, P);
					polygons.push_back(P);
					colors.push_back(SP->color);
				}
				
			} else {
			
				std::vector<CGAL_Point_3> V;
				size_t n = C->contours_size();
				V.reserve(n);

				for (size_t i = 0 ; i < n ; ++i) {
					size_t i_next = (i == n - 1 ? 0 : i + 1);
					Constraint C_curr = C->get_contour(i);
					Constraint C_next = C->get_contour(i_next);
					CGAL_Point_2 M = SP->get_intersection_point(C_curr.first, C_next.first);
					V.push_back(SP->backproject(M));
				}

				const CGAL_Plane & H = SP->plane;
				const CGAL_Vector_3 & h = H.orthogonal_vector();

				std::list<CGAL_Point_3> P(1, V[0]);
				bool forward = (CGAL::cross_product(V[1] - V[0], V[2] - V[0]) * h > 0);
				if (forward) {
					for (size_t i = 1; i < n; i++) P.push_back(V[i]);
				} else {
					for (size_t i = n - 1; i > 0; i--) P.push_back(V[i]);
				}
				polygons.push_back(P);
				colors.push_back(SP->color);
			}
		}
	}









	void Polygon_Set::get_sequence_of_3d_vertices(Polygon* polygon, const FT & t, std::list<CGAL_Point_3> & P)
	{
		Support_Plane* SP = Universe::map_of_planes[id_plane];

		const CGAL_Plane & H = SP->plane;
		const CGAL_Vector_3 & h = H.orthogonal_vector();

		std::list<Polygon_Vertex*> & vertices = polygon->vertices;
		Polygon_Vertex *v_init = *(vertices.begin()), *v_curr = v_init, *v_next = nullptr;
		Polygon_Edge* e_prev = nullptr;

		std::vector<CGAL_Point_3> V;
		size_t n = vertices.size();
		V.reserve(n);

		// Gets a list of backprojected vertices
		while (true) {
			CGAL_Point_2 M = v_curr->pt(t);
			V.push_back(SP->backproject(M));
			Polygon_Edge* e = (v_curr->e1 == e_prev ? v_curr->e2 : v_curr->e1);
			Polygon_Vertex* v_next = (e->v1 == v_curr ? e->v2 : e->v1);
			if (v_next == v_init) {
				break;
			} else {
				e_prev = e;
				v_curr = v_next;
			}
		}

		// Builds facets

		P.clear();
		P.push_back(V[0]);
		bool forward = (CGAL::cross_product(V[1] - V[0], V[2] - V[0]) * h > 0);
		if (forward) {
			for (size_t i = 1; i < n; i++) P.push_back(V[i]);
		} else {
			for (size_t i = n - 1; i > 0; i--) P.push_back(V[i]);
		}
	}


	void Polygon_Set::get_polygon_description(Polygon_Vertex_Octree* V, std::list<std::list<int> > & P, const FT & t)
	{
		double prec = 1.0 / (1 << 30) / (1 << 10);
		FT::set_relative_precision_of_to_double(prec);

		Support_Plane* SP = Universe::map_of_planes[id_plane];

		// Get the description of all polygons at time t.

		for (auto it_n = nodes.begin() ; it_n != nodes.end() ; ++it_n) {
			Polygon_Node* N = it_n->second;

			if (N->contours_are_empty()) {
			
				for (std::list<Polygon*>::const_iterator it_p = N->polygons_begin(); it_p != N->polygons_end(); it_p++) {
					Polygon* polygon = (*it_p);

					if (polygon->vertices.empty()) continue;

					std::list<CGAL_Point_3> V_node;
					get_sequence_of_3d_vertices(polygon, t, V_node);
					
					std::list<int> P_node;
					for (std::list<CGAL_Point_3>::const_iterator it_v = V_node.begin() ; it_v != V_node.end() ; ++it_v) {
						const CGAL_Point_3 & v = (*it_v);

						double vx = CGAL::to_double(v.x()), vy = CGAL::to_double(v.y()), vz = CGAL::to_double(v.z());
						CGAL_Inexact_Point_3 hint_v (vx, vy, vz);
						P_node.push_back(V->get_identifier(v, hint_v));
					}

					P.push_back(P_node);
				}

			} else {
				
				std::vector<CGAL_Point_3> V_node;
				size_t n = N->contours_size();
				V_node.reserve(n);

				for (size_t i = 0 ; i < n ; ++i) {
					size_t i_next = (i == n - 1 ? 0 : i + 1);
					Constraint C_curr = N->get_contour(i);
					Constraint C_next = N->get_contour(i_next);
					CGAL_Point_2 M = SP->get_intersection_point(C_curr.first, C_next.first);
					V_node.push_back(SP->backproject(M));
				}

				const CGAL_Plane & H = SP->plane;
				const CGAL_Vector_3 & h = H.orthogonal_vector();

				bool forward = (CGAL::cross_product(V_node[1] - V_node[0], V_node[2] - V_node[0]) * h > 0);
				int i_min, i_max, i_incr;
				if (forward) {
					i_min = 0, i_max = n, i_incr = 1;
				} else {
					i_min = n - 1, i_max = -1, i_incr = -1;
				}

				std::list<int> P_node;
				int i = i_min;
				while (i != i_max) {
					const CGAL_Point_3 & v = V_node[i];

					double vx = CGAL::to_double(v.x()), vy = CGAL::to_double(v.y()), vz = CGAL::to_double(v.z());
					CGAL_Inexact_Point_3 hint_v (vx, vy, vz);
					P_node.push_back(V->get_identifier(v, hint_v));

					i += i_incr;
				}

				P.push_back(P_node);
			}
		}
	}






	void Polygon_Set::build_graph()
	{
		// We build a graph that connects adjacent polygons together.

		// Vertices of this graph are already known : these are the different Polygon_Nodes.
		// Here, two Polygon_Nodes are linked iff there is no Polygon_Segment between them.
		// Also, we only consider couples (n1, n2) where n1->signature < n2->signature.

		for (auto it_n1 = nodes.begin() ; it_n1 != nodes.end() ; ++it_n1) {
			Polygon_Node* n1 = it_n1->second;

			size_t p = n1->contours_size();
			for (size_t i = 0 ; i < p ; ++i) {
				Constraint C = n1->get_contour(i);
				Intersection_Line* I = C.first;

				if (C.second == MINUS) {
					// Does a polygon exist on the other side of I_curr ?
					Signature S = get_adjacent_polygons_signature(n1->get_one_polygon(), I);
					auto it_n2 = nodes.find(S);
					if (it_n2 == nodes.end()) continue;

					// Gets the two constraints surrounding C_curr
					size_t i_prev = (i == 0 ? p - 1 : i - 1);
					size_t i_next = (i == p - 1 ? 0 : i + 1);
					Constraint C_prev = n1->get_contour(i_prev);
					Constraint C_next = n1->get_contour(i_next);

					// Tests the existence of a segment
					// If there is no segment, connects nodes n1 and n2
					bool nodes_connected = !I->exists_segment_adjacent_to_edge(C_prev, C, C_next);

					if (nodes_connected) {
						Polygon_Node* n2 = it_n2->second;
						n1->set_contour_as_not_on_border(i);
						n2->set_contour_as_not_on_border(I);
						edges.push_back(new Polygon_Link(n1, n2));
					}
				}
			}
		}
	}


	void Polygon_Set::get_groups(std::list<Polygon_Group*> & groups)
	{
		for (auto it_n = nodes.begin() ; it_n != nodes.end() ; ++it_n) {
			Polygon_Node* n_init = it_n->second;
			if (n_init->is_visited()) continue;

			// If the node hasn't been visited before, then we initialize a queue.
			// The queue will be progressively filled with neighbors of the current node,
			// which haven't been visited as well.

			std::list<Polygon_Node*> G;
			std::queue<Polygon_Node*> Q;
			Q.push(n_init);
			n_init->set_queued();

			while (!Q.empty()) {
				Polygon_Node* n = Q.front();
				Q.pop();
				G.push_back(n);

				for (std::list<Polygon_Link*>::const_iterator it_l = n->links_begin() ; it_l != n->links_end() ; ++it_l) {
					Polygon_Node* m = (*it_l)->other_node(n);
					if (!m->is_queued()) {
						Q.push(m);
						m->set_queued();
					}
				}

				n->set_visited();
			}

			// Appends the group G
			groups.push_back(new Polygon_Group(G));
		}
	}


	void Polygon_Set::clear_graph()
	{
		// Once all groups of connected nodes have been retrieved,
		// we can clear the structure ogf graph formed by the elements of the set.

		for (std::list<Polygon_Link*>::iterator it_e = edges.begin() ; it_e != edges.end() ; ++it_e) {
			delete (*it_e);
		}
		edges.clear();

		for (auto it_n = nodes.begin() ; it_n != nodes.end() ; ++it_n) {
			it_n->second->clear_links();
		}
	}
}