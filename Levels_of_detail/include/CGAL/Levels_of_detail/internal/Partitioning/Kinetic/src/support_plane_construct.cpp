#include "../include/support_plane.h"


namespace Skippy {
	std::vector<int> Support_Plane::box_edges_to_facets(const int i)
	{
		int _facets[12][2] = { {0, 4}, {1, 4}, {0, 5}, {1, 5}, {2, 4}, {3, 4}, {2, 5}, {3, 5}, {0, 2}, {0, 3}, {1, 2}, {1, 3} };
		std::vector<int> facets = { _facets[i][0], _facets[i][1] };
		return facets;
	}



	std::vector<int> Support_Plane::box_vertices_to_facets(const int i)
	{
		int _vertices_to_facets[8][3] = { {0, 2, 4}, {1, 2, 4}, {0, 3, 4}, {1, 3, 4}, {0, 2, 5}, {1, 2, 5}, {0, 3, 5}, {1, 3, 5} };
		std::vector<int> facets = { _vertices_to_facets[i][0], _vertices_to_facets[i][1], _vertices_to_facets[i][2] };
		return facets;
	}



	std::vector<int> Support_Plane::box_facets_to_vertices(const int i)
	{
		int _vertices[6][4] = { {0, 4, 6, 2}, {1, 5, 7, 3}, {1, 5, 4, 0}, {3, 7, 6, 2}, {1, 0, 2, 3}, {5, 4, 6, 7} };
		std::vector<int> vertices = { _vertices[i][0], _vertices[i][1], _vertices[i][2], _vertices[i][3] };
		return vertices;
	}



	std::vector<int> Support_Plane::box_facets_to_edges(const int i)
	{
		int _edges[6][4] = { {0, 8, 2, 9}, {1, 10, 3, 11}, {10, 6, 8, 4}, {11, 7, 9, 5}, {4, 0, 5, 1}, {6, 2, 7, 3} };
		std::vector<int> edges = { _edges[i][0], _edges[i][1], _edges[i][2], _edges[i][3] };
		return edges;
	}



	CGAL_Point_3 Support_Plane::intersection_plane_and_edge_of_bounding_box(const CGAL_Point_3 & pt_min, const CGAL_Point_3 & pt_max, const CGAL_Plane & P, const int i)
	{
		FT x, y, z;
		const FT & a = P.a(), &b = P.b(), &c = P.c(), &d = P.d();

		switch (i) {
		case 0: x = pt_min.x(), z = pt_min.z(); y = -(a * x + c * z + d) / b; break;
		case 1: x = pt_max.x(), z = pt_min.z(); y = -(a * x + c * z + d) / b; break;
		case 2: x = pt_min.x(), z = pt_max.z(); y = -(a * x + c * z + d) / b; break;
		case 3: x = pt_max.x(), z = pt_max.z(); y = -(a * x + c * z + d) / b; break;
		case 4: y = pt_min.y(), z = pt_min.z(); x = -(b * y + c * z + d) / a; break;
		case 5: y = pt_max.y(), z = pt_min.z(); x = -(b * y + c * z + d) / a; break;
		case 6: y = pt_min.y(), z = pt_max.z(); x = -(b * y + c * z + d) / a; break;
		case 7: y = pt_max.y(), z = pt_max.z(); x = -(b * y + c * z + d) / a; break;
		case 8: x = pt_min.x(), y = pt_min.y(); z = -(a * x + b * y + d) / c; break;
		case 9: x = pt_min.x(), y = pt_max.y(); z = -(a * x + b * y + d) / c; break;
		case 10: x = pt_max.x(), y = pt_min.y(); z = -(a * x + b * y + d) / c; break;
		case 11: x = pt_max.x(), y = pt_max.y(); z = -(a * x + b * y + d) / c; break;
		}

		return CGAL_Point_3(x, y, z);
	}



	bool Support_Plane::find_next_object_colliding_plane(const CGAL_Point_3 & pt_min,
		const CGAL_Point_3 & pt_max,
		const std::vector<CGAL_Point_3> & box_corners,
		const std::vector<std::pair<int, int> > & box_edges,
		const CGAL_Plane & P,
		const std::vector<int> & V_ind,
		const std::vector<int> & E_ind,
		const std::pair<bool, int> & previous_object,
		std::pair<bool, int> & next_object,
		CGAL_Point_3 & M)
	{
		// Given a list of vertices and edges, corresponding either to the entire bounding box
		// or to only one of its facets, we determine the next object to be crossed by P.
		// First option : one of the vertices belongs to P.
		for (int i = 0; i < V_ind.size(); i++) {
			int v_i = V_ind[i];
			if ((previous_object.first && previous_object.second != v_i) || (!previous_object.first)) {
				const CGAL_Point_3 & v = box_corners[v_i];
				if (P.a() * v.x() + P.b() * v.y() + P.c() * v.z() + P.d() == 0) {
					next_object = std::make_pair(true, v_i);
					M = v;
					return true;
				}
			}
		}

		// Second option : one edge of the facet is crossed by P
		for (int i = 0; i < E_ind.size(); i++) {
			int e_i = E_ind[i];
			if (previous_object.first || (!previous_object.first && previous_object.second != e_i)) {
				const CGAL_Point_3 & s = box_corners[box_edges[e_i].first];
				const CGAL_Point_3 & t = box_corners[box_edges[e_i].second];

				if ((P.a() * s.x() + P.b() * s.y() + P.c() * s.z() + P.d()) * (P.a() * t.x() + P.b() * t.y() + P.c() * t.z() + P.d()) < 0) {
					next_object = std::make_pair(false, e_i);

					// Determines the coordinates of the intersection point
					FT x, y, z;
					if (e_i <= 3) {
						x = box_corners[box_edges[e_i].first].x(), z = box_corners[box_edges[e_i].first].z();
						y = -(P.a() * x + P.c() * z + P.d()) / P.b();
					} else if (e_i <= 7) {
						y = box_corners[box_edges[e_i].first].y(), z = box_corners[box_edges[e_i].first].z();
						x = -(P.b() * y + P.c() * z + P.d()) / P.a();
					} else {
						x = box_corners[box_edges[e_i].first].x(), y = box_corners[box_edges[e_i].first].y();
						z = -(P.a() * x + P.b() * y + P.d()) / P.c();
					}
					M = CGAL_Point_3(x, y, z);
					return true;
				}
			}
		}

		return false;
	}



	void Support_Plane::find_next_object_colliding_plane(const CGAL_Point_3 & pt_min,
		const CGAL_Point_3 & pt_max,
		const std::vector<CGAL_Point_3> & box_corners,
		const std::vector<std::pair<int, int> > & box_edges,
		const CGAL_Plane & P,
		std::pair<bool, int> & next_object,
		CGAL_Point_3 & M)
	{
		std::vector<int> V_ind(8, -1), E_ind(12, -1);
		for (int i = 0; i < V_ind.size(); i++) V_ind[i] = i;
		for (int i = 0; i < E_ind.size(); i++) E_ind[i] = i;
		std::pair<bool, int> previous_object(false, -1);

		find_next_object_colliding_plane(pt_min, pt_max, box_corners, box_edges, P, V_ind, E_ind, previous_object, next_object, M);
	}



	void Support_Plane::construct_bounding_polygon_of_support_plane(const CGAL_Point_3 & pt_min,
		const CGAL_Point_3 & pt_max,
		const std::vector<CGAL_Point_3> & box_corners,
		const std::vector<std::pair<int, int> > & box_edges,
		const CGAL_Plane & P,
		std::list<CGAL_Point_3> & bounding_polygon)
	{
		std::vector<std::list<int> > bounding_facets;
		construct_bounding_polygon_of_support_plane(pt_min, pt_max, box_corners, box_edges, P, bounding_polygon, bounding_facets);
	}



	void Support_Plane::construct_bounding_polygon_of_support_plane(const CGAL_Point_3 & pt_min,
		const CGAL_Point_3 & pt_max,
		const std::vector<CGAL_Point_3> & box_corners,
		const std::vector<std::pair<int, int> > & box_edges,
		const CGAL_Plane & P,
		std::list<CGAL_Point_3> & bounding_polygon,
		std::vector<std::list<int> > & bounding_facets)
	{
		bounding_polygon.clear();
		bounding_facets.clear();

		// Initialization : we search for an object in B_vertices belonging to P.
		// If such an object cannot be found, we identify an edge crossed by P.
		// Each object is identify by a boolean (true iif vertex) and an index.

		CGAL_Point_3 M;
		std::pair<bool, int> init_object(true, -1), prev_object, curr_object;
		find_next_object_colliding_plane(pt_min, pt_max, box_corners, box_edges, P, init_object, M);
		bounding_polygon.push_back(M);

		// Iterates over the edges and faces of the bounding box, until reaching the initial point again
		prev_object = init_object;
		int prev_facet = -1;

		do {
			// Given a current object, we find the list of facets that are adjacent to this object
			std::vector<int> adjacent_facets;
			if (prev_object.first) {
				adjacent_facets = box_vertices_to_facets(prev_object.second);
			} else {
				adjacent_facets = box_edges_to_facets(prev_object.second);
			}

			// Now, we loop on the list of adjacent facets to find the next object.
			// For each facet, we identify the next vertex or edge of the bounding box crossed by P.
			// We must remember to avoid the facet indexed by 'previous_facet'.
			bool iteration_done = false;
			int current_facet;
			for (int f = 0; f < adjacent_facets.size(); f++) {
				current_facet = adjacent_facets[f];
				if (current_facet == prev_facet) continue;

				std::vector<int> vertices = box_facets_to_vertices(current_facet);
				std::vector<int> edges = box_facets_to_edges(current_facet);

				iteration_done = find_next_object_colliding_plane(pt_min, pt_max, box_corners, box_edges, P, vertices, edges, prev_object, curr_object, M);
				if (iteration_done) {
					// The new intersection point between P and the box should always be added, excepted when we come back to the start
					if (curr_object != init_object) {
						bounding_polygon.push_back(M);
					}

					// The facet(s) crossed for going from prev_object to curr_object should always be added
					if (curr_object.first) {
						std::list<int> facets;
						for (int g = 0; g < adjacent_facets.size(); g++) {
							if (adjacent_facets[g] != prev_facet) facets.push_back(adjacent_facets[g]);
						}
						bounding_facets.push_back(facets);
					} else {
						bounding_facets.push_back(std::list<int>(1, current_facet));
					}

					prev_object = curr_object;
					prev_facet = current_facet;
					break;
				}
			}

			// Constructs a new vertex and adds it to the polygon
			assert(iteration_done);

		} while (curr_object != init_object);
	}
}