#include "../include/partition_multiple.h"
#include "../include/partition_objects.h"
#include "../include/vars.h"
#include "../include/stats.h"


namespace Skippy
{
	Partition_Multiple::Partition_Multiple(const std::vector<CGAL_Plane> & plane_defs)
		: Partition(plane_defs)
	{
	}


	Partition_Multiple::~Partition_Multiple()
	{
	}


	void Partition_Multiple::incremental_build()
	{
		clock_t t_0 = clock();

		// Everytime a sub-volume of the global bounding box is processed,
		// we construct new elements that we add to the partition.

		// We build the facets related to the inner support planes.
		// We add vertices and edges associated to the local bounding box,
		// but do not build facets.

		// Also, we remove bivalent vertices to save memory.

		std::list<Partition_Vertex*> V;
		std::list<Partition_Edge*> E;

		build_inner_facets(V, E);
		
		build_missing_vertices_for_bounding_box(V);
		build_missing_edges_for_bounding_box(E);

		Partition::remove_bivalent_vertices(V, E);

		std::copy(E.begin(), E.end(), std::back_inserter(edges));

		KP_Stats::partition_time += double(clock() - t_0) / CLOCKS_PER_SEC;
	}


	void Partition_Multiple::finalize(const int slices)
	{
		clock_t t_0 = clock();

		// When the kinetic propagation has been performed on each of the N subvolumes,
		// we must finalize the graph by cleaning edges belonging to slicing planes.

		rebuild_edges_at_boundaries(slices);

		// We are now read to build facets on the slicing planes,
		// but ids of the edges should be contiguous, that's why
		// a preprocessing step is applied.

		make_contiguous_indices(slices);
		build_facets_for_bounding_box();
		build_facets_on_slicing_planes(slices);
		reset_indices_of_facets();

		// We build the final polyhedrons.

		Partition::remove_bivalent_vertices();
		build_polyhedrons();

		KP_Stats::final_nb_vertices = Counters::id_partition_vertex;
		KP_Stats::final_nb_edges = Counters::id_partition_edge;
		KP_Stats::final_nb_facets = Counters::id_partition_facet;
		KP_Stats::final_nb_polyhedrons = Counters::id_partition_polyhedron;

		KP_Stats::partition_time += double(clock() - t_0) / CLOCKS_PER_SEC;
		std::cout << "** Built a partition with " << polyhedrons.size() << " polyhedrons" << std::endl;
	}



	void Partition_Multiple::make_contiguous_indices(const int slices)
	{
		int n = int(planes.size());
		Counters::par_v_local_ids = std::vector<int>(n, -1);
		Counters::par_e_local_ids = std::vector<int>(n, -1);

		std::list<Partition_Vertex*> V;
		octree_vertices->get_all_vertices(V);

		for (std::list<Partition_Vertex*>::iterator it_v = V.begin() ; it_v != V.end() ; ++it_v) {
			Partition_Vertex* v = (*it_v);
			for (std::map<int, int>::iterator it_l = v->local_ids_begin() ; it_l != v->local_ids_end() ; ++it_l) {
				int p = it_l->first;
				it_l->second = ++Counters::par_v_local_ids[p];
			}
		}

		for (std::list<Partition_Edge*>::iterator it_e = edges.begin() ; it_e != edges.end() ; ++it_e) {
			Partition_Edge* e = (*it_e);
			for (std::map<int, int>::iterator it_l = e->local_ids_begin() ; it_l != e->local_ids_end() ; ++it_l) {
				int p = it_l->first;
				it_l->second = ++Counters::par_e_local_ids[p];
			}
		}
	}



	void Partition_Multiple::rebuild_edges_at_boundaries(const int slices)
	{
		int n = int(planes.size());

		// After performing a kinetic propagation in two subvolumes separated by a slicing plane,
		// and because some planes can be taken into account in one subvolume but not in the other,
		// we observe the two following phenomena in the (kind of) planar graph (V, E) associated
		// to the boundary plane of index 'id' :
		// - some paths might have been duplicated, and therefore should be simplified,
		// - some edges may cross each other and need to be shortened (for slicing planes only).

		for (int p = 0 ; p < 6 ; ++p) {
			rebuild_edges(p, planes[p], false);
		}

		for (int p = n - slices ; p < n ; ++p) {
			rebuild_edges(p, planes[p], true);
		}

		/*for (int p = n - slices ; p < n ; ++p) {
			check(p, planes[p]);
		}*/
	}


	void Partition_Multiple::rebuild_edges(const int id, const CGAL_Plane & P, bool edges_may_intersect)
	{
		// Step 1.
		// Gets all the vertices of the plane P,
		// and iterators pointing to all edges included in the plane P.
		// Then, for each vertex, we check if it is used to delimit two collinear edges,
		// included in pid and that go to the same direction.

		std::list<Partition_Vertex*> V;
		int P_or = get_orientation_of_normal_vector(P);
		octree_vertices->get_vertices_for_boundary_plane(id, P_or, P, V);

		std::map<int, std::list<Partition_Edge*>::iterator> E_map;
		for (std::list<Partition_Edge*>::iterator it_e = edges.begin() ; it_e != edges.end() ; ++it_e) {
			if ((*it_e)->get_local_id(id) != -1) E_map[(*it_e)->id] = it_e;
		}

		for (std::list<Partition_Vertex*>::iterator it_v = V.begin() ; it_v != V.end() ; ++it_v) {
			Partition_Vertex* v = (*it_v);

			std::list<std::pair<Partition_Edge*, Partition_Edge*> > E_dup;
			if (v->delimits_duplicate_paths(id, E_dup)) {

				// If the current vertex has at least one pair of duplicate edges,
				// then we exhibit the full paths and merge them.
				for (auto it_p = E_dup.begin() ; it_p != E_dup.end() ; ++it_p) {
					Partition_Edge* e1 = it_p->first;
					Partition_Edge* e2 = it_p->second;
					exhibit_and_merge_duplicate_paths(id, v, e1, e2, E_map);
				}
			}
		}

		// Step 2.
		// Detects self-crossing edges and replaces them.
		
		remove_overlapping_edges(id, P, E_map);

		if (edges_may_intersect) {
			remove_self_intersecting_edges(P, V, E_map);
		}
	}



	void Partition_Multiple::check(const int id, const CGAL_Plane & P)
	{
		std::cout << "Checks plane [" << id << "]" << std::endl;

		std::list<Partition_Edge*> E;
		std::map<int, CGAL_Segment_2> S;

		for (std::list<Partition_Edge*>::const_iterator it_e = edges.begin() ; it_e != edges.end() ; ++it_e) {
			Partition_Edge* e = (*it_e);
			if (e->get_local_id(id) != -1) {
				E.push_back(e);
				CGAL_Point_2 A = P.to_2d(e->source(true)->M);
				CGAL_Point_2 B = P.to_2d(e->target(true)->M);
				S[e->id] = CGAL_Segment_2(A, B);
			}
		}

		std::list<Partition_Edge*>::const_iterator it_e1, it_e2;
		for (it_e1 = E.begin() ; it_e1 != E.end() ; ++it_e1) {
			Partition_Edge* e1 = (*it_e1);
			const CGAL_Segment_2 & S_1 = S[e1->id];

			it_e2 = it_e1;
			++it_e2;
			while (it_e2 != E.end()) {
				Partition_Edge* e2 = (*it_e2);
				const CGAL_Segment_2 & S_2 = S[e2->id];

				CGAL::cpp11::result_of<K::Intersect_2(CGAL_Segment_2, CGAL_Segment_2)>::type object = CGAL::intersection(S_1, S_2);
				if (object) {
					if (const CGAL_Segment_2* ptr = boost::get<CGAL_Segment_2>(&*object)) {
						std::cout << "Intersection between edges : " << e1->id << ", " << e2->id << std::endl;
					}
				}

				++it_e2;
			}
		}
	}


	void Partition_Multiple::remove_overlapping_edges(const int id, const CGAL_Plane & P,
		std::map<int, std::list<Partition_Edge*>::iterator> & E_map)
	{
		bool verbose = false; // (id == 2163);

		// Step 1.
		// We loop on all edges of E_map, that we split into N sublists of edges,
		// depending on the identity of the other planes to which these edges belong.

		std::map<int, std::list<int> > E_map_sublists;
		for (auto it_e = E_map.begin(); it_e != E_map.end(); ++it_e) {
			Partition_Edge* e = (*it_e->second);

			// We suppose that each edge belongs to only one other plane
			if (e->local_ids_size() != 2) {
				throw std::logic_error("Error in Partition_Multiple::remove_overlapping_edges : unhandled case.");
			}

			for (auto it_l = e->local_ids_begin() ; it_l != e->local_ids_end() ; ++it_l) {
				if (it_l->first == id) {
					continue;
				} else {
					E_map_sublists[it_l->first].push_back(e->id);
				}
			}
		}
		
		// Step 2.
		// We loop on all edges of all edges.
		// If we detect a couple of overlapping edges (e1 e2), then e1 and e2 are shortened and replaced.

		for (std::map<int, std::list<int> >::iterator it_p = E_map_sublists.begin() ; it_p != E_map_sublists.end() ; ++it_p) {
			int id_p = it_p->first;
			if (verbose) {
				std::cout << "Line : " << id_p << std::endl;
			}

			// Step 2.1
			// Gets a 3D line, that we are going to use
			// in order to avoid the use of CGAL::intersection for CGAL_Segment_3 objects which are parallel

			CGAL_Line_3 D_p;
			bool D_p_is_set = false;
			const CGAL_Plane & H_p = planes[id_p];
			CGAL::cpp11::result_of<K::Intersect_3(CGAL_Plane, CGAL_Plane)>::type object = CGAL::intersection(P, H_p);
			if (object) {
				if (const CGAL_Line_3* ptr = boost::get<CGAL_Line_3>(&*object)) {
					D_p = *ptr;
					D_p_is_set = true;
				}
			}
			assert(D_p_is_set);

			const CGAL_Point_3 D_pt = D_p.point();
			const CGAL_Vector_3 D_dir = D_p.to_vector();

			int used_coordinate = -1;
			if (D_dir.x() != 0) {
				used_coordinate = 0;
			} else if (D_dir.y() != 0) {
				used_coordinate = 1;
			} else {
				used_coordinate = 2;
			}

			// Step 2.2
			// Loops on a list of identifiers

			std::list<int> & L_p = it_p->second;
			std::list<int>::iterator it_e1 = L_p.begin(), it_e2;

			while (it_e1 != L_p.end()) {
				int i1 = (*it_e1);
				Partition_Edge* e1 = *E_map[i1];
				Partition_Vertex *v11 = e1->source(true), *v12 = e1->target(true);
				const CGAL_Point_3 & M11 = v11->M, &M12 = v12->M;	
				bool intersection_exists = false;

				// Computes abscissas for M11 and M12 along the line D_p
				FT t_11 = 0, t_12 = 0;
				switch (used_coordinate) {
				case 0: t_11 = (M11.x() - D_pt.x()); t_12 = (M12.x() - D_pt.x()); break;
				case 1: t_11 = (M11.y() - D_pt.y()); t_12 = (M12.y() - D_pt.y()); break;
				case 2: t_11 = (M11.z() - D_pt.z()); t_12 = (M12.z() - D_pt.z()); break;
				}
				FT s1_min = (t_11 < t_12 ? t_11 : t_12), s1_max = (t_11 > t_12 ? t_11 : t_12);

				it_e2 = it_e1;
				++it_e2;

				// Step 2.2.1
				// Tests if e1 and e2 intersect
				
				while (it_e2 != L_p.end()) {
					int i2 = (*it_e2);
					Partition_Edge* e2 = *E_map[i2];
					Partition_Vertex *v21 = e2->source(true), *v22 = e2->target(true);
					const CGAL_Point_3 & M21 = v21->M, &M22 = v22->M;

					// Computes abscissas again
					FT t_21 = 0, t_22 = 0;
					switch (used_coordinate) {
					case 0: t_21 = (M21.x() - D_pt.x()); t_22 = (M22.x() - D_pt.x()); break;
					case 1: t_21 = (M21.y() - D_pt.y()); t_22 = (M22.y() - D_pt.y()); break;
					case 2: t_21 = (M21.z() - D_pt.z()); t_22 = (M22.z() - D_pt.z()); break;
					}
					FT s2_min = (t_21 < t_22 ? t_21 : t_22), s2_max = (t_21 > t_22 ? t_21 : t_22);

					// Tests intersection
					// If e1 and e2 don't intersect, fine, we continue
					// Otherwise we fill a vector of abscissas, which will later tell us
					// in which order e1 and e2 must be refined.
					intersection_exists = jinr(s1_min, s2_min, s1_max) || jinr(s1_min, s2_max, s1_max)
						|| jinr(s2_min, s1_min, s2_max) || jinr(s2_min, s1_max, s2_max)
						|| (s1_min == s2_min) || (s1_max == s2_max);

					if (!intersection_exists) {
						++it_e2;
					} else {
						if (verbose) {
							std::cout << "-- Edges : " << e1->id << " " << e2->id << std::endl;
							std::cout << "e1 = (" << v11->id << " " << v12->id << ")" << std::endl;
							std::cout << "     (" << M11 << " :: " << M12 << ")" << std::endl;
							std::cout << "     (" << t_11 << " " << t_12 << ")" << std::endl;
							std::cout << "e2 = (" << v21->id << " " << v22->id << ")" << std::endl;
							std::cout << "     (" << M21 << " :: " << M22 << ")" << std::endl;
							std::cout << "     (" << t_21 << " " << t_22 << ")" << std::endl;
						}

						// Step 2.2.2
						// Replacement of e1 and e2 by smaller edges.

						std::vector<std::tuple<FT, Partition_Vertex*, int> > V_abscissas;
						V_abscissas.push_back(std::make_tuple(t_11, v11, 1));
						V_abscissas.push_back(std::make_tuple(t_12, v12, 1));
						V_abscissas.push_back(std::make_tuple(t_21, v21, 2));
						V_abscissas.push_back(std::make_tuple(t_22, v22, 2));
						std::sort(V_abscissas.begin(), V_abscissas.end());

						if (verbose) {
							std::cout << "-- Sorted points : " << std::endl;
							for (size_t i = 0 ; i < V_abscissas.size() ; ++i) {
								std::cout << "(" << std::get<1>(V_abscissas[i]) << ", " << std::get<2>(V_abscissas[i]) << ")" << " ";
							}
							std::cout << std::endl;
						}

						// Creates three subedges

						Partition_Vertex *v_0 = std::get<1>(V_abscissas[0]), *v_1 = std::get<1>(V_abscissas[1]),
							*v_2 = std::get<1>(V_abscissas[2]), *v_3 = std::get<1>(V_abscissas[3]);

						Partition_Edge *e_01, *e_12, *e_23;
						e_01 = e_12 = e_23 = nullptr;
						if (v_0 != v_1) e_01 = new Partition_Edge(v_0, v_1);
						if (v_1 != v_2) e_12 = new Partition_Edge(v_1, v_2);
						if (v_2 != v_3) e_23 = new Partition_Edge(v_2, v_3);

						// Substitutes e1 and e2 in facets that contain these edges.
						
						std::list<Partition_Edge*> S_1, S_2;
						int a_0 = std::get<2>(V_abscissas[0]), a_1 = std::get<2>(V_abscissas[1]),
							a_2 = std::get<2>(V_abscissas[2]), a_3 = std::get<2>(V_abscissas[3]);
						if (a_0 == 1 && a_1 == 2 && a_2 == 1 && a_3 == 2) {
							// e_1 and e_2 intersect partially
							S_1.push_back(e_01); S_1.push_back(e_12);
							S_2.push_back(e_12); S_2.push_back(e_23);
						} else if (a_0 == 1 && a_1 == 2 && a_2 == 2 && a_3 == 1) {
							// e_2 is included in e_1
							S_1.push_back(e_01); S_1.push_back(e_12); S_1.push_back(e_23);
							S_2.push_back(e_12);
						} else if (a_0 == 2 && a_1 == 1 && a_2 == 1 && a_3 == 2) {
							// e_1 is included in e_2
							S_1.push_back(e_12);
							S_2.push_back(e_01); S_2.push_back(e_12); S_2.push_back(e_23);
						} else if (a_0 == 2 && a_1 == 1 && a_2 == 2 && a_3 == 1) {
							// e_1 and e_2 intersect partially
							S_1.push_back(e_12); S_1.push_back(e_23);
							S_2.push_back(e_01); S_2.push_back(e_12);
						} else {
							throw std::logic_error("Error in remove_overlapped_edges : unhandled overlap case.");
						}
						e1->substitute_in_facets(S_1);
						e2->substitute_in_facets(S_2);

						// Deletes e1 and e2
						std::map<int, std::list<Partition_Edge*>::iterator>::iterator it_r1 = E_map.find(e1->id), it_r2 = E_map.find(e2->id);
						edges.erase(it_r1->second);
						edges.erase(it_r2->second);
						E_map.erase(it_r1);
						E_map.erase(it_r2);
						delete e1;
						delete e2;

						// Inserts the three edges in the data structure
						Partition_Edge* edges_created[3] = { e_01, e_12, e_23 };
						for (int i = 0 ; i < 3 ; ++i) {
							Partition_Edge* f = edges_created[i];
							if (f == nullptr) continue;

							if (verbose) {
								std::cout << "-- Creates : " << f->id << " (" << f->source(true)->id << " " << f->target(true)->id << ")" << std::endl;
							}

							L_p.push_back(f->id);
							edges.push_back(f);
							E_map[f->id] = --edges.end();
						}

						// Finally exists the loop with it_e2 still pointing to e2 (deleted)
						break;
					}
				}

				if (intersection_exists) {
					L_p.erase(it_e2);
					it_e1 = L_p.erase(it_e1);
				} else {
					++it_e1;
				}
			}
		}
	}


	void Partition_Multiple::exhibit_and_merge_duplicate_paths(const int id, Partition_Vertex* v, Partition_Edge* e1, Partition_Edge* e2,
		std::map<int, std::list<Partition_Edge*>::iterator> & E_map)
	{
		// Provided a vertex v that delimits two duplicate paths, whose initial edges are e1 and e2,
		// we exhibit the full duplicate paths and merge them into one.

		// Step 1.
		// Exhibits the duplicate paths.

		std::vector<Partition_Edge*> E_1, E_2;
		std::vector<Partition_Vertex*> V_1, V_2;
		exhibit_duplicate_path(v, e1, E_1, V_1);
		exhibit_duplicate_path(v, e2, E_2, V_2);

		// Step 2.
		// Finds a common element to V_1 and V_2,
		// which is used to shorten lists E_1 and E_2.

		std::vector<std::tuple<Partition_Vertex*, int, FT> > V;
		Partition_Vertex* v_dest = shorten_duplicate_paths(E_1, V_1, E_2, V_2);
		sort_vertices_by_distance(v, v_dest, V_1, V_2, V);

		// Step 3.
		// Simultaneously processes elements of V_1 and V_2 to build edges between elements of V_1 and V_2.

		merge_duplicate_paths(v, E_1, E_2, V, E_map);
	}


	void Partition_Multiple::exhibit_duplicate_path(Partition_Vertex* v_init, Partition_Edge* e_init, 
		std::vector<Partition_Edge*> & E, std::vector<Partition_Vertex*> & V)
	{
		std::list<Partition_Edge*> E_0;
		std::list<Partition_Vertex*> V_0;

		Partition_Vertex* v_curr = v_init;
		Partition_Edge* e_curr = e_init, *e_next = nullptr;

		E_0.push_back(e_curr);

		do {
			Partition_Vertex* v_next = e_curr->second_vertex(v_curr);
			V_0.push_back(v_next);

			e_next = v_next->prolongate(e_curr);
			if (e_next != nullptr) {
				v_curr = v_next;
				e_curr = e_next;
				E_0.push_back(e_curr);
			}
		} while (e_next != nullptr);

		E = std::vector<Partition_Edge*>(E_0.begin(), E_0.end());
		V = std::vector<Partition_Vertex*>(V_0.begin(), V_0.end());
	}



	Partition_Vertex* Partition_Multiple::shorten_duplicate_paths(std::vector<Partition_Edge*> & E_1, std::vector<Partition_Vertex*> & V_1,
		std::vector<Partition_Edge*> & E_2, std::vector<Partition_Vertex*> & V_2)
	{
		size_t i_1 = 0, i_2 = 0;

		// Loops on elements of V_1 and V_2 till finding the common element

		Partition_Vertex* v = nullptr;

		for (i_1 = 0; i_1 < V_1.size(); ++i_1) {
			for (i_2 = 0 ; i_2 < V_2.size(); ++i_2) {
				if (V_1[i_1] == V_2[i_2]) {
					v = V_1[i_1];
					goto v_anchor;
				}
			}
		}

		if (i_1 == V_1.size()) return nullptr;

	v_anchor:
		// We have v == V_1[i_1] == V_2[i_2].
		// We remove the elements of E_1, E_2, V_1 and V_2 that correspond to elements after v.
		// Note that v itself is removed from V_1 and V_2.

		V_1.resize(i_1); E_1.resize(i_1 + 1);
		V_2.resize(i_2); E_2.resize(i_2 + 1);

		return v;
	}


	void Partition_Multiple::sort_vertices_by_distance(Partition_Vertex* v_source, Partition_Vertex* v_dest,
		std::vector<Partition_Vertex*> & V_1, std::vector<Partition_Vertex*> & V_2,
		std::vector<std::tuple<Partition_Vertex*, int, FT> > & V)
	{
		// Between v_source and v_dest, which may be nullptr,
		// there are all the elements of V_1 and V_2.
		// We insert them into V which is sorted by distance to v_source.

		size_t n_1 = V_1.size(), n_2 = V_2.size();
		V.reserve(n_1 + n_2 + 1);

		// We insert v_source as front element.
		// The second element of the tuple is 0 : it is on both paths.

		for (size_t i = 0 ; i < n_1 ; ++i) {
			Partition_Vertex* v = V_1[i];
			FT d = (V_1[i]->M - v_source->M).squared_length();
			V.push_back(std::make_tuple(v, 1, d));
		}

		for (size_t i = 0 ; i < n_2 ; ++i) {
			Partition_Vertex* v = V_2[i];
			FT d = (V_2[i]->M - v_source->M).squared_length();
			V.push_back(std::make_tuple(v, 2, d));
		}

		struct _Tuple_Comparator {
			bool operator() (const std::tuple<Partition_Vertex*, int, FT> & L, const std::tuple<Partition_Vertex*, int, FT> & R) {
				return std::get<2>(L) < std::get<2>(R);
			}
		} Tuple_Comparator;

		std::sort(V.begin(), V.end(), Tuple_Comparator);

		// If there exists a final vertex shared by E_1 and E_2,
		// we can insert it into V and its path id is the id of the vertex just before it.
		// This will avoid to recreate a vertex that already exists.
		if (v_dest != nullptr) {
			FT d = (v_dest->M - v_source->M).squared_length();
			int path_id = (V.empty() ? 1 : std::get<1>(V.back()));
			V.push_back(std::make_tuple(v_dest, path_id, d));
		}
	}


	void Partition_Multiple::merge_duplicate_paths(Partition_Vertex* v_source,
		std::vector<Partition_Edge*> & E_1, std::vector<Partition_Edge*> & E_2,
		std::vector<std::tuple<Partition_Vertex*, int, FT> > & V,
		std::map<int, std::list<Partition_Edge*>::iterator> & E_map)
	{
		// Starting from a vertex v_source,
		// we loop on all the vertices of V, sorted by increasing distance.
		
		// Declares some variables :
		// - indices of the edges currently processed on E_1 and E_2,
		// - the currently followed path,
		// - the edge of E_1 and E_2 that's going to be removed,
		// - the sequence of edges substituted to the removed edge.

		int i_1 = 0, i_2 = 0, path;
		Partition_Edge* R = nullptr;
		std::list<Partition_Edge*> S;

		// Initialization : v_source is on both paths E_1 and E_2, 
		// That's why we select the path (1 or 2) initially followed.

		if (std::get<1>(V[0]) == 1) {
			R = E_2[i_2];
			++i_2;
			path = 1;
		} else {
			R = E_1[i_1];
			++i_1;
			path = 2;
		}

		// We loop on each couple of successive vertices V[j - 1] and V[j]
		// and check if we need to switch between paths.
		
		for (size_t j = 0 ; j < V.size() ; ++j) {

			int path_next = std::get<1>(V[j]);
			if (path == path_next) {
				// We stay on the same path.
				// A edge is added to the list S, intended to replace R.
				// Finally, we iterate on E_1 or E_2.
				if (path == 1) {
					S.push_back(E_1[i_1]);
					++i_1;
				} else if (path == 2) {
					S.push_back(E_2[i_2]);
					++i_2;
				}
			}

			else {
				// We switch between paths.
				// First we create the edge V[j - 1] V[j] which doesn't exist.

				Partition_Edge* e = new Partition_Edge(j == 0 ? v_source : std::get<0>(V[j - 1]), std::get<0>(V[j]));
				
				edges.push_back(e);
				E_map[e->id] = --edges.end();

				// With this edge, we get the final sequence of edges S intended to replace R.
				// We remove any reference to R in the list of facets that contain it, and replace R by S.
				// After that, we definitely erase R.

				S.push_back(e);
				R->substitute_in_facets(S);

				std::map<int, std::list<Partition_Edge*>::iterator>::iterator it = E_map.find(R->id);
				edges.erase(it->second);
				E_map.erase(it);
				delete R;
				R = nullptr;

				// Now, we switch from one path to the other.
				// The next edge of the other path, which couldn't be selected, is going to be removed.
				// We reinitilize a list of substitute edges.

				S.clear();
				S.push_back(e);
				path = path_next;

				if (path == 1) {
					if (i_2 != E_2.size()) R = E_2[i_2];
					++i_2;
				} else {
					if (i_1 != E_1.size()) R = E_1[i_1];
					++i_1;
				}

				if (R == nullptr) return;
			}
		}

		// If v_dest != nullptr, the loop ends and one duplicate edge R is not removed.
		// We do it now.

		if (R != nullptr) {
			R->substitute_in_facets(S);

			std::map<int, std::list<Partition_Edge*>::iterator>::iterator it = E_map.find(R->id);
			edges.erase(it->second);
			E_map.erase(it);
			delete R;
		}
	}



	void Partition_Multiple::get_bounding_box_for_rtree(Partition_Edge* e, const std::map<int, CGAL_Inexact_Point_2> & hints_points_2d,
		const double & eps, Boost_Box & B)
	{
		Partition_Vertex *v_s = e->source(true), *v_t = e->target(true);
		const std::map<int, CGAL_Inexact_Point_2>::const_iterator it_s = hints_points_2d.find(v_s->id), it_t = hints_points_2d.find(v_t->id);
		const CGAL_Inexact_Point_2 & M_s = it_s->second, &M_t = it_t->second;

		const double &x_s = M_s.x(), &y_s = M_s.y(), &x_t = M_t.x(), &y_t = M_t.y();
		double x_min, x_max, y_min, y_max;
		if (x_s < x_t) {
			x_min = x_s - eps, x_max = x_t + eps;
		} else {
			x_min = x_t - eps, x_max = x_s + eps;
		}
		if (y_s < y_t) {
			y_min = y_s - eps, y_max = y_t + eps;
		} else {
			y_min = y_t - eps, y_max = y_s + eps;
		}

		B = Boost_Box(Boost_Point(x_min, y_min), Boost_Point(x_max, y_max));
	}



	void Partition_Multiple::remove_self_intersecting_edges(const CGAL_Plane & P, std::list<Partition_Vertex*> & V,
		std::map<int, std::list<Partition_Edge*>::iterator> & E_map)
	{
		// Part 1.
		// Initiliazes a map associating Partition_Vertices to 2D representations on the slicing plane P.
		
		double prec = 1.0 / (1 << 30) / (1 << 10);
		
		std::map<int, CGAL_Point_2> points_2d;
		std::map<int, CGAL_Inexact_Point_2> hints_points_2d;
		for (std::list<Partition_Vertex*>::const_iterator it_v = V.begin() ; it_v != V.end() ; ++it_v) {
			Partition_Vertex* v = (*it_v);

			CGAL_Point_2 M = P.to_2d(v->M);

			FT::set_relative_precision_of_to_double(prec);
			double x = CGAL::to_double(M.x()), y = CGAL::to_double(M.y());
			CGAL_Inexact_Point_2 hint_M(x, y);

			points_2d[v->id] = M;
			hints_points_2d[v->id] = hint_M;
		}

		// Part 2.
		// Initializes a R-tree with the list of edges of E_map, which are included in P.

		double eps = 1e-6;
		Boost_RTree rtree_edges;

		for (std::map<int, std::list<Partition_Edge*>::iterator>::iterator it_e = E_map.begin() ; it_e != E_map.end() ; ++it_e) {
			Partition_Edge* e = (*it_e->second);
			Boost_Box B;

			get_bounding_box_for_rtree(e, hints_points_2d, eps, B);
			rtree_edges.insert(std::make_pair(B, e->id));
		}

		// Part 3.
		// Loops on the list of edges and performs queries in a R-tree to find a list of edges which potentially intersect.
		// When a box B_ref intersects another box B, then we determine, using the exact kernel, if edges e1 and e2 intersect in their interiors.
		// If so we split them into smaller ones, and they are inserted at the end of E.

		std::map<int, std::list<Partition_Edge*>::iterator>::iterator it_e1 = E_map.begin();

		while (it_e1 != E_map.end()) {
			Partition_Edge* e1 = (*it_e1->second);
			Boost_Box query;

			// Gets definition and support line of e1
			Partition_Vertex *v1_1 = e1->source(true), *v1_2 = e1->target(true);
			const CGAL_Point_2 & A = points_2d[v1_1->id], & B = points_2d[v1_2->id];
			CGAL_Line_2 L_1(A, B);

			bool non_vertical_line_1 = true;
			FT xa = A.x(), xb = B.x(), ya = 0, yb = 0;
			if (xa == xb) {
				ya = A.y(), yb = B.y();
				non_vertical_line_1 = false;
			}

			// R-tree query
			get_bounding_box_for_rtree(e1, hints_points_2d, eps, query);
			std::vector<Boost_Value> candidates;
			rtree_edges.query(bgi::intersects(query), std::back_inserter(candidates));

			// Loops on the results of the query and searches for an intersection
			bool intersection_detected = false;
			std::map<int, std::list<Partition_Edge*>::iterator>::iterator it_e2 = E_map.end();

			size_t u = 0;
			while (u < candidates.size()) {

				// Gets current edge
				// Avoid e1 or an edge that has already been removed
				int e2_id = candidates[u].second;
				it_e2 = E_map.find(e2_id);
				if (it_e2 == it_e1 || it_e2 == E_map.end()) {
					++u;
					continue;
				}

				Partition_Edge* e2 = (*it_e2->second);

				// Gets definition and support line of e2
				Partition_Vertex *v2_1 = e2->source(true), *v2_2 = e2->target(true);
				const CGAL_Point_2 & C = points_2d[v2_1->id], & D = points_2d[v2_2->id];
				CGAL_Line_2 L_2(C, D);

				bool non_vertical_line_2 = true;
				FT xc = C.x(), xd = D.x(), yc = 0, yd = 0;
				if (xc == xd) {
					yc = C.y(), yd = D.y();
					non_vertical_line_2 = false;
				}

				// Computes intersection, if it exists.
				CGAL_Point_2 M;
				CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object = CGAL::intersection(L_1, L_2);
				if (object) {
					if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object)) {
						M = *ptr;
						bool intersection_detected_1, intersection_detected_2 = false;

						if (non_vertical_line_1) {
							FT x = M.x();
							intersection_detected_1 = (xa < x && x < xb) || (xb < x && x < xa);
						} else {
							FT y = M.y();
							intersection_detected_1 = (ya < y && y < yb) || (yb < y && y < ya);
						}

						if (intersection_detected_1) {
							if (non_vertical_line_2) {
								FT x = M.x();
								intersection_detected_2 = (xc < x && x < xd) || (xd < x && x < xc);
							} else {
								FT y = M.y();
								intersection_detected_2 = (yc < y && y < yd) || (yd < y && y < yc);
							}
						}

						intersection_detected = intersection_detected_1 && intersection_detected_2;
					}
				}
				
				if (!intersection_detected) {
					++u;
				} 
				
				else {
					// Creates a Partition_Vertex representing the intersection of e1 and e2.
					// We don't forget to update V and points_2d.

					std::list<int> M_def;
					Partition_Edge::definition_at_intersection(e1, e2, M_def);
					Partition_Vertex* v = new Partition_Vertex(P.to_3d(M), M_def);
					octree_vertices->add(v);

					V.push_back(v);
					points_2d[v->id] = M;

					FT::set_relative_precision_of_to_double(prec);
					double x = CGAL::to_double(M.x()), y = CGAL::to_double(M.y());
					CGAL_Inexact_Point_2 hint_M(x, y);
					hints_points_2d[v->id] = hint_M;

					// Creates four subedges

					Partition_Edge* e1_1 = new Partition_Edge(v1_1, v);
					Partition_Edge* e1_2 = new Partition_Edge(v1_2, v);
					Partition_Edge* e2_1 = new Partition_Edge(v2_1, v);
					Partition_Edge* e2_2 = new Partition_Edge(v2_2, v);

					// Substitutes e1 by (e1_1, e1_2) and e2 by (e2_1, e2_2) 
					// in facets which contain these edges.

					std::list<Partition_Edge*> S_1, S_2;
					S_1.push_back(e1_1); S_1.push_back(e1_2);
					S_2.push_back(e2_1); S_2.push_back(e2_2);
					e1->substitute_in_facets(S_1);
					e2->substitute_in_facets(S_2);

					// Deletes e1 and e2
					edges.erase(it_e1->second);
					edges.erase(it_e2->second);
					delete e1;
					delete e2;

					// Inserts the four subedges in the data structure
					Partition_Edge* edges_created[4] = { e1_1, e1_2, e2_1, e2_2 };
					for (int i = 0; i < 4; ++i) {
						Partition_Edge* f = edges_created[i];
						edges.push_back(f);
						E_map[f->id] = --edges.end();
						
						Boost_Box B_f;
						get_bounding_box_for_rtree(f, hints_points_2d, eps, B_f);
						rtree_edges.insert(std::make_pair(B_f, f->id));
					}

					// Finally breaks the loop with it_e2 still pointing to e2 (deleted)
					break;
				}
			}

			if (intersection_detected) {
				E_map.erase(it_e2);
				it_e1 = E_map.erase(it_e1);
			} else {
				++it_e1;
			}
		}
	}



#if 0
	void Partition_Multiple::remove_self_intersecting_edges(const CGAL_Plane & P, std::list<Partition_Vertex*> & V,
		std::map<int, std::list<Partition_Edge*>::iterator> & E_map)
	{
		// Initiliazes a map associating Partition_Vertices to 2D representations on the slicing plane P,
		// and a list of edges belonging to this plane.
		
		std::map<int, CGAL_Point_2> points_2d;
		for (std::list<Partition_Vertex*>::const_iterator it_v = V.begin() ; it_v != V.end() ; ++it_v) {
			Partition_Vertex* v = (*it_v);
			points_2d[v->id] = P.to_2d(v->M);
		}

		std::list<Partition_Edge*> E;
		for (std::map<int, std::list<Partition_Edge*>::iterator>::iterator it_e = E_map.begin() ; it_e != E_map.end() ; ++it_e) {
			E.push_back((*it_e->second));
		}

		// Double loop on the list of edges.
		// If we find a couple of edges e1 and e2 which intersect in their interior,
		// then we split them into smaller ones, and they are inserted at the end of E.

		std::list<Partition_Edge*>::iterator it_e1 = E.begin(), it_e2;

		while (it_e1 != E.end()) {
			Partition_Edge* e1 = (*it_e1);
			Partition_Vertex *v1_1 = e1->source(true), *v1_2 = e1->target(true);
			const CGAL_Point_2 & A = points_2d[v1_1->id], & B = points_2d[v1_2->id];
			CGAL_Line_2 L_1(A, B);

			bool intersection_detected = false;

			bool non_vertical_line_1 = true;
			FT xa = A.x(), xb = B.x(), ya = 0, yb = 0;
			if (xa == xb) {
				ya = A.y(), yb = B.y();
				non_vertical_line_1 = false;
			}

			it_e2 = it_e1;
			++it_e2;

			while (it_e2 != E.end()) {
				Partition_Edge* e2 = (*it_e2);
				Partition_Vertex *v2_1 = e2->source(true), *v2_2 = e2->target(true);
				const CGAL_Point_2 & C = points_2d[v2_1->id], & D = points_2d[v2_2->id];
				CGAL_Line_2 L_2(C, D);
				CGAL_Point_2 M;

				bool non_vertical_line_2 = true;
				FT xc = C.x(), xd = D.x(), yc = 0, yd = 0;
				if (xc == xd) {
					yc = C.y(), yd = D.y();
					non_vertical_line_2 = false;
				}

				CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object = CGAL::intersection(L_1, L_2);
				if (object) {
					if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object)) {
						M = *ptr;

						bool intersection_detected_1, intersection_detected_2;

						if (non_vertical_line_1) {
							FT x = M.x();
							intersection_detected_1 = (xa < x && x < xb) || (xb < x && x < xa);
						} else {
							FT y = M.y();
							intersection_detected_1 = (ya < y && y < yb) || (yb < y && y < ya);
						}

						if (non_vertical_line_2) {
							FT x = M.x();
							intersection_detected_2 = (xc < x && x < xd) || (xd < x && x < xc);
						} else {
							FT y = M.y();
							intersection_detected_2 = (yc < y && y < yd) || (yd < y && y < yc);
						}

						intersection_detected = intersection_detected_1 && intersection_detected_2;
					}
				}

				if (!intersection_detected) {
					++it_e2;
				} else {
					
					// Creates a Partition_Vertex representing the intersection of e1 and e2.
					// We don't forget to update V and points_2d.

					std::list<int> M_def;
					Partition_Edge::definition_at_intersection(e1, e2, M_def);
					Partition_Vertex* v = new Partition_Vertex(P.to_3d(M), M_def);
					octree_vertices->add(v);

					V.push_back(v);
					points_2d[v->id] = M;

					// Creates four subedges

					Partition_Edge* e1_1 = new Partition_Edge(v1_1, v);
					Partition_Edge* e1_2 = new Partition_Edge(v1_2, v);
					Partition_Edge* e2_1 = new Partition_Edge(v2_1, v);
					Partition_Edge* e2_2 = new Partition_Edge(v2_2, v);

					// Substitutes e1 by (e1_1, e1_2) and e2 by (e2_1, e2_2) 
					// in facets which contain these edges.

					std::list<Partition_Edge*> S_1, S_2;
					S_1.push_back(e1_1); S_1.push_back(e1_2);
					S_2.push_back(e2_1); S_2.push_back(e2_2);
					e1->substitute_in_facets(S_1);
					e2->substitute_in_facets(S_2);

					// Deletes e1 and e2

					std::map<int, std::list<Partition_Edge*>::iterator>::iterator it_r1 = E_map.find(e1->id), it_r2 = E_map.find(e2->id);
					edges.erase(it_r1->second);
					edges.erase(it_r2->second);
					E_map.erase(it_r1);
					E_map.erase(it_r2);
					delete e1;
					delete e2;

					// Inserts the four subedges in the data structure

					Partition_Edge* edges_created[4] = { e1_1, e1_2, e2_1, e2_2 };
					for (int i = 0 ; i < 4 ; ++i) {
						Partition_Edge* f = edges_created[i];
						E.push_back(f);
						edges.push_back(f);
						E_map[f->id] = --edges.end();
					}

					// Finally breaks the loop with it_e2 still pointing to e2 (deleted)
					break;
				}
			}

			if (intersection_detected) {
				E.erase(it_e2);
				it_e1 = E.erase(it_e1);
			} else {
				++it_e1;
			}
		}
	}
#endif


	void Partition_Multiple::build_facets_on_slicing_planes(const int slices)
	{
		int n = int(planes.size());

		for (int p = n - slices ; p < n ; ++p) {
			build_facets_on_boundaries(p, planes[p]);
		}
	}


	void Partition_Multiple::reset_indices_of_facets()
	{
		Counters::id_partition_facet = -1;

		for (int range = 0 ; range < 2 ; ++range) {
			size_t p_min = (range == 0 ? 6 : 0);
			size_t p_max = (range == 0 ? facets.size() : 6);

			for (size_t p = p_min; p < p_max; ++p) {
				for (std::list<Partition_Facet*>::const_iterator it_f = facets[p].begin(); it_f != facets[p].end(); ++it_f) {
					Partition_Facet* f = (*it_f);
					f->id = ++Counters::id_partition_facet;
				}
			}
		}
	}
}