#include "regularization_ordinates_quadratic.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include "trace.h"
#include "geometry.h"


typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K>    Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                    Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                      Delaunay;
typedef Delaunay::Point                                             Point;


Regularization_Ordinates_Quadratic::Regularization_Ordinates_Quadratic()
{
}


Regularization_Ordinates_Quadratic::~Regularization_Ordinates_Quadratic()
{
}


void Regularization_Ordinates_Quadratic::regularize(Kinetic_Model *model)
{
	clock_t t_begin = clock();
	Segment_Regularization_Tree* tree = model->tree;

	std::vector<double> dt_max;
	vector<Triplet<double> > v_mu, v_targets;
	set_bounds(model->segments, dt_max, model->params);

	map<double, Node_Parallel_Segments*>::iterator it_c = tree->parallel_segments.begin();
	while (it_c != tree->parallel_segments.end()) {
	
		// Accesses the cluster
		double theta = it_c->first;
		Node_Parallel_Segments* cluster = it_c->second;
		cluster->delete_colinear_nodes();

		if (cluster->parallel_segments.size() == 0) {
			it_c = tree->parallel_segments.erase(it_c);
			continue;
		}

		// Transforms coordinates to a rotated frame
		double x_min, x_max, y_min, y_max;
		transform_coordinates(theta, cluster, x_min, x_max, y_min, y_max);

		Boost_RTree rtree_segments;
		Regularization_Ordinates::build_rtree(cluster, rtree_segments);

		build_neighbors_graph(model, cluster, model->segments, rtree_segments, model->params->regp_quad_lambda, dt_max, v_mu, v_targets);
		it_c++;
	}

	// Initializes and solves a quadratic problem
    SparseMat<double> mu, targets;
    int individuals = int(model->segments.size()), variables = int(model->segments.size() + v_mu.size());
    build_sparse_matrices(individuals, v_mu, v_targets, mu, targets);

	vector<double> dt;
	vector<double> shift(individuals, 0);
	Quadratic_Problem* Q = new Quadratic_Problem (individuals, variables, model->params->regp_quad_lambda, dt_max, mu, targets, shift);
	if (Q->solve(dt)) {
		build_regularization_tree(model->segments, targets, dt, model->tree);

		for (map<double, Node_Parallel_Segments*>::iterator it_c = tree->parallel_segments.begin(); it_c != tree->parallel_segments.end() ; it_c++) {
			translate(it_c->second);
			it_c->second->parallel_segments.clear();
		}

		model->applied_regularization_ordinates = true;
#if NOT_MEASURING_PERFORMANCES
		make_layer(model);
#endif
	}
	delete Q;

	clock_t t_end = clock();
	std::cout << "** Regularization (part 2) done in " + std::to_string(float(t_end - t_begin) / CLOCKS_PER_SEC) + " s." << std::endl;
#if NOT_MEASURING_PERFORMANCES
	for (auto it_m = tree->parallel_segments.begin() ; it_m != tree->parallel_segments.end() ; it_m++) {
		Node_Parallel_Segments* node = it_m->second;
		
		for (auto it_n = node->colinear_segments.begin() ; it_n != node->colinear_segments.end() ; it_n++) {
			auto it_n_next = it_n;
			++it_n_next;
			if (it_n_next != node->colinear_segments.end()) {
				assert(fabs(it_n->first - it_n_next->first) > 1);
			}

			// assert(!it_n->second->colinear_segments.empty());
			if (!it_n->second->colinear_segments.empty()) {
				Segment* s = it_n->second->colinear_segments.front();
				double a = s->a, b = s->b, c = s->c;
			}
		}
	}
#endif
}


void Regularization_Ordinates_Quadratic::set_bounds(vector<Segment *> & segments, vector<double> & dt_max, Parameters* params)
{
	dt_max = vector<double>(segments.size(), 0);

	double (Regularization_Ordinates::*dt_function)(Segment*, Parameters*) = nullptr;
	switch (params->regp_trans_function) {
	case 0: dt_function = &Regularization_Ordinates::dt_max_const; break;
	case 1: dt_function = &Regularization_Ordinates::dt_max_lsd_width; break;
	}

	// For each segment we set its dt_max according to the selected function
	for (uint i = 0 ; i < segments.size() ; i++) {
		dt_max[i] = (this->*dt_function)(segments[i], params);
	}
}


void Regularization_Ordinates_Quadratic::build_neighbors_graph(Kinetic_Model* model, Node_Parallel_Segments* node, vector<Segment *> & all_segments, Boost_RTree & rtree_segments, double lambda, vector<double> & dt_max, vector<Triplet<double> > &mu, vector<Triplet<double> > &targets)
{
	vector<pair<Point, uint> > points;
	map<uint, Segment*> points_to_segments;

	if (model->params->rega_quad_discretize) {
		uint j = 0;
		double ds = model->params->rega_quad_discretization_step;
		for (list<Segment *>::iterator it_s = node->parallel_segments.begin() ; it_s != node->parallel_segments.end() ; it_s++) {
			Segment* s = (*it_s);
			Point2d & s_end1 = s->interEnd1;
			Point2d & s_end2 = s->interEnd2;
			Point2d & s_bar = s->barycenter;
			Vec2d s_dir = cv::normalize(Vec2d(s_end2 - s_end1));

			int subdivs = int(floor(s->length / (2 * ds)));
			for (int k = 0 ; k < subdivs ; k++) {
				points.push_back(std::make_pair(Point(s_end1.x + k * s_dir[0], s_end1.y + k * s_dir[1]), j));
				points_to_segments[j] = s;
				++j;

				points.push_back(std::make_pair(Point(s_end2.x - k * s_dir[0], s_end2.y - k * s_dir[1]), j));
				points_to_segments[j] = s;
				++j;
			}

			points.push_back(std::make_pair(Point(s_bar.x, s_bar.y), j));
			points_to_segments[j] = s;
			++j;
		}
	} else {
		uint j = 0;
		for (list<Segment *>::iterator it_s = node->parallel_segments.begin() ; it_s != node->parallel_segments.end() ; it_s++) {
			Segment* s = (*it_s);
			Point2d & s_bar = s->barycenter;
			points.push_back(std::make_pair(Point(s_bar.x, s_bar.y), j));
			points_to_segments[j] = s;
			++j;
		}
	}

	Delaunay DT;
	DT.insert(points.begin(), points.end());
	set<pair<Segment*, Segment*> > considered_potentials;

	for (Delaunay::Finite_edges_iterator it_e = DT.finite_edges_begin() ; it_e != DT.finite_edges_end() ; it_e++) {
		Delaunay::Edge e = *it_e;
		uint e_i = e.first->vertex((e.second + 1) % 3)->info();
		uint e_j = e.first->vertex((e.second + 2) % 3)->info();
		Segment* s_i = points_to_segments[e_i];
		Segment* s_j = points_to_segments[e_j];
		uint i = s_i->index;
		uint j = s_j->index;

		if (i == j) continue;

		pair<Segment*, Segment*> p_ij = (i < j ? std::make_pair(s_i, s_j) : std::make_pair(s_j, s_i));
		if (considered_potentials.find(p_ij) != considered_potentials.end()) continue;
		considered_potentials.insert(p_ij);

		double y_ij = s_i->referencing_coordinates.y - s_j->referencing_coordinates.y;
		double mu_ij = lambda;

		if (fabs(y_ij) < dt_max[i] + dt_max[j]) {
			mu.push_back(Triplet<double>(i, j, mu_ij));
			targets.push_back(Triplet<double>(i, j, y_ij));
		}
	}

/*
	for (list<Segment *>::iterator it_s1 = node->parallel_segments.begin() ; it_s1 != node->parallel_segments.end() ; it_s1++) {
		Segment* s1 = (*it_s1);
		
		int i = s1->index;
		list<Segment *> neighbors;
		Regularization_Ordinates::search_neighborhood(s1, 1000, 10, rtree_segments, all_segments, neighbors);

		for (list<Segment *>::iterator it_s2 = neighbors.begin() ; it_s2 != neighbors.end() ; it_s2++) {
			Segment* s2 = (*it_s2);
			if (s1 == s2 || s1->index > s2->index) continue;

			assert(s1->node_parallel == s2->node_parallel);
			int j = s2->index;
			double y_ij = s1->referencing_coordinates.y - s2->referencing_coordinates.y;
			double mu_ij = lambda;

			if (fabs(y_ij) < dt_max[i] + dt_max[j]) {
				mu.push_back(Triplet<double>(i, j, mu_ij));
				targets.push_back(Triplet<double>(i, j, y_ij));
			}
		}
	}
	*/
}


void Regularization_Ordinates_Quadratic::build_sparse_matrices(int individuals, vector<Triplet<double> > & v_mu, vector<Triplet<double> > & v_targets, SparseMat<double> & mu, SparseMat<double> & targets)
{
    mu = SparseMat<double>(individuals, individuals);
    mu.setFromTriplets(v_mu.begin(), v_mu.end());
    targets = SparseMat<double>(individuals, individuals);
    targets.setFromTriplets(v_targets.begin(), v_targets.end());
}


void Regularization_Ordinates_Quadratic::build_regularization_tree(vector<Segment *> & segments, SparseMat<double> & targets, vector<double> & dt, Segment_Regularization_Tree* & tree)
{
	int n = int(segments.size());
	vector<int> segments_to_groups(n, -1);
	map<int, list<int> > groups_to_segments = map<int, list<int> >();
	map<Node_Parallel_Segments*, list<int> > nodes_to_groups = map<Node_Parallel_Segments*, list<int> >();
	int g = 0;

	int p = 0;
	for (int k = 0 ; k < targets.outerSize() ; k++) {
		for (SparseMat<double>::InnerIterator it_1(targets, k) ; it_1 ; ++it_1) {
			int i = it_1.row(), j = it_1.col();
			if (fabs(dt[n + p]) < 1e-6) {
				// Then segments i and j belong to the same group of parallel segments.
				// For the moment, these groups are materialized by integers.
				if (segments_to_groups[i] == -1 && segments_to_groups[j] == -1) {
                    // Segments i and j are not assigned to any group of parallel segments
                    // So we create one with them
                    segments_to_groups[i] = segments_to_groups[j] = g;
                    groups_to_segments[g].push_back(i);
                    groups_to_segments[g].push_back(j);
					nodes_to_groups[segments[i]->node_parallel].push_back(g);
                    g++;
                } else if (segments_to_groups[i] == -1 && segments_to_groups[j] != -1) {
                    // Assigns segment i to the group of the segment j
                    int g_j = segments_to_groups[j];
                    segments_to_groups[i] = g_j;
                    groups_to_segments[g_j].push_back(i);
                } else if (segments_to_groups[i] != -1 && segments_to_groups[j] == -1) {
                    // Assigns segment j to the group of the segment i
                    int g_i = segments_to_groups[i];
                    segments_to_groups[j] = g_i;
                    groups_to_segments[g_i].push_back(j);
                } else {
                    int g_i = segments_to_groups[i];
                    int g_j = segments_to_groups[j];
                    if (g_i != g_j) {
                        // Segments i and j have been assigned to different groups, but in fact
                        // they belong to the same group. That's why we merge them.
                        for (list<int>::iterator it_2 = groups_to_segments[g_j].begin() ; it_2 != groups_to_segments[g_j].end() ; it_2++) {
                            segments_to_groups[*it_2] = g_i;
                            groups_to_segments[g_i].push_back(*it_2);
                        }
                        groups_to_segments[g_j].clear();
						// Deletes entry g_j from 'nodes_to_groups'
						list<int>::iterator it_n = nodes_to_groups[segments[i]->node_parallel].begin();
						while (it_n != nodes_to_groups[segments[i]->node_parallel].end()) {
							if ((*it_n) == g_j) {
								nodes_to_groups[segments[i]->node_parallel].erase(it_n);
								break;
							}
							++it_n;
						}
                    }
                }
			}
			++p;
		}
	}

	// We prepare the construction of the regularization tree
	double y_eps = 1;
	map<int, double> ordinates;

	for (uint i = 0; i < segments_to_groups.size() ; i++) {
		int g_i = segments_to_groups[i];
		Node_Parallel_Segments* node_i = segments[i]->node_parallel;

		if (g_i != -1) {
			if (ordinates.find(g_i) == ordinates.end()) {
				double y = segments[i]->referencing_coordinates.y + dt[i];

				// Checks if this ordinates seems to be associated to another group of segments
				int g_j = -1;
				for (map<int, double>::iterator it_m = ordinates.begin() ; it_m != ordinates.end() ; it_m++) {
					if (fabs(it_m->second - y) < y_eps) {
						// We found a value close to it_m, but does it correspond to the same group of parallel segments ?
						list<int>::iterator it_n = nodes_to_groups[node_i].begin();
						while (it_n != nodes_to_groups[node_i].end()) {
							if ((*it_n) == it_m->first) break;
							++it_n;
						}
						if (it_n != nodes_to_groups[node_i].end()) {
							g_j = it_m->first;
						}
					}
				}

				if (g_j == -1) {
					ordinates[g_i] = y;
				} else {
					// Merge groups
					for (list<int>::iterator it = groups_to_segments[g_i].begin() ; it != groups_to_segments[g_i].end() ; it++) {
						segments_to_groups[*it] = g_j;
						groups_to_segments[g_j].push_back(*it);
					}
					groups_to_segments[g_i].clear();
				}
			}
		}
	}

	// Tries to assign segments whose orientation has not been optimized thanks to the regularization process, to an existing group
	for (uint i = 0 ; i < segments_to_groups.size() ; i++) {
		int g_i = segments_to_groups[i];
		Node_Parallel_Segments* node_i = segments[i]->node_parallel;
		if (g_i == -1) {
			double y = segments[i]->referencing_coordinates.y;

			int g_j = -1;
			for (map<int, double>::iterator it_m = ordinates.begin(); it_m != ordinates.end(); it_m++) {
				double y_j = it_m->second;
				if (fabs(y_j - y) < y_eps) {
					// We found a value close to it_m, but does it correspond to the same group of parallel segments ?
					list<int>::iterator it_n = nodes_to_groups[node_i].begin();
					while (it_n != nodes_to_groups[node_i].end()) {
						if ((*it_n) == it_m->first) break;
						++it_n;
					}
					if (it_n != nodes_to_groups[node_i].end()) {
						g_j = it_m->first;
					}
				}
				if (g_j != -1) break;
			}

			if (g_j == -1) {
				g_i = ordinates.rbegin()->first + 1;
				ordinates[g_i] = y;
				nodes_to_groups[node_i].push_back(g_i);
			} else {
				g_i = g_j;
			}
			segments_to_groups[i] = g_i;
			groups_to_segments[g_i].push_back(i);
		}
	}

	// Finally builds the regularization tree.
	for (map<Node_Parallel_Segments*, list<int> >::iterator it_m = nodes_to_groups.begin() ; it_m != nodes_to_groups.end() ; it_m++) {
		Node_Parallel_Segments* node = it_m->first;
		list<int> & groups = it_m->second;
		for (list<int>::iterator it_n = groups.begin() ; it_n != groups.end() ; it_n++) {
			double y = ordinates[*it_n];
			node->create_colinear_node(y);
		}
	}

	// Assigns segments
	for (uint i = 0 ; i < segments_to_groups.size() ; i++) {
		Segment* s_i = segments[i];
		int g_i = segments_to_groups[i];
		s_i->node_parallel->assign_to_colinear_node(ordinates[g_i], s_i);
	}
}