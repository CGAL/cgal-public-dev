#include "regularization_angles_quadratic.h"
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


Regularization_Angles_Quadratic::Regularization_Angles_Quadratic()
    : Regularization_Angles()
{

}



Regularization_Angles_Quadratic::~Regularization_Angles_Quadratic()
{

}



void Regularization_Angles_Quadratic::regularize(Kinetic_Model* model)
{
	clock_t t_begin = clock();

    bool include_artificial = true;

    // Cleans the regularization tree, if another execution has been perform before
    model->tree->delete_parallel_nodes();

    // Builds a R-tree from the list of segments
    Boost_RTree rtree_segments;
    Regularization_Angles::build_rtree(model->segments, rtree_segments, include_artificial);

    // We set the angle by which each segment can be rotated
    std::vector<double> dalpha_max;
    set_bounds(model->segments, dalpha_max, model->params);

    // We build our graph of neighbors
    vector<Triplet<double> > v_mu, v_targets;
	vector<Triplet<int> > v_relations;
    build_neighbors_graph(model, model->segments, rtree_segments, model->params->rega_quad_distance, include_artificial, model->params->rega_quad_lambda, dalpha_max, v_mu, v_targets, v_relations);

    // Initializes and solves a quadratic problem
    SparseMat<double> mu, targets;
	SparseMat<int> relations;
    int individuals = int(model->segments.size()), variables = int(model->segments.size() + v_mu.size());
    build_sparse_matrices(individuals, v_mu, v_targets, v_relations, mu, targets, relations);
	
	// std::cout << individuals << " " << variables << std::endl;

    vector<double> dalpha;
	vector<double> shift(individuals, 0);
    Quadratic_Problem* Q = new Quadratic_Problem (individuals, variables, model->params->rega_quad_lambda, dalpha_max, mu, targets, shift);
    if (Q->solve(dalpha)) {

		/*
		for (size_t u = 0 ; u < dalpha.size() ; u++) {
			std::cout << dalpha[u] << " ";
			if (u % 10 == 0) std::cout << std::endl;
		} */

		// std::cout << "** rotating ..." << std::endl;
        // Rotates segments
        build_regularization_tree(model->segments, targets, relations, dalpha, model->tree, model->params->rega_epsilon);
        rotate(model->tree);

		model->rega_quad_potentials = variables - individuals;
        model->applied_regularization_angles = true;

#if NOT_MEASURING_PERFORMANCES
        make_layer(model);
#endif
    }
    delete Q;

	clock_t t_end = clock();
	std::cout << "** Regularization (part 1) done in " + std::to_string(float(t_end - t_begin) / CLOCKS_PER_SEC) + " s." << std::endl;
}



void Regularization_Angles_Quadratic::set_bounds(vector<Segment *> & segments, vector<double> & dalpha_max, Parameters *params)
{
    dalpha_max = vector<double>(segments.size(), 0);

    double (Regularization_Angles::*dalpha_function)(Segment*, Parameters*) = NULL;
    switch (params->rega_angle_function) {
	case 0: dalpha_function = &Regularization_Angles::dalpha_max_const; break;
	case 1: dalpha_function = &Regularization_Angles::dalpha_max_offset; break;
	case 2: dalpha_function = &Regularization_Angles::dalpha_max_lsd; break;
	}

    // For each segment we set its dalpha_max in a way that limits
    // the shift of both segments' ends along its initial normal axis,
    // and the longer a segment, the shorter the shift.
    for (uint i = 0 ; i < segments.size() ; i++) {
        dalpha_max[i] = (this->*dalpha_function)(segments[i], params);
    }
}



void Regularization_Angles_Quadratic::build_neighbors_graph(Kinetic_Model* model, vector<Segment *> &segments, Boost_RTree &rtree_segments, double D, bool include_artificial, double lambda, vector<double> & dalpha_max, 
	vector<Triplet<double> > &mu, vector<Triplet<double> > &targets, vector<Triplet<int> > & relations)
{
#if NOT_MEASURING_PERFORMANCES
	model->clear_line_items(model->L_ag);
#endif
	if (model->params->rega_quad_graph == 0) {
		build_neighbors_graph_euclidian(model, segments, rtree_segments, D, include_artificial, lambda, dalpha_max, mu, targets, relations);
	} else if (model->params->rega_quad_graph == 1) {
		build_neighbors_graph_delaunay(model, segments, rtree_segments, D, include_artificial, lambda, dalpha_max, mu, targets, relations);
	}
}



void Regularization_Angles_Quadratic::build_neighbors_graph_delaunay(Kinetic_Model* model, vector<Segment *> &segments, Boost_RTree &rtree_segments, double D, bool include_artificial, double lambda, vector<double> & dalpha_max,
	vector<Triplet<double> > &mu, vector<Triplet<double> > &targets, vector<Triplet<int> > & relations)
{
	vector<pair<Point, uint> > points;
	map<uint, uint> points_to_segments;

	if (model->params->rega_quad_discretize) {
		uint j = 0;
		double ds = model->params->rega_quad_discretization_step;
		for (uint i = 0 ; i < segments.size(); i++) {
			Segment* s_i = segments[i];

			Point2d & s_i_end1 = s_i->end1;
			Point2d & s_i_end2 = s_i->end2;
			Point2d & s_i_bar = s_i->barycenter;
			Vec2d & s_i_dir = s_i->direction;
			Vec2d dir = cv::normalize(s_i_dir);

			int subdivs = int(floor(s_i->length / (2 * ds)));
			for (int k = 0 ; k < subdivs ; k++) {
				points.push_back(std::make_pair(Point(s_i_end1.x + k * dir[0], s_i_end1.y + k * dir[1]), j));
				points_to_segments[j] = i;
				++j;

				points.push_back(std::make_pair(Point(s_i_end2.x - k * dir[0], s_i_end2.y - k * dir[1]), j));
				points_to_segments[j] = i;
				++j;
			}

			points.push_back(std::make_pair(Point(s_i_bar.x, s_i_bar.y), j));
			points_to_segments[j] = i;
			++j;
		}
	} else {
		for (uint i = 0; i < segments.size(); i++) {
			Point2d & bar_i = segments[i]->barycenter;
			points.push_back(std::make_pair(Point(bar_i.x, bar_i.y), i));
			points_to_segments[i] = i;
		}
	}


	Delaunay DT;
	DT.insert(points.begin(), points.end());
	set<pair<uint, uint> > considered_potentials;

	for (Delaunay::Finite_edges_iterator it_e = DT.finite_edges_begin(); it_e != DT.finite_edges_end(); it_e++) {
		Delaunay::Edge e = *it_e;
		uint e_i = e.first->vertex((e.second + 1) % 3)->info();
		uint e_j = e.first->vertex((e.second + 2) % 3)->info();
		uint i = points_to_segments[e_i];
		uint j = points_to_segments[e_j];

		if (i == j) continue;

		pair<uint, uint> p_ij = (i < j ? std::make_pair(i, j) : std::make_pair(j, i));
		if (considered_potentials.find(p_ij) != considered_potentials.end()) continue;
		considered_potentials.insert(p_ij);

		Segment* s_i = segments[i];
		Segment* s_j = segments[j];
		if (model->params->rega_quad_distance_considered && (Geometry::distance_initial_coordinates(s_i, s_j) > D)) continue;

		double mes_ij = s_i->alpha - s_j->alpha;
		double to_lower = 90 * floor(mes_ij / 90) - mes_ij;
		double to_upper = 90 * (floor(mes_ij / 90) + 1) - mes_ij;
		double t_ij = (fabs(to_lower) < fabs(to_upper) ? to_lower : to_upper);
		int r_ij;
		if (fabs(to_lower) < fabs(to_upper)) {
			r_ij = ((90 * int(floor(mes_ij / 90))) % 180 == 0 ? 0 : 1);
		} else {
			r_ij = ((90 * int(floor(mes_ij / 90) + 1)) % 180 == 0 ? 0 : 1);
		}
		double mu_ij = lambda;

		Point2d s_i_center = Point2d(jclamp(0, s_i->barycenter.x, model->I.cols - 1), jclamp(0, model->I.rows - s_i->barycenter.y, model->I.rows - 1));
		Point2d s_j_center = Point2d(jclamp(0, s_j->barycenter.x, model->I.cols - 1), jclamp(0, model->I.rows - s_j->barycenter.y, model->I.rows - 1));

		if ((r_ij == 0 && model->params->rega_quad_optimize_para) || (r_ij == 1 && model->params->rega_quad_optimize_ortho)) {
			if (fabs(t_ij) < dalpha_max[i] + dalpha_max[j]) {
				mu.push_back(Triplet<double>(i, j, mu_ij));
				targets.push_back(Triplet<double>(i, j, t_ij));
				relations.push_back(Triplet<int>(i, j, r_ij));
#if NOT_MEASURING_PERFORMANCES
				model->add_line(model->L_ag, s_i_center.x, s_i_center.y, s_j_center.x, s_j_center.y, 0, 0, 255);
#endif
			}
		}
	}
}


void Regularization_Angles_Quadratic::build_neighbors_graph_euclidian(Kinetic_Model* model, vector<Segment *> &segments, Boost_RTree &rtree_segments, double D, bool include_artificial, double lambda, vector<double> & dalpha_max,
	vector<Triplet<double> > &mu, vector<Triplet<double> > &targets, vector<Triplet<int> > & relations)
{
	for (uint i = 0; i < segments.size(); i++) {

		// Over a first phase, we search for segments that are located at less than D pixels from s_i
		Segment* s_i = segments[i];
		list<Segment *> neighbors;
		search_neighborhood(s_i, D, rtree_segments, segments, neighbors, include_artificial);

		Point2d s_i_center = Point2d(jclamp(0, s_i->barycenter.x, model->I.cols - 1), jclamp(0, model->I.rows - s_i->barycenter.y, model->I.rows - 1));

		// Now, we loop on each segment and we check if this value is close to an interesting angle
		for (list<Segment *>::iterator it_s = neighbors.begin(); it_s != neighbors.end(); it_s++) {
			Segment* s_j = (*it_s);
			if (s_i == s_j || s_i->index > s_j->index) continue;

			int j = s_j->index;
			double mes_ij = s_i->alpha - s_j->alpha;
			double to_lower = 90 * floor(mes_ij / 90) - mes_ij;
			double to_upper = 90 * (floor(mes_ij / 90) + 1) - mes_ij;

			double t_ij = (fabs(to_lower) < fabs(to_upper) ? to_lower : to_upper);
			int r_ij;
			if (fabs(to_lower) < fabs(to_upper)) {
				r_ij = ((90 * int(floor(mes_ij / 90))) % 180 == 0 ? 0 : 1);
			} else {
				r_ij = ((90 * int(floor(mes_ij / 90) + 1)) % 180 == 0 ? 0 : 1);
			}
			double mu_ij = lambda;

			if ((r_ij == 0 && model->params->rega_quad_optimize_para) || (r_ij == 1 && model->params->rega_quad_optimize_ortho)) {
				if (fabs(t_ij) < dalpha_max[i] + dalpha_max[j]) {
					mu.push_back(Triplet<double>(i, j, mu_ij));
					targets.push_back(Triplet<double>(i, j, t_ij));
					relations.push_back(Triplet<int>(i, j, r_ij));

					Point2d s_j_center = Point2d(jclamp(0, s_j->barycenter.x, model->I.cols - 1), jclamp(0, model->I.rows - s_j->barycenter.y, model->I.rows - 1));
#if NOT_MEASURING_PERFORMANCES
					model->add_line(model->L_ag, s_i_center.x, s_i_center.y, s_j_center.x, s_j_center.y, 0, 0, 255);
#endif
				}
			}
		}
	}
}



void Regularization_Angles_Quadratic::build_sparse_matrices(int individuals, vector<Triplet<double> > & v_mu, vector<Triplet<double> > & v_targets, vector<Triplet<int> > & v_relations, SparseMat<double> & mu, SparseMat<double> & targets, SparseMat<int> & relations)
{
    mu = SparseMat<double>(individuals, individuals);
    mu.setFromTriplets(v_mu.begin(), v_mu.end());
    targets = SparseMat<double>(individuals, individuals);
    targets.setFromTriplets(v_targets.begin(), v_targets.end());
	relations = SparseMat<int>(individuals, individuals);
	relations.setFromTriplets(v_relations.begin(), v_relations.end());
}



void Regularization_Angles_Quadratic::build_regularization_tree(vector<Segment *> & segments, SparseMat<double> & targets, SparseMat<int> & relations, vector<double> & dalpha, Segment_Regularization_Tree* & tree, double & theta_eps)
{
	int n = int(segments.size());
	vector<int> segments_to_groups(n, -1);
	map<int, list<int> > groups_to_segments = map<int, list<int> >();
	int g = 0, p = 0;

	for (int k = 0 ; k < targets.outerSize() ; ++k) {
		SparseMat<double>::InnerIterator it_t(targets, k);
		SparseMat<int>::InnerIterator it_r(relations, k);
		while (it_t && it_r) {
			int i = it_t.row(), j = it_t.col(), r = it_r.value();
			if (fabs(dalpha[n + p]) < 1e-6) {
				if (segments_to_groups[i] == -1 && segments_to_groups[j] == -1) {
					if (r == 0) {
						// Then segments i and j belong to the same group of parallel segments.
						// We should create a group of segments, that is initialized with these two individuals.
						segments_to_groups[i] = segments_to_groups[j] = g;
						groups_to_segments[g].push_back(i);
						groups_to_segments[g].push_back(j);
						++g;
					} else if (r == 1) {
						// The segments i and j are orthogonal.
						// We create two different groups of parallel segments.
						segments_to_groups[i] = g;
						groups_to_segments[g].push_back(i);
						segments_to_groups[j] = ++g;
						groups_to_segments[g].push_back(j);
						++g;
					}
				} else if (segments_to_groups[i] == -1 && segments_to_groups[j] != -1) {
					if (r == 0) {
						// Then segment i is parallel to j, and can be assigned to the same group
						int g_j = segments_to_groups[j];
						segments_to_groups[i] = g_j;
						groups_to_segments[g_j].push_back(i);
					} else if (r == 1) {
						// Then segment i is orthogonal to j, and we should initialize a new group with this segment.
						segments_to_groups[i] = g;
						groups_to_segments[g].push_back(i);
						++g;
					}
				} else if (segments_to_groups[i] != -1 && segments_to_groups[j] == -1) {
					// Symmetrical situation to before
					if (r == 0) {
						int g_i = segments_to_groups[i];
						segments_to_groups[j] = g_i;
						groups_to_segments[g_i].push_back(j);
					} else if (r == 1) {
						segments_to_groups[j] = g;
						groups_to_segments[g].push_back(j);
						++g;
					}
				} else {
					int g_i = segments_to_groups[i];
                    int g_j = segments_to_groups[j];
                    if (g_i != g_j) {
						if (r == 0) {
							// Segments i and j have been assigned to different groups, but in fact
							// they are parappel and belong to the same group. That's why we merge them.
							for (list<int>::iterator it_2 = groups_to_segments[g_j].begin() ; it_2 != groups_to_segments[g_j].end() ; it_2++) {
								segments_to_groups[*it_2] = g_i;
								groups_to_segments[g_i].push_back(*it_2);
							}
							groups_to_segments[g_j].clear();
						} else if (r == 1) {
							// We do nothing.
						}
                    }
				}
			}
			++p;
			++it_t;
			++it_r;
		}
	}

	/*
	map<int, double> angles_per_group;
	for (uint i = 0 ; i < segments.size() ; i++) {
		int g_i = segments_to_groups[i];
		if (g_i != -1 && angles_per_group.find(g_i) == angles_per_group.end()) {
			angles_per_group[g_i] = segments[i]->alpha + dalpha[i];
		}
	}

	vector<pair<int, int> > groups_sizes;
	for (map<int, list<int> >::iterator it_g = groups_to_segments.begin() ; it_g != groups_to_segments.end() ; it_g++) {
		groups_sizes.push_back(std::make_pair(it_g->first, it_g->second.size()));
	}
	*/

	// We prepare the construction of the regularization tree
	map<int, double> angles;

	for (uint i = 0 ; i < segments_to_groups.size() ; i++) {
		int g_i = segments_to_groups[i];
		if (g_i != -1) {
			if (angles.find(g_i) == angles.end()) {
				double theta = segments[i]->alpha + dalpha[i];
				if (theta < 0) {
					theta += 180;
				} else if (theta > 180) {
					theta -= 180;
				}

				// Checks if the angle that seems to be associated to this group of segments is not too close to another value
				int g_j = -1;
				for (map<int, double>::iterator it_m = angles.begin() ; it_m != angles.end() ; it_m++) {
					if (fabs(it_m->second - theta) < theta_eps) {
						g_j = it_m->first;
					}
				}

				if (g_j == -1) {
					angles[g_i] = theta;
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
		if (g_i == -1) {
			double alpha = segments[i]->alpha;

			int g_j = -1;
			for (map<int, double>::iterator it_m = angles.begin() ; it_m != angles.end() ; it_m++) {
				double alpha_j = it_m->second;
				for (int k = -1 ; k <= 1 ; k++) {
					if (fabs(alpha_j - alpha + k * 180) < theta_eps) {
						g_j = it_m->first;
						break;
					}
				}
				if (g_j != -1) break;
			}

			if (g_j == -1) {
				g_i = angles.rbegin()->first + 1;
				angles[g_i] = alpha;
			} else {
				g_i = g_j;
			}
			segments_to_groups[i] = g_i;
			groups_to_segments[g_i].push_back(i);
		}
	}

	// Builds the regularization tree
	for (map<int, double>::iterator it_m = angles.begin(); it_m != angles.end(); it_m++) tree->create_parallel_node(angles[it_m->first]);

	for (uint i = 0 ; i < segments_to_groups.size() ; i++) {
		// If the segments s_i is included in a group of parallel segments,
        // then it should be assigned to a leaf of the regularization tree.
        tree->assign_to_parallel_node(angles[segments_to_groups[i]], segments[i]);
	}
}
