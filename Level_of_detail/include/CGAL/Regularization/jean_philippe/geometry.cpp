#include "defs.h"
#include "geometry.h"
#include "segment_tree.h"
#include <iostream>


double Geometry::mes(Vec2d & _a, Vec2d & _b) {
    Vec2d a = normalize(_a);
    Vec2d b = normalize(_b);
    double x = a.ddot(b);
    if (x > 1.0 - 1e-7) {
        return 0.0;
    } else if (x < -1.0 + 1e-7) {
        return PI;
    } else {
        return acos(x);
    }
}


double Geometry::distance_initial_coordinates(Segment *s, Segment *t) {
    Vec2d & s_dir = s->direction;
    Vec2d & t_dir = t->direction;
    Point2d & A = s->end1;
    Point2d & B = s->end2;
    Point2d & C = t->end1;
    Point2d & D = t->end2;

    Vec2d s_nor(-s_dir[1], s_dir[0]);
    if (fabs(s_nor.ddot(t_dir)) > 1e-6) {
		// Computes the intersection
		double det = (B.x - A.x) * (C.y - D.y) - (B.y - A.y) * (C.x - D.x);
		double t_AB = ((C.x - A.x) * (C.y - D.y) - (C.y - A.y) * (C.x - D.x)) / det;
		double u_CD = ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) / det;
		if (t_AB >= 0 && t_AB <= 1 && u_CD >= 0 && u_CD <= 1) {
			// Intersection belongs to [AB] and [CD] : distance is null
			return 0;
		}
	}
	Vec2d AB = B - A, AC = C - A, AD = D - A;
	Vec2d CA = A - C, CB = B - C, CD = D - C;
	Vec2d BC = C - B, BD = D - B;
	double l_AB2 = norm(AB) * norm(AB);
	double l_CD2 = norm(CD) * norm(CD);
	double t_C = AC.ddot(AB) / l_AB2;
	double t_D = AD.ddot(AB) / l_AB2;
	double t_A = CA.ddot(CD) / l_CD2;
	double t_B = CB.ddot(CD) / l_CD2;
    double dist = MIN(MIN(norm(AC), norm(AD)), MIN(norm(BC), norm(BD)));
	if (t_C >= 0 && t_C <= 1) {
		Point2d C_prime = A + Point2d(t_C * AB);
        dist = MIN(dist, norm(C_prime - C));
	}
	if (t_D >= 0 && t_D <= 1) {
		Point2d D_prime = A + Point2d(t_D * AB);
        dist = MIN(dist, norm(D_prime - D));
	}
	if (t_A >= 0 && t_A <= 1) {
		Point2d A_prime = C + Point2d(t_A * CD);
        dist = MIN(dist, norm(A_prime - A));
	}
	if (t_B >= 0 && t_B <= 1) {
		Point2d B_prime = C + Point2d(t_B * CD);
        dist = MIN(dist, norm(B_prime - B));
	}
	return dist;
}


double Geometry::distance(cv::Point2d &p, Segment *s)
{
	return (fabs(s->a * p.x + s->b * p.y + s->c) / sqrt(s->a * s->a + s->b * s->b));
}

/*
bool Geometry::are_neighbors_euclidian(Segment *s, Segment *t, double d_lim, double & d)
{
	if (bounding_boxes_intersect(s, t, d_lim)) {
		d = distance(s, t);
		return (d <= d_lim);
	} else {
		return false;
	}
}
*/

/*
bool Geometry::bounding_boxes_intersect(Segment *s, Segment *t, double d)
{
	int d_0 = (int(ceil(d)) + 1) / 2;
	int s_xmin = s->bbox_x_min - d_0, s_xmax = s->bbox_x_max + d_0;
	int s_ymin = s->bbox_y_min - d_0, s_ymax = s->bbox_y_max + d_0;
	int t_xmin = t->bbox_x_min - d_0, t_xmax = t->bbox_x_max + d_0;
	int t_ymin = t->bbox_y_min - d_0, t_ymax = t->bbox_y_max + d_0;

	bool inter_x = ((s_xmin <= t_xmin && t_xmin <= s_xmax) || (s_xmin <= t_xmax && t_xmax <= s_xmax)) || ((t_xmin <= s_xmin && s_xmin <= t_xmax) || (t_xmin <= s_xmax && s_xmax <= t_xmax));
	if (inter_x) {
		bool inter_y = ((s_ymin <= t_ymin && t_ymin <= s_ymax) || (s_ymin <= t_ymax && t_ymax <= s_ymax)) || ((t_ymin <= s_ymin && s_ymin <= t_ymax) || (t_ymin <= s_ymax && s_ymax <= t_ymax));
		return inter_y;
	} else {
		return false;
	}
}
*/

/*
bool Geometry::bounding_boxes_intersect(Segment *s, Segment *t) {
	bool inter_x = (s->bbox_x_min <= t->bbox_x_min && t->bbox_x_min <= s->bbox_x_max) || (s->bbox_x_min <= t->bbox_x_max && t->bbox_x_max <= s->bbox_x_max);
	if (inter_x) {
		bool inter_y = (s->bbox_y_min <= t->bbox_y_min && t->bbox_y_min <= s->bbox_y_max) || (s->bbox_y_min <= t->bbox_y_max && t->bbox_y_max <= s->bbox_y_max);
		return inter_y;
	}
	return false;
}
*/



void Geometry::build_tree_when_regularization_disabled(vector<Segment *> & segments, Segment_Regularization_Tree * tree)
{
	map<int, list<Segment *> > groups;
	uint n = uint(segments.size());

	double alpha_eps = 0.25;
	double y_eps = 1;
	// Part 1. Discretization of the set of observed orientations

	// First, we associate each segment to a group of parallel segments, indexed by angle
	
	vector<double> angles;
	for (uint i = 0 ; i < n ; ++i) {
		int related_group = -1;
		double alpha = segments[i]->alpha;

		// Tries to match a segment's orientation with another angle
		for (uint j = 0 ; j < groups.size() ; ++j) {
			double alpha_j = angles[j];
			for (int k = -1 ; k <= 1 ; ++k) {
				if (fabs(alpha - alpha_j + k * 180) < alpha_eps) {
					related_group = j;
					break;
				}
			}
			if (related_group != -1) break;
		}

		if (related_group == -1) {
			// Creates a new group with alpha_i
			related_group = int(groups.size());
			angles.push_back(alpha);
		}
		groups[related_group].push_back(segments[i]);
	}

	// Builds the first level of the tree
	for (map<int, list<Segment *> >::iterator it_m = groups.begin() ; it_m != groups.end() ; ++it_m) {
		double theta = angles[it_m->first];

		Vec2d v_dir = Vec2d(cos(theta * PI / 180), sin(theta * PI / 180));
        Vec2d v_nor = Vec2d(-v_dir[1], v_dir[0]);
        double a = v_nor[0], b = v_nor[1];

		tree->create_parallel_node(theta);
		for (list<Segment *>::iterator it_s = it_m->second.begin() ; it_s != it_m->second.end() ; ++it_s) {
			Segment* s = (*it_s);

			// Computes the equation of the support line of s
			double c = -a * s->barycenter.x - b * s->barycenter.y;
            s->set_dalpha(theta - s->alpha, theta, a, b, c, v_dir);
			tree->assign_to_parallel_node(theta, s);
		}
	}

	// Part 2. Detects near-colinear segments

	for (map<double, Node_Parallel_Segments *>::iterator it_c = tree->parallel_segments.begin() ; it_c != tree->parallel_segments.end() ; ++it_c) {

		groups.clear();
		double theta = it_c->first;
		Node_Parallel_Segments* cluster = it_c->second;
		Point2d O = cluster->parallel_segments.front()->interBarycenter;

		// We transform the coordinates of the barycenters of all segments in the local frame defined by :
		// - the center of the first segment in the cluster
		// - the unit vectors I(cos(theta), sin(theta)) and J(-sin(theta), cos(theta))
		Vec2d I = Vec2d(cos(theta * PI / 180), sin(theta * PI / 180));
        Vec2d J = Vec2d(-I[1], I[0]);

		vector<double> ordinates;
		for (list<Segment *>::iterator it_s = cluster->parallel_segments.begin() ; it_s != cluster->parallel_segments.end() ; ++it_s) {
			Segment* s = (*it_s);
			double dx = s->interBarycenter.x - O.x, dy = s->interBarycenter.y - O.y;
			// double x = dx * I[0] + dy * I[1];
			double y = dx * J[0] + dy * J[1];

			// Following a process directly adapted from the first part of this function,
			// we want to merge segments whose y is too close to other values of y computed before (near-colinear segments)
			int related_group = -1;
			for (uint j = 0 ; j < groups.size() ; ++j) {
				double y_j = ordinates[j];
				if (fabs(y_j - y) < y_eps) {
					related_group = j;
					break;
				}
				if (related_group != -1) break;
			}

			if (related_group == -1) {
				related_group = groups.size();
				ordinates.push_back(y);
			}
			groups[related_group].push_back(s);
		}

		// Builds the second level of the tree
		for (map<int, list<Segment *> >::iterator it_m = groups.begin() ; it_m != groups.end() ; ++it_m) {
			
			double dt = ordinates[it_m->first];

			// Gets the longest segment
			double l_max = -FLT_MAX;
			Segment* s_longest = NULL;
			for (list<Segment *>::iterator it_s = it_m->second.begin() ; it_s != it_m->second.end() ; ++it_s) {
				if ((*it_s)->length > l_max) {
					l_max = (*it_s)->length;
					s_longest = (*it_s);
				}
			}

			cluster->create_colinear_node(dt);

			// Translates the longest segment and gets its line equation
			double dx = s_longest->interBarycenter.x - O.x, dy = s_longest->interBarycenter.y - O.y;
			double y = dx * J[0] + dy * J[1];

			s_longest->set_dt(dt - y);
			double a = s_longest->a, b = s_longest->b, c = s_longest->c;
			Vec2d dir = s_longest->finalDirection;

			// Translates other segments
			for (list<Segment *>::iterator it_s = it_m->second.begin() ; it_s != it_m->second.end() ; ++it_s) {
				Segment* s = (*it_s);
				if (s == s_longest) continue;

				dx = s->interBarycenter.x - O.x, dy = s->interBarycenter.y - O.y;
				y = dx * J[0] + dy * J[1];
				s->set_dt(dt - y, a, b, c, dir);
			}
			cluster->assign_to_colinear_node(dt, it_m->second);
		}
	}
}


#if 0
void Geometry::find_groups_of_parallel_segments(Segment_Regularization_Tree *tree, vector<vector<int> > & groups)
{
	groups.clear();
	vector<Vec2d> group_directions = vector<Vec2d>();

	// The tree passed as argument of this function already did the work, at least partially
	// Indeed, it is a tree composed of nodes of parallel segments, and a node of parallel segments has, itself, one or several nodes containing list of colinear segments
	// So, here, we just go down the tree to 
	for (map<double, Node_Parallel_Segments *>::iterator it_m1 = tree->parallel_segments.begin() ; it_m1 != tree->parallel_segments.end() ; it_m1++) {
		
		// Builds the direction vector 
		double alpha = it_m1->first;
		Vec2d alpha_dir = Vec2d(cos(alpha * PI / 180), sin(alpha * PI / 180));
		list<int> alpha_indices;

		Node_Parallel_Segments* subtree_1 = it_m1->second;

		// Explores the nodes with colinear segments
		for (map<double, Node_Colinear_Segments *>::iterator it_m2 = subtree_1->colinear_segments.begin() ; it_m2 != subtree_1->colinear_segments.end() ; it_m2++) {
			Node_Colinear_Segments* subtree_2 = it_m2->second;
			for (list<Segment *>::iterator it_s = subtree_2->colinear_segments.begin() ; it_s != subtree_2->colinear_segments.end() ; it_s++) {
				alpha_indices.push_back((*it_s)->index);
			}
		}
		// Other parallel segments of the node
		for (list<Segment *>::iterator it_s = subtree_1->other_segments.begin() ; it_s != subtree_1->other_segments.end() ; it_s++) {
			alpha_indices.push_back((*it_s)->index);
		}

		group_directions.push_back(alpha_dir);
		groups.push_back(vector<int>(alpha_indices.begin(), alpha_indices.end()));
	}
	std::cout << "find_groups_of_parallel_segments, step 1 : " << groups.size() << " groups" << std::endl;

	// The segments for which there is no parallelism relationship a priori, may be parallel with some of the other segments
	for (list<Segment *>::iterator it_s = tree->other_segments.begin() ; it_s != tree->other_segments.end() ; it_s++) {
        Vec2d s_dir = (*it_s)->finalDirection;
		Vec2d s_nor = Vec2d(-s_dir[1], s_dir[0]);

		// Checks if a group already exists
		int assigned_group = -1;
		for (int i = 0 ; i < group_directions.size() ; i++) {
			Vec2d g_dir = group_directions[i];
			if (fabs(s_nor.ddot(g_dir)) < 1e-6) {
				assigned_group = i;
				break;
			}
		}

		if (assigned_group == -1) {
			// Creates a new group, saves its normal
			groups.push_back(vector<int>(1, (*it_s)->index));
			group_directions.push_back(s_dir);
		} else {
			groups[assigned_group].push_back((*it_s)->index);
		}
	}

	// Removes groups consisted of one isolated segment
	vector< vector<int> >::iterator it = groups.begin();
	while (it != groups.end()) {
		if ((*it).size() > 1) {
			it++;
		} else {
			it = groups.erase(it);
		}
	}
	std::cout << "find_groups_of_parallel_segments, step 2 : " << groups.size() << " groups" << std::endl;
}


void Geometry::find_groups_of_parallel_segments(vector<Segment *> & segments, vector<vector<int> > & groups)
{
	int n = int(segments.size());
	groups.clear();
	vector<Vec2d> group_directions = vector<Vec2d>();

	// Associates each segment to a group of parallel segments
	for (int i = 0; i < n; i++) {
		int related_group = -1;
        Vec2d s_dir = segments[i]->finalDirection;
		Vec2d s_nor = Vec2d(-s_dir[1], s_dir[0]);

		// Loops on all existing groups to see if one corresponds
		for (int j = 0; j < groups.size(); j++) {
			Vec2d g_dir = group_directions[j];
			if (fabs(s_nor.ddot(g_dir)) < 1e-6) {
				related_group = j;
				break;
			}
		}

		if (related_group == -1) {
			// Creates a new group, saves its normal
			groups.push_back(vector<int>(1, i));
			group_directions.push_back(s_dir);
		} else {
			groups[related_group].push_back(i);
		}
	}

	// Removes groups consisted of one isolated segment
	vector< vector<int> >::iterator it = groups.begin();
	while (it != groups.end()) {
		if ((*it).size() > 1) {
			it++;
		} else {
			it = groups.erase(it);
		}
	}
}



void Geometry::find_groups_of_colinear_segments(vector<Segment *> & segments, list< pair<int, int> > & colinear_indices)
{
	colinear_indices.clear();
	int n = int(segments.size());

	// First, we are going to define groups of directions
	// Every segment, depending on its direction, will be assigned to a group
	vector<Vec2d> groups_directions;
	vector<vector<int> > groups;

	for (int i = 0; i < n; i++) {
		int related_group = -1;
        Vec2d s_dir = segments[i]->finalDirection;
		Vec2d s_nor = Vec2d(-s_dir[1], s_dir[0]);
		// Loops on all existing groups to see if one corresponds
		for (int j = 0; j < groups.size(); j++) {
			Vec2d g_dir = groups_directions[j];
			if (fabs(s_nor.ddot(g_dir)) < 1e-6) {
				related_group = j;
				break;
			}
		}
		if (related_group == -1) {
			// Creates a new group, saves its normal
			// Note : we save the index that correspond to the related segment, not to the ray
			// Indeed, every segment is duplicated to form a couple of opposite rays
			groups.push_back(vector<int>(1, i));
			groups_directions.push_back(s_dir);
		} else {
			// Assigns the segment to the related group
			groups[related_group].push_back(i);
		}
	}

	// We are now going to search for groups of colinear segments among groups of parallel segments
	for (int g = 0; g < groups.size(); g++) {
		Vec2d g_dir = groups_directions[g];
		Vec2d g_nor = Vec2d(-g_dir[1], g_dir[0]);
		bool horizontal = (fabs(g_dir[0]) < 1e-6 || fabs(g_dir[1] / g_dir[0]) < 1);
		double tol = (horizontal ? fabs(g_dir[1] / g_dir[0]) : fabs(g_dir[0] / g_dir[1]));

		// For each couple (i, j) of parallel segments, we check if they are colinear
		// To this end, we exhibit the expression of the linear function to which the respective
		// centers of the these segments belong. We already know a and b in ax + by + c = 0 : we
		// only need to compute c. If two segments have the same c, it means that they are
		// included in the same group of colinear segments.
		vector<double> subgroups_intercepts;
		vector<vector<int> > subgroups;
		for (int i = 0; i < groups[g].size(); i++) {
            Point2d O = segments[groups[g][i]]->finalBarycenter;
			double intercept = -g_nor[0] * O.x - g_nor[1] * O.y;
			if (!horizontal) {
				// If the segment is vertical, we should rather classify the segment with some
				// kind of "x-intercept" which presents less disparities
				intercept = -intercept / g_nor[0];
			}
			int related_subgroup = -1;
			// Tries to assign the intercept of the current segment to an existing group
			for (int j = 0; j < subgroups_intercepts.size(); j++) {
				double intercept_j = subgroups_intercepts[j];
				if (fabs(intercept_j - intercept) < (1 + tol)) {
					related_subgroup = j;
					break;
				}
			}
			// Decision : creates a new subgroup of colinear segments, or
			// appends the index of the current segment to a list of colinear segments
			// Warning : indices that are saved are related to the element's position in groups[g]
			if (related_subgroup == -1) {
				subgroups.push_back(vector<int>(1, i));
				subgroups_intercepts.push_back(intercept);
			} else {
				subgroups[related_subgroup].push_back(i);
			}
		}

		double y_intercept = 0;
		// For each group of colinear segments obtained
		for (int j = 0; j < subgroups.size(); j++) {
			if (horizontal) {
				y_intercept = subgroups_intercepts[j];
			} else {
				int subindex = subgroups[j][0];
                Point2d O = segments[groups[g][subindex]]->finalBarycenter;
				y_intercept = -g_nor[0] * O.x - g_nor[1] * O.y;
			}
			for (int k = 0; k < subgroups[j].size(); k++) {
				// Adjusts the coordinates of both ends of each colinear segment, so that they really
				// rest on the same affine line entirely, which is defined by the following three
				// coefficients
				int subindex_k = subgroups[j][k];
				segments[groups[g][subindex_k]]->move_to_line(g_nor[0], g_nor[1], y_intercept);
				// Adds pairs of indices of colinear segments to the dedicated list
				for (int l = k + 1; l < subgroups[j].size(); l++) {
					int subindex_l = subgroups[j][l];
					colinear_indices.push_back(std::make_pair(groups[g][subindex_k], groups[g][subindex_l]));
				}
			}
		}

		for (int j = 0; j < subgroups.size(); j++) subgroups[j].clear();
		subgroups.clear();
		subgroups_intercepts.clear();
	}

	for (int i = 0; i < groups.size(); i++) groups[i].clear();
	groups.clear();
	groups_directions.clear();
}
#endif


void Geometry::merge_overlapped_segments(Segment_Regularization_Tree *tree)
{
	for (map<double, Node_Parallel_Segments*>::iterator it_m1 = tree->parallel_segments.begin() ; it_m1 != tree->parallel_segments.end() ; it_m1++) {
		Node_Parallel_Segments* node_parallel = it_m1->second;
		for (map<double, Node_Colinear_Segments*>::iterator it_m2 = node_parallel->colinear_segments.begin() ; it_m2 != node_parallel->colinear_segments.end() ; it_m2++) {
			Node_Colinear_Segments* node_colinear = it_m2->second;
			if (node_colinear->colinear_segments.size() == 1) continue;

			bool at_least_one_merge;
			do {
				at_least_one_merge = false;
				list<Segment *>::iterator it_s1, it_s2;
				for (it_s1 = node_colinear->colinear_segments.begin() ; it_s1 != node_colinear->colinear_segments.end() ; ++it_s1) {
					Segment* s1 = (*it_s1);
					if (s1->is_disabled) {
						continue;
					}
					it_s2 = it_s1;
					while (++it_s2 != node_colinear->colinear_segments.end()) {
						Segment* s2 = (*it_s2);
						if (s2->is_disabled) {
							continue;
						}
						Point2d P, Q;
						// Tests if colinear segments s1 and s2 are overlapping
						if (are_overlapping(s1, s2, P, Q)) {
							// s1 takes the coordinates of the merged segment s1 U s2
							s1->set_final_extrema(P, Q);
#if NOT_MEASURING_PERFORMANCES
							// Update the lists of brothers
							s1->merge_brothers(s2);
							s2->clear_brothers();
#endif
							s2->disable();
							at_least_one_merge = true;
						}
					}
				}
			} while (at_least_one_merge);
		}
	}
}


#if 0
void Geometry::merge_overlapped_segments(vector<Segment *> & segments, list<pair<int, int> > & colinear_indices)
{
	Point2d P = Point2d(0, 0);
	Point2d Q = Point2d(0, 0);
    vector<int> to_disable;
	bool at_least_one_merge = false;

	do {
		at_least_one_merge = false;
		// Reminder : each element (i, j) of the list satisfies (i < j)
		list<pair<int, int> >::iterator it1 = colinear_indices.begin(), it2;
		while (it1 != colinear_indices.end()) {
			int i = it1->first, j = it1->second;
			Segment *si = segments[i], *sj = segments[j];

			// Tests if colinear segments si and sj are overlapping
			bool overlap = are_overlapping(si, sj, P, Q);
			if (overlap) {
				// Segment si takes the coordinates of the merged segment [si U sj]
                si->set_final_extrema(P, Q);
                // Update the lists of brothers
                si->merge_brothers(sj);
                sj->clear_brothers();
				// Segment sj is going to be deleted
                to_disable.push_back(j);
				// Removes further references to segment j in the list of colinear_indices
				// We initialize a new iterator that starts to explore the list from the element
				// located just after the current pair
				it2 = colinear_indices.begin();
				while (it2 != colinear_indices.end()) {
					if ((it2->first == j || it2->second == j) && it1 != it2) {
						it2 = colinear_indices.erase(it2);
					} else {
						it2++;
					}
				}
				// Finally, removes the current element properly
				it1 = colinear_indices.erase(it1);
				at_least_one_merge = true;
			} else {
				it1++;
			}
		}
	} while (at_least_one_merge);

    for (uint i = 0 ; i < to_disable.size() ; i++) {
        segments[to_disable[i]]->disable();
    }
}
#endif


void Geometry::disable_segments_outside_boundaries(vector<Segment *> & segments, int rows, int cols)
{
	uint n = uint(segments.size());
	for (uint i = 0 ; i < n ; i++) {
		Point2d & O = segments[i]->finalBarycenter;
		if (O.x <= 0 || O.x >= cols || O.y <= 0 || O.y >= rows) {
			segments[i]->disable();
		}
	}
}


void Geometry::intersection_boundary(Ray* r_i, int rows, int cols, double & t_i, Image_Boundary & r_j, double & t_j)
{
	Vec2d dir = r_i->OA;
	t_i = FLT_MAX;

	if (dir[1] > 0) {
		// Collision with the line y = rows
		double _t_i = (rows - r_i->A.y) / dir[1];
		if (_t_i < t_i) {
			t_i = _t_i;
			t_j = r_i->A.x + t_i * dir[0];
			r_j = TOP_IMAGE;
		}
	} else if (dir[1] < 0) {
		// Collision with the line y = 0
		double _t_i = (0 - r_i->A.y) / dir[1];
		if (_t_i < t_i) {
			t_i = _t_i;
			t_j = r_i->A.x + t_i * dir[0];
			r_j = BOTTOM_IMAGE;
		}
	}

	if (dir[0] < 0) {
		// Collision with the line x = 0
		double _t_i = (0 - r_i->A.x) / dir[0];
		if (_t_i < t_i) {
			t_i = _t_i;
			t_j = r_i->A.y + t_i * dir[1];
			r_j = LEFT_IMAGE;
		}
	} else if (dir[0] > 0) {
		// Collision with the line x = cols
		double _t_i = (cols - r_i->A.x) / dir[0];
		if (_t_i < t_i) {
			t_i = _t_i;
			t_j = r_i->A.y + t_i * dir[1];
			r_j = RIGHT_IMAGE;
		}
	}
}


#if 0
void Geometry::intersection_boundary(Ray *s, int rows, int cols, double & s_time, Image_Boundary & hit_boundary) {
	Vec2d dir = s->OA;
	s_time = FLT_MAX;
	if (dir[1] > 0) {
		// Intersection with y = rows
		double s_0 = (rows - s->A.y) / dir[1];
		if (s_0 < s_time) {
			s_time = s_0;
			hit_boundary = TOP_IMAGE;
		}
	} else if (dir[1] < 0) {
		// Intersection with y = 0
		double s_0 = (0 - s->A.y) / dir[1];
		if (s_0 < s_time) {
			s_time = s_0;
			hit_boundary = BOTTOM_IMAGE;
		}
	}
	if (dir[0] < 0) {
		// Intersection with x = 0
		double s_0 = (0 - s->A.x) / dir[0];
		if (s_0 < s_time) {
			s_time = s_0;
			hit_boundary = LEFT_IMAGE;
		}
	} else if (dir[0] > 0) {
		// Intersection with x = cols
		double s_0 = (cols - s->A.x) / dir[0];
		if (s_0 < s_time) {
			s_time = s_0;
			hit_boundary = RIGHT_IMAGE;
		}
	}
}
#endif


void Geometry::direct_intersection(Ray* r_i, Image_Boundary b_j, double rows, double cols, double & t_i, double & t_j)
{
	switch (b_j) {
	case TOP_IMAGE:
		t_i = (rows - r_i->A.y) / r_i->OA[1];
		t_j = r_i->A.x + t_i * r_i->OA[0];
		return;
	case BOTTOM_IMAGE:
		t_i = (0 - r_i->A.y) / r_i->OA[1];
		t_j = r_i->A.x + t_i * r_i->OA[0];
		return;
	case LEFT_IMAGE:
		t_i = (0 - r_i->A.x) / r_i->OA[0];
		t_j = r_i->A.y + t_i * r_i->OA[1];
		return;
	case RIGHT_IMAGE :
		t_i = (cols - r_i->A.x) / r_i->OA[0];
		t_j = r_i->A.y + t_i * r_i->OA[1];
		return;
	}
}

/*
void Geometry::direct_intersection(Ray* r_i, Ray* r_j, double & t_i, double & t_j) 
{
	Point2d &B = r_i->A, &D = r_j->A;
	Vec2d &AB = r_i->OA, &CD = r_j->OA;

	// Let D (t_j) = D + t_j * CD
	// Finds t_j such that a * (D + t_j * CD).x + b * (D + t_j * CD).y + c = 0
	// where a, b and c define the support line of r_i
	double a_i = r_i->parent->a, b_i = r_i->parent->b, c_i = r_i->parent->c;
	t_j = -(a_i * D.x + b_i * D.y + c_i) / (a_i * CD[0] + b_i * CD[1]);

	// Same for t_i
	double a_j = r_j->parent->a, b_j = r_j->parent->b, c_j = r_j->parent->c;
	t_i = -(a_j * B.x + b_j * B.y + c_j) / (a_j * AB[0] + b_j * AB[1]);
}
*/


void Geometry::direct_intersection(Ray* r_i, Ray* r_j, double & t_i, double & t_j) 
{
	Vec2d minusrj = -r_j->OA;

	Point2d &B = r_i->A, &D = r_j->A;
	Vec2d &AB = r_i->OA, &DC = minusrj;
	double det = AB[0] * DC[1] - AB[1] * DC[0];
	if (fabs(det) < 1e-6) {
		t_i = t_j = FLT_MAX;
	} else {
		Vec2d BD = D - B;
		t_i = (BD[0] * DC[1] - BD[1] * DC[0]) / det;
		t_j = (AB[0] * BD[1] - AB[1] * BD[0]) / det;
	}
}


bool Geometry::intersection(Ray *s, Ray *t, double max_time_s0, double max_time_s1, 
	unsigned int & s_index, unsigned int & t_index, double & s_time, double & t_time) {
	Point2d & B = s->A;
	Point2d & D = t->A;
	Vec2d & AB = s->OA;

	Vec2d minust = -t->OA;
	Vec2d & DC = minust;
	// Let B(x) = B + x * AB, D(y) = D + y * CD
	// We aim at finding x_0 > 0 and y_0 > 0 such that B(x_0) = C(y_0)
	double determinant = AB[0] * DC[1] - AB[1] * DC[0];
	// If segments are parallel, return false
	if (fabs(determinant) < 1e-6) {
		return false;
	} else {
		Vec2d BD = D - B;
		double l_AB = s->initial_length;
		double l_CD = t->initial_length;
		//double l_AB = sqrt((B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y));
		//double l_CD = sqrt((D.x - C.x) * (D.x - C.x) + (D.y - C.y) * (D.y - C.y));
		double x_i = (BD[0] * DC[1] - BD[1] * DC[0]) / determinant;
		double y_i = (AB[0] * BD[1] - AB[1] * BD[0]) / determinant;
		if (x_i >= -l_AB) {
			s_index = s->index;
			s_time = x_i;
			if (s_time > max_time_s0) return false;
		} else {
			s_index = s->index + 1;
			s_time = -2 * l_AB - x_i;
			if (s_time > max_time_s1) return false;
		}
		if (y_i >= -l_CD) {
			t_index = t->index;
			t_time = y_i;
		} else {
			t_index = t->index + 1;
			t_time = -2 * l_CD - y_i;
		}
		return true;
	}
}


bool Geometry::intersection_colinear(Ray* s, Ray* t, double & s_time, double & t_time)
{
	Point2d & B = s->A;
	Vec2d & AB = s->OA;
	Point2d & D = t->A;
	Vec2d & CD = t->OA;
	// As s and t are colinear, either AB = CD or AB = -CD
	// If AB = CD, the equation has no solution
	assert(fabs(AB[0] * CD[1] - AB[1] * CD[0]) < 1e-5);
	if (fabs(AB[0] - CD[0]) < 1e-5 && fabs(AB[1] - CD[1]) < 1e-5) {
		return false;
	} else {
		Point2d I = (B + D) / 2;
		// We determine s_0 such that I = B + s_0 * AB
		double s_0 = (fabs(AB[0]) > 1e-6 ? (I.x - B.x) / AB[0] : (I.y - B.y) / AB[1]);
		// We determine t_0 such that I = D + t_0 * CD
		double t_0 = (fabs(CD[0]) > 1e-6 ? (I.x - D.x) / CD[0] : (I.y - D.y) / CD[1]);
		// We assert s_0 and t_0 are valid values
		if (s_0 > -s->initial_length && t_0 > -t->initial_length) {
			s_time = s_0;
			t_time = t_0;
			return true;
		} else {
			return false;
		}
	}
}


bool Geometry::are_overlapping(Segment *s, Segment *t, Point2d & P, Point2d & Q)
{
	// If s and t are not parallel, they can't be overlapped
    Vec2d & AB = s->finalDirection;
    Vec2d & CD = t->finalDirection;
	if (fabs(AB[0] * CD[1] - AB[1] * CD[0]) > 1e-6) {
		return false;
	} else {
		// We determine if C or D belongs to [AB]
        Point2d & A = s->finalEnd1;
        Point2d & B = s->finalEnd2;
        Point2d & C = t->finalEnd1;
        Point2d & D = t->finalEnd2;
		double l_AB = s->length;
		double t_C, t_D;
		if (fabs(AB[0]) < 1e-6) {
			if (fabs(C.x - A.x) < 1e-3) {
				t_C = (C.y - A.y) / AB[1];
				t_D = (D.y - A.y) / AB[1];
			} else {
				return false;
			}
		} else if (fabs(AB[1]) < 1e-6) {
			if (fabs(C.y - A.y) < 1e-3) {
				t_C = (C.x - A.x) / AB[0];
				t_D = (D.x - A.x) / AB[0];
			} else {
				return false;
			}
		} else {
			double t_Cx = (C.x - A.x) / AB[0];
			double t_Cy = (C.y - A.y) / AB[1];
			double t_Dx = (D.x - A.x) / AB[0];
			double t_Dy = (D.y - A.y) / AB[1];
			if ((fabs(t_Cx - t_Cy) < 1e-3) && (fabs(t_Dx - t_Dy) < 1e-3)) {				
				t_C = t_Cx;
				t_D = t_Dx;
			} else {
				return false;
			}
		}
		if ((t_C < 0 && t_D < 0) || t_C > l_AB && t_D > l_AB) {
			return false;
		} else {
			// Now, we know are [AB] and [CD] are overlapping
			// A, B, C and D are assigned to abscissas on the line [AB]
			// The ends of the merged segment correspond to the argmin and the argmax of such abscissas
			vector<Point2d> points;
			vector<double> abscissas;
			points.push_back(A); points.push_back(B); 
			points.push_back(C); points.push_back(D);
			abscissas.push_back(0); abscissas.push_back(l_AB);
			abscissas.push_back(t_C); abscissas.push_back(t_D);
			double minimum = FLT_MAX, maximum = -FLT_MAX;
			int argmin = -1, argmax = -1;
			for (int i = 0 ; i < 4 ; i++) {
				if (abscissas[i] < minimum) {
					minimum = abscissas[i];
					argmin = i;
				}
				if (abscissas[i] > maximum) {
					maximum = abscissas[i];
					argmax = i;
				}
			}
			P = points[argmin];
			Q = points[argmax];
			return true;
		}
	}
}
