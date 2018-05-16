#pragma once
#include <utility>
#include <vector>
#include "segment_ray.h"
#include "segment_tree.h"

using std::vector;
using std::pair;
using std::make_pair;


namespace Geometry
{
    double mes(Vec2d & _a, Vec2d & _b);

    double distance_initial_coordinates(Segment *s, Segment *t);

	double distance(Point2d & p, Segment *s);

    //bool are_neighbors_euclidian(Segment *s, Segment *t, double d_lim, double & d);

	//bool bounding_boxes_intersect(Segment *s, Segment *t, double d_lim);

	//bool bounding_boxes_intersect(Segment *s, Segment *t);

	void build_tree_when_regularization_disabled(vector<Segment *> & segments, Segment_Regularization_Tree * tree);

    //void find_groups_of_parallel_segments(Segment_Regularization_Tree* tree, vector<vector<int> > & groups);

	//void find_groups_of_parallel_segments(vector<Segment *> & segments, vector<vector<int> > & groups);

	//void find_groups_of_colinear_segments(Segment_Reguarization_Tree* tree, vector<Segment *> & segments, list<pair<int, int> > & colinear_indices);

	//void find_groups_of_colinear_segments(vector<Segment *> & segments, list<pair<int, int> > & colinear_indices);

	void merge_overlapped_segments(Segment_Regularization_Tree *tree);

	void disable_segments_outside_boundaries(vector<Segment *> & segments, int rows, int cols);

	//void merge_overlapped_segments(vector<Segment *> & segments, list<pair<int, int> > & colinear_indices);

	void intersection_boundary(Ray* r_i, int rows, int cols, double & t_i, Image_Boundary & r_j, double & t_j);
	
	// void intersection_boundary(Ray *s, int rows, int cols, double & s_time, Image_Boundary & hit_boundary);

    //void intersection(Ray *s, Ray *t, int rows, int cols, double & s_time, double & t_time);

	void direct_intersection(Ray* r_i, Image_Boundary b_j, double rows, double cols, double & t_i, double & t_j);

	void direct_intersection(Ray* r_i, Ray* r_j, double & t_i, double & t_j);

	bool intersection(Ray *s, Ray *t, double max_time_s0, double max_time_s1, unsigned int & s_index, unsigned int & t_index, double & s_time, double & t_time);

	bool intersection_colinear(Ray* s, Ray* t, double & s_time, double & t_time);

	bool are_overlapping(Segment *s, Segment *t, Point2d & P, Point2d & Q);
};

bool inline in(double a, double x, double b) { return (a <= x && x <= b); }
