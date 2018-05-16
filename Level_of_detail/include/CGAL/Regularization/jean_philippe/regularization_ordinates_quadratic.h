#pragma once
#include "defs.h"
#include "kinetic_model.h"
#include "regularization_ordinates.h"
#include "quadratic_problem.h"
#include <eigen3/Eigen/SparseCore>
#include <vector>

using std::vector;


class Regularization_Ordinates_Quadratic : public Regularization_Ordinates
{
public:
	Regularization_Ordinates_Quadratic();

	virtual ~Regularization_Ordinates_Quadratic();

	void regularize(Kinetic_Model* model);

	void set_bounds(vector<Segment *> & segments, vector<double> & dt_max, Parameters* params);

	void build_neighbors_graph(Kinetic_Model* model, Node_Parallel_Segments* node, vector<Segment *> & segments, Boost_RTree & rtree_segments, double lambda, vector<double> & dt_max, vector<Triplet<double> > &mu, vector<Triplet<double> > &targets);

	void build_sparse_matrices(int individuals, vector<Triplet<double> > & v_mu, vector<Triplet<double> > & v_targets, SparseMat<double> & mu, SparseMat<double> & targets);

	void build_regularization_tree(vector<Segment *> & segments, SparseMat<double> & targets, vector<double> & dt, Segment_Regularization_Tree* & tree);
};

