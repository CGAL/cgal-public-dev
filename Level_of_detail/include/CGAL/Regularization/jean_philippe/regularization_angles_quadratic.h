#ifndef REGULARIZATION_ANGLES_QUADRATIC_H
#define REGULARIZATION_ANGLES_QUADRATIC_H

#include "defs.h"
#include "kinetic_model.h"
#include "regularization_angles.h"
#include "quadratic_problem.h"
#include <eigen3/Eigen/SparseCore>
#include <vector>
#include <map>

using std::vector;
using std::map;

class Regularization_Angles_Quadratic : public Regularization_Angles
{
public:
    Regularization_Angles_Quadratic();

    virtual ~Regularization_Angles_Quadratic();

    void regularize(Kinetic_Model* model);

    void set_bounds(vector<Segment *> & segments, vector<double> &dalpha_max, Parameters* params);

    void build_neighbors_graph(Kinetic_Model* model, vector<Segment *> & segments, Boost_RTree & rtree_segments, double D, bool include_artificial, double lambda, vector<double> &dalpha_max, 
		vector<Triplet<double> > &mu, vector<Triplet<double> > &targets, vector<Triplet<int> > & relations);

    void build_sparse_matrices(int individuals, vector<Triplet<double> > & v_mu, vector<Triplet<double> > & v_targets, vector<Triplet<int> > & v_relations, SparseMat<double> & mu, SparseMat<double> & targets, SparseMat<int> & relations);

    void build_regularization_tree(vector<Segment *> & segments, SparseMat<double> & targets, SparseMat<int> & relations, vector<double> & dalpha, Segment_Regularization_Tree* & tree, double & theta_eps);

private:
	void build_neighbors_graph_delaunay(Kinetic_Model* model, vector<Segment *> &segments, Boost_RTree &rtree_segments, double D, bool include_artificial, double lambda, vector<double> & dalpha_max,
		vector<Triplet<double> > &mu, vector<Triplet<double> > &targets, vector<Triplet<int> > & relations);

	void build_neighbors_graph_euclidian(Kinetic_Model* model, vector<Segment *> &segments, Boost_RTree &rtree_segments, double D, bool include_artificial, double lambda, vector<double> & dalpha_max,
		vector<Triplet<double> > &mu, vector<Triplet<double> > &targets, vector<Triplet<int> > & relations);
};

#endif // REGULARIZATION_ANGLES_QUADRATIC_H
