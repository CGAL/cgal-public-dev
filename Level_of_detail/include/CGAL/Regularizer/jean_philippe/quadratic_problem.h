#ifndef QUADRATIC_PROBLEM_H
#define QUADRATIC_PROBLEM_H

#include <eigen3/Eigen/SparseCore>
#include <vector>

template <typename T>
using SparseMat = Eigen::SparseMatrix<T, Eigen::RowMajor>;
using Eigen::Triplet;

using std::vector;


class Quadratic_Problem
{
public:
    Quadratic_Problem(int _individuals, int _variables, double _lambda, vector<double> & _bounds, SparseMat<double> & _mu, SparseMat<double> & _targets, vector<double> & _shift);

    virtual ~Quadratic_Problem();

protected:
    void allocate_objective_function(int** p_rowQ, int** p_colQ, double** p_dQ, double** p_c);

    void allocate_bounds(double** p_xlow, char** p_ixlow, double** p_xupp, char** p_ixupp);

    void allocate_inequality_constraints(int** p_rowC, int** p_colC, double** p_dC, double** p_clow, char** p_iclow, double** p_cupp, char** p_icupp);

    void desallocate_objective_function(int* rowQ, int* colQ, double* dQ, double* c);

    void desallocate_bounds(double* xlow, char* ixlow, double* xupp, char* ixupp);

    void desallocate_inequality_constraints(int* rowC, int* colC, double* dC, double* clow, char* iclow, double* cupp, char* icupp);

public:
    bool solve(vector<double> & solution);

public:
    int individuals;
    int variables;
    double lambda;
    vector<double> & bounds;
    SparseMat<double> & mu;
    SparseMat<double> & targets;
	vector<double> & shift;
};

#endif // QUADRATIC_PROBLEM_H
