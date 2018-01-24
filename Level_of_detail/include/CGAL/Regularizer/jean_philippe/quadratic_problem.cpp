#include "QpGenData.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "GondzioSolver.h"
#include "QpGenSparseMa27.h"
#include "quadratic_problem.h"


Quadratic_Problem::Quadratic_Problem(int _individuals, int _variables, double _lambda, vector<double> &_bounds, SparseMat<double> &_mu, SparseMat<double> &_targets, vector<double> & _shift)
    : individuals (_individuals),
      variables (_variables),
      lambda (_lambda),
      bounds (_bounds),
      mu (_mu),
      targets (_targets),
	  shift (_shift)
{

}


Quadratic_Problem::~Quadratic_Problem()
{
}


void Quadratic_Problem::allocate_objective_function(int** p_rowQ, int** p_colQ, double** p_dQ, double** p_c)
{
    // Allocates the quadratic term
    *p_rowQ = new int[variables + 1];
    *p_colQ = new int[individuals];
    *p_dQ = new double[individuals];
    int* rowQ = *p_rowQ;
    int* colQ = *p_colQ;
    double* dQ = *p_dQ;
	double M = bounds[0];

    for (int i = 0; i <= variables; i++) {
        rowQ[i] = (i < individuals ? i : individuals);
        if (i < individuals) {
            colQ[i] = i;
            //dQ[i] = 2 * (1 - lambda) / (M * M);
			dQ[i] = 100000 * 2 * (1 - lambda) / (M * M * individuals);
        }
    }

	// Allocates the linear term of the function
    *p_c = new double[variables];
    double* c = *p_c;
    for (int i = 0; i < variables; i++) {
		//double c_i = (i < individuals ? -shift[i] * dQ[i] : lambda);
		double c_i = 100000 * (i < individuals ? 0 : lambda / (4 * M * (variables - individuals)));
        c[i] = c_i;
    }
}


void Quadratic_Problem::allocate_bounds(double** p_xlow, char** p_ixlow, double** p_xupp, char** p_ixupp)
{
    // Defines the bounds for every variable
    *p_xlow = new double[variables];
    *p_xupp = new double[variables];
    *p_ixlow = new char[variables];
    *p_ixupp = new char[variables];
    double* xlow = *p_xlow;
    double* xupp = *p_xupp;
    char* ixlow = *p_ixlow;
    char* ixupp = *p_ixupp;
    for (int i = 0; i < variables; i++) {
        if (i < individuals) {
            xlow[i] = -bounds[i];
            xupp[i] = bounds[i];
            ixlow[i] = 1;
            ixupp[i] = 1;
        }
        else {
            xlow[i] = 0;
            xupp[i] = 0;
            ixlow[i] = 0;
            ixupp[i] = 0;
        }
    }
}


void Quadratic_Problem::allocate_inequality_constraints(int** p_rowC, int** p_colC, double** p_dC, double** p_clow, char** p_iclow, double** p_cupp, char** p_icupp)
{
    // Allocates the sparse matrix that corresponds to the constraints
    *p_rowC = new int[2 * mu.nonZeros() + 1];
    *p_colC = new int[6 * mu.nonZeros()];
    *p_dC = new double[6 * mu.nonZeros()];

    int* rowC = *p_rowC;
    int* colC = *p_colC;
    double* dC = *p_dC;
    for (int i = 0; i <= 2 * mu.nonZeros(); i++) rowC[i] = 3 * i;

    int p = 0;
    for (int k = 0; k < mu.outerSize(); k++) {
        for (SparseMat<double>::InnerIterator it(mu, k); it; ++it) {
            int i = it.row(), j = it.col();
            double mu_ij = it.value();
            colC[6 * p + 0] = i;
            colC[6 * p + 1] = j;
            colC[6 * p + 2] = individuals + p;
            colC[6 * p + 3] = i;
            colC[6 * p + 4] = j;
            colC[6 * p + 5] = individuals + p;
            dC[6 * p + 0] = -2 * mu_ij;
            dC[6 * p + 1] = 2 * mu_ij;
            dC[6 * p + 2] = -1;
            dC[6 * p + 3] = 2 * mu_ij;
            dC[6 * p + 4] = -2 * mu_ij;
            dC[6 * p + 5] = -1;
            p++;
        }
    }
    // Defines clow <= Cx <= cupp
    *p_clow = new double[2 * mu.nonZeros()];
    *p_cupp = new double[2 * mu.nonZeros()];
    *p_iclow = new char[2 * mu.nonZeros()];
    *p_icupp = new char[2 * mu.nonZeros()];

    double* clow = *p_clow;
    double* cupp = *p_cupp;
    char* iclow = *p_iclow;
    char* icupp = *p_icupp;
    p = 0;
    for (int k = 0; k < mu.outerSize(); k++) {
        SparseMat<double>::InnerIterator it_mu(mu, k);
        SparseMat<double>::InnerIterator it_targets(targets, k);
        while (it_mu && it_targets) {
            iclow[2 * p + 0] = 0;
            iclow[2 * p + 1] = 0;
            clow[2 * p + 0] = 0;
            clow[2 * p + 1] = 0;
            icupp[2 * p + 0] = 1;
            icupp[2 * p + 1] = 1;
            cupp[2 * p + 0] = -2 * it_mu.value() * it_targets.value();
            cupp[2 * p + 1] = 2 * it_mu.value() * it_targets.value();
            ++it_mu;
            ++it_targets;
            p++;
        }
    }
}


void Quadratic_Problem::desallocate_objective_function(int* rowQ, int* colQ, double* dQ, double* c)
{
    delete[] rowQ;
    delete[] colQ;
    delete[] dQ;
    delete[] c;
}


void Quadratic_Problem::desallocate_bounds(double* xlow, char* ixlow, double* xupp, char* ixupp)
{
    delete[] xlow;
    delete[] ixlow;
    delete[] xupp;
    delete[] ixupp;
}


void Quadratic_Problem::desallocate_inequality_constraints(int* rowC, int* colC, double* dC, double* clow, char* iclow, double* cupp, char* icupp)
{
    delete[] rowC;
    delete[] colC;
    delete[] dC;
    delete[] clow;
    delete[] iclow;
    delete[] cupp;
    delete[] icupp;
}


bool Quadratic_Problem::solve(vector<double> &solution)
{
    int *krowQ, *jcolQ;
    double *c, *dQ;
    allocate_objective_function(&krowQ, &jcolQ, &dQ, &c);

    double *xlow, *xupp;
    char *ixlow, *ixupp;
    allocate_bounds(&xlow, &ixlow, &xupp, &ixupp);

    int *krowC, *jcolC;
    double *dC, *clow, *cupp;
    char *iclow, *icupp;
    allocate_inequality_constraints(&krowC, &jcolC, &dC, &clow, &iclow, &cupp, &icupp);

    int krowA = 0;

    /*std::cout << "** Running solver HSL-MA27" << std::endl;
    std::cout << "Individuals : " << individuals << std::endl;
    std::cout << "Potentials : " << mu.nonZeros() << std::endl;
    std::cout << "Variables : " << variables << std::endl;*/

    QpGenSparseMa27* qp = new QpGenSparseMa27(variables, 0, 2 * mu.nonZeros(), individuals, 0, 6 * mu.nonZeros());
    QpGenData* prob = (QpGenData *)qp->makeData(c, krowQ, jcolQ, dQ, xlow, ixlow, xupp, ixupp, &krowA, NULL, NULL, NULL, krowC, jcolC, dC, clow, iclow, cupp, icupp);
    QpGenVars* vars = (QpGenVars *)qp->makeVariables(prob);
    QpGenResiduals* resid = (QpGenResiduals *)qp->makeResiduals(prob);
    GondzioSolver *solver = new GondzioSolver(qp, prob);

    int status = solver->solve(prob, vars, resid);
    if (status == 0) {
        solution = vector<double>(variables, 0);
        vars->x->copyIntoArray(solution.data());
    } else {
        std::string m_status;
        switch (status) {
        case 1: m_status = "NOT_FINISHED"; break;
        case 2: m_status = "MAX_ITS_EXCEEDED"; break;
        case 3: m_status = "INFEASIBLE"; break;
        case 4: m_status = "UNKNOWN"; break;
        default: m_status = "OTHER_STATUS"; break;
        }
        std::cout << "WARNING : solver status is " << m_status << " (code " << status << ")" << std::endl;
    }
    delete solver;
    delete resid;
    delete vars;
    delete prob;
    delete qp;

    desallocate_inequality_constraints(krowC, jcolC, dC, clow, iclow, cupp, icupp);
    desallocate_bounds(xlow, ixlow, xupp, ixupp);
    desallocate_objective_function(krowQ, jcolQ, dQ, c);

    return (status == 0);
}
