#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <iostream>
#include <random> // Requires C++ 11

#include <chrono>
#include <time.h>

#include <vector>
#include <algorithm>

#include <SymGEigsSolver.h>
#include <MatOp/SparseSymMatProd.h>
#include <MatOp/SparseCholesky.h>

using namespace Spectra;

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;
typedef Eigen::SparseMatrix<double> SpMatrix;

typedef std::chrono::high_resolution_clock myclock;
typedef std::chrono::microseconds TimeT;

typedef SparseSymMatProd<double> OpType;
typedef SparseCholesky<double> BOpType;


// Generate random sparse matrix
SpMatrix sprand(int size, int elem = 18, int seed = 0)
{
    SpMatrix mat(size, size);
    Vector reserve_idx = Vector::Constant(size, elem);   
    mat.reserve(reserve_idx);
    std::default_random_engine gen(seed);
    //gen.seed(seed);
    std::uniform_real_distribution<double> distr(-1.0, 1.0);

    std::mt19937 mersenne_engine(seed);
    std::uniform_int_distribution<int> dist(0, size);
    auto int_gen = std::bind(dist, mersenne_engine);

    for(int i = 0; i < size; i++){
        std::vector<int> v(elem);
        std::generate(std::begin(v), std::end(v), int_gen);

        for(int j = 0; j < elem; j++){
            mat.coeffRef(i, v[j]) = distr(gen);
        }    
    }

    return mat;
}


void gen_sparse_data(int n, SpMatrix& A, SpMatrix& B, int elem = 18, int seed = 0)
{
    // Eigen solver only uses the lower triangle of A,
    // so we don't need to make A symmetric here.
    A = sprand(n, elem, seed);
    B = A.transpose() * A;
    // To make sure B is positive definite
    for(int i = 0; i < n; i++)
        B.coeffRef(i, i) += 1.0 / elem; 
}

template <typename MatType, int SelectionRule>
unsigned int run_test(const MatType& A, const MatType& B, int k, int m)
{
    OpType op(A);
    BOpType Bop(B);
    // Make sure B is positive definite and the decompoition is successful
    assert(Bop.info() == SUCCESSFUL);

    SymGEigsSolver<double, SelectionRule, OpType, BOpType, GEIGS_CHOLESKY> eigs(&op, &Bop, k, m);
    auto start = std::chrono::steady_clock::now();
    eigs.init();
    int nconv = eigs.compute(100); // maxit = 100 to reduce running time for failed cases
    auto duration = std::chrono::duration_cast<TimeT> (std::chrono::steady_clock::now() - start);

    int niter = eigs.num_iterations();
    int nops = eigs.num_operations();

    assert(eigs.info() == SUCCESSFUL);

    Vector evals = eigs.eigenvalues();
    Matrix evecs = eigs.eigenvectors();

    Matrix resid = A.template selfadjointView<Eigen::Lower>() * evecs -
                   B.template selfadjointView<Eigen::Lower>() * evecs * evals.asDiagonal();
    const double err = resid.array().abs().maxCoeff();

    //std::cout << err << std::endl;
    assert(err < 1e-9);

    return duration.count();
}

void usage()
{
    std::cout << "Usage: ./test size k m elem seed num" << std::endl;
    std::cout << "size: an integer indicating the size of the matrix A" << std::endl;
    std::cout << "k   : an integer indicating the number of eigenvectors to be calculated" << std::endl;
    std::cout << "m   : an integer indicating the ncv (decide the converge speed)" << std::endl;
    std::cout << "elem: non-zero entries in each line" << std::endl;
    std::cout << "seed: a bool flag. If 0, the seed is fixed (0). Otherwise, it is generated randomly" << std::endl;
    std::cout << "num : the number of experiments" << std::endl;
}



int main(int argc, char** argv)
{
    // Check number of parameters
    if(argc != 7){
        usage();
        return EXIT_FAILURE;
    }

    // Initialize parameters
    int size = atoi(argv[1]), k = atoi(argv[2]), m = atoi(argv[3]), elem = atoi(argv[4]), num = atoi(argv[6]);
    bool seed = false;
    if(atoi(argv[5]) != 0) seed = true;  

    double avg_time = 0;
    double coeff = 1.0 / num;  

    for(int i = 0; i < num; i++){
        SpMatrix A, B;
        int curr_seed = 0;
        if(seed){
            curr_seed = static_cast<int>(time(NULL));
        }
        gen_sparse_data(size, A, B, elem, seed);
        unsigned int curr_time = run_test<SpMatrix, LARGEST_ALGE>(A, B, k, m);
        std::cout << curr_time << std::endl;
        avg_time += coeff * static_cast<double>(curr_time);
    }

    std::cout << "Size of matrix : " << size << std::endl;
    std::cout << "Number of tests: " << num << std::endl;
    std::cout << "Average time   : " << avg_time << " µs" << std::endl;
    
    return EXIT_SUCCESS;     
}
