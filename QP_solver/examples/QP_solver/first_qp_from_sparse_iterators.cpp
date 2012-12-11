
// example: construct a quadratic program from sparse iterator data
#include <iostream>
#include <cassert>
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>


// choose exact integral type
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif

// program and solution types
typedef typename std::vector<std::pair<int, int> > Col_sparse;
typedef CGAL::Quadratic_program_from_sparse_iterators
<Col_sparse*,                                         // for A
int*,                                                 // for b
CGAL::Const_oneset_iterator<CGAL::Comparison_result>, // for r
bool*,                                                // for fl
int*,                                                 // for l
bool*,                                                // for fu
int*,                                                 // for u
Col_sparse*,                                          // for D
int*>                                                 // for c
Program_from_sparse_iterators;


typedef CGAL::Quadratic_program_solution<ET> Solution;
typedef CGAL::Quadratic_program_options Options;



int main(const int argNr, const char **args) {
    
	// get desired level of additional logging output:
    const int verbosity = argNr < 2 ? 0 : std::atoi(args[1]);
    
    
    Options  options;
    options.set_verbosity(verbosity);
    
    Col_sparse A[] = {Col_sparse(), Col_sparse()};
    Col_sparse D[] = {Col_sparse(), Col_sparse()};
    A[0].push_back(std::make_pair(0,1));    // A_{1,1} == 1
    A[0].push_back(std::make_pair(1,-11));  // A_{2,1} == -11
    A[1].push_back(std::make_pair(0,11));   // A_{1,2} == 11
    A[1].push_back(std::make_pair(1,2));    // A_{2,2} == 2
    
    D[0].push_back(std::make_pair(0,2)); // 2D_{1,1} == 2
    D[1].push_back(std::make_pair(1,8)); // 2D_{2,2} == 8

    int   b[] = {7, 4};                         // right-hand side
    CGAL::Const_oneset_iterator<CGAL::Comparison_result>
    r(    CGAL::SMALLER);                 // constraints are "<="
    bool fl[] = {true, true};                   // both x, y are lower-bounded
    int   l[] = {0, 0};
    bool fu[] = {false, true};                  // only y is upper-bounded
    int   u[] = {0, 4};                         // x's u-entry is ignored
    int   c[] = {0, -32};
    int  c0   = 64;                             // constant term
    
    // now construct the quadratic program; the first two parameters are
    // the number of variables and the number of constraints (rows of A)
    Program_from_sparse_iterators qp (2, 2, A, b, r, fl, l, fu, u, D, c, c0);
    CGAL::print_quadratic_program (std::cout, qp);
    
    // solve the program, using ET as the exact type
    Solution s = CGAL::solve_quadratic_program(qp, ET(), options);
    
    // output solution
    std::cout << s;

    return 0;
}