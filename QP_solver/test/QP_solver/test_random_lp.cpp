#include <cstdlib>
#include <cassert>
#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/Random.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_options.h>
#include <CGAL/QP_functions.h>
#include <CGAL/Timer.h>

// choose exact floating point type
#ifndef CGAL_USE_GMP
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#else
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
#endif

// program and solution types
typedef CGAL::Quadratic_program<double> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

// random number generator
CGAL::Random rd;

// timer
CGAL::Timer timer;

CGAL::Comparison_result random_rel()
{
  int z = rd.get_int(-1,2);
  return CGAL::Comparison_result(z);
}

void statistics (const Solution& s, 
		 unsigned int& o, unsigned int& i, unsigned int& u)
{
    switch (s.status()) {
    case CGAL::QP_OPTIMAL:
      o++;
      break;
    case CGAL::QP_INFEASIBLE:
      i++;
      break;
    case CGAL::QP_UNBOUNDED:
      u++;
      break;
    default:
      CGAL_qpe_assertion(false);
    }
}

unsigned int qp_optimal = 0;
unsigned int qp_infeasible = 0;
unsigned int qp_unbounded = 0;
unsigned int nqp_optimal = 0;
unsigned int nqp_infeasible = 0;
unsigned int nqp_unbounded = 0;
unsigned int lp_optimal = 0;
unsigned int lp_infeasible = 0;
unsigned int lp_unbounded = 0;
unsigned int nlp_optimal = 0;
unsigned int nlp_infeasible = 0;
unsigned int nlp_unbounded = 0;

// parameters
int tries = 10;
int max_dim = 100; // must be >1
// int max_entry = 110; // must be >0

int main(const int argNr,const char **args) {  

  // obtain and print target density
  const double target_density = argNr > 1 ? std::atof(args[1]) : 0.5; // default target density is 0.5
  std::cout << "Target density: " << target_density << std::endl;  
  
  // choose dimensions
  int n = argNr > 2 ? std::atoi(args[2]) : rd.get_int(1,max_dim);
  int m = argNr > 3 ? std::atoi(args[3]) : rd.get_int(1,max_dim);
  
  // print seed
  std::cout << "Random seed: " << rd.get_seed() << std::endl;
  
  // counter for total time
  double total_time(0);

  // options
  CGAL::Quadratic_program_options options;
  options.set_auto_validation(true);

  // generate a set of small random qp's
  for (int i=0; i<tries; ++i) {

    // construct matrix D as C^T C, for C randomly chosen with n columns
    /*
    int k = rd.get_int (1, 2*n); // number of rows of C
    std::vector<std::vector<int> > C (k, std::vector<int>(n, 0));
    for (int j=0; j<k+n; ++j)  // sparse C
      C[rd.get_int(0, k)][rd.get_int(0,n)] = 
	rd.get_int(-max_entry, max_entry);
*/

    // now fill the program 
    Program p;
    
    
    // fill A according to target density
    std::vector<std::pair<int,int> > index_pairs;
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < m; ++k) {
        index_pairs.push_back(std::make_pair(j,k));
      }
    }
    for (int j = 0; j < static_cast<int>(target_density*m*n); ++j) {
      int index = rd.get_int(0,index_pairs.size());
      p.set_a (index_pairs[index].first, index_pairs[index].second, rd.get_double());
      index_pairs.erase(index_pairs.begin()+index);
    }
    
    
    
    
    // b, r
    for (int j = 0; j < m/2; ++j) {
      p.set_b (rd.get_int(0,m), rd.get_double());
      p.set_r (rd.get_int(0,m), random_rel());
    }
    // fl, l, fu, u
    for (int j=0; j<n; ++j) {
      double l = rd.get_double();
      double u = rd.get_double();
      if (l > u) std::swap (l, u); 
      p.set_l(j, rd.get_bool(), l);
      p.set_u(j, rd.get_bool(), u);
    }
    // D
    /*
    for (int i=0; i<n; ++i)
      for (int j=0; j<=i; ++j) {
        double entry = 0;
        for (int l=0; l<k; ++l) 
          entry += C[l][i] * C[l][j];
        p.set_d(i, j, entry);
      }
      */
    // c
    for (int j=0; j<n/2; ++j)
      p.set_c (rd.get_int(0, n), rd.get_double());
    // c0
    p.set_c0(rd.get_double());
    
    
    // write out and read back to test equality
    /*
    std::stringstream inout;
    CGAL::print_quadratic_program (inout, p);
    CGAL::Quadratic_program_from_mps<double> p2 (inout);
    assert(CGAL::QP_functions_detail::are_equal_qp (p, p2));
    */
    
    // TAG: DEBUG
    std::cout << "Trial nr: " << i << std::endl;
    std::cout << "n: " << n << ", m: " << m << std::endl;
    std::cout << "Number of elements set: "<< static_cast<int>(target_density*m*n) << std::endl;
    
    bool print = false;
    
    /*
    // solve it
    Solution s;
    timer.reset();
    timer.start();
    s = CGAL::solve_linear_program (p, ET(), options);    
    timer.stop();
    assert(s.is_valid()); 
    statistics (s, lp_optimal, lp_infeasible, lp_unbounded);
    std::cout << (s.status() == CGAL::QP_OPTIMAL ? "OPTIMAL\n" : "INFEASIBLE or UNBOUNDED\n");
    if (s.status() == CGAL::QP_OPTIMAL) print = true;
    std::cout << "Used " << timer.time() << " seconds (solve_linear_program)." << std::endl;
    total_time += timer.time();
    
    timer.reset();
    timer.start();
    s = CGAL::solve_nonnegative_linear_program (p, ET(), options);
    timer.stop();
    //CGAL::print_quadratic_program (std::cout, p);
    assert(s.is_valid());  
    statistics (s, nlp_optimal, nlp_infeasible, nlp_unbounded);
    std::cout << (s.status() == CGAL::QP_OPTIMAL ? "OPTIMAL\n" : "INFEASIBLE or UNBOUNDED\n");
    if (s.status() == CGAL::QP_OPTIMAL) print = true;
    std::cout << "Used " << timer.time() << " seconds (solve_nonnegative_linear_program)." << std::endl;
    total_time += timer.time();
    */ 
    
    //if (print) CGAL::print_quadratic_program (std::cout, p);
    if (true) CGAL::print_quadratic_program (std::cout, p);
  }
  
  // output statistics
  std::cout << "Solved " << tries 
	    << " random LP / NLP .\n"
	    << " Optimal:    " 
	    << lp_optimal << " / " 
	    << nlp_optimal << "\n"
	    << " Infeasible: "
	    << lp_infeasible << " / " 
	    << nlp_infeasible << "\n"
	    << " Unbounded:  "
	    << lp_unbounded << " / " 
	    << nlp_unbounded << std::endl;
      
  // output timer info
  std::cout << "Used " << total_time << " seconds in total." << std::endl;
  std::cout.flush();

  return 0;
}
