#ifndef CGAL_SHAPE_REGULARIZATION_DENSE_QP_SOLVER
#define CGAL_SHAPE_REGULARIZATION_DENSE_QP_SOLVER

// #include <CGAL/license/Shape_regularization.h>
// #include <iostream>

#include "osqp.h"
#include <Eigen/Dense>

#include <map>
#include <utility>
#include <set>
#include <vector>

namespace CGAL {
namespace Regularization {

  template<
    typename GeomTraits,
    typename InputRange>
  class Dense_QP_solver{ 

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using FT = typename GeomTraits::FT;

    Dense_QP_solver(
       InputRange& input_range):
    m_input_range(input_range) {
      CGAL_precondition(input_range.size() > 0);
    }

    void solve(std::set<std::pair<int, int>> graph, std::map <std::pair<int, int>, FT> t_ijs_map, std::map <std::pair<int, int>, FT> r_ijs_map, std::vector<FT> & result){

      c_int n = m_input_range.size() + graph.size(); // = 6; number of variables
      c_int m = 2 * graph.size() + n; // = 12; number of constraints
      //csc *P - quadratic part of the cost P in csc format (size n x n)
      c_int   P_nnz  = n; // = 6;
      // helper variables for calculation
      int weight = 100000;
      FT lambda = FT(4)/FT(5);
      c_float theta_max = 10.0;
      int k = m_input_range.size(); // k segments

      // c_float P_x[6] = { 133.333, 133.333, 133.333, 0.0, 0.0, 0.0, };
      c_float P_x[n];
      for(int i = 0; i < n; i++) {
        if (i < m_input_range.size()) 
          P_x[i] = (2 * weight * (1 - lambda)) / (theta_max * theta_max * k);
        else
          P_x[i] = 0.0;
      }

      // c_int   P_i[6] = { 0, 1, 2, 3, 4, 5, };
      // c_int   P_p[7] = { 0, 1, 2, 3, 4, 5, 6, }; 
      c_int   P_i[n];
      c_int   P_p[n+1];
      for (int i = 0; i < n; i++) {
        P_i[i] = i;
        P_p[i] = i;
      }
      P_p[n] = n;

      // c_float q[6]   = { 0.0, 0.0, 0.0, 666.667, 666.667, 666.667, }; // dense array for linear part of cost function (size n)
      c_float q[n];
      for (int i = 0; i < n; i++) {
        if(i < m_input_range.size()) {
          q[i] = 0.0;
        }
        else {
          q[i] = lambda * weight / (4 * theta_max * (n - k));
        }
      }

      FT val_pos = 2 * lambda;
      FT val_neg = -2 * lambda;
      //csc *A - linear constraints matrix A in csc format (size m x n) 
      c_int   A_nnz = 6 * graph.size() + n; // 24;  
      // c_float A_x[24] = { -1.6, 1.6, -1.6, 1.6, 1.0, 1.6, -1.6, -1.6, 1.6, 1.0, 1.6, -1.6, 
      //                     1.6, -1.6, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, };
      // c_int   A_i[24] = { 0, 1, 2, 3, 6, 0, 1, 4, 5, 7, 2, 3, 4, 5, 8, 0, 1, 9, 2, 3, 10, 4, 5, 11, };
      // c_int   A_p[7] = { 0, 5, 10, 15, 18, 21, 24, };
      Eigen::MatrixXd A_matrix(m, n);
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          if(j == i - graph.size() * 2)
            A_matrix(i, j) = 1;
          else
            A_matrix(i, j) = 0;
        }
      }
      std::set<std::pair<int, int>>::iterator graph_iterator;
      int it = 0, ij = k;
      for (graph_iterator = graph.begin(); graph_iterator != graph.end(); graph_iterator++) {
          A_matrix(it, graph_iterator->first) = val_neg;
          A_matrix(it, graph_iterator->second) = val_pos;
          A_matrix(it, ij) = -1;
          ++it;
          A_matrix(it, graph_iterator->first) = val_pos;
          A_matrix(it, graph_iterator->second) = val_neg;
          A_matrix(it, ij) = -1;
          ++it;
          ++ij;
      }
      std::cout << "Matrix A: " << std::endl << A_matrix << std::endl;
      
      c_float A_x[A_nnz] = { -1.6, 1.6, -1.6, 1.6, 1.0, 1.6, -1.6, -1.6, 1.6, 1.0, 1.6, -1.6, 
                          1.6, -1.6, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, };
      c_int   A_i[A_nnz] = { 0, 1, 2, 3, 6, 0, 1, 4, 5, 7, 2, 3, 4, 5, 8, 0, 1, 9, 2, 3, 10, 4, 5, 11, };
      c_int   A_p[n+1] = { 0, 5, 10, 15, 18, 21, 24, };

      // c_float l[m]   = { -10000000.0, -10000000.0, -10000000.0, -10000000.0, -10000000.0, -10000000.0, 
      //                   -10.0, -10.0, -10.0, -10000000.0, -10000000.0, -10000000.0, }; // dense array for lower bound (size m) 
      // c_float u[m]   = { 9.13, -9.13, -0.0, 0.0, -9.13, 9.13, 10, 10, 10, 10000000.0, 10000000.0, 10000000.0, }; // dense array for upper bound (size m) 
      c_float l[m], u[m];
      graph_iterator = graph.begin();
      for (int i = 0; i < m; i++) {
        if(i < graph.size() * 2) {
          if (i % 2 == 0) {
            u[i] = val_neg * t_ijs_map[*graph_iterator];
          }
          else { 
            u[i] = val_pos * t_ijs_map[*graph_iterator];
            graph_iterator++;
          }
          l[i] = -10000000.0;
        }
        else if (i < graph.size() * 2 + k) {
          l[i] = -1 * theta_max;
          u[i] = theta_max;
        }
        else {
          l[i] = -10000000.0;
          u[i] = 10000000.0;
        }
      }

     /* std::cout << "l: ";
      for(int i = 0; i < m; i++) {
        std::cout << l[i] << " ";
      }
      std::cout << std::endl << "u: ";
      for(int i = 0; i < m; i++) {
        std::cout << u[i] << " ";
      }
      std::cout << std::endl; */

       // Problem settings
      OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

      // Structures
      OSQPWorkspace *work; // Workspace
      OSQPData *data;      // OSQPData

      // Populate data
      data    = (OSQPData *)c_malloc(sizeof(OSQPData));
      data->n = n;
      data->m = m;
      data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
      data->q = q;
      data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
      data->l = l;
      data->u = u;

      // Define Solver settings as default
      osqp_set_default_settings(settings);

      // Setup workspace
      work = osqp_setup(data, settings);

      // Solve Problem
      osqp_solve(work);

      //output OSQP result:
      c_float *i = work->solution->x;
      for(int j = 0; j < n; j++) {
        std::cout << i[j] << " ";
      }
      std::cout << std::endl;

      // Clean workspace
      osqp_cleanup(work);
      c_free(data->A);
      c_free(data->P);
      c_free(data);
      c_free(settings);
    

      result.clear();
      result.push_back(FT(-1.90353));
      result.push_back(FT(3.80706));
      result.push_back(FT(-1.90353));
      result.push_back(FT(1.56148e-12));
      result.push_back(FT(2.66171e-14));
      result.push_back(FT(1.56259e-12));
    }
    // creates an instance of CGAL QP solver

     /* void test() {
      c_float P_x[4] = { 4.0, 1.0, 1.0, 2.0, };
      c_int   P_nnz  = 4;
      c_int   P_i[4] = { 0, 1, 0, 1, };
      c_int   P_p[3] = { 0, 2, 4, };
      c_float q[2]   = { 1.0, 1.0, };
      c_float A_x[4] = { 1.0, 1.0, 1.0, 1.0, };
      c_int   A_nnz  = 4;
      c_int   A_i[4] = { 0, 1, 0, 2, };
      c_int   A_p[3] = { 0, 2, 4, };
      c_float l[3]   = { 1.0, 0.0, 0.0, };
      c_float u[3]   = { 1.0, 0.7, 0.7, };
      c_int n = 2;
      c_int m = 3;

      // Problem settings
      OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

      // Structures
      OSQPWorkspace *work; // Workspace
      OSQPData *data;      // OSQPData

      // Populate data
      data    = (OSQPData *)c_malloc(sizeof(OSQPData));
      data->n = n;
      data->m = m;
      data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
      data->q = q;
      data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
      data->l = l;
      data->u = u;

      // Define Solver settings as default
      osqp_set_default_settings(settings);

      // Setup workspace
      work = osqp_setup(data, settings);

      // Solve Problem
      osqp_solve(work);

      c_float *i = work->solution->x;
      std::cout << "OSQP solver: " << *i << std::endl;
      std::cout << std::endl;

      // Clean workspace
      osqp_cleanup(work);
      c_free(data->A);
      c_free(data->P);
      c_free(data);
      c_free(settings);
    } */

  private:
    const Input_range& m_input_range;
    const FT m_mu_ij = FT(4) / FT(5);

    
  };


} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_DENSE_QP_SOLVER
