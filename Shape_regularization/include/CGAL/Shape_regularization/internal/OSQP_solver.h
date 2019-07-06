#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_OSQP_SOLVER
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_OSQP_SOLVER

// #include <CGAL/license/Shape_regularization.h>
// #include <iostream>

#include "osqp.h"
#include <Eigen/Sparse>

#include <map>
#include <utility>
#include <set>
#include <vector>

namespace CGAL {
namespace Regularization {
namespace internal {

  template<
    typename GeomTraits,
    typename InputRange>
  class OSQP_solver { 

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using FT = typename GeomTraits::FT;
    using FT_triplet = Eigen::Triplet<FT>;

    OSQP_solver(
       InputRange& input_range):
    m_input_range(input_range) {
      CGAL_precondition(input_range.size() > 0);
    }

    void solve(std::set<std::pair<std::size_t, std::size_t>> graph, std::map <std::pair<std::size_t, std::size_t>, FT> t_ijs_map, std::vector<FT> & result){

      c_int n = m_input_range.size() + graph.size(); // = 6; number of variables
      c_int m = 2 * graph.size() + n; // = 12; number of constraints
      //csc *P - quadratic part of the cost P in csc format (size n x n)
      c_int   P_nnz  = n; // = 6;
      // helper variables for calculation
      int weight = 100000;
      FT lambda = FT(4)/FT(5);
      c_float theta_max = 10.0;
      int k = m_input_range.size(); // k segments

      c_float P_x[n];
      for(int i = 0; i < n; i++) {
        if (i < m_input_range.size()) 
          P_x[i] = (2 * weight * (1 - lambda)) / (theta_max * theta_max * k);
        else
          P_x[i] = 0.0;
      }

      c_int   P_i[n];
      c_int   P_p[n+1];
      for (int i = 0; i < n; i++) {
        P_i[i] = i;
        P_p[i] = i;
      }
      P_p[n] = n;

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
      c_int   A_nnz = 6 * graph.size() + n;  
      
      std::vector<FT_triplet> vec;
      vec.reserve(A_nnz);
      std::set<std::pair<std::size_t, std::size_t>>::iterator graph_iterator;
      int it = 0, ij = k;
      for (graph_iterator = graph.begin(); graph_iterator != graph.end(); graph_iterator++) {
        vec.push_back(FT_triplet(it, graph_iterator->first, val_neg));
        vec.push_back(FT_triplet(it, graph_iterator->second, val_pos));
        vec.push_back(FT_triplet(it, ij, -1));
        ++it;
        vec.push_back(FT_triplet(it, graph_iterator->first, val_pos));
        vec.push_back(FT_triplet(it, graph_iterator->second, val_neg));
        vec.push_back(FT_triplet(it, ij, -1));
        ++it;
        ++ij;
      }
      for (int i = graph.size() * 2; i < m; i++) {
        for (int j = 0; j < n; j++) {
          if(j == i - graph.size() * 2)
            vec.push_back(FT_triplet(i, j, 1));
        }
      }

      Eigen::SparseMatrix<FT, Eigen::ColMajor> A_mat;
      A_mat.resize(m, n);
      A_mat.setFromTriplets(vec.begin(), vec.end());
      A_mat.makeCompressed();
      
      c_float A_x[A_nnz];
      it = 0;
      for (int i = 0; i < A_mat.outerSize(); ++i) {
        for (typename Eigen::SparseMatrix<FT, Eigen::ColMajor>::InnerIterator m_i(A_mat, i); m_i; ++m_i) {
          const double val = CGAL::to_double(m_i.value());
          A_x[it] = val;
          ++it;
        }
      }
      c_int A_i[A_nnz];
      it = 0;
      for (int i = 0; i < A_mat.outerSize(); ++i) {
        for (typename Eigen::SparseMatrix<FT, Eigen::ColMajor>::InnerIterator m_i(A_mat, i); m_i; ++m_i) {
          const int idx = m_i.row();
          A_i[it] = idx;
          ++it;
        }
      }
      c_int A_p[n+1];
      A_p[0] = 0;
      for (std::size_t i = 1; i <= A_mat.outerSize(); ++i) {
        const std::size_t coln = A_mat.innerVector(i-1).nonZeros();
        A_p[i] = A_p[i-1] + coln;
      }

            
      c_float l[m], u[m]; // dense arrays for lower and upper bounds (size m) 
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
      settings->eps_rel = 1.0e-15;

      // Setup workspace
      work = osqp_setup(data, settings);

      // Solve Problem
      osqp_solve(work);

      //output OSQP result:
      // c_float *i = work->solution->x;
      // for(int j = 0; j < n; j++) {
      //   std::cout << i[j] << " ";
      // }
      // std::cout << std::endl;

      result.clear();
      c_float *i = work->solution->x;
      for(int j = 0; j < n; j++) {
        result.push_back(FT(i[j]));
      }

      // Clean workspace
      osqp_cleanup(work);
      c_free(data->A);
      c_free(data->P);
      c_free(data);
      c_free(settings);
    
    }

  private:
    const Input_range& m_input_range;
    const FT m_mu_ij = FT(4) / FT(5);

    
  };


} // namespace internal
} // namespace Regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_OSQP_SOLVER
