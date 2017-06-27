#include"typedefs.h"
#include"hexextr.h"
#include<boost/numeric/ublas/matrix.hpp>

using matrix = boost::numeric::ublas::matrix;

int find_number_of_boundary_vertices(LCC_3& lcc){
  int count = 0;
  for(LCC_3::One_dart_per_cell_range<0,0>::iterator it = lcc.one_dart_per_cell<0,0>().begin(), itend = lcc.one_dart_per_cell<0,0>().end(); it!= itend; it++){
    if(is_free(it, 3)) count++;
  }
  return count;
}

void optimise_frame_field(HexExtr& h, int n){ // n is the number of smoothing iterations
  int nl = find_number_of_boundary_vertices(h.input_tet_mesh);
  int nv = (h.vertices).length(); 
  vector<vector<double>> a;
  for(int i = 0; i < n; i++){
    matrix<double> A(0, (9*nv + 2*nl + 3*nv));//vector<vector<double>> A;
    vector<double> b;
    add_smoothing_terms(h, A, b);
    add_normal_constraints(h, A, b);
    if(i>0){
      add_local_optim_constraints(h, a, A, b);
    }
    X = call_least_squares_solver(A, b);
    for(int j = 0;j<nv;j++){
      a[j] = X[9i...9i+8];
    }
  }
}
