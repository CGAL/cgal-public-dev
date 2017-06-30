#include"typedefs.h"
#include"hexextr.h"
#include<boost/numeric/ublas/matrix.hpp> 
#include<boost/numeric/ublas/vector.hpp> //is this needed?
#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>
#include <CGAL/Eigen_svd.h>
#include<cmath>
#define PI 3.14159265
typedef CGAL::Eigen_svd Svd;
#endif
typedef Svd::FT     FT;
typedef Svd::Vector Eigen_vector;
typedef Svd::Matrix Eigen_matrix;
using matrix = boost::numeric::ublas::matrix;
using vect = boost::numeric::ublas::vector;


matrix return_Rz(double gamma){
  matrix<double> Rz(9,9);
  for(int i = 0; i<9; i++){
    for(int j = 0; j<9; j++){
     Rz.set(i,j,0);
    }  
  }
  Rz(0,0) = cos(4*gamma); Rz(1,1) = cos(3*gamma); Rz(2,2) = cos(2*gamma); Rz(3,3) = cos(gamma); Rz(4,4) = 1; Rz(5,5) = cos(gamma); Rz(6,6) = cos(2*gamma); Rz(7,7) = cos(3*gamma); Rz(8,8) = cos(4*gamma);
  Rz(0,8) = sin(4*gamma); Rz(1,7) = sin(3*gamma); Rz(2,6) = sin(2*gamma); Rz(3,5) = sin(gamma); Rz(5,3) = -sin(gamma); Rz(6,2) = -sin(2*gamma); Rz(7,1) = -sin(3*gamma); Rz(8,0) = -sin(4*gamma);
  
  return Rz;
}


matrix return_Ry(double beta){
  matrix<double> Ry(9,9);
  matrix Rx_90(9,9);
  for(int i = 0; i<9; i++){
    for(int j = 0; j<9; j++){
     Rx_90.set(i,j,0);
    }  
  }
  Rx_90(0,5) = sqrt(14.0/16); Rx_90(0,7) = -sqrt(1.0/8); Rx_90(1,1) = -0.75; Rx_90(1,3) = -sqrt(7.0/16); Rx_90(2,5) = sqrt(1.0/8); Rx_90(2,7) = sqrt(7.0/8); Rx_90(3,1) = sqrt(7.0/16); Rx_90(3,3) = 0.75; Rx_90(4,4) = 3.0/8; Rx_90(4,6) = sqrt(5.0/16); Rx_90(4,8) = sqrt(35.0/64); Rx_90(5,0) = -sqrt(7.0/8); Rx_90(5,2) = -sqrt(1.0/8); Rx_90(6,4) = sqrt(5/16); Rx_90(6,6) = 0.5; Rx_90(6,8) = -sqrt(7/16); Rx_90(7,0) = sqrt(1.0/8); Rx_90(7,2) = -sqrt(7.0/8); Rx_90(8,4) = sqrt(35.0/64); Rx_90(8,6) = -sqrt(7.0/16); Rx_90(8,8) = 0.125;

  Ry = boost::numeric::ublas::prec_prod(Rx_90, return_Rz(beta));
  Ry = boost::numeric::ublas::prec_prod(Ry, boost::numeric::ublas::trans(Rx_90));
  return Ry;
}


matrix return_Rx(double alpha){
  matrix<double> Rx(9,9), Ry_90(9,9);
  Ry_90 = return_Ry(PI/2);
  Rx = boost::numeric::ublas::prec_prod(boost::numeric::ublas::trans(Ry_90), return_Rz(alpha));
  Rx = boost::numeric::ublas::prec_prod(Rx, Ry_90);
  return Rx;
}


matrix find_rotation_matrix(Vector_3 n){
  double alpha = acos((-1)*n[1]/sqrt(a*a+ b*b)); 
  double beta = acos(c/sqrt(a*a +b*b+c*c));
  double gamma = 0; //TODO: need to check this
  matrix<double> Rz(9,9), Ry(9,9), Rx(9,9);
  Rz = return_Rz(gamma);
  Ry = return_Ry(beta);
  Rx = return_Rx(alpha);
  Rx = boost::numeric::ublas::prec_prod(Rx, Ry);
  Rx = boost::numeric::ublas::prec_prod(Rx, Rz);
  return Rx;
}


void sort_vertices(HexExtr& h, vector<Vertex_handle>& vertices){
  std::sort(vertices.begin(), vertices.end(), comp);
  int i = 0;
  for(std::vector<Vertex_handle>::iterator it  = vertices.begin(), itend = vertices.end(); it!=itend; it++){
    (*it).enumeration = i; i++;
  }
}

int find_number_of_boundary_vertices(LCC_3& lcc){
  int count = 0;
  for(LCC_3::One_dart_per_cell_range<0>::iterator it = lcc.one_dart_per_cell<0>().begin(), itend = lcc.one_dart_per_cell<0>().end(); it!= itend; it++){
    if(lcc.is_free(it, 3)) count++;
  }
  return count;
}

void closest_frame(){

}

void add_smoothing_terms(HexExtr &h, vector<vector<double>>& A, vector<double>& b){
  for(std::vector<Edge_handle>::iterator it =  (h.edges).begin(), itend = (h.edges).end(); it != itend; it++){
    for(int d = 0; d < 9; d++){
      vector<double> row(9*nv+2*nl+3*nv);
      row.fill(row.begin(), row.end(), 0);
      int i = ((*it).from). enumeration, j = ((*it).from). enumeration;
      row[9*i+d] = 1;
      row[9*j+d] = -1;
      A.push_back(row);
      b.push_back(0);
   }
  }
}
  

void add_local_optim_constraints(HexExtr& h, vector<vector<double>>& a, vector<vector<double>& A, vector<double>& b){
  matrix<double> Ex(9,9), Ey(9,9), Ez(9,9);
  
  for(int i = 0; i<9; i++){
    for(int j = 0; j<9; j++){
      Ex.set(i,j,0);  Ey.set(i,j,0); Ez.set(i,j,0);
    }  
  }
  Ex(0,7)= (-1)*sqrt(2)); Ex(1,6) = (-1)*sqrt(3.5); Ex(2,5) = (-1)*sqrt(4.5); Ex(3, 4) = (-1)*sqrt(10);
  Ex(1,8) = (-1)*sqrt(2); Ex(2,7) = (-1)*sqrt(3.5); Ex(3,6) = (-1)*sqrt(4.5); 
  Ex(7,0) = sqrt(2); Ex(6,1) = sqrt(3.5); Ex(5,2) = sqrt(4.5); Ex(4,3) = sqrt(10)); 
  Ex(8,1) = sqrt(2); Ex(7,2) = sqrt(3.5); Ex(6,3) = sqrt(4.5);
  Ey(0,1) = sqrt(2); Ey(1,2) = sqrt(3.5); Ey(2,3) = sqrt(4.5); Ey(4,5) = (-1)*sqrt(10)); 
  Ey(5,6) = (-1)*sqrt(4.5); Ey(6,7) = (-1)*sqrt(3.5); Ey(7,8) = (-1)*sqrt(2); 
  Ey(1,0) =(-1)*sqrt(2); Ey(2,1) = (-1)*sqrt(3.5); Ey(3,2) = (-1)*sqrt(4.5); 
  Ey(5,4) = sqrt(10); Ey(6,5) = sqrt(4.5); Ey(7,6) = sqrt(3.5); Ey(8,7) = sqrt(2);
  Ez(0,8) = 4; Ez(1,7) = 3; Ez(2,6) = 2; Ez(3,5) = 1; Ez(5,3) = -1; Ez(6,2) = -2; Ez(7,1) = -3; Ez(8,0) = -4;
 
  for(int i = 0; i<nv; i++){
    vect ai(a[i].size());
    for(int j = 0; j<a[i].size();j++) ai[j] = (a[i])[j]; 
    vect<double> cx = boost::numeric::ublas::prec_prod(Ex, ai);
    vect<double> cy = boost::numeric::ublas::prec_prod(Ey, ai);
    vect<double> cz = boost::numeric::ublas::prec_prod(Ez, ai);
    int lambda = 100; //quadratic penalty multiplier
    for(int d = 0; d < 9; d++){
      vector<double> row(9*nv+2*nl+3*nv);
      row.fill(row.begin(), row.end(), 0);
      row[9*i + d] = lambda;
      row[9*nv + 2*nl + 3*i] = (-1)*lambda*cx(d);
      row[9*nv + 2*nl + 3*i + 1] = (-1)*lambda*cy(d);
      row[9*nv + 2*nl + 3*i + 2] = (-1)*lambda*cz(d);
      A.push_back(row);
      b.push_back(lambda*ai(d));      
    }
  }
}

void add_normal_constraints(HexExtr& h, vector<vector<double>>& A, vector<double>& b){
  for(i = 0; i < nl; i++){
    Vector_3 n = CGAL::compute_normal_of_cell_0(h.input_tet_mesh, (vertices[i].incident_dart));
    matrix<double> R(9,9);
    R = find_rotation_matrix(n); //TODO- Done
    vect<double> temp(9);
    temp(0) = 1; temp(1) = 0; temp(2) = 0; temp(3) = 0; temp(4) = 0; temp(5) = 0; temp(6) = 0; temp(7) = 0; temp(8) = 0;
    vect<double> h0 = boost::numeric::ublas::prec_prod(R, temp);
    temp(0) = 0; temp(4) = 1;
    vect<double> h4 = boost::numeric::ublas::prec_prod(R, temp);
    temp(4) = 0; temp(8) = 1;
    vect<double> h8 = boost::numeric::ublas::prec_prod(R, temp);
    int lambda = 100; //quadratic penalty multiplier
    for(int d = 0; d < 9; d++){
      vector<double> row(9*nv+2*nl+3*nv);
      row.fill(row.begin(), row.end(), 0);
      row[9*i + d] = lambda;
      row[9*nv + 2*i] = lambda*h0(d);
      row[9*nv + 2*i + 1] = lambda*h8(d);
      A.push_back(row);
      b.push_back(lambda*sqrt((double)7/12)*h4(d));
    }
  }
}

void optimise_frame_field(HexExtr& h, int n){ // n is the number of smoothing iterations
  int nl = find_number_of_boundary_vertices(h.input_tet_mesh);
  int nv = (h.vertices).length(); 
  sort_vertices(h, h.vertices); //TODO - DONE
  vector<vector<double>> a;
  for(int i = 0; i < n; i++){
    Eigen_matrix<double> A(0, (9*nv + 2*nl + 3*nv));
    vector<vector<double>> A_tobeconverted;
    Eigen_vector<double> b;
    vector<double> b_tobeconverted;
    add_smoothing_terms(h, A, b); //TODO- DONE
    add_normal_constraints(h, A_tobeconverted, b_tobeconverted);
    if(i>0){
      add_local_optim_constraints(h, a, A_tobeconverted, b_tobeconverted);
    }
    if(i == 0) vector<vector<double>> a;
    Eigen_vector X = b;
    Svd::solve(A, X); //solution gets stored in X
    a.clear();
    for(int j = 0;j<nv;j++){
      vector<double> temp;
      for(int k = 9*j; k<(9*j+9);k++){
          temp.push_back(X.at(k));
        }
      a.push_back(temp);
            
      //closest frame ? //TODO
    }
  }
}
