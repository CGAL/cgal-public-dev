// covariance_matrix_test.cpp

//----------------------------------------------------------
// Test the Covariance_matrix_3 class:
// For each input point, create a tensor either from eigenvalues 
// and eigenvectors or from the point and its normal.
//----------------------------------------------------------


// CGAL
#include <CGAL/Timer.h>
#include <CGAL/array.h>
#include <CGAL/Covariance_matrix_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/disable_warnings.h>

// C++
#include <math.h>
#include <iostream>
#include <random>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Covariance Matrix
typedef CGAL::Covariance_matrix_3<Kernel> Covariance_matrix;
typedef CGAL::cpp11::array<FT, 6>   Eigen_matrix;
typedef CGAL::cpp11::array<FT, 3>   Eigen_vector;
typedef CGAL::cpp11::array<FT, 9>   Eigen_three_vectors;

FT ut_v_u(const Vector& u, const Eigen_matrix& v)
{
  FT prod = (u.x() * u.x() * v[0] + 
          u.y() * u.y() * v[3] + 
          u.z() * u.z() * v[5] +
          u.x() * u.y() * v[1] * 2 +
          u.x() * u.z() * v[2] * 2 +
          u.y() * u.z() * v[4] * 2);
  return prod;
}


void print_matrix(const Covariance_matrix& my_conv)
{
  Eigen_matrix my_tensor = my_conv.tensors();
  Eigen_vector my_values = my_conv.eigen_values();
  Vector vmin = my_conv.eigen_vect(0);
  Vector vmid = my_conv.eigen_vect(1);
  Vector vmax = my_conv.eigen_vect(2);

  std::cout << "Initialize Covariance Matrix: " << std::endl;
  std::cout << "Tensor: [" << my_tensor[0] << ", "
                           << my_tensor[1] << ", "
                           << my_tensor[2] << ", "
                           << my_tensor[3] << ", "
                           << my_tensor[4] << ", "
                           << my_tensor[5] << "] "
                           << std::endl;
  std::cout << "Largest eigenvalue is " << my_values[2] <<
               " with eigenvector [" << vmax << "]" << std::endl;
  std::cout << "Middle eigenvalue is " << my_values[1] <<
               " with eigenvector [" << vmid << "]" << std::endl;
  std::cout << "Smallest eigenvalue is " << my_values[0] <<
               " with eigenvector [" << vmin << "]" << std::endl;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
  std::cerr << "Test the Covariance_matrix_3 class" << std::endl;

  //***************************************
  // test functions
  //***************************************

  // Constructor
  Covariance_matrix identity;
  print_matrix(identity);

  Point p(0, 0, 0);
  Vector n(0, 0, 1);
  Covariance_matrix xy_plane(p, n, 10);
  print_matrix(xy_plane);

  
  std::random_device rd;
  std::mt19937 e2(rd());
  std::uniform_real_distribution<> dist(-1, 1);

  FT normal_prod = ut_v_u(n, xy_plane.tensors());
  std::cout << "The normal product is: " << normal_prod << std::endl;

  bool flag = true;
  size_t max_iter = 100;

  for(size_t i = 0; i < max_iter; i++){
    Vector random_vect(dist(e2), dist(e2), dist(e2));
    FT new_prod = ut_v_u(random_vect, xy_plane.tensors());
    new_prod = new_prod / std::sqrt(new_prod * new_prod);
    if(new_prod > normal_prod){
      std::cout << "Find larger product " << new_prod << 
                   " by the vector [" << random_vect << "]" << std::endl;
      flag = false;
    }
  }

  if(flag == false)
    std::cout << "Error: find larger product!" << std::endl;
  else
    std::cout << "Successful!" << std::endl;

  return 0;
}
