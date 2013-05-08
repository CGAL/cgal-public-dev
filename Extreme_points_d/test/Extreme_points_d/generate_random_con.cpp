#include <vector>
#include <fstream>
#include <string>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <cassert>

#include <CGAL/config.h>
#include <CGAL/Timer.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Extreme_points_d.h>

#include <CGAL/Random.h>

//#define RANDMAX 0xFFFFFFFF
#define RANDMAX 0xFFFF

typedef CGAL::Cartesian_d<double> Kernel_d;
typedef Kernel_d::Point_d Point_d;

const std::string whitespace = "        ";

CGAL::Random custom_random;


double random_coordinate() {
  double tmp = custom_random(RANDMAX) / static_cast<double>(RANDMAX);
  return (custom_random(RANDMAX) % 2 ? -tmp : tmp);
}

int random_int_coordinate() {
  int tmp = custom_random(RANDMAX);
  return (custom_random(RANDMAX) % 2 ? -tmp : tmp);
}



void generate_random_2nonzero_con (std::ostream& out, int dimension, int no_points) {
  
  assert(dimension > 1);
  
  out << dimension << whitespace;
  out << no_points << std::endl;
  
  for (int i = 0; i < no_points; ++i) {
    int indeces[2] = {custom_random(dimension),custom_random(dimension)};
    while (indeces[1] == indeces[0]) indeces[1] = custom_random(dimension);
    if (indeces[1] < indeces[0]) std::swap(indeces[0], indeces[1]);
    //std::cout << "0: " << indeces[0] << " 1: " << indeces[1] << std::endl;
     
    for (int j = 0; j < indeces[0]; ++j) {
      out << 0.0 << whitespace;
    }
    out << random_int_coordinate() << whitespace; // at first index
    for (int j = indeces[0]+1; j < indeces[1]; ++j) {
      out << 0.0 << whitespace;
    }
    out << random_int_coordinate() << whitespace; // at second index
    for (int j = indeces[1]+1; j < dimension; ++j) {
      out << 0.0 << whitespace;
    }
    out << "\n";
    
  }
}



void generate_random_Nnonzero_con (std::ostream& out, int dimension, int no_points, int N) {
  
  assert(dimension > 1);
  assert(N <= dimension);
  
  out << dimension << whitespace;
  out << no_points << std::endl;
  
  int index;
  std::vector<int> indeces;
  
  for (int i = 0; i < no_points; ++i) {
    indeces.clear();
    indeces.push_back(-1); // boundary index
    indeces.push_back(custom_random(dimension));
    for (int i = 1; i < N; ++i) {
      do {
        index = custom_random(dimension);
      } while (std::find(indeces.begin(), indeces.end(), index) != indeces.end());
      indeces.push_back(index);
    }
    std::sort(indeces.begin()+1, indeces.end());
    for (int k = 0; k < N; ++k) {
      for (int j = indeces[k]+1; j < indeces[k+1]; ++j) {
        out << 0.0 << whitespace;
      }
      out << random_int_coordinate() << whitespace; // produce non-zero coordinate
    }
    for (int j = indeces[N]+1; j < dimension; ++j) {
      out << 0.0 << whitespace;
    }
    out << "\n";  
  }
}


void generate_random_con (std::ostream& out, int dimension, int no_points, double coord_probability) {

  assert(coord_probability >= 0.0 && coord_probability <= 1.0);
  
  out << dimension << whitespace;
  out << no_points << std::endl;
  
  for (int i = 0; i < no_points; ++i) {
    for (int j = 0; j < dimension; ++j) {
      double dice = custom_random(RANDMAX) / static_cast<double>(RANDMAX);
      out << (dice < coord_probability ? random_coordinate() : 0.0) << whitespace;
    }
    out << "\n";
  }
}


int main(int argc, char **argv) {

  if (argc < 5) {
    std::cout << "%%% Usage: generate_random_con d n p s\n"
              << "%%% where d is the dimension, n is the number of points,\n"
              << "%%%      p is the number of non-zero entries per point,\n"
              << "%%%      and s is the seed for the random number generator.\n"
              << "%%%      Defaults are: d=5, n=100, p=1, and s=time(0)." << ".\n"
              << "%%%      Random numbers are generated between " << (-RANDMAX) << " and " << RANDMAX << std::endl;
  }

  const int   dimension = argc < 2 ? 5 : std::atoi(argv[1]);
  const int   no_points = argc < 3 ? 100 : std::atoi(argv[2]);
  const double coord_probability = argc < 4 ? 1.0 : std::atof(argv[3]); // probability of non-zero entry
  const int rand_seed = argc < 5 ? time(0) : std::atoi(argv[4]);
  
  custom_random = CGAL::Random(rand_seed);
  

  //generate_random_con(std::cout, dimension, no_points, coord_probability);
  //generate_random_2nonzero_con(std::cout, dimension, no_points);
  generate_random_Nnonzero_con(std::cout, dimension, no_points, static_cast<int>(coord_probability));
  
  
  return 0;
}
