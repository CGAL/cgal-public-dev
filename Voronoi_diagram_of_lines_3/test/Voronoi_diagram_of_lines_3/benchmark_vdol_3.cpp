
// TODO REMOVE THIS MACRO
#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1 

// ----------------------------------------------------
// includes for VOL_3
#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <CGAL/macros.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Voronoi_diagram_of_lines_3.h> 

#include <CGAL/VDOL_3/io.h> 

// TODO:
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
struct Linear_kernel: public CGAL::Cartesian<CGAL::Arithmetic_kernel::Rational> {};

typedef CGAL::Voronoi_diagram_of_lines_3<Linear_kernel> VDOL_3;

typedef VDOL_3::Line_3  Line_3;
typedef VDOL_3::Point_3 Point_3;
typedef VDOL_3::FT      FT; 

template<typename Line, typename Point> inline 
CGAL::Comparison_result compare_squared_distance(const Line& l1, const Line& l2, const Point& p){
  return CGAL::compare(
      CGAL::squared_distance(l1,p),
      CGAL::squared_distance(l2,p));
}

int main(int argc, char **argv)
{
  CGAL::Timer total_timer; 
  total_timer.start(); 
  
  Linear_kernel linear_kernel;
  
  //CGAL::set_pretty_mode(std::cout);
  //CGAL::set_pretty_mode(std::cerr);
  
  // Get the name of the input file from the command line, or use the default
  // last.dat file if no command-line parameters are given.
  const char * filename = (argc >= 2) ? argv[1] : "last.dat"; 
  if(argc > 1){
    CGAL::VDOL_3::copy_file_into(filename,"last.dat");
  }

  std::cout << "Run test with file: " << filename << std::endl;

  std::vector<Line_3> lines;
  CGAL::VDOL_3::read_lines_3(linear_kernel,filename,std::back_inserter(lines));
  if(lines.size()==0) 
    std::cerr << "Failed to read file " << filename << std::endl;
  
  std::cout << "Number of lines (incl base line): "<< lines.size() << std::endl;
  

  VDOL_3 vdol_3(lines.begin(),lines.end());
  
 
  total_timer.stop();
  std::cout << "Total time spent: " <<  total_timer.time() << std::endl;
  return 0;
}
