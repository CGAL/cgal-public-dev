
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

template <typename Linear_kernel, typename OutputIterator>
OutputIterator read_lines(const Linear_kernel& linear_kernel, const char *filename, OutputIterator oi){

  typedef typename Linear_kernel::Point_3 Point_3; 
  typedef typename Linear_kernel::Line_3  Line_3; 
  typedef typename Linear_kernel::FT      FT; 
  
  typename Linear_kernel::Construct_point_3 cpoint_3 = linear_kernel.construct_point_3_object(); 
  typename Linear_kernel::Construct_line_3  cline_3  = linear_kernel.construct_line_3_object(); 
  
  // Open the input file.
  std::ifstream     in_file (filename);
  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << filename << "!" << std::endl;
    return oi;
  }
  // Read the lines from the file
  // The input file format should be (all coordinate values are integers):
  // <n>                                       // number of lines.
  // <a1_x> <a1_y> <a1_z> <b1_x> <b1_y> <b1_z> // line #1.
  // <a2_x> <a2_y> <a2_z> <b2_x> <b2_y> <b2_z> // line #2.
  //   :      :       :      :
  // <an_x> <an_y> <an_z> <bn_x> <bn_y> <bn_z> // line #n.
  
  // read number of lines   
  unsigned int n; in_file >> n;

  for (int k = 0; k < n; ++k) {
    int a,b,c,d,e,f;
    in_file >> a >> b >> c>> d >> e >> f;
    *oi++ = cline_3(cpoint_3(FT(a),FT(b),FT(c)),cpoint_3(FT(d),FT(e),FT(f)));
  }
  in_file.close();
}



template<typename Line, typename Point> inline 
CGAL::Comparison_result compare_squared_distance(const Line& l1, const Line& l2, const Point& p){
  return CGAL::compare(
      CGAL::squared_distance(l1,p),
      CGAL::squared_distance(l2,p));
}

int main(int argc, char **argv)
{

  Linear_kernel linear_kernel;
  Linear_kernel::Construct_point_3 cpoint_3 = linear_kernel.construct_point_3_object(); 

  
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
  
  assert(lines.size()>=2);
  
  for(int i = -10; i<10; i+=1){ 
    for(int j = -10; j<10; j+=1){ 
      for(int k = -10; k<10; k+=1){
        Point_3 p = cpoint_3(FT(i),FT(j),FT(k));
        int min_l = compare_squared_distance(lines[0],lines[1],p)==CGAL::SMALLER?0:1;
        for(int l = 0; l < lines.size(); l++){
          if(compare_squared_distance(lines[l],lines[min_l],p) == CGAL::SMALLER){
            min_l=l;
          }
        }
        std::vector<int> result_1; 
        for(int l = 0; l < lines.size(); l++){
          if(compare_squared_distance(lines[min_l],lines[l],p)==CGAL::EQUAL){
            result_1.push_back(l);
          }
        }  
        std::vector<int> result_2;
        vdol_3.locate(p,std::back_inserter(result_2));
           
//         for(int l =0; l < result_1.size();l++){
//           std::cout << result_1[l] << " ";
//         }
//         std::cout<<std::endl;
//         for(int l =0; l < result_2.size();l++){
//           std::cout << result_2[l] << " ";
//         }
//         std::cout<<std::endl;
//         std::cout<<std::endl;
        assert(result_1 == result_2);
      }
    }
  }

  return 0;
}
