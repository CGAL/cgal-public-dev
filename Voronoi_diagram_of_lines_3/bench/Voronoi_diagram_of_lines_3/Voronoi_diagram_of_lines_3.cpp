
// TODO REMOVE THIS MACRO
#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1 

// #define CGAL_SL_VERBOSE 1
// #define CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE 1
// ----------------------------------------------------

// includes for VOL_3
#include <CGAL/basic.h>
#include <CGAL/Timer.h>
#include <CGAL/macros.h>
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Voronoi_diagram_of_lines_3.h> 
#include <CGAL/VDOL_3/io.h> 
#include <CGAL/iterator.h> 


// TODO:
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
struct Linear_kernel: public CGAL::Cartesian<CGAL::Arithmetic_kernel::Rational> {};

typedef CGAL::Voronoi_diagram_of_lines_3<Linear_kernel> VDOL_3;
typedef VDOL_3::Line_3  Line_3;
typedef VDOL_3::Point_3 Point_3;
typedef VDOL_3::FT      FT; 
typedef VDOL_3::RT      RT;
typedef VDOL_3::Poly_int_3   Poly_int_3;
typedef VDOL_3::Poly_int_2   Poly_int_2;
typedef VDOL_3::Poly_int_1   Poly_int_1;

template <class Integer> 
int get_bits(Integer& x){
  int count = 1;
  x=CGAL::abs(x);
  while(x!=0){
    x = CGAL::div(x,2);
    count++;
  }
  return count;
}

template<class PolyContainer>
void report_bits(const PolyContainer& container, double& min, double& ave, double& max){
  typedef typename PolyContainer::value_type Polynomial; 
  typedef CGAL::Polynomial_traits_d<Polynomial> PT; 
  typename PT::Monomial_representation monom_rep;
  typedef typename PT::Innermost_coefficient_type ICoeff; 
  typedef std::pair<CGAL::Exponent_vector,ICoeff> Monom; 
  std::vector<Monom> monoms; 
  BOOST_FOREACH(Polynomial poly, container){
    monom_rep(poly,std::back_inserter(monoms)); 
  }
  BOOST_FOREACH(Monom monom, monoms){
    double bits = get_bits(monom.second);
    min = (min==0)?bits:(std::min)(min,bits);
    max = (std::max)(max,bits);
    ave += bits;      
  }
  ave/=monoms.size(); 
  return; 
}


void compute_input_bits(
    const char *filename, 
    double& min_bits, double& ave_bits, double& max_bits){
  
  max_bits = min_bits = ave_bits = 0;
  
  std::ifstream     in_file (filename);
  unsigned int n; in_file >> n;
  for (unsigned int k = 0; k < n; ++k) {
    for (int i = 0; i<6;i++){
      FT a;
      in_file >> a;
      a = CGAL::abs(a); 
      int count = 0;  
      while(a >= 1){
        a /= 2; 
        count++;
      }
      max_bits = (std::max)(max_bits,double(count));   
      min_bits = (!min_bits)?count:(std::min)(min_bits,double(count)); 
      ave_bits += count; 
    }    
  }
  ave_bits /= 6*n; 
  in_file.close();
}

template <typename Linear_kernel, typename OutputIterator>
OutputIterator read_lines(const Linear_kernel& linear_kernel, const char *filename, OutputIterator oi){

  CGAL::set_pretty_mode(std::cout);

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

  for (unsigned int k = 0; k < n; ++k) {
    FT a,b,c,d,e,f;
    in_file >> a >> b >> c>> d >> e >> f;
    *oi++ = cline_3(cpoint_3(FT(a),FT(b),FT(c)),cpoint_3(FT(d),FT(e),FT(f)));
  }
  in_file.close();
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

  std::cerr << "Run test with file: " << filename << std::endl;

  std::vector<Line_3> lines;
  read_lines(linear_kernel,filename,std::back_inserter(lines));
  if(lines.size()==0) 
    std::cerr << "Failed to read file " << filename << std::endl;
  
  std::cerr << "Number of lines (incl base line): "<< lines.size() << std::endl;
  

  double line_bits_min(0),line_bits_ave(0),line_bits_max(0); 
  compute_input_bits(filename,line_bits_min,line_bits_ave,line_bits_max);
  

  CGAL::Real_timer timer; 
  timer.start(); 
  VDOL_3 vdol_3(lines.begin(),lines.end());
  timer.stop(); 

  
  std::cerr << "Total number of lines was:   " <<  lines.size() << std::endl; 
  std::cerr << "construct vdol_3 total time: " <<  timer.time() << std::endl; 
  
  std::pair<int,int> number_of_vertices = vdol_3.number_of_vertices(); 
  double bs_bits_min(0) ,bs_bits_ave(0) ,bs_bits_max(0); 

  std::set<Poly_int_3> bisectors; 
  vdol_3.report_bisectors(CGAL::inserter(bisectors));
  report_bits(bisectors,bs_bits_min,bs_bits_ave,bs_bits_max);

  double tbs_bits_min(0) ,tbs_bits_ave(0) ,tbs_bits_max(0); 
  std::set<Poly_int_3> trans_bisectors;
  vdol_3.report_transformed_bisectors(CGAL::inserter(trans_bisectors));
  report_bits(trans_bisectors,tbs_bits_min,tbs_bits_ave,tbs_bits_max); 
    
  double res_bits_min(0),res_bits_ave(0),res_bits_max(0); 
  std::set<Poly_int_1> resultants;
  vdol_3.report_resultants(CGAL::inserter(resultants));
  report_bits(resultants,res_bits_min,res_bits_ave,res_bits_max); 

  std::cout << "BENCH_RESULT" <<  " " 
            << filename << " " 
            << lines.size() << " "
            << timer.time() << "  "  
            
            << number_of_vertices.first << " " 
            << number_of_vertices.second << "  "
            
            << line_bits_min << " " 
            << line_bits_ave << " " 
            << line_bits_max << "  "

            << bs_bits_min << " " 
            << bs_bits_ave << " " 
            << bs_bits_max << "  "

            << tbs_bits_min << " " 
            << tbs_bits_ave << " " 
            << tbs_bits_max << "  "

            << res_bits_min << " " 
            << res_bits_ave << " " 
            << res_bits_max << "  "
    
            << std::endl; 
  
  return 0;
}
