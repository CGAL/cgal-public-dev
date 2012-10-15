#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1 
#define CGAL_VDOL_USE_TIMER 1
#define CGAL_VDOL_USE_COUNTER 1

// ----------------------------------------------------
// includes for VOL_3
#include <CGAL/basic.h>
#include <CGAL/Timer.h>

CGAL::Timer timer_total, timer_vdol_total, timer_vdol_hierarchy_total, timer_query;


#ifdef CGAL_VDOL_USE_COUNTER
long count_cell_contains = 0; 
long count_catch_frame =0; 
long count_all_vertex =0;
long count_extra_vertex=0; 
#endif 


#include <CGAL/macros.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Voronoi_diagram_of_lines_3.h> 
#include <CGAL/Voronoi_diagram_of_lines_hierarchy_3.h> 


// TODO:
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
struct Linear_kernel: public CGAL::Cartesian<CGAL::Arithmetic_kernel::Rational> {};

typedef CGAL::Voronoi_diagram_of_lines_3<Linear_kernel> VDOL_3;
typedef CGAL::Voronoi_diagram_of_lines_hierarchy_3<VDOL_3> VDOLH_3;

typedef VDOLH_3::Line_3  Line_3;
typedef VDOLH_3::Point_3 Point_3;
typedef VDOLH_3::FT      FT; 

// // get arithmetic kernel 
// typedef CGAL::Arithmetic_kernel AK;
// typedef AK::Integer Integer;
// typedef AK::Rational Rational;

// // define linear kernel 
// #if 1
// typedef CGAL::Cartesian<Rational> EK_3;
// struct K: public CGAL::Exact_predicates_inexact_constructions_kernel {};
// #else
// typedef CGAL::Cartesian<Rational> EK_3;
// typedef EK_3 K;
// // struct K: public CGAL::Exact_predicates_exact_constructions_kernel {};
// #endif

// typedef EK_3::Line_3                                    ELine_3;
// typedef EK_3::Point_3                                   EPoint_3;

// // define algebraic kerenel 
// typedef CGAL::Algebraic_kernel_d_1_generator<Integer>     GEN_AK_1;
// typedef GEN_AK_1::Default_algebraic_kernel_1            AK_1;
// typedef CGAL::Algebraic_curve_kernel_2<AK_1>            AK_2;

// // define curved kernel 
// typedef CGAL::Curved_kernel_via_analysis_2< AK_2 >      CK_2;

// // define traits for lower envelop diagram 
// typedef CGAL::VDOL_3::
// Single_voronoi_cell_envelope_traits_3<EK_3, CK_2>       SVCET_3;
 

typedef VDOL_3::VCOL_3  VCOL_3;
typedef VCOL_3::Envelope_diagram_2 Envelope_diagram_2;
typedef Envelope_diagram_2::Vertex_const_iterator Vertex_const_iterator;
typedef Envelope_diagram_2::Surface_const_iterator Surface_const_iterator;

void copy_file_into(const char *infilename, const char *outfilename){
  std::ifstream     in_file (infilename);
  std::ofstream     out_file (outfilename);
  std::cerr << "save " << infilename << " to " << outfilename << "... " << std::flush;
  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << infilename << "!" << std::endl;
    return; 
  }
  if (! out_file.is_open()) {
    std::cerr << "Failed to open " << outfilename << "!" << std::endl;
    return; 
  }
  out_file << in_file.rdbuf();
  std::cerr << "  done "<< std::endl;
}


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
// ---------------------------------------- 
// the voronoi cell of one line 
  
  //CGAL::set_pretty_mode(std::cout);
  //CGAL::set_pretty_mode(std::cerr);
  
  // Get the name of the input file from the command line, or use the default
  // last.dat file if no command-line parameters are given.
  const char * filename = (argc >= 2) ? argv[1] : "last.dat"; 
  if(argc > 1){
    copy_file_into(filename,"last.dat");
  }
  
  int curr_arg=2, hierarchy_ratio, no_queries, no_runs;
  
  hierarchy_ratio=(argc<3)?3:atoi(argv[curr_arg]);
  curr_arg++; 
  no_queries=(argc<4)?2000:atoi(argv[curr_arg]);
  curr_arg++; 
  no_runs=(argc<5)?20:atoi(argv[curr_arg]);
  curr_arg++;
  
  

  std::cerr << "Run test with file: " << filename << std::endl;
  std::cerr << "hierarchy_ratio   : " << hierarchy_ratio << std::endl;
  std::cerr << "no_queries   :      " << no_queries << std::endl;
  std::cerr << "no_runs   :         " << no_runs << std::endl;
  

  std::vector<Line_3> lines;
  read_lines(linear_kernel,filename,std::back_inserter(lines));
  if(lines.size()==0) 
    std::cerr << "Failed to read file " << filename << std::endl;
  
  std::cerr << "Number of lines (incl base line): "<< lines.size() << std::endl;
   
//   EK_3::Equal_3 equal_3;
//   for(int i = 0; i <lines.size();i++){
//     for(int j = i+1; j <lines.size();j++){
//       if(CGAL::do_intersect(lines[i],lines[j])){
//         //std::cout <<CGAL::VDOL_3::construct_bisector_3(lines[i],lines[j])<< std::endl;
//         std::cout <<"input contains equal lines: "<< i << " " << j << std::endl;
//         return 0; 
//       }
//     }
//   }
  
  int hierarchy_size; 
  for(int RUN = 0; RUN < no_runs; RUN++){

    std::random_shuffle(lines.begin(),lines.end());    
    timer_total.start();
    VDOLH_3 vdolh_3(lines.begin(),lines.end(),hierarchy_ratio);
    timer_total.stop();
    
    assert(lines.size()>=2);
    
    CGAL::Random rand;
    for(int i = 0; i<no_queries; i++){ 
      Point_3 p(rand.get_int(0,1000),rand.get_int(0,1000),rand.get_int(0,1000));
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
      timer_query.start();
      vdolh_3.locate(p,std::back_inserter(result_2));
      timer_query.stop();
      assert(result_1 == result_2);
    }
    
    
    // statistics 
    typedef VDOLH_3::VDOL_3::Cell_iterator Cell_iterator;
    for( Cell_iterator cit = vdolh_3.vdol_3().cells().begin();
         cit != vdolh_3.vdol_3().cells().end();
         cit++){
      count_extra_vertex+=cit->number_of_vertical_asymptotes();
      count_all_vertex+=cit->number_of_vertical_asymptotes();
      
      for(Vertex_const_iterator vit = cit->envelope_diagram_2(true).vertices_begin();
          vit != cit->envelope_diagram_2(true).vertices_end();
          vit++){
        count_all_vertex++;
        std::set<long> ids;
        for(Surface_const_iterator sit = vit->surfaces_begin();
            sit != vit->surfaces_end();
            sit++){
          ids.insert(sit->id());
        }
        if(ids.size()<3) 
          count_extra_vertex++;
      }
      
      for(Vertex_const_iterator vit = cit->envelope_diagram_2(false).vertices_begin();
          vit != cit->envelope_diagram_2(false).vertices_end();
        vit++){
        count_all_vertex++;
        count_extra_vertex+=cit->number_of_vertical_asymptotes();
        count_all_vertex+=cit->number_of_vertical_asymptotes();
        std::set<long> ids;
        for(Surface_const_iterator sit = vit->surfaces_begin();
            sit != vit->surfaces_end();
            sit++){
          ids.insert(sit->id());
        }
        if(ids.size()<3) 
          count_extra_vertex++;
      }
    }
    hierarchy_size = vdolh_3.hierarchy().size();
  }

  std::cerr << " ILE"
            << " LINES " 
            << " HRATIO "
            << " SVDOLH " 
            << " TTOTAL " 
            << " TVODL " 
            << " TVDOLH " 
            << " TQUERY "
            << " CCONTAIN " 
            << " NCATCH " 
            << " CAVERTEX " 
            << " CEVERTEX "  
            << std::endl;
#if 0
  std::cout 
    << filename << " " 
    << lines.size() << " " 
    << hierarchy_ratio << " "
    << hierarchy_size << " "
    << timer_total.time() << " " 
    << timer_vdol_total.time() << " " 
    << timer_vdol_hierarchy_total.time() << " " 
    << timer_query.time() << " " 
    << count_cell_contains << " "
    << count_catch_frame << " "
    << count_all_vertex << " "
    << count_extra_vertex << " "
    << std::endl;
#else  
  std::cout.precision(2);
  std::cout.setf(std::ios::fixed,std::ios::floatfield); 
  std::cout << double(count_cell_contains) / (no_queries * no_runs)  << " " << std::flush;
#endif 

  return 0;
}
