
#define CGAL_LEDA_VERSION 620
//#define CGAL_LAZY_KERNEL_DEBUG 1 // comment out if not needed 

/* Set the required traits */
#define USE_CONIC_TRAITS 1
#define USE_RATIONAL_ARC_TRAITS 0
#define USE_SQRT_TRAITS 0
#define USE_LAZY 0
#define USE_LEDA 0

#include <iostream>
#include <fstream>
#include <cstdio>
#include <boost/variant.hpp>
#include <string>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <vector>
#include <algorithm>
#include <CGAL/Timer.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Lazy_kernel.h>

#include "color_print.h"

//#include <Arrangement_general_functions.h>
#define LTS_DRAW_ARR 0
#define ARR_ON_SUR_DEBUG 0
#define LINES_DEBUG 0
#define ISOLATED_POINTS_DEBUG 0
#define BOUNDED_SEGMENTS_VECTOR_DEBUG 0
#define OBSERVER_PRINTS 0
#define OBJ_ON_ARR_DEBUG 0
#define CGAL_DEBUG_OUTPUT 0

#define PRINT_OUTPUT 0
#define CGAL_IDENTIFICATION_XY CGAL_X_MINUS_11_Y_7

#define LTS_WITH_SEGMENTS 1

#if USE_CONIC_TRAITS
#include <CGAL/Arr_conic_traits_2.h>
#endif
#if USE_RATIONAL_ARC_TRAITS
#include <CGAL/Arr_rational_function_traits_2.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#endif

#if USE_SQRT_TRAITS
#include <CGAL/Arr_rational_function_traits_2.h>
#include <CGAL/Algebraic_kernel_2_1.h>
#include <CGAL/Lazy_exact_nt.h>
#endif

#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>


#include <CGAL/CORE_algebraic_number_traits.h>

#if USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_rational.h>
#include <CGAL/leda_real.h>
#endif


/* Lines_through_segments includes */
#include <CGAL/IO/Lines_through_segments_io_stream.h>
#include <CGAL/Lines_through_segments_3.h>
#include <CGAL/Lines_through_segments_traits_3.h>

#if USE_LEDA
typedef leda_real                                       _Algebraic;
typedef leda_rational                                   _Rational;
typedef leda_integer                                    _Integer;
#else
typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Algebraic                            _Algebraic;
typedef Nt_traits::Rational                             _Rational;
typedef CORE::BigInt                                    _Integer;
#endif

#if  USE_LAZY
typedef CGAL::Lazy_exact_nt<_Algebraic>                 Algebraic;
typedef CGAL::Lazy_exact_nt<_Rational>                  Rational;
typedef CGAL::Lazy_exact_nt<_Integer>                   Integer;

#if 1
// seems faster 
typedef CGAL::Cartesian<CGAL::Lazy_exact_nt<_Rational> >  Rational_kernel;
#else
// but I leave this here for now

// class Rational_kernel
//   : public CGAL::internal::Static_filters<
//       CGAL::Type_equality_wrapper<
//              CGAL::Lazy_kernel_base< CGAL::Simple_cartesian<_Rational>, CGAL::Simple_cartesian<CGAL::Interval_nt_advanced>,
// 	                       CGAL::Cartesian_converter< CGAL::Simple_cartesian<_Rational>, CGAL::Simple_cartesian<CGAL::Interval_nt_advanced> >, Rational_kernel>,
//              Rational_kernel >, false>
// {};

class Rational_kernel
  : public CGAL::Type_equality_wrapper<
             CGAL::Lazy_kernel_base< CGAL::Simple_cartesian<_Rational>, CGAL::Simple_cartesian<CGAL::Interval_nt_advanced>,
	                       CGAL::Cartesian_converter< CGAL::Simple_cartesian<_Rational>, CGAL::Simple_cartesian<CGAL::Interval_nt_advanced> >, Rational_kernel>,
             Rational_kernel >
{};
#endif 

#else
typedef _Algebraic                                      Algebraic;
typedef _Rational                                       Rational;
typedef _Integer                                        Integer; 
typedef CGAL::Cartesian<Rational>                       Rational_kernel;

#endif




typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Segment_3<Rational_kernel>                Rational_segment_3;
typedef CGAL::Line_3<Rational_kernel>                   Rational_line_3;
typedef CGAL::Point_3<Rational_kernel>                  Rational_point_3;
   
#if LTS_WITH_SEGMENTS
typedef boost::true_type With_segments;
#else
typedef boost::false_type With_segments;
#endif

#if USE_CONIC_TRAITS
typedef CGAL::Arr_conic_traits_2<Rational_kernel, Alg_kernel, Nt_traits>    Conic_traits_arr_on_plane_2;
#endif
#if USE_RATIONAL_ARC_TRAITS
typedef CGAL::Algebraic_kernel_d_1<Integer>	  AK1;
typedef CGAL::Arr_rational_function_traits_2<AK1>                  Rational_arc_arr_traits_arr_on_plane_2;
#endif
#if USE_SQRT_TRAITS
typedef CGAL::Algebraic_kernel_2_1<Rational >	  AK1;
typedef CGAL::Arr_rational_function_traits_2<AK1>                  Rational_arc_arr_traits_arr_on_plane_2;
#endif
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Rational_kernel>          Traits_arr_on_sphere_2;
#if USE_CONIC_TRAITS
typedef CGAL::Lines_through_segments_traits_3<Alg_kernel,
                                              Rational_kernel,
                                              Conic_traits_arr_on_plane_2,
                                              Traits_arr_on_sphere_2> Lines_through_segs_traits_using_conic_2;
#endif
#if (USE_RATIONAL_ARC_TRAITS || USE_SQRT_TRAITS)
typedef CGAL::Lines_through_segments_traits_3<Alg_kernel,
                                              Rational_kernel,
                                              Rational_arc_arr_traits_arr_on_plane_2,
                                              Traits_arr_on_sphere_2> Lines_through_segs_traits_using_rational_arc_2;
#endif
using namespace std;

void ReportError(std::string str, int line, std::string file)
{
  std::cout << change_color(CGAL_RED,str) << std::endl;
  std::cout << change_color(CGAL_RED,"File = ", file, "\tLine = ", line) << std::endl;
  exit(EXIT_FAILURE);
}

/* Parse line at the inputfile */
#define NUM_OF_COORDINATES 6
#define LINE_SIZE 400

template <class Rational>
static void ParseInputFileLine(Rational coordinates[NUM_OF_COORDINATES],char line[LINE_SIZE])
{
  int pos = 0,next_pos = 0;
  string temp_str;
  int number_index = 0;
  string line_str(line);
      
  while (number_index < 6)
    {
      while ((next_pos = line_str.find(" ",pos)) != string::npos)
        {
          temp_str = line_str.substr(pos,next_pos-pos);
#if 0
          std::cout << change_color(CGAL_RED,"NOTE the parser can't parse fractions") << std::endl;
          double temp = atof(temp_str.c_str());
          coordinates[number_index] = Rational(temp);
#else
          coordinates[number_index] = Rational(strtod (temp_str.c_str(),NULL));
#endif
         
          number_index++;
          pos = next_pos + 1;
        }

      temp_str = line_str.substr(pos,line_str.size()-pos);
#if 0
      double temp = atof(temp_str.c_str());

      coordinates[number_index] = Rational(temp);
#else
      coordinates[number_index] = Rational(strtod (temp_str.c_str(),NULL));
#endif
      number_index++;
    }
   
  if (number_index != 6)
    {
      ReportError("Bad number format",__LINE__,__FILE__);		
    }

}


/*************************************************************
 * The following function parses an input file into array of lines.
 *
 * Input:
 *      input_file_name - A file with list of lines. each line is represented
 *                        by 6 coordinates, the first 3 of the first point and the last three
 *                        of another point. The line passes through these two points.
 *                        The first line at the file holds the number of lines at the file.
 *       
 * Output:
 *      lines -           Vector of lines.
 *
 *
 **************************************************************/

template <class Rational_segment_3,class Rational,class Rational_point_3>
int ReadInputFileFlip(char* input_file_name,vector<Rational_segment_3> &lines,
    int &expected_num_of_lines,
    int &expected_num_of_arr_curves,
    int &expected_num_of_arr_polygons,
    int &expected_num_of_point_3,
    int &expected_num_of_arr_arcs,
    int &expected_num_of_overlap_seg_3,
    bool flip,
    int &expected_num_of_lines_flip,
    int &expected_num_of_arr_curves_flip,
    int &expected_num_of_arr_polygons_flip,
    int &expected_num_of_point_3_flip,
    int &expected_num_of_arr_arcs_flip,
    int &expected_num_of_overlap_seg_3_flip)
{
  int ii = 0;
  int num_of_lines = 0;
  char line[LINE_SIZE];
  Rational coordinates[NUM_OF_COORDINATES];
  cout<<input_file_name<<endl;

  ifstream input_file(input_file_name);
   
  if (input_file.fail())
    {
      ReportError("Error opening file",__LINE__,__FILE__);
    }
  input_file.getline(line,LINE_SIZE);
   
   
  if (sscanf(line,"%d",&num_of_lines) != 1)
    {
      ReportError("Number of lines at file not specified.",__LINE__,__FILE__);
    }
   
  /* Read the file */
  while (!input_file.eof() && ii < num_of_lines)
    {
      input_file.getline(line,sizeof(line));
      ParseInputFileLine(coordinates,line);
      ii++;
      lines.push_back(Rational_segment_3(
                          Rational_point_3(coordinates[0],coordinates[1],coordinates[2]),
                          Rational_point_3(coordinates[3],coordinates[4],coordinates[5])));
    }
   
  if (ii != num_of_lines)
    {
      ReportError("Two many lines in file",__LINE__,__FILE__);
    }

  input_file.getline(line,sizeof(line));
  if (strncmp("OUTPUT:",line,7) != 0)
    {
      ReportError("OUTPUT is missing",__LINE__,__FILE__);
    }
   
  input_file.getline(line,sizeof(line));
  if (sscanf(line,"%d",&expected_num_of_lines) != 1)
    {
      ReportError("Number of output lines at file not specified.",__LINE__,__FILE__);
    }

  input_file.getline(line,sizeof(line));
  if (sscanf(line,"%d",&expected_num_of_arr_curves) != 1)
    {
      ReportError("Number of output curves at file not specified.",__LINE__,__FILE__);
    }

  input_file.getline(line,sizeof(line));
  if (sscanf(line,"%d",&expected_num_of_arr_polygons) != 1)
    {
      ReportError("Number of output polygons at file not specified.",__LINE__,__FILE__);
    }

  input_file.getline(line,sizeof(line));
  if (sscanf(line,"%d",&expected_num_of_point_3) != 1)
    {
      ReportError("Number of output point 3 at file not specified.",__LINE__,__FILE__);
    }

  input_file.getline(line,sizeof(line));
  if (sscanf(line,"%d",&expected_num_of_arr_arcs) != 1)
    {
      ReportError("Number of output arcs at file not specified.",__LINE__,__FILE__);
    }

  input_file.getline(line,sizeof(line));
  if (sscanf(line,"%d",&expected_num_of_overlap_seg_3) != 1)
    {
      ReportError("Number of output overlap segs at file not specified.",__LINE__,__FILE__);
    }
  input_file.getline(line,sizeof(line));
  if (flip && strncmp("AFTER FLIP:",line,11) == 0)
    {
      input_file.getline(line,sizeof(line));
      if (sscanf(line,"%d",&expected_num_of_lines_flip) != 1)
        {
          ReportError("Number of output lines at file not specified.",__LINE__,__FILE__);
        }
      
      input_file.getline(line,sizeof(line));
      if (sscanf(line,"%d",&expected_num_of_arr_curves_flip) != 1)
        {
          ReportError("Number of output curves at file not specified.",__LINE__,__FILE__);
        }

      input_file.getline(line,sizeof(line));
      if (sscanf(line,"%d",&expected_num_of_arr_polygons_flip) != 1)
        {
          ReportError("Number of output polygons at file not specified.",__LINE__,__FILE__);
        }
      
      input_file.getline(line,sizeof(line));
      if (sscanf(line,"%d",&expected_num_of_point_3_flip) != 1)
        {
          ReportError("Number of output point 3 at file not specified.",__LINE__,__FILE__);
        }
      
      input_file.getline(line,sizeof(line));
      if (sscanf(line,"%d",&expected_num_of_arr_arcs_flip) != 1)
        {
          ReportError("Number of output arcs at file not specified.",__LINE__,__FILE__);
        }
      
      input_file.getline(line,sizeof(line));
      if (sscanf(line,"%d",&expected_num_of_overlap_seg_3_flip) != 1)
        {
          ReportError("Number of output overlap segs at file not specified.",__LINE__,__FILE__);
        }
    }
   
  input_file.close();
  return num_of_lines;
}

template <class Rational_segment_3,class Rational,class Rational_point_3>
int ReadInputFile(char* input_file_name,vector<Rational_segment_3> &lines,
    int &expected_num_of_lines,
    int &expected_num_of_arr_curves,
    int &expected_num_of_arr_polygons,
    int &expected_num_of_point_3,
    int &expected_num_of_arr_arcs,
    int &expected_num_of_overlap_seg_3)
{
  int unused;
  return ReadInputFileFlip
    <Rational_segment_3,Rational,Rational_point_3>(input_file_name,lines,
        expected_num_of_lines,
        expected_num_of_arr_curves,
        expected_num_of_arr_polygons,
        expected_num_of_point_3,
        expected_num_of_arr_arcs,
        expected_num_of_overlap_seg_3,
        false,
        unused,
        unused,
        unused,
        unused,
        unused,
        unused);
}

template <typename Lines_through_segs>
bool get_all_common_lines(
    Lines_through_segs& line_through_segs,
    vector<Rational_segment_3> &lines,
    int &expected_num_of_lines,
    int &expected_num_of_arr_curves,
    int &expected_num_of_arr_polygons,
    int &expected_num_of_point_3,
    int &expected_num_of_arr_arcs,
    int &expected_num_of_overlap_seg_3)
{
  int num_of_lines = 0;
  int num_of_arr_curves = 0;
  int num_of_arr_polygons = 0;
  int num_of_points = 0;
  int num_of_arr_arcs = 0;
  int num_of_overlap_seg_3 = 0;
#if LTS_WITH_SEGMENTS
  typedef typename Lines_through_segs::Transversal_with_segments LTS_output_obj;
#else
  typedef typename Lines_through_segs::Transversal LTS_output_obj;
#endif
  
  list<LTS_output_obj> output_list;
  
  CGAL::Timer t;
  t.reset();
  t.start();

  int count =0; 
  int runs = 1;
  while (count ++ < runs ){
    output_list.clear();
    line_through_segs(lines.begin(),lines.end(),
                      std::front_inserter(output_list));
  }
  t.stop();
    
  typename list<LTS_output_obj>::iterator it_output_list;

  {
    typename std::vector< boost::shared_ptr<typename Lines_through_segs::Arr_on_plane> >::iterator
      b, e;
    boost::tie(b, e) = line_through_segs.planar_arrangements();
    typename std::vector< boost::shared_ptr<typename Lines_through_segs::Arr_on_sphere> >::iterator
      bs, es;
    boost::tie(bs, es) = line_through_segs.spherical_arrangements();
  }

#if PRINT_OUTPUT
  std::cout << t.time()/runs << std::endl;
  std::cout << "list size = " << output_list.size() << std::endl;
  cout << "OUTPUT:" << endl;
#endif
   
  typename Lines_through_segs::Mapped_general_polygon_2   *polygon_obj;
  typename Lines_through_segs::Mapped_x_monotone_curve_2  *curve_obj;
  typename Lines_through_segs::Mapped_point_2             *arr_point_2_obj;
  typename Lines_through_segs::Line_3                     *line_obj;
  typename Lines_through_segs::Through_point_3            *point_obj;
  typename Lines_through_segs::Through_point_3_segment_3  *arc_obj;
  typename Lines_through_segs::Through_segment_3          *seg_obj;

  typename Lines_through_segs::Mapped_2                   *mapped_obj;
  typename Lines_through_segs::Through_3                  *through_obj;
   
  typename Lines_through_segs::Transversal                   transversal;
  typename Lines_through_segs::Transversal_with_segments     *transversalw;
   
  for (it_output_list = output_list.begin(); it_output_list != output_list.end(); it_output_list++)
    {
#if LTS_WITH_SEGMENTS
#if PRINT_OUTPUT
      cout << "S1 = " << *it_output_list->second[0] << std::endl;
      cout << "S2 = " << *it_output_list->second[1] << std::endl;
      cout << "S3 = " << *it_output_list->second[2] << std::endl;
      cout << "S4 = " << *it_output_list->second[3] << std::endl;
      cout << it_output_list->first << std::endl;
#endif
      transversalw = (&(*it_output_list));
      transversal = transversalw->first;
#else
#if PRINT_OUTPUT
      cout << (*it_output_list) << endl;
#endif
      transversal = (*it_output_list);
#endif
      
      if ((line_obj = boost::get < typename Lines_through_segs::Line_3 >(&transversal)))
        {
          num_of_lines++;
        }
      else if ((mapped_obj = boost::get<typename Lines_through_segs::Mapped_2>(&transversal)))
        {
          typename Lines_through_segs::Mapped_2::Mapped_line_3 line = mapped_obj->line();
          typename Lines_through_segs::Mapped_transversal mapped_transversal = mapped_obj->mapped_transversal();
          if ((curve_obj = boost::get<typename Lines_through_segs::Mapped_x_monotone_curve_2>(&mapped_transversal)))
            {
              num_of_arr_curves++;
            }
          else if((polygon_obj = boost::get<typename Lines_through_segs::Mapped_general_polygon_2> (&mapped_transversal)))
            {
              num_of_arr_polygons++;
            }
          else if ((arr_point_2_obj = boost::get<typename Lines_through_segs::Mapped_point_2> (&mapped_transversal)))
            {
              num_of_lines++;
            }
          else
            {
              ReportError("Unexpected Error - invalid mapped obj value",__LINE__,__FILE__);
            }
        }
      else if ((through_obj = boost::get<typename Lines_through_segs::Through_3>(&transversal)))
        {
          typename Lines_through_segs::Through_transversal through_transversal = through_obj->through_transversal();
          if ((arc_obj = boost::get<typename Lines_through_segs::Through_point_3_segment_3> (&through_transversal)))
            {
              num_of_arr_arcs++;
            }
          else if ((seg_obj = boost::get<typename Lines_through_segs::Through_segment_3> (&through_transversal)))
            {
              num_of_overlap_seg_3++;
            }
          else if ((point_obj = boost::get<typename Lines_through_segs::Through_point_3> (&through_transversal)))
            {
              num_of_points++;
            }
          else
            {
              ReportError("Unexpected Error - Invalid through obj type",__LINE__,__FILE__);
            }
        }
      else
        {
          ReportError("Unexpected Error - Invalied output obj type",__LINE__,__FILE__);
        }
    }
   
#if PRINT_OUTPUT
  cout<<endl<<endl<<endl;
#endif
//    std::cout << "num_of_lines = " << num_of_lines << std::endl;
//    std::cout << "num_of_arr_polygons = " << num_of_arr_polygons << std::endl;
//    std::cout << "num_of_arr_curves = " << num_of_arr_curves << std::endl;
//   std::cout << "num_of_arr_arcs = " << num_of_arr_arcs << std::endl;
   
  if (expected_num_of_lines == -1)
    {
      expected_num_of_lines = num_of_lines;
      expected_num_of_arr_curves = num_of_arr_curves;
      expected_num_of_arr_polygons = num_of_arr_polygons;
      expected_num_of_arr_arcs = num_of_arr_arcs;
      expected_num_of_overlap_seg_3 = num_of_overlap_seg_3;
      expected_num_of_point_3 = num_of_points;
    }
  else
    {
      if (expected_num_of_lines != num_of_lines)
        {
          std::cout << "num_of_lines = " << num_of_lines << std::endl;
          return false;
        }
      
      if (expected_num_of_arr_curves != num_of_arr_curves)
        {
          return false;
        }
      
      if (expected_num_of_arr_polygons != num_of_arr_polygons)
        {
          return false;
        }

      if (expected_num_of_arr_arcs != num_of_arr_arcs)
        {
          return false;
        }

      if (expected_num_of_overlap_seg_3 != num_of_overlap_seg_3)
        {
          return false;
        }

      if (expected_num_of_point_3 != num_of_points)
        {
          return false;
        }
      
    }
  return true;
}




static int iteration = -1;

template <typename Lines_through_segs>
void run_all_permutations(
    Lines_through_segs& line_through_segs, 
    vector<Rational_segment_3> &lines,
    int *lines_index,
    unsigned int len ,
    unsigned int min, 
    bool print,
    int &expected_num_of_lines,
    int &expected_num_of_arr_curves,
    int &expected_num_of_arr_polygons,
    int &expected_num_of_point_3,
    int &expected_num_of_arr_arcs,
    int &expected_num_of_overlap_seg_3)
{
  if(min == len || iteration > 2000)
    {
      return;
    }
          
  if (print)
    {
      iteration++;
         
      if (!get_all_common_lines(line_through_segs,
              lines,
              expected_num_of_lines,
              expected_num_of_arr_curves,
              expected_num_of_arr_polygons,
              expected_num_of_point_3,
              expected_num_of_arr_arcs,
              expected_num_of_overlap_seg_3))
        {
           
          cout << change_color(CGAL_RED,"Iteration Error: ",iteration,"(");
          for (unsigned int ii = 0; ii < len-1; ii++)
            {
              cout << change_color(CGAL_RED,lines_index[ii],", ");
            }
          cout << change_color(CGAL_RED,lines_index[len-1],")") << endl;
          exit(EXIT_FAILURE);
        }
    }
   
  for (unsigned int ii = min; ii < len;ii++)
    {
      swap(lines.at(min),lines.at(ii));
      swap(lines_index[min],lines_index[ii]);
         
      run_all_permutations(line_through_segs,
          lines, lines_index, 
          len ,min + 1, ii != min,
          expected_num_of_lines,
          expected_num_of_arr_curves,
          expected_num_of_arr_polygons,
          expected_num_of_point_3,
          expected_num_of_arr_arcs,
          expected_num_of_overlap_seg_3);
         
      swap(lines.at(min),lines.at(ii));
      swap(lines_index[min],lines_index[ii]);
         
    }
}


template <typename Lines_through_segs>
void flip(
    Lines_through_segs& line_through_segs, 
    vector<Rational_segment_3> &lines,
    int &expected_num_of_lines,
    int &expected_num_of_arr_curves,
    int &expected_num_of_arr_polygons,
    int &expected_num_of_point_3,
    int &expected_num_of_arr_arcs,
    int &expected_num_of_overlap_seg_3,
    int &expected_num_of_lines_flip,
    int &expected_num_of_arr_curves_flip,
    int &expected_num_of_arr_polygons_flip,
    int &expected_num_of_point_3_flip,
    int &expected_num_of_arr_arcs_flip,
    int &expected_num_of_overlap_seg_3_flip)

{
  if (expected_num_of_lines_flip < 0)
    {
      ReportError("Unexpected error",__LINE__,__FILE__);
    }
   
#if LINES_DEBUG   
  cout <<"  get_all_common_lines iteration = " << 0 << endl;
#endif
  if (!get_all_common_lines(line_through_segs,
          lines,
          expected_num_of_lines,
          expected_num_of_arr_curves,
          expected_num_of_arr_polygons,
          expected_num_of_point_3,
          expected_num_of_arr_arcs,
          expected_num_of_overlap_seg_3))
    {
      ReportError("Error iteration 0",__LINE__,__FILE__);
    }
   

  unsigned int len = lines.size()/2;
  for (unsigned int ii = 0; ii < len;ii++)
    {
      swap(lines.at(ii),lines.at((lines.size()-1)-ii));
    }
   
  if (!get_all_common_lines(line_through_segs,
          lines,
          expected_num_of_lines_flip,
          expected_num_of_arr_curves_flip,
          expected_num_of_arr_polygons_flip,
          expected_num_of_point_3_flip,
          expected_num_of_arr_arcs_flip,
          expected_num_of_overlap_seg_3_flip))
    {
      ReportError("Error iteration 1",__LINE__,__FILE__);
    }
}

void create_random_input(int num_of_lines)
{
  std::ofstream outputfile;
  outputfile.open("input.txt");
  outputfile << num_of_lines << std::endl;
    
  for (int ii = 0; ii < num_of_lines; ii++)
    {
      for(int jj = 0;jj<5;jj++)
        {
          outputfile << std::rand() % 10000 << " ";
        }
      outputfile << std::rand() % 10000 << std::endl;
    }
}

int main (int argc,char **args)
{
    int expected_num_of_lines = -1;
    int expected_num_of_arr_curves = -1;
    int expected_num_of_arr_polygons = -1;
    int expected_num_of_point_3 = -1;
    int expected_num_of_arr_arcs = -1;
    int expected_num_of_overlap_seg_3 = -1;

    int expected_num_of_lines_flip = -1;
    int expected_num_of_arr_curves_flip = -1;
    int expected_num_of_arr_polygons_flip = -1;
    int expected_num_of_point_3_flip = -1;
    int expected_num_of_arr_arcs_flip = -1;
    int expected_num_of_overlap_seg_3_flip = -1;

    vector<Rational_segment_3> lines;
    int arr[40];
    for (int ii = 0; ii < 40; ii++)
      {
        arr[ii] = ii;
      }
      
    if (argc < 3)
      {
        ReportError("File name not specified\n",__LINE__,__FILE__);
      }
      
    Alg_kernel alg_kernel;
    Rational_kernel rat_kernel;
   
#if USE_CONIC_TRAITS      
    CGAL::Lines_through_segments_3<
    Lines_through_segs_traits_using_conic_2,With_segments> 
      line_through_segs_use_conic(alg_kernel,rat_kernel);
#endif
#if (USE_RATIONAL_ARC_TRAITS || USE_SQRT_TRAITS)
    CGAL::Lines_through_segments_3<
    Lines_through_segs_traits_using_rational_arc_2,With_segments> 
      line_through_segs_use_rat_arc(alg_kernel,rat_kernel);
#endif

    if (argc >= 2 && strncmp(args[2],"run_perm",8) == 0)
      {
        ReadInputFile<Rational_segment_3,Rational,Rational_point_3>(args[1],lines,
            expected_num_of_lines,
            expected_num_of_arr_curves,
            expected_num_of_arr_polygons,
            expected_num_of_point_3,
            expected_num_of_arr_arcs,
            expected_num_of_overlap_seg_3);


#if USE_CONIC_TRAITS
        run_all_permutations(line_through_segs_use_conic,
            lines, arr, lines.size(), 0,
            true,
            expected_num_of_lines,
            expected_num_of_arr_curves,
            expected_num_of_arr_polygons,
            expected_num_of_point_3,
            expected_num_of_arr_arcs,
            expected_num_of_overlap_seg_3);
#endif
#if (USE_RATIONAL_ARC_TRAITS || USE_SQRT_TRAITS)
        run_all_permutations(line_through_segs_use_rat_arc,
            lines, arr, lines.size(), 0,
            true,
            expected_num_of_lines,
            expected_num_of_arr_curves,
            expected_num_of_arr_polygons,
            expected_num_of_point_3,
            expected_num_of_arr_arcs,
            expected_num_of_overlap_seg_3);
#endif

      }
    else if (argc >= 2 && strncmp(args[2],"flip",4) == 0)
      {
        ReadInputFileFlip<Rational_segment_3,Rational,Rational_point_3>(args[1],lines,
            expected_num_of_lines,
            expected_num_of_arr_curves,
            expected_num_of_arr_polygons,
            expected_num_of_point_3,
            expected_num_of_arr_arcs,
            expected_num_of_overlap_seg_3,
            true,
            expected_num_of_lines_flip,
            expected_num_of_arr_curves_flip,
            expected_num_of_arr_polygons_flip,
            expected_num_of_point_3_flip,
            expected_num_of_arr_arcs_flip,
            expected_num_of_overlap_seg_3_flip);

#if USE_CONIC_TRAITS
        flip(line_through_segs_use_conic,lines,
            expected_num_of_lines,
            expected_num_of_arr_curves,
            expected_num_of_arr_polygons,
            expected_num_of_point_3,
            expected_num_of_arr_arcs,
            expected_num_of_overlap_seg_3,
            expected_num_of_lines_flip,
            expected_num_of_arr_curves_flip,
            expected_num_of_arr_polygons_flip,
            expected_num_of_point_3_flip,
            expected_num_of_arr_arcs_flip,
            expected_num_of_overlap_seg_3_flip);
#endif
#if (USE_RATIONAL_ARC_TRAITS || USE_SQRT_TRAITS)
        flip(line_through_segs_use_rat_arc,lines,
            expected_num_of_lines,
            expected_num_of_arr_curves,
            expected_num_of_arr_polygons,
            expected_num_of_point_3,
            expected_num_of_arr_arcs,
            expected_num_of_overlap_seg_3,
            expected_num_of_lines_flip,
            expected_num_of_arr_curves_flip,
            expected_num_of_arr_polygons_flip,
            expected_num_of_point_3_flip,
            expected_num_of_arr_arcs_flip,
            expected_num_of_overlap_seg_3_flip);
#endif
      }
    else
      {
        ReadInputFile<Rational_segment_3,Rational,Rational_point_3>(args[1],lines,
            expected_num_of_lines,
            expected_num_of_arr_curves,
            expected_num_of_arr_polygons,
            expected_num_of_point_3,
            expected_num_of_arr_arcs,
            expected_num_of_overlap_seg_3);
  
        cout <<"get_all_common_lines iteration = " << iteration << endl;
#if USE_CONIC_TRAITS
        if(!get_all_common_lines(
             line_through_segs_use_conic,
             lines,
             expected_num_of_lines,
             expected_num_of_arr_curves,
             expected_num_of_arr_polygons,
             expected_num_of_point_3,
             expected_num_of_arr_arcs,
             expected_num_of_overlap_seg_3))
          return EXIT_FAILURE;
#endif
        // std::cout << "Num of lines = " << expected_num_of_lines << std::endl;
        // std::cout << "Num of curves = " << expected_num_of_arr_curves << std::endl;
        // std::cout << "Num of polygons = " << expected_num_of_arr_polygons << std::endl;
        // std::cout << "Num of overlap segs = " << expected_num_of_overlap_seg_3 << std::endl;
        // std::cout << "Num of points = " << expected_num_of_point_3 << std::endl;
        // std::cout << "Num of arcs = " << expected_num_of_arr_arcs << std::endl;

#if (USE_RATIONAL_ARC_TRAITS || USE_SQRT_TRAITS)
        expected_num_of_lines = -1;
        expected_num_of_arr_curves = -1;
        expected_num_of_arr_polygons = -1;
        expected_num_of_point_3 = -1;
        expected_num_of_arr_arcs = -1;
        expected_num_of_overlap_seg_3 = -1;
        if(!get_all_common_lines(
             line_through_segs_use_rat_arc,
             lines,
             expected_num_of_lines,
             expected_num_of_arr_curves,
             expected_num_of_arr_polygons,
             expected_num_of_point_3,
             expected_num_of_arr_arcs,
             expected_num_of_overlap_seg_3))
          return EXIT_FAILURE;
#endif
     
        // std::cout << "Num of lines = " << expected_num_of_lines << std::endl;
        // std::cout << "Num of curves = " << expected_num_of_arr_curves << std::endl;
        // std::cout << "Num of polygons = " << expected_num_of_arr_polygons << std::endl;
        // std::cout << "Num of overlap segs = " << expected_num_of_overlap_seg_3 << std::endl;
        // std::cout << "Num of points = " << expected_num_of_point_3 << std::endl;
        // std::cout << "Num of arcs = " << expected_num_of_arr_arcs << std::endl;

      }
      
    cout << "Program finished successfully" << endl;
 
    return EXIT_SUCCESS;
}

// {0*x^2 + 0*y^2 + -690*xy + 650*x + 480*y + -433} : (0.6435,0.409195) --cw--> (0.666154,0)
// 1 1 1 5 5 5
// {0*x^2 + 0*y^2 + -690*xy + 650*x + 480*y + -433} : (0.538914,0.764738) --cw--> (0.666154,0)
// 3 3 3 7 7 7
// before_create_edge
// c = {0*x^2 + 0*y^2 + -690*xy + 650*x + 480*y + -433} : (0.538914,0.764738) --cw--> (0.6435,0.409195)
// c.data() = 0x11050a8
// m_last_inserted_segment = 3 3 3 7 7 7
// after_create_edge
// curve = ({0*x^2 + 0*y^2 + -690*xy + 650*x + 480*y + -433} : (0.538914,0.764738) --cw--> (0.6435,0.409195))
// before_create_edge
// c = {0*x^2 + 0*y^2 + -690*xy + 650*x + 480*y + -433} : (0.6435,0.409195) --cw--> (0.666154,0)
// c.data() = 0x11050a8
// m_last_inserted_segment = 3 3 3 7 7 7
// after_create_edge
// curve = ({0*x^2 + 0*y^2 + -690*xy + 650*x + 480*y + -433} : (0.6435,0.409195) --cw--> (0.666154,0))
