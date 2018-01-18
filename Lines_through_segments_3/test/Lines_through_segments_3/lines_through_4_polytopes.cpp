
#define CGAL_LEDA_VERSION 620

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <boost/variant.hpp>
#include <string>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <vector>
#include <algorithm>
#include <CGAL/Timer.h>

#define USE_LEDA 0
/* Set the required traits */
#define USE_CONIC_TRAITS 1
#define USE_RATIONAL_ARC_TRAITS 0

#include "color_print.h"
// #include "Arrangement_general_functions.h"

#define LTS_DRAW_ARR 1
#define ARR_ON_SUR_DEBUG 0
#define LINES_DEBUG 0
#define ISOLATED_POINTS_DEBUG 0
#define BOUNDED_SEGMENTS_VECTOR_DEBUG 0
#define OBSERVER_PRINTS 0
#define OBJ_ON_ARR_DEBUG 0
#define CGAL_DEBUG_OUTPUT 1
#define CGAL_LTS_MEASURE_TIME 0
#define LTS_POLY_IMPL_DEBUG 0
#define LTS_POLY_CON_COMP_DEBUG 0
#define LTS_POLY_OVERLAY_DEBUG 0

#define PRINT_OUTPUT 0
#define CGAL_IDENTIFICATION_XY CGAL_X_MINUS_11_Y_7

#define LTS_WITH_SEGMENTS 0

#if USE_CONIC_TRAITS
#include <CGAL/Arr_conic_traits_2.h>
#endif
#if USE_RATIONAL_ARC_TRAITS
#include <CGAL/Arr_rational_arc_traits_d_1.h>
#include <CGAL/Arr_vertical_segment_traits.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#endif

#include <CGAL/Arrangement_2.h>

#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>

#include <CGAL/CORE_algebraic_number_traits.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

/* Lines_through_segments includes */
#include <CGAL/IO/Lines_through_segments_io_stream.h>
#include <CGAL/Lines_through_polytopes_3/Lines_through_polytopes_3.h>
#include <CGAL/Lines_through_segments_traits_3.h>


typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Algebraic                            Algebraic;
typedef Nt_traits::Rational                             Rational;

typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Cartesian<Rational>                       Rational_kernel;
typedef CGAL::Cartesian<double>                         Double_kernel;
typedef CGAL::Segment_3<Rational_kernel>                Rational_segment_3;
typedef CGAL::Polyhedron_3<Double_kernel>               Double_polyhedron_3;
typedef CGAL::Polyhedron_3<Rational_kernel>             Rational_polyhedron_3;
typedef CGAL::Line_3<Rational_kernel>                   Rational_line_3;
typedef CGAL::Point_3<Rational_kernel>                  Rational_point_3;
typedef CGAL::Point_2<Rational_kernel>                  Rational_point_2;
typedef CGAL::Segment_2<Rational_kernel>                Rational_segment_2;

typedef CGAL::Point_3<Alg_kernel>                 Alg_point_3;
typedef CGAL::Point_2<Alg_kernel>                 Alg_point_2;
typedef CGAL::Line_3<Alg_kernel>                  Alg_line_3;
   
   
#if USE_CONIC_TRAITS
typedef CGAL::Arr_conic_traits_2<Rational_kernel, Alg_kernel, Nt_traits>    Conic_traits_arr_on_plane_2;
#endif
#if USE_RATIONAL_ARC_TRAITS
typedef CGAL::Arr_rational_arc_traits_d_1<Rational_kernel>                  Traits_d_1;
typedef CGAL::Arr_traits_with_vertical_segments <Traits_d_1>                Rational_arc_arr_traits_arr_on_plane_2;
#endif
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Alg_kernel>               Traits_arr_on_sphere_2;
#if USE_CONIC_TRAITS
typedef CGAL::Lines_through_segments_traits_3<Alg_kernel,
                                              Rational_kernel,
                                              Conic_traits_arr_on_plane_2,
                                              Traits_arr_on_sphere_2> Lines_through_segs_traits_using_conic_2;
#endif
#if USE_RATIONAL_ARC_TRAITS
typedef CGAL::Lines_through_segments_traits_3<Alg_kernel,
                                              Rational_kernel,
                                              Rational_arc_arr_traits_arr_on_plane_2,
                                              Traits_arr_on_sphere_2> Lines_through_segs_traits_using_rational_arc_2;
#endif
using namespace std;


#define BOLD 1
#define DARK 2
#define UNDERLINE 4
#define BLINK 5
#define NEGATIVE 7
#define BLACK 30
#define RED 31
#define GREEN 32
#define YELLOW 33
#define BLUE 34
#define MAGNETA 35
#define CYAN 36
#define WHITE 37


#if 1
std::string changeColorRec(std::ostringstream &o, int color)
{
   o << "\033[1;" << color << "m" << "" << "\033[0m";
   return o.str();
}

template <typename T, typename... Args>
std::string changeColorRec(std::ostringstream &o, int color, T &value, Args&... args)
{
   o << "\033[1;" << color << "m" << value << "\033[0m";
   return changeColorRec(o, color, args...);
}

template<typename T, typename... Args>
std::string changeColor(int color, T &value,Args&... args)
{
   std::ostringstream o;
   o << "\033[1;" << color << "m" << value << "\033[0m";
   return changeColorRec(o, color, args...);
}
#else
template<typename T, typename... Args>
std::string changeColor(int color, T &value,Args&... args)
{
   std::ostringstream o;
   o << "\033[1;" << color << "m" << value << "\033[0m";
   return changeColorRec(o, color, args...);
}

#endif

void ReportError(std::string str, int line, std::string file)
{
#if 1
   std::cout << changeColor(RED,str) << std::endl;
   std::cout << changeColor(RED,"File = ", file, "\tLine = ", line) << std::endl;
#else
   std::cout << str << std::endl;
   std::cout << "File = " << file << "\tLine = " << line << std::endl;
#endif
   exit(0);
}

/* Parse line at the inputfile */
#define LINE_SIZE 400

/*************************************************************
 * The following function parses an input file into array of lines.
 *
 * Input:
 *      input_file_name - A file with list of off files of polyhedrons.
 *       
 * Output:
 *      polytopes -           Vector of polytops.
 *
 *
 **************************************************************/
template <typename Polyhedron_3>
int ReadInputFile(char* input_file_name,vector<Polyhedron_3> &polytopes)
{
   int ii = 0;
   int num_of_lines = 0;
   char line[LINE_SIZE];
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
      Polyhedron_3 P;
      input_file.getline(line,sizeof(line));
      std::cout << line << std::endl;
      
      ifstream input_polyhedron(line);
      if (input_polyhedron.fail())
      {
         ReportError("Error opening file",__LINE__,__FILE__);
      }
      input_polyhedron >> P;
      polytopes.push_back(P);
      std::cout << polytopes[ii] << std::endl;
      
      ii++;
      input_polyhedron.close();
   }
   
   if (ii != num_of_lines)
   {
      ReportError("Two many lines in file",__LINE__,__FILE__);
   }

   input_file.close();
   return num_of_lines;
}

template<typename LTS_output_obj>
class My_front_insert_iterator
{
   typedef list<LTS_output_obj> My_container;
protected:
   My_container* container;
      
public:
      
   /* The explicit function specifier controls unwanted implicit type conversions. 
      It can only be used in declarations of constructors within a class declaration. 
      For example, except for the default constructor, 
      the constructors in the following class are converting constructors.*/
   explicit My_front_insert_iterator(My_container& x)
   {
      container = &x;
   }
      
   My_front_insert_iterator& operator= (typename My_container::const_reference value)
   { 
      container->push_front(value);
      return *this;
   }

   My_front_insert_iterator& operator* ()
   { 
      return *this;
   }

   My_front_insert_iterator& operator++ (int)
   { 
      return *this;
   }

   My_front_insert_iterator& operator++ ()
   { 
      return *this;
   }
};


template <typename Lines_through_poly>
bool get_all_common_lines(
   Lines_through_poly& line_through_poly,
   vector<Rational_polyhedron_3> &polytopes)
{
   int num_of_lines = 0;
   int num_of_arr_curves = 0;
   int num_of_arr_polygons = 0;
   int num_of_points = 0;
   int num_of_arr_arcs = 0;
   int num_of_overlap_seg_3 = 0;

   typedef typename Lines_through_poly::Transversal LTS_output_obj;
   
   list<LTS_output_obj> output_list;
   
#if CGAL_LTS_MEASURE_TIME
   CGAL::Timer t;
   t.reset();
   t.start();
// for (int ii = 0;ii<1;ii++)
// {
//    output_list.clear();
#endif
   line_through_poly(
      polytopes.begin(),
      polytopes.end(),
      My_front_insert_iterator<LTS_output_obj>(output_list));
#if CGAL_LTS_MEASURE_TIME
// }
   t.stop();
   std::cout << t.time() << std::endl;
#endif
   
   typename list<LTS_output_obj>::iterator it_output_list;
#if PRINT_OUTPUT
   cout<<"OUTPUT:" << endl;
#endif
   
   typename Lines_through_poly::Mapped_general_polygon_2   *polygon_obj;
   typename Lines_through_poly::Mapped_x_monotone_curve_2  *curve_obj;
   typename Lines_through_poly::Mapped_point_2             *arr_point_2_obj;
   typename Lines_through_poly::Line_3                     *line_obj;
   typename Lines_through_poly::Through_point_3            *point_obj;
   typename Lines_through_poly::Through_point_3_segment_3  *arc_obj;
   typename Lines_through_poly::Through_segment_3          *seg_obj;

   typename Lines_through_poly::Mapped_2                   *mapped_obj;
   typename Lines_through_poly::Through_3                  *through_obj;
   
   typename Lines_through_poly::Transversal                   transversal;
   typename Lines_through_poly::Transversal_with_segments     *transversalw;
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
      
      if ((line_obj = 
           boost::get < typename Lines_through_poly::Line_3 > (&transversal)))
      {
         num_of_lines++;
      } else if ((mapped_obj = boost::get<typename Lines_through_poly::Mapped_2>(&transversal)))
      {
         typename Lines_through_poly::Mapped_2::Mapped_line_3 line = mapped_obj->line();
         typename Lines_through_poly::Mapped_transversal mapped_transversal = mapped_obj->mapped_transversal();
         if ((curve_obj = 
             boost::get<typename Lines_through_poly::Mapped_x_monotone_curve_2> 
              (&mapped_transversal)))
         {
            num_of_arr_curves++;
         }
         else if (polygon_obj = 
                  boost::get<typename Lines_through_poly::Mapped_general_polygon_2> (&mapped_transversal))
         {
            num_of_arr_polygons++;
         }
         else if (arr_point_2_obj = 
                  boost::get<typename Lines_through_poly::Mapped_point_2> (&mapped_transversal))
         {
            num_of_lines++;
         }
         else
         {
            ReportError("Unexpected Error - invalid mapped obj value",__LINE__,__FILE__);
         }
      }
      else if (through_obj = boost::get<typename Lines_through_poly::Through_3> 
               (&transversal))
      {
         typename Lines_through_poly::Through_transversal through_transversal = through_obj->through_transversal();
         if (arc_obj = 
             boost::get<typename Lines_through_poly::Through_point_3_segment_3> 
             (&through_transversal))
         {
            num_of_arr_arcs++;
         }
         else if (seg_obj = 
                  boost::get<typename Lines_through_poly::Through_segment_3> (&through_transversal))
         {
            num_of_overlap_seg_3++;
         }
         else if (point_obj = 
                  boost::get<typename Lines_through_poly::Through_point_3> (&through_transversal))
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
   
   return true;
}


int main (int argc,char **args)
{
   vector<Rational_polyhedron_3> polytopes;
   int arr[40];
   for (int ii = 0; ii < 40; ii++)
   {
      arr[ii] = ii;
   }
      
   if (argc < 3)
   {
      ReportError("File name not specified\n",__LINE__,__FILE__);
   }
   
   ReadInputFile(args[1],polytopes);
#if 0
   typedef Double_polyhedron_3::Edge_iterator   Edge_iterator;
   for (Edge_iterator e_p1 = polytopes[0].edges_begin(); 
        e_p1 != polytopes[0].edges_end();
        ++e_p1)
   {
      Rational_point_3 p1(
         Rational(e_p1->vertex()->point().x()),
         Rational(e_p1->vertex()->point().y()),
         Rational(e_p1->vertex()->point().z()));
      
      std::cout << "Edge\ne_p1 = "
                << e_p1->vertex()->point().x()
                << " " 
                << e_p1->vertex()->point().y()
                << " " 
                << e_p1->vertex()->point().z()
                << " " 
                << e_p1->opposite()->vertex()->point().x()
                << " " 
                << e_p1->opposite()->vertex()->point().y()
                << " " 
                << e_p1->opposite()->vertex()->point().z()
                << std::endl;
      std::cout << "p1 = " << CGAL::to_double(p1.x()) << std::endl;
   }
#endif

   Alg_kernel alg_kernel;
   Rational_kernel rat_kernel;

#if USE_CONIC_TRAITS      
   CGAL::Lines_through_polytopes_3<Lines_through_segs_traits_using_conic_2> 
      line_through_polys_use_conic(alg_kernel,rat_kernel);
   get_all_common_lines(
      line_through_polys_use_conic,
      polytopes);

#endif
#if USE_RATIONAL_ARC_TRAITS
   CGAL::Lines_through_polytopes_3<Lines_through_segs_traits_using_rational_arc_2> 
      line_through_polys_use_rat_arc(alg_kernel,rat_kernel);
   get_all_common_lines(
      line_through_polys_use_rat_arc,
      polytopes);

#endif

   
   cout << "Program finished successfully" << endl;
   
   return 0;
}

