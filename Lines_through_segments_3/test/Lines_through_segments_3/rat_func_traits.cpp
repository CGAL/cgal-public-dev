
#include <iostream>
#include <fstream>
#include <stdio.h>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <vector>

#include <CGAL/Arr_rational_function_traits_2.h>
#include <CGAL/Algebraic_kernel_d_1.h>

#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arrangement_on_surface_with_history_2.h>

#include <CGAL/CORE_algebraic_number_traits.h>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;

typedef CGAL::Cartesian<Rational>                       Rational_kernel;
typedef CORE::BigInt                                    Integer;   
typedef CGAL::Algebraic_kernel_d_1<Integer>	   AK1;   
typedef CGAL::Arr_rational_function_traits_2<AK1>                  Rational_arc_arr_traits_arr_on_plane_2;


using namespace std;
typedef Rational_arc_arr_traits_arr_on_plane_2::Curve_2 Curve_2;

Curve_2 construct_curve(double start, double end,
                     const char* num, const char* den)
{
   //typedef CORE::Polynomial<Rational>                    PNT_1;//
   typedef Rational_arc_arr_traits_arr_on_plane_2::Polynomial_1 PNT_1;
   typedef Rational_arc_arr_traits_arr_on_plane_2::Algebraic_real_1 NT;

   Rational_arc_arr_traits_arr_on_plane_2 traits_2;
   std::stringstream ss_num(num);
   std::stringstream ss_den(den);
   PNT_1 f_num;
   ss_num >> f_num;

   PNT_1 f_den;
   ss_den >> f_den;
   Rational rat_start(start);
   Rational rat_end(end);
   NT nt_start(rat_start);
   NT nt_end(rat_end);
   return traits_2.construct_curve_2_object()(f_num, f_den, nt_start, nt_end);
}

int main ()
{
   typedef CGAL::Arrangement_with_history_2<Rational_arc_arr_traits_arr_on_plane_2> Arrangement_2;
   typedef Rational_arc_arr_traits_arr_on_plane_2::Point_2 Point_2;
   std::list<Curve_2> arcs;
   Arrangement_2 arr;
// NEW ARRANGEMENT S1 = 5679 9630 620 6094 3900 5064
// S2 = 9710 1901 5718 4780 9163 3250

// Add Element to Arr = 9187 7768 4788 5516 8274 5945
// cv = _f = P[1(0,-36915044052873212461927352624251288)(1,25794859833267407254694964792965724)]
// _f = P[1(0,-37065654280263952377341959104571664)(1,9274260669763038094357576406681912)]
// y = (25794859833267407254694964792965724*x - 36915044052873212461927352624251288) / (9274260669763038094357576406681912*x - 37065654280263952377341959104571664) on [0.374422, 0.798057]

   Curve_2 cv = construct_curve((double)0.374422, (double)0.798057,
                        "P[1(0,-36915044052873212461927352624251288)(1,25794859833267407254694964792965724)]",
                        "P[1(0,-37065654280263952377341959104571664)(1,9274260669763038094357576406681912)]");
   std::cout << cv << std::endl;
   arcs.push_back(cv);
   insert (arr, arcs.begin(), arcs.end());
   std::cout << "Arr on plane arrangement size:" 
             << "   V = " << arr.number_of_vertices()
             << ",  E = " << arr.number_of_edges()
             << ",  F = " << arr.number_of_faces() << std::endl;

   arcs.clear();
   // Arr on plane arrangement size:   V = 2,  E = 1,  F = 1

// Add Element to Arr = 8968 3698 7217 3958 9143 1221
// cv = _f = P[1(0,5754200987313160936360708869680838400)(1,-11498247355274517569649510388690080000)]
// _f = P[1(0,9162519012574528914051000323161356800)(1,-17851487279937433953457938595684480000)]
// y = (-11498247355274517569649510388690080000*x + 5754200987313160936360708869680838400) / (-17851487279937433953457938595684480000*x + 9162519012574528914051000323161356800) on [0, 0.439515]
   cv = construct_curve((double)0, (double)0.439515,
                        "P[1(0,5754200987313160936360708869680838400)(1,-11498247355274517569649510388690080000)]",
                        "P[1(0,9162519012574528914051000323161356800)(1,-17851487279937433953457938595684480000)]");
   std::cout << cv << std::endl;
   arcs.push_back(cv);
// cv = _f = P[1(0,5754200987313160936360708869680838400)(1,-11498247355274517569649510388690080000)]
// _f = P[1(0,9162519012574528914051000323161356800)(1,-17851487279937433953457938595684480000)]
// y = (-11498247355274517569649510388690080000*x + 5754200987313160936360708869680838400) / (-17851487279937433953457938595684480000*x + 9162519012574528914051000323161356800) on [0.495669, 0.500442]
   cv = construct_curve((double)0.495669, (double)0.500442,
                        "P[1(0,5754200987313160936360708869680838400)(1,-11498247355274517569649510388690080000)]",
                        "P[1(0,9162519012574528914051000323161356800)(1,-17851487279937433953457938595684480000)]");
   std::cout << cv << std::endl;
   arcs.push_back(cv);
// cv = _f = P[1(0,5754200987313160936360708869680838400)(1,-11498247355274517569649510388690080000)]
// _f = P[1(0,9162519012574528914051000323161356800)(1,-17851487279937433953457938595684480000)]
// y = (-11498247355274517569649510388690080000*x + 5754200987313160936360708869680838400) / (-17851487279937433953457938595684480000*x + 9162519012574528914051000323161356800) on [0.536469, 1]
   cv = construct_curve((double)0.536469, (double)1,
                        "P[1(0,5754200987313160936360708869680838400)(1,-11498247355274517569649510388690080000)]",
                        "P[1(0,9162519012574528914051000323161356800)(1,-17851487279937433953457938595684480000)]");
   std::cout << cv << std::endl;
   arcs.push_back(cv);
   insert (arr, arcs.begin(), arcs.end());
   std::cout << "Arr on plane arrangement size:" 
             << "   V = " << arr.number_of_vertices()
             << ",  E = " << arr.number_of_edges()
             << ",  F = " << arr.number_of_faces() << std::endl;

   arcs.clear();

// Arr on plane arrangement size:   V = 21,  E = 14,  F = 1

// Add Element to Arr = 9058 6042 4290 9648 6028 1538
// cv = _f = P[1(0,-340592195143155629078557092473733120)(1,3989154696577795775073806545059840)]
// _f = P[1(0,-261512520774634962813299797335736320)(1,-528543499654802052134583887305113600)]
// y = (3989154696577795775073806545059840*x - 340592195143155629078557092473733120) / (-528543499654802052134583887305113600*x - 261512520774634962813299797335736320) on [0.998485, 1]
   cv = construct_curve((double)0.998485, (double)1,
                        "P[1(0,-340592195143155629078557092473733120)(1,3989154696577795775073806545059840)]",
                        "P[1(0,-261512520774634962813299797335736320)(1,-528543499654802052134583887305113600)]");
   std::cout << cv << std::endl;
   arcs.push_back(cv);
   insert (arr, arcs.begin(), arcs.end());
   std::cout << "Arr on plane arrangement size:" 
             << "   V = " << arr.number_of_vertices()
             << ",  E = " << arr.number_of_edges()
             << ",  F = " << arr.number_of_faces() << std::endl;

// Arr on plane arrangement size:   V = 6,  E = 3,  F = 1

   Arrangement_2::Edge_iterator   eit;
   for (eit = arr.edges_begin(); 
        eit != arr.edges_end();
        ++eit)
   {
      std::cout << "arr_on_plane.number_of_originating_curves(eit) = " <<
         arr.number_of_originating_curves(eit) << std::endl;
   }
   
   cout << "Program finished successfully" << endl;

   return 0;
}
