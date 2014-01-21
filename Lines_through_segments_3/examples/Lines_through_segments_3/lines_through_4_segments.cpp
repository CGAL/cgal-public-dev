#define CGAL_IDENTIFICATION_XY CGAL_X_MINUS_11_Y_7

#include <boost/variant.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/IO/Lines_through_segments_io_stream.h>
#include <CGAL/Lines_through_segments_3.h>
#include <CGAL/Lines_through_segments_traits_3.h>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Algebraic                            Algebraic;
typedef Nt_traits::Rational                             Rational;

typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef Rat_kernel::Line_3                              Rat_line_3;
typedef Rat_kernel::Segment_3                           Rat_segment_3;
typedef Rat_kernel::Point_3                             Rat_point_3;
typedef Rat_kernel::FT                                  Rat_number_type;

typedef boost::true_type With_segments;
typedef CGAL::Lines_through_segments_traits_3<Alg_kernel, Rat_kernel>
                                                        Traits_3;

typedef CGAL::Lines_through_segments_3<Traits_3,With_segments>
  Lines_through_segments_3;
typedef Lines_through_segments_3::Transversal_with_segments
  Transversal_with_segments;

int main()
{
  Alg_kernel alg_kernel;
  Rat_kernel rat_kernel;
  Rat_number_type _5_3 = Rat_number_type(5) / Rat_number_type(3);
  std::vector<Rat_segment_3> segments(4);
   
  std::cout << "INPUT 1:" << std::endl;
  std::cout << "--------" << std::endl;
   
  /* Input1: Two pairs of intersecting segments. */
  /* Output1: A line that passes through the 2 intersection points.   */
  segments[0] = Rat_segment_3(Rat_point_3(-2, 12, 0), Rat_point_3(11, 1, 0));
  segments[1] = Rat_segment_3(Rat_point_3(-2, -5, 0), Rat_point_3(11, 3, 0));
  segments[2] = Rat_segment_3(Rat_point_3(0, 2, 8), Rat_point_3(0, 11, 1));
  segments[3] = Rat_segment_3(Rat_point_3(0, 2, -1), Rat_point_3(0, 11, 3));

  std::list<Transversal_with_segments> output_list1;
  Lines_through_segments_3 lines_through_segs(alg_kernel,
                                              rat_kernel);
  lines_through_segs(segments.begin(), 
                     segments.end(),
                     std::back_inserter(output_list1), true);
  
  copy(output_list1.begin(), output_list1.end(),
       std::ostream_iterator<Transversal_with_segments>(std::cout, "\n"));

  typedef Lines_through_segments_3::Arr_on_plane Arr_on_plane;
  typedef Lines_through_segments_3::Arr_on_sphere Arr_on_sphere;
  {
    Lines_through_segments_3::Arr_on_plane_iterator b, e;
    for(boost::tie(b, e) = lines_through_segs.planar_arrangements(); b != e; ++b) {
      std::cout << "Planar Arrangment" << std::endl;
      for(Arr_on_plane::Halfedge_iterator it = (*b)->halfedges_begin(); it != (*b)->halfedges_end(); ++it) { 
        std::cout << "  Halfedge: " << it->curve() << std::endl;
        for(Arr_on_plane::Halfedge::const_iterator it2 = it->segs_begin(); it2 != it->segs_end(); ++it2) { 
          std::cout << "    Segment: " << **it2 << std::endl;
        }
      }
    }
  }

  {
    Lines_through_segments_3::Arr_on_sphere_iterator b, e;
    for(boost::tie(b, e) = lines_through_segs.spherical_arrangements(); b != e; ++b) {
      std::cout << "Spherical Arrangment" << std::endl;
      for(Arr_on_sphere::Halfedge_iterator it = (*b)->halfedges_begin(); it != (*b)->halfedges_end(); ++it) { 
        std::cout << "  Halfedge: " << it->curve() << std::endl;
        for(Arr_on_sphere::Halfedge::const_iterator it2 = it->segs_begin(); it2 != it->segs_end(); ++it2) { 
          std::cout << "    Segment: " << **it2 << std::endl;
        }
      }
    }
  }

  std::cout << "INPUT 2:" << std::endl;
  std::cout << "--------" << std::endl;
   
  /* Input2: Triangle with a perpendicular segment at the middle.
   * Output2: Three lines, each through the perpendicular segment and vertex
   * of the triangle
   */
  segments[0] = Rat_segment_3(Rat_point_3(2, 1, 5), Rat_point_3(4, 3, 7));
  segments[1] = Rat_segment_3(Rat_point_3(4, 3, 7), Rat_point_3(6, 10, 12));
  segments[2] = Rat_segment_3(Rat_point_3(6, 10, 12), Rat_point_3(2, 1, 5));
  segments[3] = Rat_segment_3(Rat_point_3(8, 11, -2), Rat_point_3(0, - _5_3, 18));

  std::list<Transversal_with_segments> output_list2;
  lines_through_segs(segments.begin(), segments.end(),
                     std::back_inserter(output_list2), true);

  std::cout << output_list2.size() << std::endl;

  typedef Lines_through_segments_3::Line_3      Line_3;
  typedef Lines_through_segments_3::Mapped_2    Mapped_2;
  typedef Lines_through_segments_3::Through_3   Through_3;

  Line_3* line_obj;
  Mapped_2* mapped_obj;
  Through_3* through_obj;
  
  std::list<Transversal_with_segments>::iterator it;
  for (it = output_list2.begin(); it != output_list2.end(); ++it) {
    Lines_through_segments_3::Transversal transversal = it->first;
    if (line_obj = boost::get<Line_3>(&transversal)) {
      std::cout << "line = " << *line_obj << std::endl;
    }
    else if (mapped_obj = boost::get<Mapped_2>(&transversal)) {
      typedef Lines_through_segments_3::Mapped_general_polygon_2
        Mapped_general_polygon_2;
      typedef Lines_through_segments_3::Mapped_x_monotone_curve_2
        Mapped_x_monotone_curve_2;
      typedef Lines_through_segments_3::Mapped_point_2
        Mapped_point_2;

      Mapped_general_polygon_2* polygon_obj;
      Mapped_x_monotone_curve_2* curve_obj;
      Mapped_point_2* point_obj;
      Lines_through_segments_3::Mapped_transversal mapped_transversal =
        mapped_obj->mapped_transversal();
      // if (rat_point_obj =
      //     boost::get<Mapped_rat_point_2>(&mapped_transversal))
      // {
      //    Rat_line_3 line = mapped_obj->rational_line();
      //    std::cout << "Mapped_rat_point_2 = " << *rat_point_obj << ", ";
      //    std::cout << "Line = " << line << std::endl;
      // }
      // else
      {
         Mapped_2::Mapped_line_3 line = mapped_obj->line();
         if (curve_obj =
             boost::get<Mapped_x_monotone_curve_2>(&mapped_transversal))
         {
            std::cout << "Mapped_x_monotone_curve_2 = " << *curve_obj << ", ";
         }
         else if (polygon_obj = 
                  boost::get<Mapped_general_polygon_2>(&mapped_transversal))
         {
            std::cout << "Mapped_general_polygon_2 = " << *polygon_obj << ", ";
         }
         else if (point_obj = boost::get<Mapped_point_2>(&mapped_transversal))
         {
            std::cout << "Mapped_point_2 = " << *point_obj<< ", ";
         }
         std::cout << "Line = " << line << std::endl;
      }
      
    }
    else if (through_obj = boost::get<Through_3>(&transversal)) {
      typedef Lines_through_segments_3::Through_point_3
        Through_point_3;
      typedef Lines_through_segments_3::Through_point_3_segment_3
        Through_point_3_segment_3;
      typedef Lines_through_segments_3::Through_segment_3
        Through_segment_3;

      Through_point_3* point_obj;
      Through_point_3_segment_3* arc_obj;
      Through_segment_3* seg_obj;

      Lines_through_segments_3::Through_transversal through_transversal =
        through_obj->through_transversal();
      if (arc_obj = boost::get<Through_point_3_segment_3>(&through_transversal))
      {
        std::cout << "Through_point_3_segment_3 = (" << arc_obj->first
                  << "," << arc_obj->second << ")" << std::endl;
      }
      else if (seg_obj = boost::get<Through_segment_3>(&through_transversal))
      {
        std::cout << "Through_segment_3 = " << *seg_obj << std::endl;
      }
      else if (point_obj = boost::get<Through_point_3>(&through_transversal))
      {
        std::cout << "Through_point_3 = " << *point_obj << std::endl;
      }
    }
  }
   
  return 0;
}
