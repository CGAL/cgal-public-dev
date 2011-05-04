// Offset_types_2.h

#ifndef Offset_types_2_h
#define Offset_types_2_h

#include <CGAL/approximated_offset_2.h>
#include <CGAL/offset_polygon_2.h>

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/General_polygon_2.h>
#include <iostream>

namespace CGAL {


/*! \classf
 * A class containing common offset polygon types.
 */

template <typename CPolygon_2>
class Offset_types_2
{
public:

  Offset_types_2(void)
  {
  }

  ~Offset_types_2(void)
  {
  }

  /// input polygon types
  typedef CPolygon_2 Polygon_2;
  typedef typename Polygon_2::Traits  Polygon_traits_2;
  typedef typename Polygon_traits_2::FT       Input_rational;
  typedef typename CGAL::Polygon_with_holes_2<Polygon_traits_2>  Polygon_with_holes_2;
  typedef Polygon_with_holes_2 Kgon_sum_polygon_2;
  typedef typename Polygon_traits_2::Point_2   Point_2;
  typedef typename Polygon_traits_2::Segment_2   Segment_2;
  typedef typename Polygon_traits_2::Direction_2   Direction_2;
  typedef typename CGAL::Aff_transformation_2<Polygon_traits_2>  Transformation_2;

  /// approximate (rational) traits 
  typedef typename CGAL::Gps_circle_segment_traits_2<Polygon_traits_2>  Rat_traits_2;
  typedef typename Rat_traits_2::Polygon_2                    Approximate_polygon_2;
  typedef typename Rat_traits_2::Polygon_with_holes_2         Approximate_offset_polygon_2;
  typedef typename std::list<Approximate_polygon_2> Approximate_polygon_list_2;

  /// exact traits 

  typedef typename CGAL::CORE_algebraic_number_traits     Nt_traits;
  typedef typename Nt_traits::Rational                    Exact_rational;
  typedef typename Nt_traits::Algebraic                   Exact_algebraic;

// instead of
//typedef CGAL::Cartesian<Rational>              Rat_kernel;
//typedef CGAL::Cartesian<Algebraic>             Alg_kernel;
//typedef CGAL::Arr_conic_traits_2<Rat_kernel,
//                                 Alg_kernel,
//                                 Nt_traits>    Conic_traits_2;
// workaround for VC++
  struct Exact_rat_kernel : public CGAL::Cartesian<Exact_rational> {};
  struct Exact_alg_kernel : public CGAL::Cartesian<Exact_algebraic> {};
  struct Conic_traits_2 : public CGAL::Arr_conic_traits_2<Exact_rat_kernel,
                                   Exact_alg_kernel,
              Nt_traits> {};

  typedef typename CGAL::Gps_traits_2<Conic_traits_2>     Exact_traits_2;


  typedef CGAL::Polygon_2<Exact_rat_kernel>   Rat_polygon_2;
  typedef CGAL::Polygon_with_holes_2<Exact_rat_kernel>   Rat_polygon_with_holes_2;

  typedef typename Exact_traits_2::Polygon_2                Exact_polygon_2;
  typedef typename Exact_traits_2::Polygon_with_holes_2     Exact_offset_polygon_2;

  typedef typename std::list<Exact_polygon_2> Exact_polygon_list_2;
  typedef typename std::list<Exact_offset_polygon_2> Exact_offset_polygon_list_2;

private:


};

template<typename OT, typename IT>
OT type_cast(const IT& it)
{
//  return OT(it);
  std::stringstream string_str;
  string_str << it;
  CORE::BigRat ot;
  string_str >> ot;

//  LOG_DEBUG << ot << " = type_cast " << it << std::endl;

  return ot;
}
/*
template<>
CORE::BigRat type_cast<CORE::BigRat>(const CGAL::Gmpq& it)
{
  std::stringstream string_str;
  string_str << it;
  CORE::BigRat ot;
  string_str >> ot;

  LOG_DEBUG << ot << " = type_cast " << it << std::endl;

  return ot;
}
*/

template <typename PolygonType1, typename PolygonType2>
void copy(const PolygonType1& p1, PolygonType2& p2)
{
    typedef typename PolygonType1::Vertex_const_iterator Vertex_const_iterator1;
    typedef typename PolygonType1::Traits::Point_2 Point1_2;
    typedef typename PolygonType2::Traits::Point_2 Point2_2;
    typedef typename PolygonType1::Traits::FT FT1;
    typedef typename PolygonType2::Traits::FT FT2;

    p2.clear();

    Vertex_const_iterator1 vi = p1.vertices_begin();
    Vertex_const_iterator1 vi_end = p1.vertices_end();

    for(; vi != vi_end; ++vi)
    {
        //Point1_2 point1 = *vi->point();
        p2.push_back(Point2_2(type_cast<FT2>(vi->x()), type_cast<FT2>(vi->y())));
    }
}

} // namespace CGAL


#endif // Offset_types_2_h
