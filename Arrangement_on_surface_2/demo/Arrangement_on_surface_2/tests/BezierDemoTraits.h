#ifndef BEZIER_DEMO_TRAITS_H
#define BEZIER_DEMO_TRAITS_H
#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#error "CGAL needs CORE."
#endif

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include "ArrangementGraphicsItem.h"
#include "ArrangementDemoTab.h"
#include "DcelTypes.h"

struct BezierDemoTraits
{
  typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
  typedef Nt_traits::Rational                             Rational;
  typedef Nt_traits::Algebraic                            Algebraic;

  typedef CGAL::Cartesian<Rational>                       Rat_kernel;
  typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;

  typedef Rat_kernel::Point_2                             Rat_point_2;
  typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                          Bezier_traits_2;
  typedef Bezier_traits_2::Curve_2                        Bezier_curve_2;
  typedef Dcel<Bezier_traits_2>                           Bezier_dcel;
  typedef CGAL::Arrangement_with_history_2<Bezier_traits_2, Bezier_dcel>
    Bezier_arrangement_2;

  typedef CGAL::Qt::ArrangementGraphicsItem< Bezier_arrangement_2 >
    BezierArrangementGraphicsItem;

  typedef ArrangementDemoTab< Bezier_arrangement_2 > BezierTab;

  // exports
  typedef Bezier_traits_2 ArrTraitsType;
  typedef Bezier_arrangement_2 ArrangementType;
};
#endif //BEZIER_DEMO_TRAITS_H
