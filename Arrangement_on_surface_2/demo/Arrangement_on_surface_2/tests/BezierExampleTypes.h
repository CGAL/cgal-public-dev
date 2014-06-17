#ifndef BEZIER_EXAMPLE_TYPES_H
#define BEZIER_EXAMPLE_TYPES_H
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
#include "DcelTypes.h"

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
//typedef Nt_traits::Rational                             NT;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef Rat_kernel::Point_2                             Rat_point_2;
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                        Traits_2;
typedef Traits_2::Curve_2                               Bezier_curve_2;
typedef Dcel<Traits_2>                                  Bezier_dcel;
typedef CGAL::Arrangement_with_history_2<Traits_2, Bezier_dcel>      Arrangement_2;

typedef CGAL::Qt::ArrangementGraphicsItem< Arrangement_2 > BezierArrangementGraphicsItem;

#endif // BEZIER_EXAMPLE_TYPES_H
