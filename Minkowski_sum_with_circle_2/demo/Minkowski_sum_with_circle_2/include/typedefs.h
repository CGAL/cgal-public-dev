#ifndef POLYGON_OFFSET_2_TYPEDEFS_H
#define POLYGON_OFFSET_2_TYPEDEFS_H

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
//#include <CGAL/minkowski_sum_2.h>
//#include <CGAL/Circle_approximation_2.h>
#include <CGAL/Minkowski_sum_with_circle_2.h>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/PolygonGraphicsItem.h>
#include <CGAL/Qt/PolygonWithHolesGraphicsItem.h>
#include <CGAL/Qt/GeneralPolygonGraphicsItem.h>
#include <CGAL/Qt/GeneralPolygonWithHolesGraphicsItem.h>
#include <CGAL/Qt/SlopedPolygonGraphicsItem.h> 
#include <CGAL/Qt/LineRangeGraphicsItem.h>
#include <CGAL/Qt/CircleRangeGraphicsItem.h>
#include <CGAL/Qt/SegmentRangeGraphicsItem.h>
#include <CGAL/Qt/StraightSkeletonGraphicsItem.h>


// #include <CGAL/Arithmetic_kernel.h>
// #define
// CGAL_USE_GMP
// CGAL_USE_CORE
// CGAL_USE_LEDA
//
// typedef CGAL::Arithmetic_kernel::Rational Rational;
// typedef CGAL::Arithmetic_kernel::Integer Integer;
// typedef Rational NT;

// polygon definitions
typedef CGAL::Gmpq                              NT;
//typedef CORE::BigRat                              NT;
typedef CGAL::Cartesian<NT>                     Kernel;
typedef Kernel::Line_2                          Line_2;
typedef Kernel::Circle_2                        Circle_2;

typedef CGAL::Polygon_2<Kernel>                       Polygon_2;
typedef CGAL::Offset_types_2<Polygon_2>  Types_2;

typedef Types_2::Polygon_with_holes_2      Polygon_with_holes_2;
typedef Types_2::Point_2                   Point_2;
typedef Types_2::Segment_2                 Segment_2;

typedef Types_2::Kgon_sum_polygon_2	Kgon_sum_polygon_2;
typedef Types_2::Approximate_offset_polygon_2	Approximate_offset_polygon_2;
typedef Types_2::Exact_offset_polygon_2	Exact_offset_polygon_2;

typedef Types_2::Approximate_polygon_2	Approximate_polygon_2;
typedef Types_2::Exact_polygon_2	Exact_polygon_2;

typedef Types_2::Conic_traits_2 Conic_traits_2;

typedef Types_2::Exact_polygon_list_2 Exact_polygon_list_2;
typedef Types_2::Exact_offset_polygon_list_2 Exact_offset_polygon_list_2;

typedef CGAL::Minkowski_sum_with_circle_2<Polygon_2>  Minsum_with_circle_2;

typedef Minsum_with_circle_2::Circle_app_2            Circle_app_2;
typedef Minsum_with_circle_2::Exact_segment_list_2 Exact_segment_list_2;
typedef Minsum_with_circle_2::Straight_skeleton_ptr_2 Straight_skeleton_ptr_2;
typedef Minsum_with_circle_2::Straight_skeleton_2 Straight_skeleton_2;

typedef CGAL::Qt::GraphicsViewPolylineInput< Kernel >                GVPolylineInput;

typedef CGAL::Qt::GraphicsItem                                       GI;
typedef CGAL::Qt::LineRangeGraphicsItem< Kernel >                    LinesGI;
typedef CGAL::Qt::CircleRangeGraphicsItem< Kernel >                  CirclesGI;
typedef CGAL::Qt::SegmentRangeGraphicsItem< Conic_traits_2 >       SegmentsGI;

typedef CGAL::Qt::PolygonGraphicsItem<Polygon_2>                      PolygonGI;
typedef CGAL::Qt::SlopedPolygonGraphicsItem<Polygon_2, Circle_app_2>  SlopedPolygonGI;
typedef CGAL::Qt::PolygonWithHolesGraphicsItem<Kgon_sum_polygon_2>    HoledPolygonGI;

typedef CGAL::Qt::GeneralPolygonGraphicsItem<Approximate_polygon_2>   AppPolygonGI;
typedef CGAL::Qt::GeneralPolygonGraphicsItem<Exact_polygon_2>         ExactPolygonGI;
typedef CGAL::Qt::StraightSkeletonGraphicsItem<Straight_skeleton_2>      StraightSkeletonGI;

typedef CGAL::Qt::GeneralPolygonWithHolesGraphicsItem<Approximate_offset_polygon_2>        ApproximateHoledPolygonGI;
typedef CGAL::Qt::GeneralPolygonWithHolesGraphicsItem<Exact_offset_polygon_2>        ExactHoledPolygonGI;


// straight skeleton
typedef CGAL::Straight_skeleton_2< Kernel > Ss_2;
typedef boost::shared_ptr<Ss_2> SsPtr_2;

#endif // POLYGON_OFFSET_2_TYPEDEFS_H
