#ifndef SIMULATE_H
#define SIMULATE_H

/*****************************************************************************/

#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/GraphicsViewNavigation.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>


/*****************************************************************************/

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};

/*****************************************************************************/


// 3D triangulation types.
typedef CGAL::Triangulation_3<K>                        Triangulation_3;
typedef CGAL::Delaunay_triangulation_3<K>               Delaunay_3;
typedef Triangulation_3::Point                          Point_3;
typedef Delaunay_3::Cell                                Cell_3;
typedef Cell_3::Cell_handle                             Cell_handle_3;
typedef Cell_3::Cell_handle                             Cell_handle;
typedef CGAL::Creator_uniform_2<double,Point_3>         Creator_3;

// 2D triangulation types.
typedef CGAL::Triangulation_2<K>                        Triangulation_2;
typedef CGAL::Delaunay_triangulation_2<K>               Delaunay_2;
typedef Triangulation_2::Point                          Point_2;
typedef Delaunay_2::Face                                Face_2;
typedef Delaunay_2::Line_face_circulator                Line_face_circulator;
typedef Face_2::Face_handle                             Face_handle_2;
typedef CGAL::Creator_uniform_2<double,Point_2>         Creator_2;
typedef CGAL::Qt::TriangulationGraphicsItem<Delaunay_2> QTriangulationGraphics;
typedef Delaunay_2::Locate_type                         Locate_type;


/*****************************************************************************/
#endif
