// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,

// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF describingIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Kabir Kedia<kabirkedia0111@gmail.com>
//             Ronnie Gandhi<ronniegandhi19999@gmail.com>
//             Apurva Bhatt <response2apurva@gmail.com>
//             Efi Fogel <efifogel@gmail.com>

//TODO:- add error handling

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <list>
#include <limits>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <boost/shared_ptr.hpp>

#include <QApplication>
#include <qmessagebox.h>
// #include <QDesktopWidget>

#include <QMainWindow>
#include <QGraphicsScene>
#include <QToolButton>
#include <QActionGroup>
#include <QSize>
#include <QPen>
#include <QtGui>
#include <QString>
#include <QFileDialog>
#include <QInputDialog>
#include <QWheelEvent>
#include <QPainter>
#include <QSlider>
#include <QGraphicsLineItem>
#include <QProgressBar>
#include <QMessageBox>
#include <QGraphicsPathItem>
#include <QGraphicsView>
#include <QGraphicsPixmapItem>
#include <QKeyEvent>
#include <QTimer>

#include <CGAL/auto_link/Qt.h>
#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>

#include <CGAL/basic.h>
#include <CGAL/Point_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Timer.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/iterator.h>
#include <CGAL/assertions_behaviour.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/General_polygon_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Small_side_angle_bisector_decomposition_2.h>
#include <CGAL/approximated_offset_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/squared_distance_2.h>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpq.h>
#else

#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#endif

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/Qt/utility.h>
#include <CGAL/IO/Dxf_bsop_reader.h>
#include <CGAL/Qt/GraphicsViewNavigation.h>

#include "QT5/MinkowskiSum.h"   //for future supports
#include "QT5/Circular_polygons.h"
#include "QT5/Linear_polygons.h"
#include "QT5/BezierCurves.h"
#include "QT5/Polyline_curves.h"
#include "QT5/Graphics_view_circular_polygon_input.h"
#include "QT5/Graphics_view_linear_polygon_input.h"
#include "QT5/Graphics_view_polyline_input.h"
#include "QT5/GraphicsViewBezierPolygonInput.h"
#include "QT5/Graphics_view_minkowski_input.h"
#include "QT5/General_polygon_2.h"
#include "QT5/General_polygon_set_2.h"
#include "QT5/General_polygon_set_on_surface_2.h"
#include "QT5/Gps_circle_segment_traits_2.h"
#include "QT5/Gps_segment_traits_2.h"
#include "QT5/Gps_traits_2.h"
#include "QT5/connect_holes.h"
#include "QT5/Polygon_set_2.h"

#include "ui_Boolean_set_operations_2.h"

#include "Typedefs.h"

// Forward Declarations
void show_warning(std::string aS);

void show_error(std::string aS);

void error(std::string aS);

void error_handler(char const *what, char const *expr, char const *file,
                   int line, char const *msg);


// Circular_polygon linearPart_2_circ(Circular_Linear_polygon const& pgn);
// Circular_polygon_with_holes linearPart_2_circ(Circular_Linear_polygon_with_holes const& pwh);

// Typedefs
typedef CGAL::Qt::Circular_set_graphics_item <Circular_polygon_set,
Circular_traits> Circular_GI;
typedef CGAL::Qt::Linear_set_graphics_item <Linear_polygon_set,
Linear_traits> Linear_GI;
typedef CGAL::Qt::Polyline_set_graphics_item <Polyline_polygon_set,
Polyline_traits> Polyline_GI;


//typedef CGAL::Qt::Bezier_set_graphics_item<Bezier_polygon_set,
//                                         Bezier_traits>      Bezier_GI;


// Functions to show errors

//show warnings but proceed
void show_warning(std::string aS) { QMessageBox::warning(NULL, "Warning", QString(aS.c_str())); }

//show errors in a message box and don't proceed
void show_error(std::string aS) { QMessageBox::critical(NULL, "Critical Error", QString(aS.c_str())); }

//show errors in a message box and then crash.
// Needs to be replaced by show error or a reset button
void error(std::string aS) {
  show_error(aS);
  throw std::runtime_error(aS);
}

void error_handler(char const *what, char const *expr, char const *file,
                   int line, char const *msg) {
  if (expr != "first != last" && expr != "_inc_to_right != cv._inc_to_right" && line != 1684 &&
      expr != "is_new == true") {
    std::ostringstream ss;

    ss << "CGAL error: " << what << " violation!" << std::endl
       << "Expr: " << expr << std::endl
       << "File: " << file << std::endl
       << "Line: " << line << std::endl;
    if (msg != 0 || msg != "" || msg != NULL) {
      ss << "Explanation:" << msg << std::endl;
    }
    if (expr == "is_simple_2(first, last, traits)") {
      ss << "\nCompleting Polygon from this point!!!" << std::endl;
    }
    if (msg == "The outer boundary self intersects at edges.") {
      ss << "\nDeleting Polygon causing Error!!!" << std::endl;
    }

    error(ss.str());
  }
}

enum {
  BLUE_GROUP, RED_GROUP, BLACK_GROUP, BROWN_GROUP, YELLOW_GROUP,
  MAGENTA_GROUP, AQUA_GROUP, RESULT_GROUP
};

//A way to maintain 5 category of polygons namely linear,circular,bezier, polyline(to be added later)
//enum generates errors so, we will use LINEAR_TYPE = 1 ,CIRCULAR_TYPE = 2 ,BEZIER_TYPE = 3 ,POLYLINE_TYPE = 4 and CONIC_TYPE = 5(Not this year)

//drawing tools for different colors
QPen sPens[] = {
  QPen(QColor(255, 0, 0), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),   //blue
  QPen(QColor(0, 0, 0), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),     //red
  QPen(QColor(0, 0, 255), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),   //black
  QPen(QColor(210, 105, 30), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),//brown
  QPen(QColor(255, 255, 0), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin), //yellow
  QPen(QColor(255, 0, 255), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin), //magenta
  QPen(QColor(0, 255, 255), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin), //aqua
  QPen(QColor(0, 0, 0, 0), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin),   //Result
};

QBrush sBrushes[] = {
  QBrush(QColor(255, 0, 0, 75)),           //blue
  QBrush(QColor(0, 0, 0, 75)),             //red
  QBrush(QColor(0, 0, 255, 75)),           //black
  QBrush(QColor(210, 105, 30, 75)),        //brown
  QBrush(QColor(255, 255, 0, 75)),         //yellow
  QBrush(QColor(255, 0, 255, 75)),         //magenta
  QBrush(QColor(0, 255, 255, 75)),         //aqua
  QBrush(QColor(0, 0, 0, 0)),              //Result
};

//A base call for rep class
struct Rep_base {
  virtual ~Rep_base() {}

  virtual int type() const = 0;

  virtual CGAL::Qt::GraphicsItem *gi() const = 0;

  virtual CGAL::Qt::GraphicsItem *gi() = 0;

  virtual void set_pen(QPen const &aPen) = 0;

  virtual void set_brush(QBrush const &aBrush) = 0;

  virtual QRectF bounding_rect() const { return gi()->boundingRect(); }

  virtual bool is_empty() const = 0;

  virtual void clear() = 0;

  virtual void complement() = 0;

  virtual void assign(Rep_base const &aOther) = 0;

  virtual void intersect(Rep_base const &aOther) = 0;

  virtual void join(Rep_base const &aOther) = 0;

  virtual void difference(Rep_base const &aOther) = 0;

  virtual void symmetric_difference(Rep_base const &aOther) = 0;
};

//Class for initializing Rep class.
//Used for linear and circular
template<typename GI_, typename Set_, typename Gps_traits>
class Rep : public Rep_base {
public:
  typedef GI_ GI;
  typedef Set_ Set;

  typedef Rep<GI, Set, Gps_traits> Self;

  Rep() { m_GI = new GI(&m_set, m_traits); }

  Set const &set() const { return m_set; }

  Set &set() { return m_set; }

  virtual CGAL::Qt::GraphicsItem *gi() const { return m_GI; }

  virtual CGAL::Qt::GraphicsItem *gi() { return m_GI; }

  virtual void set_pen(QPen const &aPen) { m_GI->setPen(aPen); }

  virtual void set_brush(QBrush const &aBrush) { m_GI->setBrush(aBrush); }

  virtual bool is_empty() const { return m_set.is_empty(); }

  //virtual functions for error handling of every operation
  virtual void clear() {
    try {
      m_set.clear();
    }
    catch (...) {
      show_error("Exception thrown during boolean operation clear");
    }
  }

  virtual void complement() {
    try {
      m_set.complement();
    }
    catch (...) {
      show_error("Exception thrown during boolean operation complement");
    }
  }

  virtual void assign(Rep_base const &aOther) {
    try {
      m_set = cast(aOther).m_set;
    }
    catch (...) {
      show_error("Exception thrown during boolean operation assign");
    }
  }

  virtual void intersect(Rep_base const &aOther) {
    try {
      m_set.intersection(cast(aOther).m_set);
    }
    catch (...) {
      show_error("Exception thrown during boolean operation intersect");
    }
  }

  virtual void join(Rep_base const &aOther) {
    try {
      m_set.join(cast(aOther).m_set);
    }
    catch (...) {
      show_error("Exception thrown during boolean operation union");
    }
  }

  virtual void difference(Rep_base const &aOther) {
    try {
      m_set.difference(cast(aOther).m_set);
    }
    catch (...) {
      show_error("Exception thrown during boolean operation difference");
    }
  }

  virtual void symmetric_difference(Rep_base const &aOther) {
    try {
      m_set.symmetric_difference(cast(aOther).m_set);
    }
    catch (...) {
      show_error("Exception thrown during boolean operation symmetric difference");
    }
  }

  static Self const &cast(Rep_base const &aOther) { return dynamic_cast<Self const &>(aOther); }

  static Self &cast(Rep_base &aOther) { return dynamic_cast<Self &>(aOther); }

private:
  //For maintaining all drawing operations
  GI *m_GI;
  //Storage for all polygons of one type. It is used as a base to perform all boolean operations
  Set m_set;
protected:
  //pass it
  Gps_traits m_traits;
};

//Implementing Bezier's rep class
template<class GI_, class Set_>
class Rep_o : public Rep_base {
public:
  typedef GI_ GI;
  typedef Set_ Set;

  typedef Rep_o<GI, Set> Self;

  Rep_o() { m_GI = new GI(&m_set); }

  Set const &set() const { return m_set; }

  Set &set() { return m_set; }

  virtual CGAL::Qt::GraphicsItem *gi() const { return m_GI; }

  virtual CGAL::Qt::GraphicsItem *gi() { return m_GI; }

  virtual void set_pen(QPen const &aPen) { m_GI->setPen(aPen); }

  virtual void set_brush(QBrush const &aBrush) { m_GI->setBrush(aBrush); }

  virtual bool is_empty() const { return m_set.is_empty(); }

  virtual void clear() {
    try {
      m_set.clear();
    }
    catch (...) {
      show_error("Exception thrown during boolean operation clear");
    }
  }

  virtual void complement() {
    try {
      m_set.complement();
    }
    catch (...) {
      show_error("Exception thrown during boolean operation complement");
    }
  }

  virtual void assign(Rep_base const &aOther) {
    try {
      m_set = cast(aOther).m_set;
    }
    catch (...) {
      show_error("Exception thrown during boolean operation assign");
    }
  }

  virtual void intersect(Rep_base const &aOther) {
    try {
      m_set.intersection(cast(aOther).m_set);
    }
    catch (...) {
      show_error("Exception thrown during boolean operation intersect");
    }
  }

  virtual void join(Rep_base const &aOther) {
    try {
      m_set.join(cast(aOther).m_set);
    }
    catch (...) {
      show_error("Exception thrown during boolean operation union");
    }
  }

  virtual void difference(Rep_base const &aOther) {
    try {
      m_set.difference(cast(aOther).m_set);
    }
    catch (...) {
      show_error("Exception thrown during boolean operation difference");
    }
  }

  virtual void symmetric_difference(Rep_base const &aOther) {
    try {
      m_set.symmetric_difference(cast(aOther).m_set);
    }
    catch (...) {
      show_error("Exception thrown during boolean operation symmetric difference");
    }
  }

  static Self const &cast(Rep_base const &aOther) { return dynamic_cast<Self const &>(aOther); }

  static Self &cast(Rep_base &aOther) { return dynamic_cast<Self &>(aOther); }

private:
  GI *m_GI;
  Set m_set;
};

//A class for connecting linear polygon GUI and this file
class Linear_rep : public Rep<Linear_GI, Linear_polygon_set, Linear_traits> {
  typedef Rep<Linear_GI, Linear_polygon_set, Linear_traits> Base;
public:
  Linear_rep() : Base() {}

  virtual int type() const { return 1; }
};

//A class for connecting circular polygon GUI and this file
class Circular_rep : public Rep<Circular_GI, Circular_polygon_set, Circular_traits> {
  typedef Rep<Circular_GI, Circular_polygon_set, Circular_traits> Base;
public:
  Circular_rep() : Base() {}

  virtual int type() const { return 2; }
};

//A class for connecting polyline polygon GUI and this file
class Polyline_rep : public Rep<Polyline_GI, Polyline_polygon_set, Polyline_traits> {
  typedef Rep<Polyline_GI, Polyline_polygon_set, Polyline_traits> Base;
public:
  Polyline_rep() : Base() {}

  virtual int type() const { return 4; }
};

//A class for connecting bezier polygon GUI and this file
class Bezier_rep : public Rep_o<Bezier_GI, Bezier_polygon_set> {
  typedef Rep_o<Bezier_GI, Bezier_polygon_set> Base;
public:
  Bezier_rep() : Base() {}

  virtual int type() const { return 3; }
};

// a container which deletes an object when last shared_ptr gets deleted or re-initiated
class Curve_set {
  typedef boost::shared_ptr <Rep_base> Rep_ptr;

public:
  //constructor
  Curve_set(int aType, QPen aPen, QBrush aBrush) : m_pen(aPen), m_brush(aBrush) { reset_type(aType); }

  void reset_type(int aType) {
    //setting shared_ptr for respective polygon
    if (aType == 1) m_rep = Rep_ptr(new Linear_rep());
    else if (aType == 2) m_rep = Rep_ptr(new Circular_rep());
    else if (aType == 3) m_rep = Rep_ptr(new Bezier_rep());
    else if (aType == 4) m_rep = Rep_ptr(new Polyline_rep());

    //setting pen and brush
    m_rep->set_pen(m_pen);
    m_rep->set_brush(m_brush);
  }

  CGAL::Qt::GraphicsItem const *gi() const { return m_rep->gi(); }

  CGAL::Qt::GraphicsItem *gi() { return m_rep->gi(); }

  QRectF bounding_rect() const { return m_rep->bounding_rect(); }

  bool is_empty() const { return !m_rep || m_rep->is_empty(); }

  void clear() { m_rep->clear(); }

  //boolean operations
  void complement() { m_rep->complement(); }

  void assign(Curve_set const &aOther) {
    if (is_linear() && aOther.is_linear())
      get_linear_rep()->assign(*aOther.get_linear_rep());
    else if (is_circular() && aOther.is_circular())
      get_circular_rep()->assign(*aOther.get_circular_rep());
    else if (is_bezier() && aOther.is_bezier())
      get_bezier_rep()->assign(*aOther.get_bezier_rep());
    else if (is_polyline() && aOther.is_polyline())
      get_polyline_rep()->assign(*aOther.get_polyline_rep());
  }

  void intersect(Curve_set const &aOther) {
    if (is_linear() && aOther.is_linear())
      get_linear_rep()->intersect(*aOther.get_linear_rep());
    else if (is_circular() && aOther.is_circular())
      get_circular_rep()->intersect(*aOther.get_circular_rep());
    else if (is_bezier() && aOther.is_bezier())
      get_bezier_rep()->intersect(*aOther.get_bezier_rep());
    else if (is_polyline() && aOther.is_polyline())
      get_polyline_rep()->intersect(*aOther.get_polyline_rep());
  }

  void join(Curve_set const &aOther) {
    if (is_linear() && aOther.is_linear())
      get_linear_rep()->join(*aOther.get_linear_rep());
    else if (is_circular() && aOther.is_circular())
      get_circular_rep()->join(*aOther.get_circular_rep());
    else if (is_bezier() && aOther.is_bezier())
      get_bezier_rep()->join(*aOther.get_bezier_rep());
    else if (is_polyline() && aOther.is_polyline())
      get_polyline_rep()->join(*aOther.get_polyline_rep());
  }

  void difference(Curve_set const &aOther) {
    if (is_linear() && aOther.is_linear())
      get_linear_rep()->difference(*aOther.get_linear_rep());
    else if (is_circular() && aOther.is_circular())
      get_circular_rep()->difference(*aOther.get_circular_rep());
    else if (is_bezier() && aOther.is_bezier())
      get_bezier_rep()->difference(*aOther.get_bezier_rep());
    else if (is_polyline() && aOther.is_polyline())
      get_polyline_rep()->difference(*aOther.get_polyline_rep());
  }

  void symmetric_difference(Curve_set const &aOther) {
    if (is_linear() && aOther.is_linear())
      get_linear_rep()->symmetric_difference(*aOther.get_linear_rep());
    else if (is_circular() && aOther.is_circular())
      get_circular_rep()->symmetric_difference(*aOther.get_circular_rep());
    else if (is_bezier() && aOther.is_bezier())
      get_bezier_rep()->symmetric_difference(*aOther.get_bezier_rep());
    else if (is_polyline() && aOther.is_polyline())
      get_polyline_rep()->symmetric_difference(*aOther.get_polyline_rep());
  }

  void change_brush_color(QBrush aBrush) {
    m_brush = aBrush;
  }

  void change_pen_color(QPen aPen) {
    m_pen = aPen;
  }

  //Minkowski Sum handled separately
  //see its need keep it for now
  const Rep_base &rep() const { return *m_rep; }

  Rep_base &rep() { return *m_rep; }

  bool is_linear() const { return m_rep->type() == 1; }

  bool is_circular() const { return m_rep->type() == 2; }

  bool is_bezier() const { return m_rep->type() == 3; }

  bool is_polyline() const { return m_rep->type() == 4; }

  //to get rep for linear polygons
  const Linear_rep *get_linear_rep() const { return dynamic_cast<Linear_rep const *>(boost::get_pointer(m_rep)); }

  Linear_rep *get_linear_rep() { return dynamic_cast<Linear_rep *>(boost::get_pointer(m_rep)); }

  //to get Linear_polygon_set
  const Linear_polygon_set &linear() const { return get_linear_rep()->set(); }

  Linear_polygon_set &linear() { return get_linear_rep()->set(); }

  //to get rep for circular polygons
  const Circular_rep *get_circular_rep() const { return dynamic_cast<Circular_rep const *>(boost::get_pointer(m_rep)); }

  Circular_rep *get_circular_rep() { return dynamic_cast<Circular_rep *>(boost::get_pointer(m_rep)); }

  //to get Circular_polygon_set
  const Circular_polygon_set &circular() const { return get_circular_rep()->set(); }

  Circular_polygon_set &circular() { return get_circular_rep()->set(); }

  //get rep for bezier polygons
  Bezier_rep const *get_bezier_rep() const { return dynamic_cast<Bezier_rep const *>( boost::get_pointer(m_rep)); }

  Bezier_rep *get_bezier_rep() { return dynamic_cast<Bezier_rep *        >( boost::get_pointer(m_rep)); }

  //to get Bezier_polygon_set
  const Bezier_polygon_set &bezier() const { return get_bezier_rep()->set(); }

  Bezier_polygon_set &bezier() { return get_bezier_rep()->set(); }

  //to get rep for polylines
  const Polyline_rep *get_polyline_rep() const {
    return dynamic_cast<Polyline_rep const *>( boost::get_pointer(m_rep));
  }

  Polyline_rep *get_polyline_rep() { return dynamic_cast<Polyline_rep *>( boost::get_pointer(m_rep)); }

  //to get Polyline_polygon_set
  const Polyline_polygon_set &polyline() const { return get_polyline_rep()->set(); }

  Polyline_polygon_set &polyline() { return get_polyline_rep()->set(); }

public:
  //drawing tools
  QPen m_pen;
  QBrush m_brush;
  // a container which deletes an object when last shared_ptr gets deleted or
  // re-initiated
  boost::shared_ptr <Rep_base> m_rep;
};

typedef std::vector <Curve_set> Curve_set_container;
typedef Curve_set_container::const_iterator Curve_set_const_iterator;
typedef Curve_set_container::iterator Curve_set_iterator;

// This container is for returning different color polygons while performing an operation
class State_current {
  /* COMPLEMENT_OP = 0
  // INTERSECTION_OP = 1
  // UNION_OP = 2
  // DIFFERENCE_OP = 3
  // SYMMETRIC_DIFFERENCE_OP = 4
  // MINKOWSKI_SUM_OP = 5
  // COPY_OP = 6
  // MOVE_OP = 7
  // CLEAR_OP = 9
  // DELETE_ALL_OP = 10
  // START_OP = 11*/

public:
  State_current(size_t operation_name, size_t aType);

public:
  Curve_set &set(size_t aGroup) { return m_curve_sets[aGroup]; }
  /*gets which group is currently active now
    size_t active_group() const { return m_color_active; }
    size_t complement_group() const {return m_color_complement; } //see if needed
  */

  //sets the current active group
  Curve_set &active_set(size_t m_color_active) { return set(m_color_active); }

  Curve_set &move_set(size_t m_color_move) { return set(m_color_move); }

  //setting curve
  Curve_set &blue_set() { return set(BLUE_GROUP); }

  Curve_set &red_set() { return set(RED_GROUP); }

  Curve_set &black_set() { return set(BLACK_GROUP); }

  Curve_set &brown_set() { return set(BROWN_GROUP); }

  Curve_set &yellow_set() { return set(YELLOW_GROUP); }

  Curve_set &magenta_set() { return set(MAGENTA_GROUP); }

  Curve_set &aqua_set() { return set(AQUA_GROUP); }

  Curve_set &result_set() { return set(RESULT_GROUP); }

  void ToggleView(size_t aGROUP, bool a_check) {
    if (a_check) set(aGROUP).gi()->show();
    else set(aGROUP).gi()->hide();
  }

  //returns linear containers
  Linear_region_source_container &blue_linear_sources() { return m_blue_linear_sources; }

  Linear_region_source_container &red_linear_sources() { return m_red_linear_sources; }

  Linear_region_source_container &black_linear_sources() { return m_black_linear_sources; }

  Linear_region_source_container &brown_linear_sources() { return m_brown_linear_sources; }

  Linear_region_source_container &yellow_linear_sources() { return m_yellow_linear_sources; }

  Linear_region_source_container &magenta_linear_sources() { return m_magenta_linear_sources; }

  Linear_region_source_container &aqua_linear_sources() { return m_aqua_linear_sources; }

  Linear_region_source_container &result_linear_sources() { return m_result_linear_sources; }

  //returns circular containers
  Circular_region_source_container &blue_circular_sources() { return m_blue_circular_sources; }

  Circular_region_source_container &red_circular_sources() { return m_red_circular_sources; }

  Circular_region_source_container &black_circular_sources() { return m_black_circular_sources; }

  Circular_region_source_container &brown_circular_sources() { return m_brown_circular_sources; }

  Circular_region_source_container &yellow_circular_sources() { return m_yellow_circular_sources; }

  Circular_region_source_container &magenta_circular_sources() { return m_magenta_circular_sources; }

  Circular_region_source_container &aqua_circular_sources() { return m_aqua_circular_sources; }

  Circular_region_source_container &result_circular_sources() { return m_result_circular_sources; }

  //returns bezier containers
  Bezier_region_source_container &blue_bezier_sources() { return m_blue_bezier_sources; }

  Bezier_region_source_container &red_bezier_sources() { return m_red_bezier_sources; }

  Bezier_region_source_container &black_bezier_sources() { return m_black_bezier_sources; }

  Bezier_region_source_container &brown_bezier_sources() { return m_brown_bezier_sources; }

  Bezier_region_source_container &yellow_bezier_sources() { return m_yellow_bezier_sources; }

  Bezier_region_source_container &magenta_bezier_sources() { return m_magenta_bezier_sources; }

  Bezier_region_source_container &aqua_bezier_sources() { return m_aqua_bezier_sources; }

  Bezier_region_source_container &result_bezier_sources() { return m_result_bezier_sources; }

  //returns polyline containers
  Polyline_region_source_container &blue_polyline_sources() { return m_blue_polyline_sources; }

  Polyline_region_source_container &red_polyline_sources() { return m_red_polyline_sources; }

  Polyline_region_source_container &black_polyline_sources() { return m_black_polyline_sources; }

  Polyline_region_source_container &brown_polyline_sources() { return m_brown_polyline_sources; }

  Polyline_region_source_container &yellow_polyline_sources() { return m_yellow_polyline_sources; }

  Polyline_region_source_container &magenta_polyline_sources() { return m_magenta_polyline_sources; }

  Polyline_region_source_container &aqua_polyline_sources() { return m_aqua_polyline_sources; }

  Polyline_region_source_container &result_polyline_sources() { return m_result_polyline_sources; }

  //choosing colors in case of linear polygons
  Linear_region_source_container &active_linear_sources(size_t m_color_active) {
    switch (m_color_active) {
     case 0:return m_blue_linear_sources;
     case 1:return m_red_linear_sources;
     case 2:return m_black_linear_sources;
     case 3:return m_brown_linear_sources;
     case 4:return m_yellow_linear_sources;
     case 5:return m_magenta_linear_sources;
     case 6:return m_aqua_linear_sources;
     case 7:return m_result_linear_sources;

     default:break;
    }
    CGAL_warning_msg(true, "Should not reach here!");
    return m_blue_linear_sources;
  }

  //choosing colors in case of circular polygons
  Circular_region_source_container &active_circular_sources(size_t m_color_active) {
    switch (m_color_active) {
     case 0:return m_blue_circular_sources;
     case 1:return m_red_circular_sources;
     case 2:return m_black_circular_sources;
     case 3:return m_brown_circular_sources;
     case 4:return m_yellow_circular_sources;
     case 5:return m_magenta_circular_sources;
     case 6:return m_aqua_circular_sources;
     case 7:return m_result_circular_sources;

     default:break;
    }
    CGAL_warning_msg(true, "Should not reach here!");
    return m_blue_circular_sources;
  }

  //choosing colors in case of bezier polygons
  Bezier_region_source_container &active_bezier_sources(size_t m_color_active) {
    switch (m_color_active) {
     case 0:return m_blue_bezier_sources;
     case 1:return m_red_bezier_sources;
     case 2:return m_black_bezier_sources;
     case 3:return m_brown_bezier_sources;
     case 4:return m_yellow_bezier_sources;
     case 5:return m_magenta_bezier_sources;
     case 6:return m_aqua_bezier_sources;
     case 7:return m_result_bezier_sources;

     default:break;
    }
    CGAL_warning_msg(true, "Should not reach here!");
    return m_blue_bezier_sources;
  }

  //choosing colors in case of polyline polygons
  Polyline_region_source_container &active_polyline_sources(size_t m_color_active) {
    switch (m_color_active) {
     case 0:return m_blue_polyline_sources;
     case 1:return m_red_polyline_sources;
     case 2:return m_black_polyline_sources;
     case 3:return m_brown_polyline_sources;
     case 4:return m_yellow_polyline_sources;
     case 5:return m_magenta_polyline_sources;
     case 6:return m_aqua_polyline_sources;
     case 7:return m_result_polyline_sources;

     default:break;
    }
    CGAL_warning_msg(true, "Should not reach here!");
    return m_blue_polyline_sources;
  }

public:
  Curve_set_container m_curve_sets;
  //container for curves
  Linear_region_source_container m_blue_linear_sources;
  Linear_region_source_container m_red_linear_sources;
  Linear_region_source_container m_black_linear_sources;
  Linear_region_source_container m_brown_linear_sources;
  Linear_region_source_container m_yellow_linear_sources;
  Linear_region_source_container m_magenta_linear_sources;
  Linear_region_source_container m_aqua_linear_sources;
  Linear_region_source_container m_result_linear_sources;

  Circular_region_source_container m_blue_circular_sources;
  Circular_region_source_container m_red_circular_sources;
  Circular_region_source_container m_black_circular_sources;
  Circular_region_source_container m_brown_circular_sources;
  Circular_region_source_container m_yellow_circular_sources;
  Circular_region_source_container m_magenta_circular_sources;
  Circular_region_source_container m_aqua_circular_sources;
  Circular_region_source_container m_result_circular_sources;

  Bezier_region_source_container m_blue_bezier_sources;
  Bezier_region_source_container m_red_bezier_sources;
  Bezier_region_source_container m_black_bezier_sources;
  Bezier_region_source_container m_brown_bezier_sources;
  Bezier_region_source_container m_yellow_bezier_sources;
  Bezier_region_source_container m_magenta_bezier_sources;
  Bezier_region_source_container m_aqua_bezier_sources;
  Bezier_region_source_container m_result_bezier_sources;

  Polyline_region_source_container m_blue_polyline_sources;
  Polyline_region_source_container m_red_polyline_sources;
  Polyline_region_source_container m_black_polyline_sources;
  Polyline_region_source_container m_brown_polyline_sources;
  Polyline_region_source_container m_yellow_polyline_sources;
  Polyline_region_source_container m_magenta_polyline_sources;
  Polyline_region_source_container m_aqua_polyline_sources;
  Polyline_region_source_container m_result_polyline_sources;

  //size_t m_state_num;
  size_t m_operation;
  size_t m_type;
};

//defining the State_current constructor
State_current::State_current(size_t operation_name, size_t aType) : m_operation(operation_name), m_type(aType) {
  //for setting Curve_set of aGroup type an int representing a set of polygon of a specific type

  //default setups with Linear
  m_curve_sets.push_back(Curve_set(m_type, sPens[BLUE_GROUP], sBrushes[BLUE_GROUP]));
  m_curve_sets.push_back(Curve_set(m_type, sPens[RED_GROUP], sBrushes[RED_GROUP]));
  m_curve_sets.push_back(Curve_set(m_type, sPens[BLACK_GROUP], sBrushes[BLACK_GROUP]));
  m_curve_sets.push_back(Curve_set(m_type, sPens[BROWN_GROUP], sBrushes[BROWN_GROUP]));
  m_curve_sets.push_back(Curve_set(m_type, sPens[YELLOW_GROUP], sBrushes[YELLOW_GROUP]));
  m_curve_sets.push_back(Curve_set(m_type, sPens[MAGENTA_GROUP], sBrushes[MAGENTA_GROUP]));
  m_curve_sets.push_back(Curve_set(m_type, sPens[AQUA_GROUP], sBrushes[AQUA_GROUP]));
  m_curve_sets.push_back(Curve_set(m_type, sPens[RESULT_GROUP], sBrushes[RESULT_GROUP]));
}

typedef std::vector <State_current> state_container;
typedef state_container::const_iterator state_const_iterator;
typedef state_container::iterator state_iterator;

//The main page for displaying the window and all the icons
class MainWindow : public CGAL::Qt::DemosMainWindow,
                   public Ui::Boolean_set_operations_2 {
  Q_OBJECT// removing it gives error ui not declared

private:
  QGraphicsScene m_scene;
  //keep it intact for now check it out
  bool m_circular_active;
  bool m_bezier_active;
  bool m_polyline_active;

  //to determine which type is currently active
  size_t m_color_active;
  size_t m_color_move;
  bool empty_warn;
  size_t m_color_visible;
  size_t m_color_cm;
  size_t m_color_complement;
  size_t m_color_copy;
  size_t m_color_move_ext;
  bool m_blue_int;
  bool m_red_int;
  bool m_black_int;
  bool m_brown_int;
  bool m_yellow_int;
  bool m_magenta_int;
  bool m_aqua_int;
  bool m_grid;
  bool m_pan;
  bool m_reset;

  bool m_blue_union;
  bool m_red_union;
  bool m_black_union;
  bool m_brown_union;
  bool m_yellow_union;
  bool m_magenta_union;
  bool m_aqua_union;

  bool m_blue_sym_diff;
  bool m_red_sym_diff;
  bool m_black_sym_diff;
  bool m_brown_sym_diff;
  bool m_yellow_sym_diff;
  bool m_magenta_sym_diff;
  bool m_aqua_sym_diff;

  bool m_blue_mink;
  bool m_red_mink;
  bool m_black_mink;
  bool m_brown_mink;
  bool m_yellow_mink;
  bool m_magenta_mink;
  bool m_aqua_mink;

  bool m_visible_black;
  bool m_visible_brown;
  bool m_visible_yellow;
  bool m_visible_magenta;
  bool m_visible_aqua;

  bool m_clear_blue;
  bool m_clear_red;
  bool m_clear_black;
  bool m_clear_brown;
  bool m_clear_yellow;
  bool m_clear_magenta;
  bool m_clear_aqua;

  bool minkowski_sum_operated;
  bool m_disjoint;

  size_t m_state_num;

  QGraphicsLineItem *xAxis = new QGraphicsLineItem();
  QGraphicsLineItem *yAxis = new QGraphicsLineItem();

  //typedefs of classes used to draw circular and linear polygon

  state_container states_stack;

  CGAL::Qt::Graphics_view_linear_polygon_input <Kernel> *m_linear_input;
  CGAL::Qt::Graphics_view_circular_polygon_input <Kernel> *m_circular_input;
  CGAL::Qt::GraphicsViewBezierPolygonInput <Bezier_traits> *m_bezier_input;
  CGAL::Qt::Graphics_view_polyline_input <Kernel> *m_polyline_input;

public:
  MainWindow();

private:
  void zoomToFit();

protected
slots:
  //for file handling
  void open(QString
            filename);

  //defining QT slots for GUI integration
public
slots:
  //for all the operations such as union, complement, difference

  void processInput(CGAL::Object o);

  void on_actionUndo_triggered();

  void on_actionNew_triggered();

  void on_actionAxis_triggered();

  void on_actionRecenter_triggered();

  void on_actionComplementH_toggled(bool aChecked);

  void on_actionUnionH_toggled(bool aChecked);

  void on_actionIntersectionH_toggled(bool aChecked);

  void on_actionDifferenceH_toggled(bool aChecked);

  void on_actionSymmetric_DifferenceH_toggled(bool aChecked);

  void on_actionMinkowski_SumH_toggled(bool aChecked);

  void on_actionInsertLinear_toggled(bool aChecked);

  void on_actionInsertCircular_toggled(bool aChecked);

  void on_actionInsertBezier_toggled(bool aChecked);

  void on_actionInsertPolyline_toggled(bool aChecked);

  void on_showColorBucket_toggled(bool aChecked);

  void on_sceneDockWidget_visibilityChanged(bool visible);

  void on_infoDockWidget_visibilityChanged(bool visible);

  void on_consoleDockWidget_visibilityChanged(bool visible);

  void on_showConsole_toggled(bool aChecked);

  void on_showInfo_toggled(bool aChecked);

  void on_actionOpenLinear_triggered();

  void on_actionOpenDXF_triggered();

  void on_actionOpenBezier_triggered();

  void on_actionOpenPolyline_triggered();

  void wheelEvent(QWheelEvent *event);

  void on_actionSaveCurrentBucket_triggered();

  void on_actionAddColor_triggered();

  void on_actionMinusColor_triggered();

  void on_actionCopyH_toggled(bool aChecked);

  void on_actionMoveH_toggled(bool aChecked);

  void show_not_empty_warning();

  //showing colors when clicked on the checkbox
  void on_showBlue_toggled(bool a_check);

  void on_showRed_toggled(bool a_check);

  void on_showBlack_toggled(bool a_check);

  void on_showBrown_toggled(bool a_check);

  void on_showYellow_toggled(bool a_check);

  void on_showMagenta_toggled(bool a_check);

  void on_showAqua_toggled(bool a_check);

  void on_showResult_toggled(bool a_check);

  void on_VisibleHeader_toggled(bool a_check);

  void on_showClear_toggled(bool a_check);

  void on_copyBlue_toggled(bool aCheck);

  void on_copyRed_toggled(bool aCheck);

  void on_copyBlack_toggled(bool aCheck);

  void on_copyBrown_toggled(bool aCheck);

  void on_copyYellow_toggled(bool aCheck);

  void on_copyMagenta_toggled(bool aCheck);

  void on_copyAqua_toggled(bool aCheck);

  void on_moveBlue_toggled(bool aCheck);

  void on_moveRed_toggled(bool aCheck);

  void on_moveBlack_toggled(bool aCheck);

  void on_moveBrown_toggled(bool aCheck);

  void on_moveYellow_toggled(bool aCheck);

  void on_moveMagenta_toggled(bool aCheck);

  void on_moveAqua_toggled(bool aCheck);

  void on_showBlueComp_toggled(bool aCheck);

  void on_showRedComp_toggled(bool aCheck);

  void on_showBlackComp_toggled(bool aCheck);

  void on_showBrownComp_toggled(bool aCheck);

  void on_showYellowComp_toggled(bool aCheck);

  void on_showMagentaComp_toggled(bool aCheck);

  void on_showAquaComp_toggled(bool aCheck);

  void on_showBlueInt_toggled(bool aCheck);

  void on_showRedInt_toggled(bool aCheck);

  void on_showBlackInt_toggled(bool aCheck);

  void on_showBrownInt_toggled(bool aCheck);

  void on_showYellowInt_toggled(bool aCheck);

  void on_showMagentaInt_toggled(bool aCheck);

  void on_showAquaInt_toggled(bool aCheck);

  //union
  void on_showBlueUnion_toggled(bool aCheck);

  void on_showRedUnion_toggled(bool aCheck);

  void on_showBlackUnion_toggled(bool aCheck);

  void on_showBrownUnion_toggled(bool aCheck);

  void on_showYellowUnion_toggled(bool aCheck);

  void on_showMagentaUnion_toggled(bool aCheck);

  void on_showAquaUnion_toggled(bool aCheck);

  void on_showBlueDiff_toggled(bool aCheck);

  void on_showRedDiff_toggled(bool aCheck);

  void on_showBlackDiff_toggled(bool aCheck);

  void on_showBrownDiff_toggled(bool aCheck);

  void on_showYellowDiff_toggled(bool aCheck);

  void on_showMagentaDiff_toggled(bool aCheck);

  void on_showAquaDiff_toggled(bool aCheck);

  //difference
  void on_showBlueSym_Diff_toggled(bool aCheck);

  void on_showRedSym_Diff_toggled(bool aCheck);

  void on_showBlackSym_Diff_toggled(bool aCheck);

  void on_showBrownSym_Diff_toggled(bool aCheck);

  void on_showYellowSym_Diff_toggled(bool aCheck);

  void on_showMagentaSym_Diff_toggled(bool aCheck);

  void on_showAquaSym_Diff_toggled(bool aCheck);

  //minkowski sum
  void on_showBlueMink_Sum_toggled(bool aCheck);

  void on_showRedMink_Sum_toggled(bool aCheck);

  void on_showBlackMink_Sum_toggled(bool aCheck);

  void on_showBrownMink_Sum_toggled(bool aCheck);

  void on_showYellowMink_Sum_toggled(bool aCheck);

  void on_showMagentaMink_Sum_toggled(bool aCheck);

  void on_showAquaMink_Sum_toggled(bool aCheck);

  void on_drawBlue_toggled(bool a_check);

  void on_drawRed_toggled(bool a_check);

  void on_drawBlack_toggled(bool a_check);

  void on_drawBrown_toggled(bool a_check);

  void on_drawYellow_toggled(bool a_check);

  void on_drawMagenta_toggled(bool a_check);

  void on_drawAqua_toggled(bool a_check);

  //about demo,cgal
  void on_aboutBolop_triggered();

  void on_aboutCGAL_triggered();

  void on_aboutMinkop_triggered();

  void on_aboutDemo_triggered();

  void resizeEvent(QResizeEvent *event);

  void exception_handler();

  void on_actionDeleteResult();

  void on_actionDeleteAll_triggered();

  void on_actionClearH_toggled(bool aChecked);

  void on_actionPAN_triggered();

  Polygon_with_holes_2 getMinkInputPolygon(size_t color);

  //for reading polygons from a file
  bool read_linear(QString aFileName, Linear_polygon_set &rSet,
                   Linear_region_source_container &rSources);

  bool read_circular(QString aFileName, Circular_polygon_set &rSet,
                     Circular_region_source_container &rSources);

  bool read_polyline(QString aFilename, Polyline_polygon_set &rSet,
                     Polyline_region_source_container &rSources);

  bool read_bezier(QString aFileName);

  Bezier_curve read_bezier_curve(std::istream &is, bool aDoubleFormat);

signals:
  void changed();

private:
  void get_new_state(size_t operation_name) {
    //operation_name
    // COMPLEMENT_OP = 0
    // INTERSECTION_OP = 1
    // UNION_OP = 2
    // DIFFERENCE_OP = 3
    // SYMMETRIC_DIFFERENCE_OP = 4
    // MINKOWSKI_SUM_OP = 5
    // RESET_OP = 6
    // COPY_OP = 7
    // MOVE_OP = 8
    if (operation_name == 12) {
      if (!m_circular_active && !m_bezier_active && !m_polyline_active) {
        State_current m_state = State_current(operation_name, 1);
        states_stack.push_back(m_state);
        states_stack.back();

        for (auto si = states_stack.back().m_curve_sets.begin();
             si != states_stack.back().m_curve_sets.end(); ++si)
        { link_GI(si->gi()); }
      } else if (!m_bezier_active && !m_polyline_active) {
        State_current m_state = State_current(operation_name, 2);
        states_stack.push_back(m_state);
        states_stack.back();

        for (auto si = states_stack.back().m_curve_sets.begin();
             si != states_stack.back().m_curve_sets.end(); ++si)
        { link_GI(si->gi()); }
      } else if (!m_polyline_active) {
        State_current m_state = State_current(operation_name, 3);
        states_stack.push_back(m_state);
        states_stack.back();

        for (auto si = states_stack.back().m_curve_sets.begin();
             si != states_stack.back().m_curve_sets.end(); ++si)
        { link_GI(si->gi()); }
      }
      else {
        State_current m_state = State_current(operation_name, 4);
        states_stack.push_back(m_state);
        states_stack.back();

        for (auto si = states_stack.back().m_curve_sets.begin();
             si != states_stack.back().m_curve_sets.end(); ++si)
        { link_GI(si->gi()); }
      }
    }
    else {
      if (!m_circular_active && !m_bezier_active && !m_polyline_active) {
        State_current m_state_x = State_current(operation_name, 1);

        for (int i = 0; i < 8; i++) {
          m_state_x.m_curve_sets[i].assign(states_stack.back().m_curve_sets[i]);
        }

        modelChanged();

        for (auto si = states_stack.back().m_curve_sets.begin();
             si != states_stack.back().m_curve_sets.end(); ++si) { unlink_GI(si->gi()); }

        modelChanged();

        states_stack.push_back(m_state_x);

        for (auto si = states_stack.back().m_curve_sets.begin();
             si != states_stack.back().m_curve_sets.end(); ++si) { link_GI(si->gi()); }

        modelChanged();

      } else if (!m_bezier_active && !m_polyline_active) {
        State_current m_state_x = State_current(operation_name, 2);

        for (int i = 0; i < 8; i++) {
          m_state_x.m_curve_sets[i].assign(states_stack.back().m_curve_sets[i]);
        }

        modelChanged();

        for (auto si = states_stack.back().m_curve_sets.begin();
             si != states_stack.back().m_curve_sets.end(); ++si) { unlink_GI(si->gi()); }

        modelChanged();

        states_stack.push_back(m_state_x);

        for (auto si = states_stack.back().m_curve_sets.begin();
             si != states_stack.back().m_curve_sets.end(); ++si) { link_GI(si->gi()); }

        modelChanged();
      } else if (!m_polyline_active) {
        State_current m_state_x = State_current(operation_name, 3);

        for (int i = 0; i < 8; i++) {
          m_state_x.m_curve_sets[i].assign(states_stack.back().m_curve_sets[i]);
        }

        modelChanged();

        for (auto si = states_stack.back().m_curve_sets.begin();
             si != states_stack.back().m_curve_sets.end(); ++si) { unlink_GI(si->gi()); }

        modelChanged();

        states_stack.push_back(m_state_x);

        for (auto si = states_stack.back().m_curve_sets.begin();
             si != states_stack.back().m_curve_sets.end(); ++si) { link_GI(si->gi()); }

        modelChanged();

      } else {
        State_current m_state_x = State_current(operation_name, 4);

        for (int i = 0; i < 8; i++) {
          m_state_x.m_curve_sets[i].assign(states_stack.back().m_curve_sets[i]);
        }

        modelChanged();

        for (auto si = states_stack.back().m_curve_sets.begin();
             si != states_stack.back().m_curve_sets.end(); ++si) { unlink_GI(si->gi()); }

        modelChanged();

        states_stack.push_back(m_state_x);

        for (auto si = states_stack.back().m_curve_sets.begin();
             si != states_stack.back().m_curve_sets.end(); ++si) { link_GI(si->gi()); }

        modelChanged();
      }
    }

    m_state_num = states_stack.size();
    actionUndo->setEnabled(true);

    // for confirming the num of each state during debugging
    //show_warning("State num: "+to_string(states_stack.size()));
    modelChanged();
  }

  void modelChanged() { emit(changed()); }

  //warning message for user
  bool ask_user_yesno(const char* aTitle, const char* aQuestion) {
    return QMessageBox::warning(this, aTitle, QString(aQuestion),
                                "&Yes", "&No", QString(), 1, 1) == 0;
  }

  //warning message for user
  void ask_user_ok(std::string aTitle, std::string aQuestion) {
    QMessageBox::warning(this, QString(aTitle.c_str()), QString(aQuestion.c_str()));
  }

  //functions for setting the color of polygon
  void SetViewBlue(bool a_check) { showBlue->setChecked(a_check); }

  void SetViewRed(bool a_check) { showRed->setChecked(a_check); }

  void SetViewBlack(bool a_check) { showBlack->setChecked(a_check); }

  void SetViewBrown(bool a_check) { showBrown->setChecked(a_check); }

  void SetViewYellow(bool a_check) { showYellow->setChecked(a_check); }

  void SetViewMagenta(bool a_check) { showMagenta->setChecked(a_check); }

  void SetViewAqua(bool a_check) { showAqua->setChecked(a_check); }

  void SetViewResult(bool a_check) { showResult->setChecked(a_check); }

  //changes the set of polygons of a specific type

  void link_GI(CGAL::Qt::GraphicsItem *aGI) {
    QObject::connect(this, SIGNAL(changed()), aGI, SLOT(modelChanged()));
    m_scene.addItem(aGI);
  }

  void unlink_GI(CGAL::Qt::GraphicsItem *aGI) {
    m_scene.removeItem(aGI);
    QObject::disconnect(this, SIGNAL(changed()), aGI, SLOT(modelChanged()));
  }

  void switch_set_type(Curve_set &aSet, int aType);

  void switch_sets_type(int aType);

  bool ensure_circular_mode();

  bool ensure_bezier_mode();

  bool ensure_linear_mode();

  bool ensure_polyline_mode();

  //bool ensure_mink_mode();

};

//constructor for MainWindow
MainWindow::MainWindow() :
  DemosMainWindow(),
  m_bezier_active(false), //default
  m_circular_active(false), //default
  m_polyline_active(false), //default
  m_color_active(0), //default
  m_color_move(1111), //default
  m_color_copy(1111),//default
  m_color_move_ext(111),//default
  m_color_cm(1111), //default
  m_color_visible(7), //default
  m_color_complement(0), //default
  m_blue_int(true), //default
  minkowski_sum_operated(false), //default
  m_red_int(true), //default
  m_black_int(false), //default
  m_brown_int(false), //default
  m_yellow_int(false), //default
  m_magenta_int(false), //default
  m_aqua_int(false), //default
  m_blue_union(true), //default
  m_red_union(true), //default
  m_black_union(false), //default
  m_brown_union(false), //default
  m_yellow_union(false), //default
  m_magenta_union(false), //default
  m_aqua_union(false), //default
  m_blue_sym_diff(true), //default
  m_red_sym_diff(true), //default
  m_black_sym_diff(false), //default
  m_brown_sym_diff(false), //default
  m_yellow_sym_diff(false), //default
  m_magenta_sym_diff(false), //default
  m_aqua_sym_diff(false), //default
  m_blue_mink(true), //default
  m_red_mink(true), //default
  m_black_mink(false), //default
  m_brown_mink(false), //default
  m_yellow_mink(false), //default
  m_magenta_mink(false), //default
  m_aqua_mink(false), //default
  m_visible_black(true), //default
  m_visible_brown(false), //default
  m_visible_yellow(false), //default
  m_visible_magenta(false), //default
  m_visible_aqua(false), //default
  empty_warn(true), // default
  m_disjoint(false), //default
  m_pan(false),
  m_grid(false),
  m_reset(false),
  m_state_num(0) {
  CGAL::set_error_handler(error_handler);
  CGAL::set_warning_handler(error_handler);
  CGAL::set_error_behaviour(CGAL::CONTINUE);

  setupUi(this);

  //
  // Setup the m_scene and the view
  //
  m_scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  m_scene.setSceneRect(-320, -210, 640, 420);
  this->graphicsView->setScene(&m_scene);
  this->graphicsView->setMouseTracking(true);
  this->graphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  this->graphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);


  //this->on_actionInsertLinear_triggered();

  // Turn the vertical axis upside down
  this->graphicsView->scale(2.5, -2.5);

  //adding basic setups

  //todo:need to be given finishing touch
  //The navigation adds zooming and translation functionality to the
  //QGraphicsView.
  this->addNavigation(this->graphicsView);

  //setting the menus
  this->setupStatusBar();
  this->setupOptionsMenu();

  //adding files to select out of polygon
  this->addRecentFiles(this->menuFile, this->actionClose);


  /*  operation_name
  // COMPLEMENT_OP = 0
  // INTERSECTION_OP = 1
  // UNION_OP = 2
  // DIFFERENCE_OP = 3
  // SYMMETRIC_DIFFERENCE_OP = 4
  // MINKOWSKI_SUM_OP = 5
  // RESET_OP = 6
  // COPY_OP = 7
  // MOVE_OP = 8
  // CLEAR_OP = 9
  // DELETEALL_OP = 10
  // START_OP = 11
  */

  //initializing and pushing first state in States_stack
  get_new_state(12);
  actionUndo->setEnabled(false);

  this->setWindowState(Qt::WindowMaximized);

  //! \todo what is it?
  // actionUse_OpenGL->setVisible(false);

  //initializing classes to draw respective polygons using mouse
  m_bezier_input = new CGAL::Qt::GraphicsViewBezierPolygonInput<Bezier_traits>(this, &m_scene);
  m_linear_input = new CGAL::Qt::Graphics_view_linear_polygon_input<Kernel>(this, &m_scene);
  m_circular_input = new CGAL::Qt::Graphics_view_circular_polygon_input<Kernel>(this, &m_scene);
  m_polyline_input = new CGAL::Qt::Graphics_view_polyline_input<Kernel>(this, &m_scene);


  //connecting GUI and the code base
  QObject::connect(m_linear_input, SIGNAL(generate(CGAL::Object)), this,
                   SLOT(processInput(CGAL::Object)));
  //QObject::connect(m_linear_input, SIGNAL(generate(std::variant<Linear_polygon>)), this,
  //                 SLOT(processInput(std::variant<Linear_polygon>)));
  QObject::connect(m_circular_input, SIGNAL(generate(CGAL::Object)), this,
                   SLOT(processInput(CGAL::Object)));
  QObject::connect(m_bezier_input, SIGNAL(generate(CGAL::Object)), this,
                   SLOT(processInput(CGAL::Object)));
  QObject::connect(m_polyline_input, SIGNAL(generate(CGAL::Object)), this,
                   SLOT(processInput(CGAL::Object)));

  m_scene.installEventFilter(m_linear_input);
  QObject::connect(this->actionQuit, SIGNAL(triggered()), this,
                   SLOT(close()));

  //for file handling
  QObject::connect(this, SIGNAL(openRecentFile(QString)), this,
                   SLOT(open(QString)));

  //drawing polygons of various colors
  QObject::connect(drawBlue, SIGNAL(toggled(bool)), this,
                   SLOT(on_drawBlue_toggled(bool)));
  QObject::connect(drawRed, SIGNAL(toggled(bool)), this,
                   SLOT(on_drawRed_toggled(bool)));
  QObject::connect(drawBlack, SIGNAL(toggled(bool)), this,
                   SLOT(on_drawBlack_toggled(bool)));
  QObject::connect(drawBrown, SIGNAL(toggled(bool)), this,
                   SLOT(on_drawBrown_toggled(bool)));
  QObject::connect(drawYellow, SIGNAL(toggled(bool)), this,
                   SLOT(on_drawYellow_toggled(bool)));
  QObject::connect(drawMagenta, SIGNAL(toggled(bool)), this,
                   SLOT(on_drawMagenta_toggled(bool)));
  QObject::connect(drawAqua, SIGNAL(toggled(bool)), this,
                   SLOT(on_drawAqua_toggled(bool)));

  QObject::connect(showBlue, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlue_toggled(bool)));
  QObject::connect(showRed, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRed_toggled(bool)));
  QObject::connect(showBlack, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlack_toggled(bool)));
  QObject::connect(showBrown, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrown_toggled(bool)));
  QObject::connect(showYellow, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellow_toggled(bool)));
  QObject::connect(showMagenta, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagenta_toggled(bool)));
  QObject::connect(showAqua, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAqua_toggled(bool)));
  QObject::connect(showResult, SIGNAL(toggled(bool)), this,
                   SLOT(on_showResult_toggled(bool)));

  /* QObject::connect(sceneDockWidget, SIGNAL(visibilityChanged(bool)), this,
  //                  SLOT(on_sceneDockWidget_visibilityChanged(bool)));
  // QObject::connect(infoDockWidget, SIGNAL(visibilityChanged(bool)), this,
  //                  SLOT(on_infoDockWidget_visibilityChanged(bool)));
  // QObject::connect(consoleDockWidget, SIGNAL(visibilityChanged(bool)), this,
  SLOT(on_consoleDockWidget_visibilityChanged(bool)));*/

  //showing all the side panels
  //todo::implement the console and the information channel quickly
  QObject::connect(showColorBucket, SIGNAL(toggled(bool)), this,
                   SLOT(on_showColorBucket_toggled(bool)));
  QObject::connect(showConsole, SIGNAL(toggled(bool)), this,
                   SLOT(on_showConsole_toggled(bool)));
  QObject::connect(showInfo, SIGNAL(toggled(bool)), this,
                   SLOT(on_showInfo_toggled(bool)));

  QObject::connect(showClear, SIGNAL(toggled(bool)), this,
                   SLOT(on_showClear_toggled(bool)));
  QObject::connect(VisibleHeader, SIGNAL(toggled(bool)), this,
                   SLOT(on_VisibleHeader_toggled(bool)));

  //copy paste for all polygons(not a boolean operation)
  QObject::connect(copyBlue, SIGNAL(toggled(bool)), this,
                   SLOT(on_copyBlue_toggled(bool)));
  QObject::connect(copyRed, SIGNAL(toggled(bool)), this,
                   SLOT(on_copyRed_toggled(bool)));
  QObject::connect(copyBlack, SIGNAL(toggled(bool)), this,
                   SLOT(on_copyBlack_toggled(bool)));
  QObject::connect(copyBrown, SIGNAL(toggled(bool)), this,
                   SLOT(on_copyBrown_toggled(bool)));
  QObject::connect(copyYellow, SIGNAL(toggled(bool)), this,
                   SLOT(on_copyYellow_toggled(bool)));
  QObject::connect(copyMagenta, SIGNAL(toggled(bool)), this,
                   SLOT(on_copyMagenta_toggled(bool)));
  QObject::connect(copyAqua, SIGNAL(toggled(bool)), this,
                   SLOT(on_copyAqua_toggled(bool)));


  //move button usage (not an operation)
  QObject::connect(moveBlue, SIGNAL(toggled(bool)), this,
                   SLOT(on_moveBlue_toggled(bool)));
  QObject::connect(moveRed, SIGNAL(toggled(bool)), this,
                   SLOT(on_moveRed_toggled(bool)));
  QObject::connect(moveBlack, SIGNAL(toggled(bool)), this,
                   SLOT(on_moveBlack_toggled(bool)));
  QObject::connect(moveBrown, SIGNAL(toggled(bool)), this,
                   SLOT(on_moveBrown_toggled(bool)));
  QObject::connect(moveYellow, SIGNAL(toggled(bool)), this,
                   SLOT(on_moveYellow_toggled(bool)));
  QObject::connect(moveMagenta, SIGNAL(toggled(bool)), this,
                   SLOT(on_moveMagenta_toggled(bool)));
  QObject::connect(moveAqua, SIGNAL(toggled(bool)), this,
                   SLOT(on_moveAqua_toggled(bool)));


  //complement operation
  QObject::connect(showBlueComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlueComp_toggled(bool)));
  QObject::connect(showRedComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRedComp_toggled(bool)));
  QObject::connect(showBlackComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlackComp_toggled(bool)));
  QObject::connect(showBrownComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrownComp_toggled(bool)));
  QObject::connect(showYellowComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellowComp_toggled(bool)));
  QObject::connect(showMagentaComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagentaComp_toggled(bool)));
  QObject::connect(showAquaComp, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAquaComp_toggled(bool)));

  //intersection operation
  QObject::connect(showBlueInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlueInt_toggled(bool)));
  QObject::connect(showRedInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRedInt_toggled(bool)));
  QObject::connect(showBlackInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlackInt_toggled(bool)));
  QObject::connect(showBrownInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrownInt_toggled(bool)));
  QObject::connect(showYellowInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellowInt_toggled(bool)));
  QObject::connect(showMagentaInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagentaInt_toggled(bool)));
  QObject::connect(showAquaInt, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAquaInt_toggled(bool)));


  //union operation
  QObject::connect(showBlueUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlueUnion_toggled(bool)));
  QObject::connect(showRedUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRedUnion_toggled(bool)));
  QObject::connect(showBlackUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlackUnion_toggled(bool)));
  QObject::connect(showBrownUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrownUnion_toggled(bool)));
  QObject::connect(showYellowUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellowUnion_toggled(bool)));
  QObject::connect(showMagentaUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagentaUnion_toggled(bool)));
  QObject::connect(showAquaUnion, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAquaUnion_toggled(bool)));

  //difference operation
  QObject::connect(showBlueDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlueDiff_toggled(bool)));
  QObject::connect(showRedDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRedDiff_toggled(bool)));
  QObject::connect(showBlackDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlackDiff_toggled(bool)));
  QObject::connect(showBrownDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrownDiff_toggled(bool)));
  QObject::connect(showYellowDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellowDiff_toggled(bool)));
  QObject::connect(showMagentaDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagentaDiff_toggled(bool)));
  QObject::connect(showAquaDiff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAquaDiff_toggled(bool)));


  //sym_difference operation
  QObject::connect(showBlueSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlueSym_Diff_toggled(bool)));
  QObject::connect(showRedSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRedSym_Diff_toggled(bool)));
  QObject::connect(showBlackSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlackSym_Diff_toggled(bool)));
  QObject::connect(showBrownSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrownSym_Diff_toggled(bool)));
  QObject::connect(showYellowSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellowSym_Diff_toggled(bool)));
  QObject::connect(showMagentaSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagentaSym_Diff_toggled(bool)));
  QObject::connect(showAquaSym_Diff, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAquaSym_Diff_toggled(bool)));


  //minkowski_sum
  QObject::connect(showBlueMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlueMink_Sum_toggled(bool)));
  QObject::connect(showRedMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showRedMink_Sum_toggled(bool)));
  QObject::connect(showBlackMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBlackMink_Sum_toggled(bool)));
  QObject::connect(showBrownMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showBrownMink_Sum_toggled(bool)));
  QObject::connect(showYellowMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showYellowMink_Sum_toggled(bool)));
  QObject::connect(showMagentaMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showMagentaMink_Sum_toggled(bool)));
  QObject::connect(showAquaMink_Sum, SIGNAL(toggled(bool)), this,
                   SLOT(on_showAquaMink_Sum_toggled(bool)));

  //clear
  this->graphicsView->setTransformationAnchor(QGraphicsView::AnchorViewCenter);

  // axis
  QPen *dashedLine{new QPen(QBrush(Qt::black), 1, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin)};
  // dashedLine->setCosmetic(true);
  dashedLine->setWidth(0);
  xAxis->setLine(-155000000000000, 0, 155000000000000, 0);
  xAxis->setPen(*dashedLine);
  yAxis->setLine(0, -100000000000000, 0, 100000000000000);
  yAxis->setPen(*dashedLine);

  /* m_scene.addLine(0,-10000000,0,10000000, QPen(Qt::black));
  // m_scene.addLine(-15500000,0,15500000,0, QPen(Qt::black));

  // // Add the vertical lines first, paint them red
  // for (int x=-15500000; x<=15500000; x+=50)
  //     m_scene.addLine(x,-10000000,x,10000000, QPen(Qt::black));

  // // Now add the horizontal lines, paint them green
  // for (int y=-10000000; y<=10000000; y+=50)
  /     m_scene.addLine(-15500000,y,15500000,y, QPen(Qt::black));
  */

  m_linear_input->mOngoingPieceGI->setPen(sPens[0]);
  m_linear_input->mLinearGI->setPen(sPens[0]);
  m_linear_input->mHandleGI->setPen(sPens[0]);
  actionOpenLinear->setEnabled(true);
  actionOpenDXF->setEnabled(false);
  actionOpenBezier->setEnabled(false);
  actionOpenPolyline->setEnabled(false);

  states_stack.back().blue_set().clear();
  states_stack.back().red_set().clear();
  states_stack.back().black_set().clear();
  states_stack.back().brown_set().clear();
  states_stack.back().yellow_set().clear();
  states_stack.back().magenta_set().clear();
  states_stack.back().aqua_set().clear();
}

//started defining Main Window class functions
void MainWindow::on_showBlue_toggled(bool a_check) {
  states_stack.back().ToggleView(BLUE_GROUP, a_check);
  if (a_check) {
    if (m_color_visible == 0) {
      if (!m_circular_active && !m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_linear_input);
      } else if (!m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_circular_input);
      } else if (!m_polyline_active) {
        m_scene.installEventFilter(m_bezier_input);
      } else {
        m_scene.installEventFilter(m_polyline_input);
      }
    }
    m_color_visible++;
    drawBlue->setEnabled(true);
  } else {
    m_color_visible--;
    if (showRed->isChecked()) drawRed->setChecked(true);
    else if (showBlack->isChecked()) drawBlack->setChecked(true);
    else if (showBrown->isChecked()) drawBrown->setChecked(true);
    else if (showYellow->isChecked()) drawYellow->setChecked(true);
    else if (showMagenta->isChecked()) drawMagenta->setChecked(true);
    else if (showAqua->isChecked()) drawAqua->setChecked(true);
    else {
      m_scene.removeEventFilter(m_linear_input);
      m_scene.removeEventFilter(m_bezier_input);
      m_scene.removeEventFilter(m_circular_input);
      m_scene.removeEventFilter(m_polyline_input);
      m_linear_input->Reset();
      m_circular_input->Reset();
      m_bezier_input->Reset();
      m_polyline_input->Reset();
    }
    drawBlue->setEnabled(false);
  }
  if (m_color_visible == 7) VisibleHeader->setChecked(true);
  //else VisibleHeader->setChecked(false);
}

void MainWindow::on_showRed_toggled(bool a_check) {
  states_stack.back().ToggleView(RED_GROUP, a_check);
  if (a_check) {
    if (m_color_visible == 0) {
      if (!m_circular_active && !m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_linear_input);
      } else if (!m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_circular_input);
      } else if (!m_bezier_active) {
        m_scene.installEventFilter(m_bezier_input);
      } else {
        m_scene.installEventFilter(m_polyline_input);
      }
    }
    m_color_visible++;
    drawRed->setEnabled(true);
  } else {
    m_color_visible--;
    if (showBlue->isChecked()) drawBlue->setChecked(true);
    else if (showBlack->isChecked()) drawBlack->setChecked(true);
    else if (showBrown->isChecked()) drawBrown->setChecked(true);
    else if (showYellow->isChecked()) drawYellow->setChecked(true);
    else if (showMagenta->isChecked()) drawMagenta->setChecked(true);
    else if (showAqua->isChecked()) drawAqua->setChecked(true);
    else {
      m_scene.removeEventFilter(m_linear_input);
      m_scene.removeEventFilter(m_bezier_input);
      m_scene.removeEventFilter(m_circular_input);
      m_scene.removeEventFilter(m_polyline_input);
      m_linear_input->Reset();
      m_circular_input->Reset();
      m_bezier_input->Reset();
      m_polyline_input->Reset();
    }
    drawRed->setEnabled(false);
  }
  if (m_color_visible == 7) VisibleHeader->setChecked(true);
  // else VisibleHeader->setChecked(false);
}

void MainWindow::on_showBlack_toggled(bool a_check) {
  states_stack.back().ToggleView(BLACK_GROUP, a_check);
  if (a_check) {
    if (m_color_visible == 0) {
      if (!m_circular_active && !m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_linear_input);
      } else if (!m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_circular_input);
      } else if (!m_polyline_active) {
        m_scene.installEventFilter(m_bezier_input);
      } else {
        m_scene.installEventFilter(m_polyline_input);
      }
    }
    m_color_visible++;
    drawBlack->setEnabled(true);
  } else {
    m_color_visible--;
    if (showBlue->isChecked()) drawBlue->setChecked(true);
    else if (showRed->isChecked()) drawRed->setChecked(true);
    else if (showBrown->isChecked()) drawBrown->setChecked(true);
    else if (showYellow->isChecked()) drawYellow->setChecked(true);
    else if (showMagenta->isChecked()) drawMagenta->setChecked(true);
    else if (showAqua->isChecked()) drawAqua->setChecked(true);
    else {
      m_scene.removeEventFilter(m_linear_input);
      m_scene.removeEventFilter(m_bezier_input);
      m_scene.removeEventFilter(m_circular_input);
      m_scene.removeEventFilter(m_polyline_input);

      m_linear_input->Reset();
      m_circular_input->Reset();
      m_bezier_input->Reset();
      m_polyline_input->Reset();
    }
    drawBlack->setEnabled(false);
  }
  if (m_color_visible == 7) VisibleHeader->setChecked(true);
  //else VisibleHeader->setChecked(false);
}

void MainWindow::on_showBrown_toggled(bool a_check) {
  states_stack.back().ToggleView(BROWN_GROUP, a_check);
  if (a_check) {
    if (m_color_visible == 0) {
      if (!m_circular_active && !m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_linear_input);
      } else if (!m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_circular_input);
      } else if (m_polyline_active) {
        m_scene.installEventFilter(m_bezier_input);
      } else {
        m_scene.installEventFilter(m_polyline_input);
      }
    }
    m_color_visible++;
    drawBrown->setEnabled(true);
  } else {
    m_color_visible--;
    if (showBlue->isChecked()) drawBlue->setChecked(true);
    else if (showRed->isChecked()) drawRed->setChecked(true);
    else if (showBlack->isChecked()) drawBlack->setChecked(true);
    else if (showYellow->isChecked()) drawYellow->setChecked(true);
    else if (showMagenta->isChecked()) drawMagenta->setChecked(true);
    else if (showAqua->isChecked()) drawAqua->setChecked(true);
    else {
      m_scene.removeEventFilter(m_linear_input);
      m_scene.removeEventFilter(m_bezier_input);
      m_scene.removeEventFilter(m_circular_input);
      m_scene.removeEventFilter(m_polyline_input);
      m_linear_input->Reset();
      m_circular_input->Reset();
      m_bezier_input->Reset();
      m_polyline_input->Reset();
    }
    drawBrown->setEnabled(false);
  }
  if (m_color_visible == 7) VisibleHeader->setChecked(true);
  //else VisibleHeader->setChecked(false);
}

void MainWindow::on_showYellow_toggled(bool a_check) {
  states_stack.back().ToggleView(YELLOW_GROUP, a_check);
  if (a_check) {
    if (m_color_visible == 0) {
      if (!m_circular_active && !m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_linear_input);
      }
      else if (!m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_circular_input);
      }
      else if (!m_polyline_active) {
        m_scene.installEventFilter(m_bezier_input);
      }
      else {
        m_scene.installEventFilter(m_polyline_input);
      }
    }
    m_color_visible++;
    drawYellow->setEnabled(true);
  } else {
    m_color_visible--;
    if (showBlue->isChecked()) drawBlue->setChecked(true);
    else if (showRed->isChecked()) drawRed->setChecked(true);
    else if (showBlack->isChecked()) drawBlack->setChecked(true);
    else if (showBrown->isChecked()) drawBrown->setChecked(true);
    else if (showMagenta->isChecked()) drawMagenta->setChecked(true);
    else if (showAqua->isChecked()) drawAqua->setChecked(true);
    else {
      m_scene.removeEventFilter(m_linear_input);
      m_scene.removeEventFilter(m_bezier_input);
      m_scene.removeEventFilter(m_circular_input);
      m_scene.removeEventFilter(m_polyline_input);
      m_linear_input->Reset();
      m_circular_input->Reset();
      m_bezier_input->Reset();
      m_polyline_input->Reset();
    }
    drawYellow->setEnabled(false);
  }
  if (m_color_visible == 7) VisibleHeader->setChecked(true);
  //else VisibleHeader->setChecked(false);
}

void MainWindow::on_showMagenta_toggled(bool a_check) {
  states_stack.back().ToggleView(MAGENTA_GROUP, a_check);
  if (a_check) {
    if (m_color_visible == 0) {
      if (!m_circular_active && !m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_linear_input);
      } else if (!m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_circular_input);
      } else if (m_polyline_active) {
        m_scene.installEventFilter(m_bezier_input);
      } else {
        m_scene.installEventFilter(m_polyline_input);
      }
    }
    m_color_visible++;
    drawMagenta->setEnabled(true);
  } else {
    m_color_visible--;
    if (showBlue->isChecked()) drawBlue->setChecked(true);
    else if (showRed->isChecked()) drawRed->setChecked(true);
    else if (showBlack->isChecked()) drawBlack->setChecked(true);
    else if (showBrown->isChecked()) drawBrown->setChecked(true);
    else if (showYellow->isChecked()) drawYellow->setChecked(true);
    else if (showAqua->isChecked()) drawAqua->setChecked(true);
    else {
      m_scene.removeEventFilter(m_linear_input);
      m_scene.removeEventFilter(m_bezier_input);
      m_scene.removeEventFilter(m_circular_input);
      m_scene.removeEventFilter(m_polyline_input);
      m_linear_input->Reset();
      m_circular_input->Reset();
      m_bezier_input->Reset();
      m_polyline_input->Reset();
    }
    drawMagenta->setEnabled(false);
  }
  if (m_color_visible == 7) VisibleHeader->setChecked(true);
  //else VisibleHeader->setChecked(false);
}

void MainWindow::on_showAqua_toggled(bool a_check) {
  states_stack.back().ToggleView(AQUA_GROUP, a_check);
  if (a_check) {
    if (m_color_visible == 0) {
      if (!m_circular_active && !m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_linear_input);
      } else if (!m_bezier_active && !m_polyline_active) {
        m_scene.installEventFilter(m_circular_input);
      } else if (!m_polyline_active) {
        m_scene.installEventFilter(m_bezier_input);
      } else {
        m_scene.installEventFilter(m_polyline_input);
      }
    }
    m_color_visible++;
    drawAqua->setEnabled(true);
  } else {
    m_color_visible--;
    if (showBlue->isChecked()) drawBlue->setChecked(true);
    else if (showRed->isChecked()) drawRed->setChecked(true);
    else if (showBlack->isChecked()) drawBlack->setChecked(true);
    else if (showBrown->isChecked()) drawBrown->setChecked(true);
    else if (showYellow->isChecked()) drawYellow->setChecked(true);
    else if (showMagenta->isChecked()) drawMagenta->setChecked(true);
    else {
      m_scene.removeEventFilter(m_linear_input);
      m_scene.removeEventFilter(m_bezier_input);
      m_scene.removeEventFilter(m_circular_input);
      m_scene.removeEventFilter(m_polyline_input);
      m_linear_input->Reset();
      m_circular_input->Reset();
      m_bezier_input->Reset();
      m_polyline_input->Reset();
    }
    drawAqua->setEnabled(false);
  }
  if (m_color_visible == 7) VisibleHeader->setChecked(true);
  //else VisibleHeader->setChecked(false);
}

void MainWindow::on_showResult_toggled(bool a_check) { states_stack.back().ToggleView(RESULT_GROUP, a_check); }

void MainWindow::on_VisibleHeader_toggled(bool a_check) {
  if (a_check) {
    showRed->setChecked(true);
    showBlue->setChecked(true);
    showBlack->setChecked(true);
    showBrown->setChecked(true);
    showYellow->setChecked(true);
    showMagenta->setChecked(true);
    showAqua->setChecked(true);
    drawBlue->setChecked(true);
  } else {
    showRed->setChecked(false);
    showBlue->setChecked(false);
    showBlack->setChecked(false);
    showBrown->setChecked(false);
    showYellow->setChecked(false);
    showMagenta->setChecked(false);
    showAqua->setChecked(false);
  }
}

//clear button
void MainWindow::on_showClear_toggled(bool a_check) {
  if (a_check) {
    clearRed->setChecked(true);
    clearBlue->setChecked(true);
    clearBrown->setChecked(true);
    clearYellow->setChecked(true);
    clearMagenta->setChecked(true);
    clearAqua->setChecked(true);
    clearBlack->setChecked(true);
  } else {
    clearRed->setChecked(false);
    clearBlue->setChecked(false);
    clearBrown->setChecked(false);
    clearYellow->setChecked(false);
    clearMagenta->setChecked(false);
    clearAqua->setChecked(false);
    clearBlack->setChecked(false);
  }
}

void MainWindow::on_sceneDockWidget_visibilityChanged(bool a_check) {
  if (a_check) showColorBucket->setChecked(true);
  else if (!this->windowState().testFlag(Qt::WindowMinimized)) {
    //sceneDockWidget->setVisible(false);
    showColorBucket->setChecked(false);
  }
}

void MainWindow::on_infoDockWidget_visibilityChanged(bool a_check) {
  if (a_check) showInfo->setChecked(true);
  else if (!this->windowState().testFlag(Qt::WindowMinimized)) {
    //infoDockWidget->setVisible(false);
    showInfo->setChecked(false);
  }
}

void MainWindow::on_consoleDockWidget_visibilityChanged(bool a_check) {
  if (a_check) showConsole->setChecked(true);
  else if (!this->windowState().testFlag(Qt::WindowMinimized)) {
    //consoleDockWidget->setVisible(false);
    showConsole->setChecked(false);
  }
}

void MainWindow::on_showColorBucket_toggled(bool a_check) {
  if (a_check) sceneDockWidget->setVisible(true);
  else sceneDockWidget->setVisible(false);
}

void MainWindow::on_showConsole_toggled(bool a_check) {
  if (a_check) consoleDockWidget->setVisible(true);
  else consoleDockWidget->setVisible(false);
}

void MainWindow::on_showInfo_toggled(bool a_check) {
  if (a_check) infoDockWidget->setVisible(true);
  else infoDockWidget->setVisible(false);
}

//complement operation button clicked
void MainWindow::on_showBlueComp_toggled(bool aCheck) {
  if (aCheck) {
    m_color_complement = 0;
    showRedComp->setChecked(false);
    showBlackComp->setChecked(false);
    showBrownComp->setChecked(false);
    showYellowComp->setChecked(false);
    showMagentaComp->setChecked(false);
    showAquaComp->setChecked(false);
  } else {
    showBlueComp->setChecked(false);
  }
}

void MainWindow::on_showRedComp_toggled(bool aCheck) {
  if (aCheck) {
    m_color_complement = 1;
    showBlueComp->setChecked(false);
    showBlackComp->setChecked(false);
    showBrownComp->setChecked(false);
    showYellowComp->setChecked(false);
    showMagentaComp->setChecked(false);
    showAquaComp->setChecked(false);
  } else {
    showRedComp->setChecked(false);
  }
}

void MainWindow::on_showBlackComp_toggled(bool aCheck) {
  if (aCheck) {
    m_color_complement = 2;
    showRedComp->setChecked(false);
    showBlueComp->setChecked(false);
    showBrownComp->setChecked(false);
    showYellowComp->setChecked(false);
    showMagentaComp->setChecked(false);
    showAquaComp->setChecked(false);
  }
  else {
    showBlackComp->setChecked(false);
  }
}

void MainWindow::on_showBrownComp_toggled(bool aCheck) {
  if (aCheck) {
    m_color_complement = 3;
    showRedComp->setChecked(false);
    showBlackComp->setChecked(false);
    showBlueComp->setChecked(false);
    showYellowComp->setChecked(false);
    showMagentaComp->setChecked(false);
    showAquaComp->setChecked(false);
  } else {
    showBrownComp->setChecked(false);
  }
}

void MainWindow::on_showYellowComp_toggled(bool aCheck) {
  if (aCheck) {
    m_color_complement = 4;
    showRedComp->setChecked(false);
    showBlackComp->setChecked(false);
    showBrownComp->setChecked(false);
    showBlueComp->setChecked(false);
    showMagentaComp->setChecked(false);
    showAquaComp->setChecked(false);
  } else {
    showYellowComp->setChecked(false);
  }
}

void MainWindow::on_showMagentaComp_toggled(bool aCheck) {
  if (aCheck) {
    m_color_complement = 5;
    showRedComp->setChecked(false);
    showBlackComp->setChecked(false);
    showBrownComp->setChecked(false);
    showYellowComp->setChecked(false);
    showBlueComp->setChecked(false);
    showAquaComp->setChecked(false);
  } else {
    showMagentaComp->setChecked(false);
  }
}

void MainWindow::on_showAquaComp_toggled(bool aCheck) {
  if (aCheck) {
    m_color_complement = 6;
    showRedComp->setChecked(false);
    showBlackComp->setChecked(false);
    showBrownComp->setChecked(false);
    showYellowComp->setChecked(false);
    showMagentaComp->setChecked(false);
    showBlueComp->setChecked(false);
  }
  else {
    showAquaComp->setChecked(false);
  }
}

//copy button clicked
void MainWindow::on_copyBlue_toggled(bool aCheck) {
  if (aCheck) {
    m_color_copy = 0;
    copyRed->setChecked(false);
    copyBlack->setChecked(false);
    copyBrown->setChecked(false);
    copyYellow->setChecked(false);
    copyMagenta->setChecked(false);
    copyAqua->setChecked(false);
  }
  else {
    copyBlue->setChecked(false);
  }
}

void MainWindow::on_copyRed_toggled(bool aCheck) {
  if (aCheck) {
    m_color_copy = 1;
    copyBlue->setChecked(false);
    copyBlack->setChecked(false);
    copyBrown->setChecked(false);
    copyYellow->setChecked(false);
    copyMagenta->setChecked(false);
    copyAqua->setChecked(false);
  }
  else {
    copyRed->setChecked(false);
  }
}

void MainWindow::on_copyBlack_toggled(bool aCheck) {
  if (aCheck) {
    m_color_copy = 2;
    copyRed->setChecked(false);
    copyBlue->setChecked(false);
    copyBrown->setChecked(false);
    copyYellow->setChecked(false);
    copyMagenta->setChecked(false);
    copyAqua->setChecked(false);
  }
  else {
    copyBlack->setChecked(false);
  }
}

void MainWindow::on_copyBrown_toggled(bool aCheck) {
  if (aCheck) {
    m_color_copy = 3;
    copyRed->setChecked(false);
    copyBlack->setChecked(false);
    copyBlue->setChecked(false);
    copyYellow->setChecked(false);
    copyMagenta->setChecked(false);
    copyAqua->setChecked(false);
  }
  else {
    copyBrown->setChecked(false);
  }
}

void MainWindow::on_copyYellow_toggled(bool aCheck) {
  if (aCheck) {
    m_color_copy = 4;
    copyRed->setChecked(false);
    copyBlack->setChecked(false);
    copyBrown->setChecked(false);
    copyBlue->setChecked(false);
    copyMagenta->setChecked(false);
    copyAqua->setChecked(false);
  } else {
    copyYellow->setChecked(false);
  }
}

void MainWindow::on_copyMagenta_toggled(bool aCheck) {
  if (aCheck) {
    m_color_copy = 5;
    copyRed->setChecked(false);
    copyBlack->setChecked(false);
    copyBrown->setChecked(false);
    copyYellow->setChecked(false);
    copyBlue->setChecked(false);
    copyAqua->setChecked(false);
  } else {
    copyMagenta->setChecked(false);
  }
}

void MainWindow::on_copyAqua_toggled(bool aCheck) {
  if (aCheck) {
    m_color_copy = 6;
    copyRed->setChecked(false);
    copyBlack->setChecked(false);
    copyBrown->setChecked(false);
    copyYellow->setChecked(false);
    copyMagenta->setChecked(false);
    copyBlue->setChecked(false);
  }
  else {
    copyAqua->setChecked(false);
  }
}

//move function clicked
void MainWindow::on_moveBlue_toggled(bool aCheck) {
  if (aCheck) {
    m_color_move_ext = 0;
    moveRed->setChecked(false);
    moveBlack->setChecked(false);
    moveBrown->setChecked(false);
    moveYellow->setChecked(false);
    moveMagenta->setChecked(false);
    moveAqua->setChecked(false);
  } else {
    moveBlue->setChecked(false);
  }
}

void MainWindow::on_moveRed_toggled(bool aCheck) {
  if (aCheck) {
    m_color_move_ext = 1;
    moveBlue->setChecked(false);
    moveBlack->setChecked(false);
    moveBrown->setChecked(false);
    moveYellow->setChecked(false);
    moveMagenta->setChecked(false);
    moveAqua->setChecked(false);
  } else {
    moveRed->setChecked(false);
  }
}

void MainWindow::on_moveBlack_toggled(bool aCheck) {
  if (aCheck) {
    m_color_move_ext = 2;
    moveRed->setChecked(false);
    moveBlue->setChecked(false);
    moveBrown->setChecked(false);
    moveYellow->setChecked(false);
    moveMagenta->setChecked(false);
    moveAqua->setChecked(false);
  }
  else {
    moveBlack->setChecked(false);
  }
}

void MainWindow::on_moveBrown_toggled(bool aCheck) {
  if (aCheck) {
    m_color_move_ext = 3;
    moveRed->setChecked(false);
    moveBlack->setChecked(false);
    moveBlue->setChecked(false);
    moveYellow->setChecked(false);
    moveMagenta->setChecked(false);
    moveAqua->setChecked(false);
  }
  else {
    copyBrown->setChecked(false);
  }
}

void MainWindow::on_moveYellow_toggled(bool aCheck) {
  if (aCheck) {
    m_color_move_ext = 4;
    moveRed->setChecked(false);
    moveBlack->setChecked(false);
    moveBrown->setChecked(false);
    moveBlue->setChecked(false);
    moveMagenta->setChecked(false);
    moveAqua->setChecked(false);
  }
  else {
    moveYellow->setChecked(false);
  }
}

void MainWindow::on_moveMagenta_toggled(bool aCheck) {
  if (aCheck) {
    m_color_move_ext = 5;
    moveRed->setChecked(false);
    moveBlack->setChecked(false);
    moveBrown->setChecked(false);
    moveYellow->setChecked(false);
    moveBlue->setChecked(false);
    moveAqua->setChecked(false);
  } else {
    moveMagenta->setChecked(false);
  }
}

void MainWindow::on_moveAqua_toggled(bool aCheck) {
  if (aCheck) {
    m_color_move_ext = 6;
    moveRed->setChecked(false);
    moveBlack->setChecked(false);
    moveBrown->setChecked(false);
    moveYellow->setChecked(false);
    moveMagenta->setChecked(false);
    moveBlue->setChecked(false);
  } else {
    moveAqua->setChecked(false);
  }
}

//show button checkbox clicked
void MainWindow::on_showBlueInt_toggled(bool aCheck) {
  if (aCheck) {
    m_blue_int = true;
  } else {
    m_blue_int = false;
  }
}

void MainWindow::on_showRedInt_toggled(bool aCheck) {
  if (aCheck) {
    m_red_int = true;
  } else {
    m_red_int = false;
  }

}

void MainWindow::on_showBlackInt_toggled(bool aCheck) {
  if (aCheck) {
    m_black_int = true;
  } else {
    m_black_int = false;
  }

}

void MainWindow::on_showBrownInt_toggled(bool aCheck) {
  if (aCheck) {
    m_brown_int = true;
  } else {
    m_brown_int = false;
  }

}

void MainWindow::on_showYellowInt_toggled(bool aCheck) {

  if (aCheck) {
    m_yellow_int = true;
  } else {
    m_yellow_int = false;
  }
}

void MainWindow::on_showMagentaInt_toggled(bool aCheck) {
  if (aCheck) {
    m_magenta_int = true;
  } else {
    m_magenta_int = false;
  }

}

void MainWindow::on_showAquaInt_toggled(bool aCheck) {
  if (aCheck) {
    m_aqua_int = true;
  } else {
    m_aqua_int = false;
  }
}

//union operation button clicked
void MainWindow::on_showBlueUnion_toggled(bool aCheck) {
  if (aCheck) {
    m_blue_union = true;
  } else {
    m_blue_union = false;
  }
}

void MainWindow::on_showRedUnion_toggled(bool aCheck) {
  if (aCheck) {
    m_red_union = true;
  } else {
    m_red_union = false;
  }

}

void MainWindow::on_showBlackUnion_toggled(bool aCheck) {

  if (aCheck) {
    m_black_union = true;
  } else {
    m_black_union = false;
  }
}

void MainWindow::on_showBrownUnion_toggled(bool aCheck) {

  if (aCheck) {
    m_brown_union = true;
  } else {
    m_brown_union = false;
  }
}

void MainWindow::on_showYellowUnion_toggled(bool aCheck) {

  if (aCheck) {
    m_yellow_union = true;
  } else {
    m_yellow_union = false;
  }
}

void MainWindow::on_showMagentaUnion_toggled(bool aCheck) {
  if (aCheck) {
    m_magenta_union = true;
  } else {
    m_magenta_union = false;
  }
}

void MainWindow::on_showAquaUnion_toggled(bool aCheck) {
  if (aCheck) {
    m_aqua_union = true;
  } else {
    m_aqua_union = false;
  }
}

//difference button clicked
void MainWindow::on_showBlueDiff_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showRedDiff->isChecked()) count++;
    if (showBlackDiff->isChecked()) count++;
    if (showBrownDiff->isChecked()) count++;
    if (showYellowDiff->isChecked()) count++;
    if (showMagentaDiff->isChecked()) count++;
    if (showAquaDiff->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionDifferenceH->setEnabled(false);
    } else {
      actionDifferenceH->setEnabled(true);
    }
  } else {
    size_t count = 0;
    if (showRedDiff->isChecked()) count++;
    if (showBlackDiff->isChecked()) count++;
    if (showBrownDiff->isChecked()) count++;
    if (showYellowDiff->isChecked()) count++;
    if (showMagentaDiff->isChecked()) count++;
    if (showAquaDiff->isChecked()) count++;

    if (count == 2) {
      actionDifferenceH->setEnabled(true);
    } else {
      actionDifferenceH->setEnabled(false);
    }
  }
}

void MainWindow::on_showRedDiff_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showAquaDiff->isChecked()) count++;
    if (showBlackDiff->isChecked()) count++;
    if (showBrownDiff->isChecked()) count++;
    if (showYellowDiff->isChecked()) count++;
    if (showMagentaDiff->isChecked()) count++;
    if (showBlueDiff->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionDifferenceH->setEnabled(false);
    } else {
      actionDifferenceH->setEnabled(true);
    }
  } else {
    size_t count = 0;
    if (showAquaDiff->isChecked()) count++;
    if (showBlackDiff->isChecked()) count++;
    if (showBrownDiff->isChecked()) count++;
    if (showYellowDiff->isChecked()) count++;
    if (showMagentaDiff->isChecked()) count++;
    if (showBlueDiff->isChecked()) count++;

    if (count == 2) {
      actionDifferenceH->setEnabled(true);
    } else {
      actionDifferenceH->setEnabled(false);
    }
  }
}

void MainWindow::on_showBlackDiff_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showRedDiff->isChecked()) count++;
    if (showAquaDiff->isChecked()) count++;
    if (showBrownDiff->isChecked()) count++;
    if (showYellowDiff->isChecked()) count++;
    if (showMagentaDiff->isChecked()) count++;
    if (showBlueDiff->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionDifferenceH->setEnabled(false);
    } else {
      actionDifferenceH->setEnabled(true);
    }
  } else {
    size_t count = 0;
    if (showRedDiff->isChecked()) count++;
    if (showAquaDiff->isChecked()) count++;
    if (showBrownDiff->isChecked()) count++;
    if (showYellowDiff->isChecked()) count++;
    if (showMagentaDiff->isChecked()) count++;
    if (showBlueDiff->isChecked()) count++;

    if (count == 2) {
      actionDifferenceH->setEnabled(true);
    } else {
      actionDifferenceH->setEnabled(false);
    }
  }
}

void MainWindow::on_showBrownDiff_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showRedDiff->isChecked()) count++;
    if (showBlackDiff->isChecked()) count++;
    if (showAquaDiff->isChecked()) count++;
    if (showYellowDiff->isChecked()) count++;
    if (showMagentaDiff->isChecked()) count++;
    if (showBlueDiff->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionDifferenceH->setEnabled(false);
    } else {
      actionDifferenceH->setEnabled(true);
    }
  } else {
    size_t count = 0;
    if (showRedDiff->isChecked()) count++;
    if (showBlackDiff->isChecked()) count++;
    if (showAquaDiff->isChecked()) count++;
    if (showYellowDiff->isChecked()) count++;
    if (showMagentaDiff->isChecked()) count++;
    if (showBlueDiff->isChecked()) count++;

    if (count == 2) {
      actionDifferenceH->setEnabled(true);
    } else {
      actionDifferenceH->setEnabled(false);
    }
  }
}

void MainWindow::on_showYellowDiff_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showRedDiff->isChecked()) count++;
    if (showBlackDiff->isChecked()) count++;
    if (showBrownDiff->isChecked()) count++;
    if (showAquaDiff->isChecked()) count++;
    if (showMagentaDiff->isChecked()) count++;
    if (showBlueDiff->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionDifferenceH->setEnabled(false);
    } else {
      actionDifferenceH->setEnabled(true);
    }
  } else {
    size_t count = 0;
    if (showRedDiff->isChecked()) count++;
    if (showBlackDiff->isChecked()) count++;
    if (showBrownDiff->isChecked()) count++;
    if (showAquaDiff->isChecked()) count++;
    if (showMagentaDiff->isChecked()) count++;
    if (showBlueDiff->isChecked()) count++;

    if (count == 2) {
      actionDifferenceH->setEnabled(true);
    } else {
      actionDifferenceH->setEnabled(false);
    }
  }
}

void MainWindow::on_showMagentaDiff_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showRedDiff->isChecked()) count++;
    if (showBlackDiff->isChecked()) count++;
    if (showBrownDiff->isChecked()) count++;
    if (showYellowDiff->isChecked()) count++;
    if (showAquaDiff->isChecked()) count++;
    if (showBlueDiff->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionDifferenceH->setEnabled(false);
    } else {
      actionDifferenceH->setEnabled(true);
    }
  } else {
    size_t count = 0;
    if (showRedDiff->isChecked()) count++;
    if (showBlackDiff->isChecked()) count++;
    if (showBrownDiff->isChecked()) count++;
    if (showYellowDiff->isChecked()) count++;
    if (showAquaDiff->isChecked()) count++;
    if (showBlueDiff->isChecked()) count++;

    if (count == 2) {
      actionDifferenceH->setEnabled(true);
    } else {
      actionDifferenceH->setEnabled(false);
    }
  }
}

void MainWindow::on_showAquaDiff_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showRedDiff->isChecked()) count++;
    if (showBlackDiff->isChecked()) count++;
    if (showBrownDiff->isChecked()) count++;
    if (showYellowDiff->isChecked()) count++;
    if (showMagentaDiff->isChecked()) count++;
    if (showBlueDiff->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionDifferenceH->setEnabled(false);
    } else {
      actionDifferenceH->setEnabled(true);
    }
  }
  else {
    size_t count = 0;
    if (showRedDiff->isChecked()) count++;
    if (showBlackDiff->isChecked()) count++;
    if (showBrownDiff->isChecked()) count++;
    if (showYellowDiff->isChecked()) count++;
    if (showMagentaDiff->isChecked()) count++;
    if (showBlueDiff->isChecked()) count++;

    if (count == 2) {
      actionDifferenceH->setEnabled(true);
    } else {
      actionDifferenceH->setEnabled(false);
    }
  }
}

//symmetric difference button clicked
void MainWindow::on_showBlueSym_Diff_toggled(bool aCheck) {
  if (aCheck) {
    m_blue_sym_diff = true;
  } else {
    m_blue_sym_diff = false;
  }
}

void MainWindow::on_showRedSym_Diff_toggled(bool aCheck) {

  if (aCheck) {
    m_red_sym_diff = true;
  } else {
    m_red_sym_diff = false;
  }
}

void MainWindow::on_showBlackSym_Diff_toggled(bool aCheck) {
  if (aCheck) {
    m_black_sym_diff = true;
  } else {
    m_black_sym_diff = false;
  }
}

void MainWindow::on_showBrownSym_Diff_toggled(bool aCheck) {
  if (aCheck) {
    m_brown_sym_diff = true;
  } else {
    m_brown_sym_diff = false;
  }
}

void MainWindow::on_showYellowSym_Diff_toggled(bool aCheck) {
  if (aCheck) {
    m_yellow_sym_diff = true;
  } else {
    m_yellow_sym_diff = false;
  }
}

void MainWindow::on_showMagentaSym_Diff_toggled(bool aCheck) {
  if (aCheck) {
    m_magenta_sym_diff = true;
  } else {
    m_magenta_sym_diff = false;
  }
}

void MainWindow::on_showAquaSym_Diff_toggled(bool aCheck) {
  if (aCheck) {
    m_aqua_sym_diff = true;
  } else {
    m_aqua_sym_diff = false;
  }
}

//Minkowski button clicked
void MainWindow::on_showBlueMink_Sum_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showRedMink_Sum->isChecked()) count++;
    if (showBlackMink_Sum->isChecked()) count++;
    if (showBrownMink_Sum->isChecked()) count++;
    if (showYellowMink_Sum->isChecked()) count++;
    if (showMagentaMink_Sum->isChecked()) count++;
    if (showAquaMink_Sum->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionMinkowski_SumH->setEnabled(false);
    } else {
      actionMinkowski_SumH->setEnabled(true);
    }
  } else {
    size_t count = 0;
    if (showRedMink_Sum->isChecked()) count++;
    if (showBlackMink_Sum->isChecked()) count++;
    if (showBrownMink_Sum->isChecked()) count++;
    if (showYellowMink_Sum->isChecked()) count++;
    if (showMagentaMink_Sum->isChecked()) count++;
    if (showAquaMink_Sum->isChecked()) count++;

    if (count == 2) {
      actionMinkowski_SumH->setEnabled(true);
    } else {
      actionMinkowski_SumH->setEnabled(false);
    }
  }
}

void MainWindow::on_showRedMink_Sum_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showAquaMink_Sum->isChecked()) count++;
    if (showBlackMink_Sum->isChecked()) count++;
    if (showBrownMink_Sum->isChecked()) count++;
    if (showYellowMink_Sum->isChecked()) count++;
    if (showMagentaMink_Sum->isChecked()) count++;
    if (showBlueMink_Sum->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionMinkowski_SumH->setEnabled(false);
    } else {
      actionMinkowski_SumH->setEnabled(true);
    }
  } else {
    size_t count = 0;
    if (showAquaMink_Sum->isChecked()) count++;
    if (showBlackMink_Sum->isChecked()) count++;
    if (showBrownMink_Sum->isChecked()) count++;
    if (showYellowMink_Sum->isChecked()) count++;
    if (showMagentaMink_Sum->isChecked()) count++;
    if (showBlueMink_Sum->isChecked()) count++;

    if (count == 2) {
      actionMinkowski_SumH->setEnabled(true);
    } else {
      actionMinkowski_SumH->setEnabled(false);
    }
  }
}

void MainWindow::on_showBlackMink_Sum_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showRedMink_Sum->isChecked()) count++;
    if (showAquaMink_Sum->isChecked()) count++;
    if (showBrownMink_Sum->isChecked()) count++;
    if (showYellowMink_Sum->isChecked()) count++;
    if (showMagentaMink_Sum->isChecked()) count++;
    if (showBlueMink_Sum->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionMinkowski_SumH->setEnabled(false);
    } else {
      actionMinkowski_SumH->setEnabled(true);
    }
  } else {
    size_t count = 0;
    if (showRedMink_Sum->isChecked()) count++;
    if (showAquaMink_Sum->isChecked()) count++;
    if (showBrownMink_Sum->isChecked()) count++;
    if (showYellowMink_Sum->isChecked()) count++;
    if (showMagentaMink_Sum->isChecked()) count++;
    if (showBlueMink_Sum->isChecked()) count++;

    if (count == 2) {
      actionMinkowski_SumH->setEnabled(true);
    } else {
      actionMinkowski_SumH->setEnabled(false);
    }
  }
}

void MainWindow::on_showBrownMink_Sum_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showRedMink_Sum->isChecked()) count++;
    if (showBlackMink_Sum->isChecked()) count++;
    if (showAquaMink_Sum->isChecked()) count++;
    if (showYellowMink_Sum->isChecked()) count++;
    if (showMagentaMink_Sum->isChecked()) count++;
    if (showBlueMink_Sum->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionMinkowski_SumH->setEnabled(false);
    } else {
      actionMinkowski_SumH->setEnabled(true);
    }
  } else {
    size_t count = 0;
    if (showRedMink_Sum->isChecked()) count++;
    if (showBlackMink_Sum->isChecked()) count++;
    if (showAquaMink_Sum->isChecked()) count++;
    if (showYellowMink_Sum->isChecked()) count++;
    if (showMagentaMink_Sum->isChecked()) count++;
    if (showBlueMink_Sum->isChecked()) count++;

    if (count == 2) {
      actionMinkowski_SumH->setEnabled(true);
    } else {
      actionMinkowski_SumH->setEnabled(false);
    }
  }
}

void MainWindow::on_showYellowMink_Sum_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showRedMink_Sum->isChecked()) count++;
    if (showBlackMink_Sum->isChecked()) count++;
    if (showBrownMink_Sum->isChecked()) count++;
    if (showAquaMink_Sum->isChecked()) count++;
    if (showMagentaMink_Sum->isChecked()) count++;
    if (showBlueMink_Sum->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionMinkowski_SumH->setEnabled(false);
    } else {
      actionMinkowski_SumH->setEnabled(true);
    }
  } else {
    size_t count = 0;
    if (showRedMink_Sum->isChecked()) count++;
    if (showBlackMink_Sum->isChecked()) count++;
    if (showBrownMink_Sum->isChecked()) count++;
    if (showAquaMink_Sum->isChecked()) count++;
    if (showMagentaMink_Sum->isChecked()) count++;
    if (showBlueMink_Sum->isChecked()) count++;

    if (count == 2) {
      actionMinkowski_SumH->setEnabled(true);
    } else {
      actionMinkowski_SumH->setEnabled(false);
    }
  }
}

void MainWindow::on_showMagentaMink_Sum_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showRedMink_Sum->isChecked()) count++;
    if (showBlackMink_Sum->isChecked()) count++;
    if (showBrownMink_Sum->isChecked()) count++;
    if (showYellowMink_Sum->isChecked()) count++;
    if (showAquaMink_Sum->isChecked()) count++;
    if (showBlueMink_Sum->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionMinkowski_SumH->setEnabled(false);
    } else {
      actionMinkowski_SumH->setEnabled(true);
    }
  } else {
    size_t count = 0;
    if (showRedMink_Sum->isChecked()) count++;
    if (showBlackMink_Sum->isChecked()) count++;
    if (showBrownMink_Sum->isChecked()) count++;
    if (showYellowMink_Sum->isChecked()) count++;
    if (showAquaMink_Sum->isChecked()) count++;
    if (showBlueMink_Sum->isChecked()) count++;

    if (count == 2) {
      actionMinkowski_SumH->setEnabled(true);
    } else {
      actionMinkowski_SumH->setEnabled(false);
    }
  }
}

void MainWindow::on_showAquaMink_Sum_toggled(bool aCheck) {
  if (aCheck) {
    size_t count = 0;
    if (showRedMink_Sum->isChecked()) count++;
    if (showBlackMink_Sum->isChecked()) count++;
    if (showBrownMink_Sum->isChecked()) count++;
    if (showYellowMink_Sum->isChecked()) count++;
    if (showMagentaMink_Sum->isChecked()) count++;
    if (showBlueMink_Sum->isChecked()) count++;

    if (count >= 2 or count == 0) {
      actionMinkowski_SumH->setEnabled(false);
    } else {
      actionMinkowski_SumH->setEnabled(true);
    }
  } else {
    size_t count = 0;
    if (showRedMink_Sum->isChecked()) count++;
    if (showBlackMink_Sum->isChecked()) count++;
    if (showBrownMink_Sum->isChecked()) count++;
    if (showYellowMink_Sum->isChecked()) count++;
    if (showMagentaMink_Sum->isChecked()) count++;
    if (showBlueMink_Sum->isChecked()) count++;

    if (count == 2) {
      actionMinkowski_SumH->setEnabled(true);
    } else {
      actionMinkowski_SumH->setEnabled(false);
    }
  }
}

//////////###################### Add Color Plus Button ###########################////////////
//done in GSoC 2020

//add color button clicked
void MainWindow::on_actionAddColor_triggered() {
  //these colors are not present by default
  if (!m_visible_brown) {
    m_visible_brown = true;
    showBrown->setVisible(true);
    drawBrown->setVisible(true);
    clearBrown->setVisible(true);
    showBrownComp->setVisible(true);
    showBrownDiff->setVisible(true);
    showBrownUnion->setVisible(true);
    showBrownInt->setVisible(true);
    showBrownSym_Diff->setVisible(true);
    showBrownMink_Sum->setVisible(true);
    copyBrown->setVisible(true);
    moveBrown->setVisible(true);
    showBrownLabel->setVisible(true);

    line10->setVisible(true);

    line1->setGeometry(QRect(115, 0, 7, 185));
    line2->setGeometry(QRect(155, 0, 7, 185));
    line3->setGeometry(QRect(200, 0, 7, 185));
    line4->setGeometry(QRect(245, 0, 7, 185));
    line5->setGeometry(QRect(290, 0, 7, 185));
    line6->setGeometry(QRect(335, 0, 7, 185));
    line06->setGeometry(QRect(425, 0, 7, 185));
    line061->setGeometry(QRect(470, 0, 7, 185));
    line062->setGeometry(QRect(515, 0, 7, 185));
    line063->setGeometry(QRect(560, 0, 7, 185));
    line0007->setGeometry(QRect(380, 0, 7, 185));
    line0006->setGeometry(QRect(70, 0, 7, 185));

    actionMinusColor->setEnabled(true);
    actionMinusColor->setText("Remove Brown");
  } else if (!m_visible_yellow) {
    m_visible_yellow = true;
    showYellow->setVisible(true);
    drawYellow->setVisible(true);
    clearYellow->setVisible(true);
    showYellowComp->setVisible(true);
    showYellowDiff->setVisible(true);
    showYellowUnion->setVisible(true);
    showYellowInt->setVisible(true);
    showYellowSym_Diff->setVisible(true);
    showYellowMink_Sum->setVisible(true);
    copyYellow->setVisible(true);
    moveYellow->setVisible(true);
    showYellowLabel->setVisible(true);

    line9->setVisible(true);

    line1->setGeometry(QRect(115, 0, 7, 215));
    line2->setGeometry(QRect(155, 0, 7, 215));
    line3->setGeometry(QRect(200, 0, 7, 215));
    line4->setGeometry(QRect(245, 0, 7, 215));
    line5->setGeometry(QRect(290, 0, 7, 215));
    line6->setGeometry(QRect(335, 0, 7, 215));
    line06->setGeometry(QRect(425, 0, 7, 215));
    line061->setGeometry(QRect(470, 0, 7, 215));
    line062->setGeometry(QRect(515, 0, 7, 215));
    line063->setGeometry(QRect(560, 0, 7, 215));
    line0007->setGeometry(QRect(380, 0, 7, 215));
    line0006->setGeometry(QRect(70, 0, 7, 215));

    actionMinusColor->setText("Remove Yellow");
  } else if (!m_visible_magenta) {
    m_visible_magenta = true;
    showMagenta->setVisible(true);
    drawMagenta->setVisible(true);
    clearMagenta->setVisible(true);
    showMagentaComp->setVisible(true);
    showMagentaDiff->setVisible(true);
    showMagentaUnion->setVisible(true);
    showMagentaInt->setVisible(true);
    showMagentaSym_Diff->setVisible(true);
    showMagentaMink_Sum->setVisible(true);
    copyMagenta->setVisible(true);
    moveMagenta->setVisible(true);
    showMagentaLabel->setVisible(true);

    line8->setVisible(true);

    line1->setGeometry(QRect(115, 0, 7, 245));
    line2->setGeometry(QRect(155, 0, 7, 245));
    line3->setGeometry(QRect(200, 0, 7, 245));
    line4->setGeometry(QRect(245, 0, 7, 245));
    line5->setGeometry(QRect(290, 0, 7, 245));
    line6->setGeometry(QRect(335, 0, 7, 245));
    line06->setGeometry(QRect(425, 0, 7, 245));
    line061->setGeometry(QRect(470, 0, 7, 245));
    line062->setGeometry(QRect(515, 0, 7, 245));
    line063->setGeometry(QRect(560, 0, 7, 245));
    line0007->setGeometry(QRect(380, 0, 7, 245));
    line0006->setGeometry(QRect(70, 0, 7, 245));

    actionMinusColor->setText("Remove Magenta");
  } else if (!m_visible_aqua) {
    m_visible_aqua = true;
    showAqua->setVisible(true);
    drawAqua->setVisible(true);
    clearAqua->setVisible(true);
    showAquaComp->setVisible(true);
    showAquaDiff->setVisible(true);
    showAquaUnion->setVisible(true);
    showAquaInt->setVisible(true);
    showAquaSym_Diff->setVisible(true);
    showAquaMink_Sum->setVisible(true);
    copyAqua->setVisible(true);
    moveAqua->setVisible(true);
    showAquaLabel->setVisible(true);

    line7->setVisible(true);

    line1->setGeometry(QRect(115, 0, 7, 275));
    line2->setGeometry(QRect(155, 0, 7, 275));
    line3->setGeometry(QRect(200, 0, 7, 275));
    line4->setGeometry(QRect(245, 0, 7, 275));
    line5->setGeometry(QRect(290, 0, 7, 275));
    line6->setGeometry(QRect(335, 0, 7, 275));
    line06->setGeometry(QRect(425, 0, 7, 275));
    line061->setGeometry(QRect(470, 0, 7, 275));
    line062->setGeometry(QRect(515, 0, 7, 275));
    line063->setGeometry(QRect(560, 0, 7, 275));
    line0007->setGeometry(QRect(380, 0, 7, 275));
    line0006->setGeometry(QRect(70, 0, 7, 275));

    actionAddColor->setEnabled(false);
    actionAddColor->setText("Maximum Bucket Limit Reached");
    actionMinusColor->setText("Remove Aqua");
  }

  modelChanged();
}

//////#########Remove Color Minus
//minus color button clicked
void MainWindow::on_actionMinusColor_triggered() {
  if (!(m_visible_yellow || m_visible_magenta || m_visible_aqua)) {
    states_stack.back().m_curve_sets[3].clear();
    //states_stack.back().m_curve_sets[7].clear();
    m_color_active = 0;

    states_stack.back().brown_circular_sources().clear();
    states_stack.back().brown_linear_sources().clear();
    states_stack.back().brown_bezier_sources().clear();
    states_stack.back().brown_polyline_sources().clear();

    m_visible_brown = false;
    showBrown->setVisible(false);
    drawBrown->setVisible(false);
    clearBrown->setVisible(false);
    showBrownComp->setVisible(false);
    showBrownDiff->setVisible(false);
    showBrownUnion->setVisible(false);
    showBrownInt->setVisible(false);
    showBrownSym_Diff->setVisible(false);
    showBrownMink_Sum->setVisible(false);
    copyBrown->setVisible(false);
    moveBrown->setVisible(false);
    showBrownLabel->setVisible(false);

    drawBlue->setChecked(true);

    line10->setVisible(false);

    line1->setGeometry(QRect(115, 0, 7, 155));
    line2->setGeometry(QRect(155, 0, 7, 155));
    line3->setGeometry(QRect(200, 0, 7, 155));
    line4->setGeometry(QRect(245, 0, 7, 155));
    line5->setGeometry(QRect(290, 0, 7, 155));
    line6->setGeometry(QRect(335, 0, 7, 155));
    line06->setGeometry(QRect(425, 0, 7, 155));
    line061->setGeometry(QRect(470, 0, 7, 155));
    line062->setGeometry(QRect(515, 0, 7, 155));
    line063->setGeometry(QRect(560, 0, 7, 155));
    line0007->setGeometry(QRect(380, 0, 7, 155));
    line0006->setGeometry(QRect(70, 0, 7, 155));

    actionMinusColor->setText("Bucket Removal Not Allowed");
    actionMinusColor->setEnabled(false);
  } else if (!(m_visible_magenta || m_visible_aqua)) {
    states_stack.back().m_curve_sets[4].clear();
    m_color_active = 0;

    states_stack.back().yellow_circular_sources().clear();
    states_stack.back().yellow_linear_sources().clear();
    states_stack.back().yellow_bezier_sources().clear();
    states_stack.back().yellow_polyline_sources().clear();

    m_visible_yellow = false;
    showYellow->setVisible(false);
    drawYellow->setVisible(false);
    clearYellow->setVisible(false);
    showYellowComp->setVisible(false);
    showYellowDiff->setVisible(false);
    showYellowUnion->setVisible(false);
    showYellowInt->setVisible(false);
    showYellowSym_Diff->setVisible(false);
    showYellowMink_Sum->setVisible(false);
    copyYellow->setVisible(false);
    moveYellow->setVisible(false);
    showYellowLabel->setVisible(false);

    drawBlue->setChecked(true);

    line9->setVisible(false);

    line1->setGeometry(QRect(115, 0, 7, 185));
    line2->setGeometry(QRect(155, 0, 7, 185));
    line3->setGeometry(QRect(200, 0, 7, 185));
    line4->setGeometry(QRect(245, 0, 7, 185));
    line5->setGeometry(QRect(290, 0, 7, 185));
    line6->setGeometry(QRect(335, 0, 7, 185));
    line06->setGeometry(QRect(425, 0, 7, 185));
    line061->setGeometry(QRect(470, 0, 7, 185));
    line062->setGeometry(QRect(515, 0, 7, 185));
    line063->setGeometry(QRect(560, 0, 7, 185));
    line0007->setGeometry(QRect(380, 0, 7, 185));
    line0006->setGeometry(QRect(70, 0, 7, 185));

    actionMinusColor->setText("Remove Brown");
  } else if (!m_visible_aqua) {
    states_stack.back().m_curve_sets[5].clear();
    m_color_active = 0;

    states_stack.back().magenta_circular_sources().clear();
    states_stack.back().magenta_linear_sources().clear();
    states_stack.back().magenta_bezier_sources().clear();
    states_stack.back().magenta_polyline_sources().clear();

    m_visible_magenta = false;
    showMagenta->setVisible(false);
    drawMagenta->setVisible(false);
    clearMagenta->setVisible(false);
    showMagentaComp->setVisible(false);
    showMagentaDiff->setVisible(false);
    showMagentaUnion->setVisible(false);
    showMagentaInt->setVisible(false);
    showMagentaSym_Diff->setVisible(false);
    showMagentaMink_Sum->setVisible(false);
    copyMagenta->setVisible(false);
    moveMagenta->setVisible(false);
    showMagentaLabel->setVisible(false);

    drawBlue->setChecked(true);

    line8->setVisible(false);

    line1->setGeometry(QRect(115, 0, 7, 215));
    line2->setGeometry(QRect(155, 0, 7, 215));
    line3->setGeometry(QRect(200, 0, 7, 215));
    line4->setGeometry(QRect(245, 0, 7, 215));
    line5->setGeometry(QRect(290, 0, 7, 215));
    line6->setGeometry(QRect(335, 0, 7, 215));
    line06->setGeometry(QRect(425, 0, 7, 215));
    line061->setGeometry(QRect(470, 0, 7, 215));
    line062->setGeometry(QRect(515, 0, 7, 215));
    line063->setGeometry(QRect(560, 0, 7, 215));
    line0007->setGeometry(QRect(380, 0, 7, 215));
    line0006->setGeometry(QRect(70, 0, 7, 215));

    actionMinusColor->setText("Remove Yellow");
  } else {
    states_stack.back().m_curve_sets[6].clear();
    m_color_active = 0;

    states_stack.back().aqua_circular_sources().clear();
    states_stack.back().aqua_linear_sources().clear();
    states_stack.back().aqua_bezier_sources().clear();
    states_stack.back().aqua_polyline_sources().clear();

    m_visible_aqua = false;
    showAqua->setVisible(false);
    drawAqua->setVisible(false);
    clearAqua->setVisible(false);
    showAquaComp->setVisible(false);
    showAquaDiff->setVisible(false);
    showAquaUnion->setVisible(false);
    showAquaInt->setVisible(false);
    showAquaSym_Diff->setVisible(false);
    showAquaMink_Sum->setVisible(false);
    copyAqua->setVisible(false);
    moveAqua->setVisible(false);
    showAquaLabel->setVisible(false);

    drawBlue->setChecked(true);

    line7->setVisible(false);

    line1->setGeometry(QRect(115, 0, 7, 245));
    line2->setGeometry(QRect(155, 0, 7, 245));
    line3->setGeometry(QRect(200, 0, 7, 245));
    line4->setGeometry(QRect(245, 0, 7, 245));
    line5->setGeometry(QRect(290, 0, 7, 245));
    line6->setGeometry(QRect(335, 0, 7, 245));
    line06->setGeometry(QRect(425, 0, 7, 245));
    line061->setGeometry(QRect(470, 0, 7, 245));
    line062->setGeometry(QRect(515, 0, 7, 245));
    line063->setGeometry(QRect(560, 0, 7, 245));
    line0007->setGeometry(QRect(380, 0, 7, 245));
    line0006->setGeometry(QRect(70, 0, 7, 245));

    actionAddColor->setEnabled(true);
    actionAddColor->setText("Add a Bucket");
    actionMinusColor->setText("Remove Magenta");
  }

  modelChanged();
}

void MainWindow::on_actionUndo_triggered() {
  if (m_state_num > 1) {
    /*operation_name
    // COMPLEMENT_OP = 0
    // INTERSECTION_OP = 1
    // UNION_OP = 2
    // DIFFERENCE_OP = 3
    // SYMMETRIC_DIFFERENCE_OP = 4
    // MINKOWSKI_SUM_OP = 5
    // RESET_OP = 6
    // COPY_OP = 7
    // MOVE_OP = 8
    // CLEAR_OP = 9
    // DELETEALL_OP = 10
    // DRAWPOL_OP = 11
    START_OP = 12
    */

    //unlink GI for last state
    for (auto si = states_stack.back().m_curve_sets.begin();
         si != states_stack.back().m_curve_sets.end(); ++si) { unlink_GI(si->gi()); }

    switch (states_stack.back().m_operation) {
     case 0:if (actionComplementH->isChecked()) actionComplementH->setChecked(false);
      break;
     case 1:if (actionIntersectionH->isChecked()) actionIntersectionH->setChecked(false);
      break;
     case 2:if (actionUnionH->isChecked()) actionUnionH->setChecked(false);
      break;
     case 3:if (actionDifferenceH->isChecked()) actionDifferenceH->setChecked(false);
      break;
     case 4:if (actionSymmetric_DifferenceH->isChecked()) actionSymmetric_DifferenceH->setChecked(false);
      break;
     case 5:if (actionMinkowski_SumH->isChecked()) actionMinkowski_SumH->setChecked(false);
      break;
     case 10:if (actionComplementH->isChecked()) actionComplementH->setChecked(false);
      break;
     case 11:if (actionComplementH->isChecked()) actionComplementH->setChecked(false);
      break;
    }

    modelChanged();

    // pop last state
    states_stack.pop_back();

    for (auto si = states_stack.back().m_curve_sets.begin();
         si != states_stack.back().m_curve_sets.end(); ++si) { link_GI(si->gi()); }

    if (states_stack.size() == 1) {
      actionUndo->setEnabled(false);
    }

    m_state_num = states_stack.size();

    // for confirming the num of each state during debugging
    // show_error("State num: "+to_string(states_stack.size()));
    modelChanged();
  } else {
    show_error("Nothing to Undo!!!");
  }
}

void MainWindow::on_actionNew_triggered() {

  if (m_state_num > 1) get_new_state(6);

  for (Curve_set_iterator si = states_stack.back().m_curve_sets.begin();
       si != states_stack.back().m_curve_sets.end(); ++si)
    si->clear();

  m_scene.setSceneRect(-320, -210, 640, 420);
  this->graphicsView->setScene(&m_scene);
  this->graphicsView->fitInView(m_scene.sceneRect());
  this->graphicsView->scale(2.5, 2.5);

  states_stack.back().result_set().clear();
  states_stack.back().blue_set().clear();
  states_stack.back().black_set().clear();
  states_stack.back().brown_set().clear();
  states_stack.back().red_set().clear();
  states_stack.back().yellow_set().clear();
  states_stack.back().magenta_set().clear();
  states_stack.back().aqua_set().clear();

  states_stack.back().blue_circular_sources().clear();
  states_stack.back().red_circular_sources().clear();
  states_stack.back().black_circular_sources().clear();
  states_stack.back().brown_circular_sources().clear();
  states_stack.back().yellow_circular_sources().clear();
  states_stack.back().magenta_circular_sources().clear();
  states_stack.back().aqua_circular_sources().clear();
  states_stack.back().result_circular_sources().clear();

  states_stack.back().blue_linear_sources().clear();
  states_stack.back().red_linear_sources().clear();
  states_stack.back().black_linear_sources().clear();
  states_stack.back().brown_linear_sources().clear();
  states_stack.back().yellow_linear_sources().clear();
  states_stack.back().magenta_linear_sources().clear();
  states_stack.back().result_linear_sources().clear();

  states_stack.back().blue_bezier_sources().clear();
  states_stack.back().red_bezier_sources().clear();
  states_stack.back().black_bezier_sources().clear();
  states_stack.back().brown_bezier_sources().clear();
  states_stack.back().yellow_bezier_sources().clear();
  states_stack.back().magenta_bezier_sources().clear();
  states_stack.back().result_bezier_sources().clear();

  states_stack.back().blue_polyline_sources().clear();
  states_stack.back().red_polyline_sources().clear();
  states_stack.back().black_polyline_sources().clear();
  states_stack.back().brown_polyline_sources().clear();
  states_stack.back().yellow_polyline_sources().clear();
  states_stack.back().magenta_polyline_sources().clear();
  states_stack.back().result_polyline_sources().clear();

  SetViewBlue(true);
  SetViewRed(true);
  SetViewBlack(true);
  SetViewBrown(true);
  SetViewYellow(true);
  SetViewMagenta(true);
  SetViewAqua(true);
  SetViewResult(true);

  actionAddColor->setEnabled(true);
  actionMinusColor->setEnabled(false);

  m_disjoint = false;

  m_linear_input->Reset();
  m_circular_input->Reset();
  m_bezier_input->Reset();
  m_polyline_input->Reset();

  m_pan = false;

  m_color_active = 0;
  m_color_move = 1111;
  m_color_cm = 1111;
  m_color_visible = 7;

  m_linear_input->mOngoingPieceGI->setPen(sPens[0]);
  m_linear_input->mLinearGI->setPen(sPens[0]);
  m_linear_input->mHandleGI->setPen(sPens[0]);

  sceneDockWidget->setVisible(true);
  consoleDockWidget->setVisible(true);
  infoDockWidget->setVisible(true);

  showColorBucket->setChecked(true);
  showConsole->setChecked(true);
  showInfo->setChecked(true);

  actionCopyH->setChecked(false);
  actionMoveH->setChecked(false);
  //actionAxis->setChecked(false);

  actionComplementH->setChecked(false);
  actionUnionH->setChecked(false);
  actionIntersectionH->setChecked(false);
  actionDifferenceH->setChecked(false);
  actionSymmetric_DifferenceH->setChecked(false);
  actionMinkowski_SumH->setChecked(false);

  m_reset = true;
  if (!m_circular_active && !m_bezier_active && !m_polyline_active) {
    actionInsertLinear->setChecked(true);
    actionInsertCircular->setChecked(false);
    actionInsertBezier->setChecked(false);
    actionInsertPolyline->setChecked(false);
  } else if (!m_bezier_active && !m_polyline_active) {
    actionInsertLinear->setChecked(false);
    actionInsertCircular->setChecked(true);
    actionInsertBezier->setChecked(false);
    actionInsertPolyline->setChecked(false);
  } else if (!m_polyline_active) {
    actionInsertLinear->setChecked(false);
    actionInsertCircular->setChecked(false);
    actionInsertBezier->setChecked(true);
    actionInsertPolyline->setChecked(false);
  } else {
    actionInsertLinear->setChecked(false);
    actionInsertCircular->setChecked(false);
    actionInsertBezier->setChecked(false);
    actionInsertPolyline->setChecked(true);
  }

  m_reset = false;

  // actionMinkowski_SumH -> setEnabled(true);
  actionDifferenceH->setEnabled(true);

  m_color_copy = 111; //default
  m_color_move_ext = 111; //default
  m_color_complement = 0; //default
  m_blue_int = true; //default
  m_red_int = true; //default
  m_black_int = false; //default
  m_brown_int = false; //default
  m_yellow_int = false; //default
  m_magenta_int = false; //default
  m_aqua_int = false; //default

  m_blue_union = true; //default
  m_red_union = true; //default
  m_black_union = false; //default
  m_brown_union = false; //default
  m_yellow_union = false; //default
  m_magenta_union = false; //default
  m_aqua_union = false; //default

  m_blue_sym_diff = true; //default
  m_red_sym_diff = true; //default
  m_black_sym_diff = false; //default
  m_brown_sym_diff = false; //default
  m_yellow_sym_diff = false; //default
  m_magenta_sym_diff = false; //default
  m_aqua_sym_diff = false; //default

  m_blue_mink = true; //default
  m_red_mink = true; //default
  m_black_mink = false; //default
  m_brown_mink = false; //default
  m_yellow_mink = false; //default
  m_magenta_mink = false; //default
  m_aqua_mink = false; //default

  drawBlue->setChecked(true);

  showBlue->setChecked(true);
  showRed->setChecked(true);
  showBlack->setChecked(true);
  showBrown->setChecked(true);
  showYellow->setChecked(true);
  showMagenta->setChecked(true);
  showAqua->setChecked(true);
  showResult->setChecked(false);

  clearBlue->setChecked(false);
  clearRed->setChecked(false);
  clearBlack->setChecked(false);
  clearBrown->setChecked(false);
  clearYellow->setChecked(false);
  clearMagenta->setChecked(false);
  clearAqua->setChecked(false);

  clearBrown->setVisible(false);
  clearYellow->setVisible(false);
  clearMagenta->setVisible(false);
  clearAqua->setVisible(false);

  showBrownLabel->setVisible(false);
  showYellowLabel->setVisible(false);
  showMagentaLabel->setVisible(false);
  showAquaLabel->setVisible(false);

  showBlueComp->setChecked(true);
  showRedComp->setChecked(false);
  showBlackComp->setChecked(false);
  showBrownComp->setChecked(false);
  showYellowComp->setChecked(false);
  showMagentaComp->setChecked(false);
  showAquaComp->setChecked(false);

  copyBlue->setChecked(false);
  copyRed->setChecked(false);
  copyBlack->setChecked(false);
  copyBrown->setChecked(false);
  copyYellow->setChecked(false);
  copyMagenta->setChecked(false);
  copyAqua->setChecked(false);

  moveBlue->setChecked(false);
  moveRed->setChecked(false);
  moveBlack->setChecked(false);
  moveBrown->setChecked(false);
  moveYellow->setChecked(false);
  moveMagenta->setChecked(false);
  moveAqua->setChecked(false);

  showBlackInt->setChecked(false);
  showBrownInt->setChecked(false);
  showYellowInt->setChecked(false);
  showMagentaInt->setChecked(false);
  showAquaInt->setChecked(false);
  showBlueInt->setChecked(true);
  showRedInt->setChecked(true);

  showBlackUnion->setChecked(false);
  showBrownUnion->setChecked(false);
  showYellowUnion->setChecked(false);
  showMagentaUnion->setChecked(false);
  showAquaUnion->setChecked(false);
  showBlueUnion->setChecked(true);
  showRedUnion->setChecked(true);

  showBlackDiff->setChecked(false);
  showBrownDiff->setChecked(false);
  showYellowDiff->setChecked(false);
  showMagentaDiff->setChecked(false);
  showAquaDiff->setChecked(false);
  showBlueDiff->setChecked(true);
  showRedDiff->setChecked(true);

  showBlackSym_Diff->setChecked(false);
  showBrownSym_Diff->setChecked(false);
  showYellowSym_Diff->setChecked(false);
  showMagentaSym_Diff->setChecked(false);
  showAquaSym_Diff->setChecked(false);
  showBlueSym_Diff->setChecked(true);
  showRedSym_Diff->setChecked(true);

  showBlackMink_Sum->setChecked(false);
  showBrownMink_Sum->setChecked(false);
  showYellowMink_Sum->setChecked(false);
  showMagentaMink_Sum->setChecked(false);
  showAquaMink_Sum->setChecked(false);
  showBlueMink_Sum->setChecked(true);
  showRedMink_Sum->setChecked(true);

  m_visible_black = true;
  showBlack->setVisible(true);
  drawBlack->setVisible(true);
  showBlackComp->setVisible(true);
  showBlackDiff->setVisible(true);
  showBlackUnion->setVisible(true);
  showBlackInt->setVisible(true);
  showBlackSym_Diff->setVisible(true);
  showBlackMink_Sum->setVisible(true);
  copyBlack->setVisible(true);
  moveBlack->setVisible(true);

  m_visible_brown = false;
  showBrown->setVisible(false);
  drawBrown->setVisible(false);
  showBrownComp->setVisible(false);
  showBrownDiff->setVisible(false);
  showBrownUnion->setVisible(false);
  showBrownInt->setVisible(false);
  showBrownSym_Diff->setVisible(false);
  showBrownMink_Sum->setVisible(false);
  copyBrown->setVisible(false);
  moveBrown->setVisible(false);

  m_visible_yellow = false;
  showYellow->setVisible(false);
  drawYellow->setVisible(false);
  showYellowComp->setVisible(false);
  showYellowDiff->setVisible(false);
  showYellowUnion->setVisible(false);
  showYellowInt->setVisible(false);
  showYellowSym_Diff->setVisible(false);
  showYellowMink_Sum->setVisible(false);
  copyYellow->setVisible(false);
  moveYellow->setVisible(false);

  m_visible_magenta = false;
  showMagenta->setVisible(false);
  drawMagenta->setVisible(false);
  showMagentaComp->setVisible(false);
  showMagentaDiff->setVisible(false);
  showMagentaUnion->setVisible(false);
  showMagentaInt->setVisible(false);
  showMagentaSym_Diff->setVisible(false);
  showMagentaMink_Sum->setVisible(false);
  copyMagenta->setVisible(false);
  moveMagenta->setVisible(false);

  m_visible_aqua = false;
  showAqua->setVisible(false);
  drawAqua->setVisible(false);
  showAquaComp->setVisible(false);
  showAquaDiff->setVisible(false);
  showAquaUnion->setVisible(false);
  showAquaInt->setVisible(false);
  showAquaSym_Diff->setVisible(false);
  showAquaMink_Sum->setVisible(false);
  copyAqua->setVisible(false);
  moveAqua->setVisible(false);

  actionAddColor->setEnabled(true);
  actionAddColor->setText("Add a Bucket");
  actionMinusColor->setText("Bucket Removal Not Allowed");
  VisibleHeader->setChecked(true);
  showClear->setChecked(false);

  line7->setVisible(false);
  line8->setVisible(false);
  line9->setVisible(false);
  line10->setVisible(false);
  line11->setVisible(true);

  line1->setGeometry(QRect(115, 0, 7, 155));
  line2->setGeometry(QRect(155, 0, 7, 155));
  line3->setGeometry(QRect(200, 0, 7, 155));
  line4->setGeometry(QRect(245, 0, 7, 155));
  line5->setGeometry(QRect(290, 0, 7, 155));
  line6->setGeometry(QRect(335, 0, 7, 155));
  line06->setGeometry(QRect(425, 0, 7, 155));
  line061->setGeometry(QRect(470, 0, 7, 155));
  line062->setGeometry(QRect(515, 0, 7, 155));
  line063->setGeometry(QRect(560, 0, 7, 155));
  line0007->setGeometry(QRect(380, 0, 7, 155));
  line0006->setGeometry(QRect(70, 0, 7, 155));

  zoomToFit();
  modelChanged();
}

void MainWindow::on_actionDeleteAll_triggered() {
  bool lDone = false;
  //bool lProceed = states_stack.back().result_set().is_empty() ? true : ask_user_yesno("Store result","All polygons will be deleted\n continue anyway?\n");

  //ask: what is the use of if(true)
  if (true) {
    actionComplementH->setChecked(false);
    actionIntersectionH->setChecked(false);
    actionUnionH->setChecked(false);
    actionDifferenceH->setChecked(false);
    actionSymmetric_DifferenceH->setChecked(false);
    actionMinkowski_SumH->setChecked(false);
    actionComplementH->setChecked(false);
    actionComplementH->setChecked(false);

    get_new_state(10);

    /*operation_name
    // COMPLEMENT_OP = 0
    // INTERSECTION_OP = 1
    // UNION_OP = 2
    // DIFFERENCE_OP = 3
    // SYMMETRIC_DIFFERENCE_OP = 4
    // MINKOWSKI_SUM_OP = 5
    // RESET_OP = 6
    // COPY_OP = 7
    // MOVE_OP = 8
    // CLEAR_OP = 9
    // DELETEALL_OP = 10
    // START_OP = 11*/

    states_stack.back().blue_set().clear();
    states_stack.back().blue_circular_sources().clear();
    states_stack.back().blue_bezier_sources().clear();
    states_stack.back().blue_linear_sources().clear();
    states_stack.back().blue_polyline_sources().clear();
    states_stack.back().red_set().clear();
    states_stack.back().red_circular_sources().clear();
    states_stack.back().red_bezier_sources().clear();
    states_stack.back().red_linear_sources().clear();
    states_stack.back().red_polyline_sources().clear();
    states_stack.back().black_set().clear();
    states_stack.back().black_circular_sources().clear();
    states_stack.back().black_bezier_sources().clear();
    states_stack.back().black_linear_sources().clear();
    states_stack.back().black_polyline_sources().clear();
    states_stack.back().brown_set().clear();
    states_stack.back().brown_circular_sources().clear();
    states_stack.back().brown_bezier_sources().clear();
    states_stack.back().brown_linear_sources().clear();
    states_stack.back().brown_polyline_sources().clear();
    states_stack.back().yellow_set().clear();
    states_stack.back().yellow_circular_sources().clear();
    states_stack.back().yellow_bezier_sources().clear();
    states_stack.back().brown_linear_sources().clear();
    states_stack.back().brown_polyline_sources().clear();
    states_stack.back().magenta_set().clear();
    states_stack.back().magenta_circular_sources().clear();
    states_stack.back().magenta_bezier_sources().clear();
    states_stack.back().magenta_linear_sources().clear();
    states_stack.back().magenta_polyline_sources().clear();
    states_stack.back().aqua_set().clear();
    states_stack.back().aqua_circular_sources().clear();
    states_stack.back().aqua_bezier_sources().clear();
    states_stack.back().aqua_linear_sources().clear();
    states_stack.back().aqua_polyline_sources().clear();

    on_actionDeleteResult();
  }
  lDone = true;
  if (lDone) modelChanged();
}

void MainWindow::on_actionDeleteResult() {
  bool lDone = false;
  bool lProceed;

  /*if(states_stack.back().result_set().is_empty()||states_stack.back().result_set().is_empty()||states_stack.back().result_set().is_empty()||states_stack.back().result_set().is_empty()||states_stack.back().result_set().is_empty()||states_stack.back().result_set().is_empty()||states_stack.back().result_set().is_empty())
    {
    bool lProceed = ask_user_yesno("Store result","Result will be deleted\n continue anyway?\n");
    }*/

  if (true) {
    states_stack.back().result_set().clear();
    states_stack.back().result_circular_sources().clear();
    states_stack.back().result_bezier_sources().clear();
    states_stack.back().result_linear_sources().clear();
    states_stack.back().result_polyline_sources().clear();
    lDone = true;
  }
  if (lDone) modelChanged();
}

void MainWindow::on_actionClearH_toggled(bool aChecked) {
  if (aChecked) {
    bool lDone = false;
    bool lProceed;

    actionComplementH->setChecked(false);
    actionIntersectionH->setChecked(false);
    actionUnionH->setChecked(false);
    actionDifferenceH->setChecked(false);
    actionSymmetric_DifferenceH->setChecked(false);
    actionMinkowski_SumH->setChecked(false);
    actionComplementH->setChecked(false);
    actionComplementH->setChecked(false);

    get_new_state(9);

    /*operation_name
    // COMPLEMENT_OP = 0
    // INTERSECTION_OP = 1
    // UNION_OP = 2
    // DIFFERENCE_OP = 3
    // SYMMETRIC_dIFFERENCE_OP = 4
    // MINKOWSKI_SUM_OP = 5
    // RESET_OP = 6
    // COPY_OP = 7
    // MOVE_OP = 8
    // CLEAR_OP = 9
    // DELETEALL_OP = 10
    // START_OP = 11*/

    if (clearBlue->isChecked()) {
      states_stack.back().blue_set().clear();
      states_stack.back().blue_linear_sources().clear();
      states_stack.back().blue_circular_sources().clear();
      states_stack.back().blue_bezier_sources().clear();
      states_stack.back().blue_polyline_sources().clear();
    }
    if (clearBlack->isChecked()) {
      states_stack.back().black_set().clear();
      states_stack.back().black_linear_sources().clear();
      states_stack.back().black_circular_sources().clear();
      states_stack.back().black_bezier_sources().clear();
      states_stack.back().black_polyline_sources().clear();
    }
    if (clearRed->isChecked()) {
      states_stack.back().red_set().clear();
      states_stack.back().red_linear_sources().clear();
      states_stack.back().red_circular_sources().clear();
      states_stack.back().red_bezier_sources().clear();
      states_stack.back().red_polyline_sources().clear();
    }
    if (clearBrown->isChecked()) {
      states_stack.back().brown_set().clear();
      states_stack.back().brown_linear_sources().clear();
      states_stack.back().brown_circular_sources().clear();
      states_stack.back().brown_bezier_sources().clear();
      states_stack.back().brown_polyline_sources().clear();
    }
    if (clearYellow->isChecked()) {
      states_stack.back().yellow_set().clear();
      states_stack.back().yellow_linear_sources().clear();
      states_stack.back().yellow_circular_sources().clear();
      states_stack.back().yellow_bezier_sources().clear();
      states_stack.back().yellow_polyline_sources().clear();
    }
    if (clearMagenta->isChecked()) {
      states_stack.back().magenta_set().clear();
      states_stack.back().magenta_linear_sources().clear();
      states_stack.back().magenta_circular_sources().clear();
      states_stack.back().magenta_bezier_sources().clear();
      states_stack.back().magenta_polyline_sources().clear();
    }
    if (clearAqua->isChecked()) {
      states_stack.back().aqua_set().clear();
      states_stack.back().aqua_linear_sources().clear();
      states_stack.back().aqua_circular_sources().clear();
      states_stack.back().aqua_bezier_sources().clear();
      states_stack.back().aqua_polyline_sources().clear();
    }

    if (!(clearBlue->isChecked()) && !(clearRed->isChecked()) && !(clearBlack->isChecked()) &&
        !(clearBrown->isChecked()) && !(clearYellow->isChecked()) && !(clearMagenta->isChecked()) &&
        !(clearAqua->isChecked()))
      show_error("Please select color bucket to be cleaned!!!");

    actionClearH->setChecked(false);
    lDone = true;

    if (lDone) modelChanged();
  }
}

//draw polygons
//todo::redefine it's parameters i.e remove bool from params
void MainWindow::on_drawBlue_toggled(bool /* a_check */) {
  m_color_active = 0;

  if (states_stack.back().blue_set().is_linear()) {
    m_linear_input->mOngoingPieceGI->setPen(sPens[0]);
    m_linear_input->mLinearGI->setPen(sPens[0]);
    m_linear_input->mHandleGI->setPen(sPens[0]);
  } else if (states_stack.back().blue_set().is_circular()) {
    m_circular_input->mOngoingPieceGI->setPen(sPens[0]);
    m_circular_input->mCircularGI->setPen(sPens[0]);
    m_circular_input->mHandleGI->setPen(sPens[0]);
  } else if (states_stack.back().blue_set().is_bezier()) {
    m_bezier_input->mOngoingPieceGI->setPen(sPens[0]);
    m_bezier_input->mBezierGI->setPen(sPens[0]);
    m_bezier_input->mHandle0GI->setPen(sPens[0]);
    m_bezier_input->mHandle1GI->setPen(sPens[0]);
  } else {
    m_polyline_input->mOngoingPieceGI->setPen(sPens[0]);
    m_polyline_input->mPolylineGI->setPen(sPens[0]);
    m_polyline_input->mHandleGI->setPen(sPens[0]);
  }

}

void MainWindow::on_drawRed_toggled(bool /* a_check */) {
  m_color_active = 1;

  if (states_stack.back().red_set().is_linear()) {
    m_linear_input->mOngoingPieceGI->setPen(sPens[1]);
    m_linear_input->mLinearGI->setPen(sPens[1]);
    m_linear_input->mHandleGI->setPen(sPens[1]);
  } else if (states_stack.back().red_set().is_circular()) {
    m_circular_input->mOngoingPieceGI->setPen(sPens[1]);
    m_circular_input->mCircularGI->setPen(sPens[1]);
    m_circular_input->mHandleGI->setPen(sPens[1]);
  } else if (states_stack.back().red_set().is_bezier()) {
    m_bezier_input->mOngoingPieceGI->setPen(sPens[1]);
    m_bezier_input->mBezierGI->setPen(sPens[1]);
    m_bezier_input->mHandle0GI->setPen(sPens[1]);
    m_bezier_input->mHandle1GI->setPen(sPens[1]);
  } else {
    m_polyline_input->mOngoingPieceGI->setPen(sPens[1]);
    m_polyline_input->mPolylineGI->setPen(sPens[1]);
    m_polyline_input->mHandleGI->setPen(sPens[1]);
  }

}

void MainWindow::on_drawBlack_toggled(bool /* a_check */) {
  m_color_active = 2;

  if (states_stack.back().black_set().is_linear()) {
    m_linear_input->mOngoingPieceGI->setPen(sPens[2]);
    m_linear_input->mLinearGI->setPen(sPens[2]);
    m_linear_input->mHandleGI->setPen(sPens[2]);
  } else if (states_stack.back().black_set().is_circular()) {
    m_circular_input->mOngoingPieceGI->setPen(sPens[2]);
    m_circular_input->mCircularGI->setPen(sPens[2]);
    m_circular_input->mHandleGI->setPen(sPens[2]);
  } else if (states_stack.back().black_set().is_bezier()) {
    m_bezier_input->mOngoingPieceGI->setPen(sPens[2]);
    m_bezier_input->mBezierGI->setPen(sPens[2]);
    m_bezier_input->mHandle0GI->setPen(sPens[2]);
    m_bezier_input->mHandle1GI->setPen(sPens[2]);
  } else {
    m_polyline_input->mOngoingPieceGI->setPen(sPens[2]);
    m_polyline_input->mPolylineGI->setPen(sPens[2]);
    m_polyline_input->mHandleGI->setPen(sPens[2]);
  }

}

void MainWindow::on_drawBrown_toggled(bool /* a_check */) {
  m_color_active = 3;

  if (states_stack.back().brown_set().is_linear()) {
    m_linear_input->mOngoingPieceGI->setPen(sPens[3]);
    m_linear_input->mLinearGI->setPen(sPens[3]);
    m_linear_input->mHandleGI->setPen(sPens[3]);
  } else if (states_stack.back().brown_set().is_circular()) {
    m_circular_input->mOngoingPieceGI->setPen(sPens[3]);
    m_circular_input->mCircularGI->setPen(sPens[3]);
    m_circular_input->mHandleGI->setPen(sPens[3]);
  } else if (states_stack.back().brown_set().is_bezier()) {
    m_bezier_input->mOngoingPieceGI->setPen(sPens[3]);
    m_bezier_input->mBezierGI->setPen(sPens[3]);
    m_bezier_input->mHandle0GI->setPen(sPens[3]);
    m_bezier_input->mHandle1GI->setPen(sPens[3]);
  } else {
    m_polyline_input->mOngoingPieceGI->setPen(sPens[3]);
    m_polyline_input->mPolylineGI->setPen(sPens[3]);
    m_polyline_input->mHandleGI->setPen(sPens[3]);
  }

}

void MainWindow::on_drawYellow_toggled(bool /* a_check */) {
  m_color_active = 4;

  if (states_stack.back().yellow_set().is_linear()) {
    m_linear_input->mOngoingPieceGI->setPen(sPens[4]);
    m_linear_input->mLinearGI->setPen(sPens[4]);
    m_linear_input->mHandleGI->setPen(sPens[4]);
  } else if (states_stack.back().yellow_set().is_circular()) {
    m_circular_input->mOngoingPieceGI->setPen(sPens[4]);
    m_circular_input->mCircularGI->setPen(sPens[4]);
    m_circular_input->mHandleGI->setPen(sPens[4]);
  } else if (states_stack.back().yellow_set().is_bezier()) {
    m_bezier_input->mOngoingPieceGI->setPen(sPens[4]);
    m_bezier_input->mBezierGI->setPen(sPens[4]);
    m_bezier_input->mHandle0GI->setPen(sPens[4]);
    m_bezier_input->mHandle1GI->setPen(sPens[4]);
  } else {
    m_polyline_input->mOngoingPieceGI->setPen(sPens[4]);
    m_polyline_input->mPolylineGI->setPen(sPens[4]);
    m_polyline_input->mHandleGI->setPen(sPens[4]);
  }

}

void MainWindow::on_drawMagenta_toggled(bool /* a_check */) {
  m_color_active = 5;

  if (states_stack.back().magenta_set().is_linear()) {
    m_linear_input->mOngoingPieceGI->setPen(sPens[5]);
    m_linear_input->mLinearGI->setPen(sPens[5]);
    m_linear_input->mHandleGI->setPen(sPens[5]);
  } else if (states_stack.back().magenta_set().is_circular()) {
    m_circular_input->mOngoingPieceGI->setPen(sPens[5]);
    m_circular_input->mCircularGI->setPen(sPens[5]);
    m_circular_input->mHandleGI->setPen(sPens[5]);
  } else if (states_stack.back().magenta_set().is_bezier()) {
    m_bezier_input->mOngoingPieceGI->setPen(sPens[5]);
    m_bezier_input->mBezierGI->setPen(sPens[5]);
    m_bezier_input->mHandle0GI->setPen(sPens[5]);
    m_bezier_input->mHandle1GI->setPen(sPens[5]);
  } else {
    m_polyline_input->mOngoingPieceGI->setPen(sPens[5]);
    m_polyline_input->mPolylineGI->setPen(sPens[5]);
    m_polyline_input->mHandleGI->setPen(sPens[5]);
  }

}

void MainWindow::on_drawAqua_toggled(bool /* a_check */) {
  m_color_active = 6;

  if (states_stack.back().aqua_set().is_linear()) {
    m_linear_input->mOngoingPieceGI->setPen(sPens[6]);
    m_linear_input->mLinearGI->setPen(sPens[6]);
    m_linear_input->mHandleGI->setPen(sPens[6]);
  } else if (states_stack.back().aqua_set().is_circular()) {
    m_circular_input->mOngoingPieceGI->setPen(sPens[6]);
    m_circular_input->mCircularGI->setPen(sPens[6]);
    m_circular_input->mHandleGI->setPen(sPens[6]);
  } else if (states_stack.back().aqua_set().is_bezier()) {
    m_bezier_input->mOngoingPieceGI->setPen(sPens[6]);
    m_bezier_input->mBezierGI->setPen(sPens[6]);
    m_bezier_input->mHandle0GI->setPen(sPens[6]);
    m_bezier_input->mHandle1GI->setPen(sPens[6]);
  } else {
    m_polyline_input->mOngoingPieceGI->setPen(sPens[6]);
    m_polyline_input->mPolylineGI->setPen(sPens[6]);
    m_polyline_input->mHandleGI->setPen(sPens[6]);
  }

}

//extra utilities
void MainWindow::on_actionRecenter_triggered() { zoomToFit(); }

void MainWindow::on_actionOpenLinear_triggered() {
  open(QFileDialog::getOpenFileName(this,
                                    tr("Open Linear Polygon"), "./data",
                                    tr("Linear Curve files (*.lps)")));
}

void MainWindow::on_actionOpenDXF_triggered() {
  open(QFileDialog::getOpenFileName(this,
                                    tr("Open Line segments and circular arcs"), "./data",
                                    tr("Circular curve files (*.cps)")));
}

void MainWindow::on_actionOpenBezier_triggered() {
  open(QFileDialog::getOpenFileName(this,
                                    tr("Open Bezier Polygon"), "./data",
                                    tr("Bezier Curve files (*.bps)")));
}

void MainWindow::on_actionOpenPolyline_triggered() {
  open(QFileDialog::getOpenFileName(this,
                                    tr("Open Polyline Polygon"), "./data",
                                    tr("Polyline Curve files (*.mps)")));
}

bool MainWindow::read_linear(QString aFileName, Linear_polygon_set& rSet,
                             Linear_region_source_container& rSources) {
  bool rOK = false;

  std::ifstream in_file(qPrintable(aFileName));
  if (in_file) {
    try {
      typedef Linear_traits Gps_traits;
      typedef typename Gps_traits::X_monotone_curve_2 Linear_X_monotone_curve;
      typedef typename Gps_traits::Curve_2 Linear_curve;
      typedef std::vector <Linear_curve> Linear_curve_vector;
      typedef typename Gps_traits::X_monotone_curve_2 Linear_X_monotone_curve;
      typedef typename Gps_traits::Polygon_2 Linear_polygon;
      typedef typename Gps_traits::Point_2 Linear_point;
      typedef typename Kernel::Point_2 Point;
      typedef typename Kernel::FT FT;
      Linear_curve_vector mLinearPolygonPieces;

      std::string format;
      std::getline(in_file, format);

      bool lDoubleFormat = (format.length() >= 6 && format.substr(0, 6) == "DOUBLE");

      //number of linear polygon with holes
      unsigned int n_regions;
      in_file >> n_regions;

      for (unsigned int r = 0; r < n_regions; ++r) {
        //number of linear polygon including holes
        unsigned int n_boundaries;
        in_file >> n_boundaries;

        std::vector <Linear_polygon> Linear_polygons;

        for (unsigned int r = 0; r < n_boundaries; ++r) {
          // number of points in linear polygon
          unsigned int n_points;
          in_file >> n_points;
          std::vector <Point> points;

          std::vector <Linear_X_monotone_curve> xcvs;
          Gps_traits traits;
          auto make_x_monotone = traits.make_x_monotone_2_object();
          using Pnt = typename Gps_traits::Point_2;
          using Xcv = typename Gps_traits::X_monotone_curve_2;
          using Make_x_monotone_result = std::variant<Pnt, Xcv>;

          for (unsigned int j = 0; j < n_points; ++j) {
            double x, y;
            in_file >> x >> y;
            points.push_back(Point(x, y));
            if (j > 0) {
              mLinearPolygonPieces.push_back(Linear_curve(points[points.size() - 2],
                                                          points[points.size() - 1]));
            }
          }
          points.clear();

          for (auto it = mLinearPolygonPieces.begin();
               it != mLinearPolygonPieces.end(); ++it) {
            std::vector<Make_x_monotone_result> x_objs;

            make_x_monotone(*it, std::back_inserter(x_objs));
            for (auto i = 0; i < x_objs.size(); ++i) {
              auto *xcv = std::get_if<Linear_X_monotone_curve>(&x_objs[i]);
              CGAL_assertion(xcv != nullptr);
              xcvs.push_back(*xcv);
            }
          }

          if (xcvs.size() > 0) {
            Linear_point const &first_point = xcvs.front().source();
            Linear_point const &last_point = xcvs.back().target();
            FT fxs = first_point.x();
            FT fys = first_point.y();
            FT lxs = last_point.x();
            FT lys = last_point.y();
            xcvs.push_back(Linear_X_monotone_curve(Point(lxs, lys), Point(fxs, fys)));
            //for (auto &i:xcvs)
            //cout<<i<<endl;

            Linear_polygon lp(xcvs.begin(), xcvs.end());
            CGAL::Orientation orient = lp.orientation();
            if (orient == CGAL::CLOCKWISE) {
              lp.reverse_orientation();
            }
            Linear_polygon_with_holes lCPWH(lp);

            if (r == 0) {
              get_new_state(11);
              states_stack.back().active_set(m_color_active).linear().join(lCPWH);
              states_stack.back().active_linear_sources(m_color_active).push_back(lCPWH);
            }
            else {
              states_stack.back().result_set().clear();
              states_stack.back().result_linear_sources().clear();
              states_stack.back().result_set().linear().join(lCPWH);
              states_stack.back().result_linear_sources().push_back(lCPWH);
              states_stack.back().active_set(m_color_active).difference(states_stack.back().result_set());
              states_stack.back().result_set().clear();
              states_stack.back().result_linear_sources().clear();
            }
            mLinearPolygonPieces.clear();
          }
        }
        rOK = true;
      }
    }
    catch (...) {
      show_warning("Exception ocurred during reading of linear polygon set.");
    }
  }
  return rOK;
}

bool MainWindow::read_circular(QString aFileName, Circular_polygon_set& rSet,
                               Circular_region_source_container& rSources) {
  bool rOK = false;

  std::ifstream in_file(qPrintable(aFileName));
  if (in_file) {
    try {
      typedef CGAL::Gps_circle_segment_traits_2 <Kernel> Gps_traits;
      typedef typename Gps_traits::Curve_2 Circular_curve;
      typedef typename Gps_traits::X_monotone_curve_2 Circular_X_monotone_curve;
      typedef typename Gps_traits::Polygon_2 Circular_polygon;
      typedef typename Circular_polygon::Point_2 Arc_point;
      typedef typename Kernel::FT FT;
      typedef typename Kernel::Vector_2 Vector;
      typedef typename Kernel::Point_2 Point;
      typedef typename Kernel::Circle_2 Circle;

      typedef std::vector <Circular_curve> Circular_curve_vector;

      Circular_curve_vector mCircularPolygonPieces;

      std::string format;
      std::getline(in_file, format);

      bool lDoubleFormat = (format.length() >= 6 && format.substr(0, 6) == "DOUBLE");

      //number of linear polygon with holes
      unsigned int n_regions;
      in_file >> n_regions;

      for (unsigned int r = 0; r < n_regions; ++r) {
        //number of linear polygon including holes
        unsigned int n_boundaries;
        in_file >> n_boundaries;

        std::vector <Circular_polygon> Circular_polygons;

        for (unsigned int r = 0; r < n_boundaries; ++r) {
          // number of points in linear polygon
          unsigned int n_points;
          in_file >> n_points;
          std::vector <Point> points;

          std::vector <Circular_X_monotone_curve> xcvs;
          Gps_traits traits;
          auto make_x_monotone = traits.make_x_monotone_2_object();
          using Pnt = typename Gps_traits::Point_2;
          using Xcv = typename Gps_traits::X_monotone_curve_2;
          using Make_x_monotone_result = std::variant<Pnt, Xcv>;

          for (unsigned int j = 0; j < n_points; ++j) {
            unsigned int bulge;
            in_file >> bulge;
            if (bulge) {
              double sx, sy, mx, my, tx, ty;

              in_file >> sx >> sy;
              in_file >> mx >> my;
              in_file >> tx >> ty;

              mCircularPolygonPieces.push_back(Circular_curve(Point(sx, sy),
                                                              Point(mx, my),
                                                              Point(tx, ty)));
            }
            else  //3576.382857 3576.382064
            {
              double sx, sy, tx, ty;

              in_file >> sx >> sy;
              in_file >> tx >> ty;

              mCircularPolygonPieces.push_back(Circular_curve(Point(sx, sy), Point(tx, ty)));
            }
          }

          if (mCircularPolygonPieces.size() > 0) {
            Gps_traits traits;
            auto make_x_monotone = traits.make_x_monotone_2_object();;

            std::vector <Circular_X_monotone_curve> xcvs;
            for (auto it = mCircularPolygonPieces.begin(); it != mCircularPolygonPieces.end(); ++it) {
              std::vector<Make_x_monotone_result> x_objs;
              make_x_monotone(*it, std::back_inserter(x_objs));

              auto *xcv = std::get_if<Circular_X_monotone_curve>(&x_objs[0]);
              //if ((*xcv).is_linear() && mCircularPolygonPieces.size() == 1) return;
              for (auto i = 0; i < x_objs.size(); ++i) {
                auto *xcv = std::get_if<Circular_X_monotone_curve>(&x_objs[i]);
                CGAL_assertion(xcv != nullptr);
                xcvs.push_back(*xcv);
              }
            }
            if (xcvs.size() > 0) {
              Circular_polygon cp(xcvs.begin(), xcvs.end());

              CGAL::Orientation orient = cp.orientation();

              if (orient == CGAL::CLOCKWISE) {
                cp.reverse_orientation();
              }

              Circular_polygon_with_holes lCPWH(cp);

              if (r == 0) {
                // processInput(CGAL::make_object(cp));
                get_new_state(11);
                states_stack.back().active_set(m_color_active).circular().join(lCPWH);
                states_stack.back().active_circular_sources(m_color_active).push_back(lCPWH);
              }
              else {
                states_stack.back().result_set().clear();
                states_stack.back().result_circular_sources().clear();
                states_stack.back().result_set().circular().join(lCPWH);
                states_stack.back().result_circular_sources().push_back(lCPWH);
                states_stack.back().active_set(m_color_active).difference(states_stack.back().result_set());
                states_stack.back().result_set().clear();
                states_stack.back().result_circular_sources().clear();
              }
            }
            mCircularPolygonPieces.clear();
          }
        }
        rOK = true;
      }
    }
    catch (...) {
      show_warning("Exception occured during reading of circular polygon set.");
    }
  }
  return rOK;
}

bool MainWindow::read_polyline(QString aFileName, Polyline_polygon_set& rSet,
                               Polyline_region_source_container& rSources) {
  bool rOK = false;

  std::ifstream in_file(qPrintable(aFileName));
  if (in_file) {
    try {
      typedef Polyline_traits Gps_traits;
      typedef typename Gps_traits::X_monotone_curve_2 Polyline_X_monotone_curve;
      typedef typename Gps_traits::Curve_2 Polyline_curve;
      typedef typename Gps_traits::X_monotone_curve_2 Polyline_X_monotone_curve;
      typedef typename Gps_traits::Polygon_2 Polyline_polygon;
      typedef typename Gps_traits::Point_2 Polyline_point;
      typedef typename Kernel::Point_2 Point;
      typedef typename Kernel::FT FT;
      typedef std::vector <Polyline_curve> Polyline_curve_vector;

      typedef CGAL::Arr_segment_traits_2 <Kernel> Segment_traits_2;
      typedef CGAL::Gps_traits_2 <CGAL::Arr_polycurve_traits_2<Segment_traits_2>> Base;
      typedef typename Base::X_monotone_subcurve_2 X_monotone_subcurve_2;

      Polyline_curve_vector mPolylinePolygonPieces;

      std::string format;
      std::getline(in_file, format);

      bool lDoubleFormat = (format.length() >= 6 && format.substr(0, 6) == "DOUBLE");

      unsigned int n_regions;
      in_file >> n_regions;

      for (unsigned int r = 0; r < n_regions; ++r) {

        unsigned int n_boundaries;
        in_file >> n_boundaries;

        std::vector <Polyline_polygon> Polyline_polygons;

        for (unsigned int r = 0; r < n_boundaries; ++r) {

          unsigned int n_points;
          in_file >> n_points;
          std::vector<Point> points;

          std::vector<Polyline_X_monotone_curve> xcvs;
          Gps_traits traits;
          auto make_x_monotone = traits.make_x_monotone_2_object();
          using Pnt = typename Gps_traits::Point_2;
          using Xcv = typename Gps_traits::X_monotone_curve_2;
          using Make_x_monotone_result = std::variant<Pnt, Xcv>;

          for (unsigned int j = 0; j < n_points; ++j) {
            double x, y;
            in_file >> x >> y;
            points.push_back(Point(x, y));
            int pt_pos=0;
            if (points.size() > 1) {
              Gps_traits x;
              X_monotone_subcurve_2 seg = x.subcurve_traits_2()->
                construct_x_monotone_curve_2_object()(points[points.size()-2],
                                                      points[points.size()-1]);
              //cout<<points[points.size()-2]<<" "<<points[points.size() - 1]<<endl;

              mPolylinePolygonPieces.push_back(Polyline_curve(seg));
            }
          }
          points.clear();

          for (auto it = mPolylinePolygonPieces.begin();
               it != mPolylinePolygonPieces.end(); ++it) {
            std::vector<Make_x_monotone_result> x_objs;

            make_x_monotone(*it, std::back_inserter(x_objs));
            for (auto i = 0; i < x_objs.size(); ++i) {
              auto* xcv = std::get_if<Polyline_X_monotone_curve>(&x_objs[i]);
              CGAL_assertion(xcv != nullptr);
              xcvs.push_back(*xcv);
            }
          }
          if (xcvs.size() > 0) {
            Polyline_point const &first_point = xcvs.front()[0].source();
            Polyline_point const &last_point = xcvs.back()[xcvs.back().number_of_subcurves() - 1].target();
            FT fxs = first_point.x();
            FT fys = first_point.y();
            FT lxs = last_point.x();
            FT lys = last_point.y();

            Gps_traits x;
            X_monotone_subcurve_2 seg = x.subcurve_traits_2()->
              construct_x_monotone_curve_2_object()(Point(lxs, lys), Point(fxs, fys));

            xcvs.push_back(Polyline_X_monotone_curve(seg));
            //for(auto &i:xcvs)
            //cout<<i<<endl;

            Polyline_polygon pp(xcvs.begin(), xcvs.end());
            CGAL::Orientation orient = pp.orientation();
            Polyline_polygon_with_holes PPWH(pp);
            if (orient == CGAL::CLOCKWISE) {
              pp.reverse_orientation();
            }
            if (r == 0) {
              get_new_state(11);
              states_stack.back().active_set(m_color_active).polyline().join(PPWH);
              states_stack.back().active_polyline_sources(m_color_active).push_back(PPWH);
            }
            else {
              states_stack.back().result_set().clear();
              states_stack.back().result_polyline_sources().clear();
              states_stack.back().result_set().polyline().join(PPWH);
              states_stack.back().result_polyline_sources().push_back(PPWH);
              states_stack.back().active_set(m_color_active).difference(states_stack.back().result_set());
              states_stack.back().result_set().clear();
              states_stack.back().result_polyline_sources().clear();
            }
            mPolylinePolygonPieces.clear();
          }
        }
        rOK = true;
      }
      //cout << "end of read polyline" << endl<<endl;
    }
    catch (...) {
      show_warning("Exception ocurred during reading of linear polygon set.");
    }
  }
  return rOK;
}

bool MainWindow::read_bezier(QString aFileName) {

  bool rOK = false;

  std::ifstream in_file(qPrintable(aFileName));

  if (in_file) {
    try {
      std::string format;
      std::getline(in_file, format);

      bool lDoubleFormat = (format.length() >= 6 && format.substr(0, 6) == "DOUBLE");

      // Read the number of bezier polygon with holes
      unsigned int n_regions;
      in_file >> n_regions;

      for (unsigned int r = 0; r < n_regions; ++r) {
        Bezier_polygon_vector bezier_polygons;
        Bezier_region_source br_source;
        // Read the number of bezier curves.
        unsigned int n_boundaries;
        in_file >> n_boundaries;

        for (unsigned int b = 0; b < n_boundaries; ++b) {
          Bezier_boundary_source bb_source;

          // Read the number of bezier curves.
          unsigned int n_curves;
          in_file >> n_curves;

          // Read the curves one by one, and construct the general polygon these
          // curve form (the outer boundary and the holes inside it).

          std::vector <Bezier_X_monotone_curve> xcvs;

          for (unsigned int k = 0; k < n_curves; ++k) {
            // Read the current curve and subdivide it into x-monotone subcurves.

            using Pnt = typename Bezier_traits::Point_2;
            using Xcv = typename Bezier_traits::X_monotone_curve_2;
            using Make_x_monotone_result = std::variant<Pnt, Xcv>;
            Bezier_X_monotone_curve xcv;
            Bezier_traits traits;
            std::vector<Make_x_monotone_result> x_objs;
            auto make_x_monotone = traits.make_x_monotone_2_object();

            Bezier_curve b = read_bezier_curve(in_file, lDoubleFormat);
            if (b.number_of_control_points() >= 2) {
              bb_source.push_back(b);
              make_x_monotone(b, std::back_inserter(x_objs));

              for (auto i = 0; i < x_objs.size(); ++i) {
                auto *xcv = std::get_if<Bezier_X_monotone_curve>(&x_objs[i]);
                CGAL_assertion(xcv != nullptr);
                xcvs.push_back(*xcv);
              }
            }
          }

          Bezier_polygon pgn(xcvs.begin(), xcvs.end());

          CGAL::Orientation orient = pgn.orientation();
          //TRACE( "  Orientation: " << orient ) ;
          if (orient == CGAL::CLOCKWISE) pgn.reverse_orientation();

          if (b == 0) {
            get_new_state(11);
            states_stack.back().active_set(m_color_active).bezier().join(Bezier_polygon_with_holes(pgn));
            Bezier_region_source br;
            br.push_back(bb_source);
            states_stack.back().active_bezier_sources(m_color_active).push_back(br);
          } else {
            states_stack.back().result_set().clear();
            states_stack.back().result_bezier_sources().clear();
            states_stack.back().result_set().bezier().join(Bezier_polygon_with_holes(pgn));
            Bezier_region_source br;
            br.push_back(bb_source);
            states_stack.back().result_bezier_sources().push_back(br);
            states_stack.back().active_set(m_color_active).difference(states_stack.back().result_set());
            states_stack.back().result_set().clear();
            states_stack.back().result_bezier_sources().clear();
          }
        }
        rOK = true;
      }
    }
    catch (const std::exception &e) {
      // std::string s = e.what();
      // show_error(s);
      on_actionUndo_triggered();
    }
  }
  return rOK;
}

Bezier_curve MainWindow::read_bezier_curve(std::istream& is,
                                           bool aDoubleFormat) {
  // Read the number of control points.
  unsigned int n;
  is >> n;

  // Read the control points.
  std::vector <Bezier_rat_point> ctrl_pts;

  for (unsigned int k = 0; k < n; k++) {
    Bezier_rat_point p;
    if (aDoubleFormat) {
      double x, y;
      is >> x >> y;
      Bezier_rational rx(static_cast<int> (1000 * x + 0.5), 1000);
      Bezier_rational ry(static_cast<int> (1000 * y + 0.5), 1000);
      p = Bezier_rat_point(rx, ry);
    }
    else {
      is >> p;
    }

    if (k == 0 || ctrl_pts[k - 1] != p) {
      ctrl_pts.push_back(p);
    }
  }
  std::vector <Bezier_rat_point> ctrl_pts2;
  typedef std::vector<Bezier_rat_point>::const_iterator cp_const_iterator;
  cp_const_iterator beg = ctrl_pts.begin();
  cp_const_iterator end = ctrl_pts.end();
  cp_const_iterator last = end - 1;
  ctrl_pts2.push_back(*beg);
  if (ctrl_pts.size() > 2) {
    cp_const_iterator curr = beg;
    cp_const_iterator next1 = curr + 1;
    cp_const_iterator next2 = next1 + 1;
    do {
      CGAL::Orientation lOrient = orientation(*curr, *next1, *next2);
      if (lOrient != CGAL::COLLINEAR)
        ctrl_pts2.push_back(*next1);
      ++curr;
      ++next1;
      ++next2;
    } while (next2 != end);
  }
  ctrl_pts2.push_back(*last);
  return Bezier_curve(ctrl_pts2.begin(), ctrl_pts2.end());
}

//save linear polygon
bool save_linear(QString aFileName, Linear_polygon_set &rSet) {
  bool rOK = false;

  Linear_polygon_set lps;
  typedef typename Kernel::Point_2 Point;
  std::list <Polygon_2> pgns;

  std::ofstream out_file(qPrintable(aFileName));
  if (out_file) {
    out_file << "DOUBLE" << std::endl;
    std::vector <Linear_polygon_with_holes> lpwh_container;
    rSet.polygons_with_holes(std::back_inserter(lpwh_container));

    out_file << lpwh_container.size() << std::endl;

    for (auto rit = lpwh_container.begin();
         rit != lpwh_container.end(); ++rit) {
      Linear_polygon_with_holes lpwh = *rit;

      out_file << " " << (1 + lpwh.number_of_holes()) << std::endl;

      int cc = lpwh.outer_boundary().size();
      int lc = cc - 1;
      out_file << "  " << cc << std::endl;
      int i = 0;

      for (auto cit = lpwh.outer_boundary().curves_begin();
           cit != lpwh.outer_boundary().curves_end(); ++cit, ++i) {
        auto pt = cit->source();
        out_file << "   " << CGAL::to_double(pt.x()) << " " << CGAL::to_double(pt.y()) << std::endl;
      }

      for (auto hit = lpwh.holes_begin(); hit != lpwh.holes_end(); ++hit) {
        int cc = hit->size();
        int lc = cc - 1;
        out_file << "  " << cc << std::endl;
        int i = 0;

        for (auto cit = hit->curves_begin(); cit != hit->curves_end();
             ++cit, ++i) {
          auto pt = cit->source();
          out_file << "   " << CGAL::to_double(pt.x()) << " "
                   << CGAL::to_double(pt.y()) << std::endl;
        }
      }
      rOK = true;
    }
  }
  return rOK;
}

bool save_circular(QString aFileName, Circular_polygon_set &rSet)
{
  bool rOK = false;
  std::ofstream out_file(qPrintable(aFileName));

  if (out_file) {
    out_file << "DOUBLE" << std::endl;
    std::vector <Circular_polygon_with_holes> cpwh_container;
    rSet.polygons_with_holes(std::back_inserter(cpwh_container));
    out_file << cpwh_container.size() << std::endl;

    for (auto rit = cpwh_container.begin();
         rit != cpwh_container.end(); ++rit) {
      Circular_polygon_with_holes cpwh = *rit;
      out_file << " " << (1 + cpwh.number_of_holes()) << std::endl;
      int cc = cpwh.outer_boundary().size();
      int lc = cc - 1;
      out_file << "  " << cc << std::endl;
      int i = 0;

      for (auto cit = cpwh.outer_boundary().curves_begin();
           cit != cpwh.outer_boundary().curves_end(); ++cit, ++i) {
        auto pt = cit->source();
        if (cit->is_circular()) {
          out_file << "   1" << std::endl;
          out_file << "    " << CGAL::to_double(pt.x()) << " "
                   << CGAL::to_double(pt.y()) << std::endl;

          double r, mx, my, tx, ty, cx, cy, sx, sy, st, tnx, tny, frmst;
          r = sqrt(CGAL::to_double(cit->supporting_circle().squared_radius())); //radius r
          sx = CGAL::to_double(cit->source().x()); // source x()
          sy = CGAL::to_double(cit->source().y()); // source y()
          tx = CGAL::to_double(cit->target().x()); // target x()
          ty = CGAL::to_double(cit->target().y()); // target y()
          cx = CGAL::to_double(cit->supporting_circle().center().x()); // center x()
          cy = CGAL::to_double(cit->supporting_circle().center().y()); // center y()
          st = sqrt((sx - tx) * (sx - tx) + (sy - ty) * (sy - ty)); //Source to target euclides distance ST
          tnx = tx - cx; // origin shifting to center
          tny = ty - cy; // origin shifting to center
          frmst = sqrt((4 * r * r) - (st * st)); // intermediate value

          if (cit->supporting_circle().orientation() == CGAL::CLOCKWISE) {
            mx = (((frmst * tnx) - (st * tny)) / (2 * r)) + cx;
            my = (((st * tnx) + (frmst * tny)) / (2 * r)) + cy;
          } else {
            mx = (((tnx * frmst) + (tny * st)) / (2 * r)) + cx;
            my = (((tny * frmst) - (st * tnx)) / (2 * r)) + cy;
          }

          out_file << "    " << mx << " " << my << std::endl;
        }
        else {
          out_file << "   0" << std::endl;
          out_file << "    " << CGAL::to_double(pt.x()) << " "
                   << CGAL::to_double(pt.y()) << std::endl;
        }
        pt = cit->target();
        out_file << "    " << CGAL::to_double(pt.x()) << " "
                 << CGAL::to_double(pt.y()) << std::endl;
      }

      for (auto hit = cpwh.holes_begin(); hit != cpwh.holes_end(); ++hit) {
        int cc = hit->size();
        int lc = cc - 1;
        out_file << "  " << cc << std::endl;
        int i = 0;

        for (auto cit = hit->curves_begin(); cit != hit->curves_end();
             ++cit, ++i) {
          auto pt = cit->source();
          if (cit->is_circular()) {
            out_file << "   1" << std::endl;
            out_file << "    " << CGAL::to_double(pt.x()) << " "
                     << CGAL::to_double(pt.y()) << std::endl;

            double r, mx, my, tx, ty, cx, cy, sx, sy, st, tnx, tny, frmst;
            r = sqrt(CGAL::to_double(cit->supporting_circle().squared_radius())); //radius r
            sx = CGAL::to_double(cit->source().x()); // source x()
            sy = CGAL::to_double(cit->source().y()); // source y()
            tx = CGAL::to_double(cit->target().x()); // target x()
            ty = CGAL::to_double(cit->target().y()); // target y()
            cx = CGAL::to_double(cit->supporting_circle().center().x()); // center x()
            cy = CGAL::to_double(cit->supporting_circle().center().y()); // center y()
            st = sqrt((sx - tx) * (sx - tx)
                      + (sy - ty) * (sy - ty)); //Source to target euclides distance ST
            tnx = tx - cx; // origin shifting to center
            tny = ty - cy; // origin shifting to center
            frmst = sqrt((4 * r * r) - (st * st)); // intermediate value

            if (cit->supporting_circle().orientation() == CGAL::CLOCKWISE) {
              mx = (((frmst * tnx) - (st * tny)) / (2 * r)) + cx;
              my = (((st * tnx) + (frmst * tny)) / (2 * r)) + cy;
            } else {
              mx = (((tnx * frmst) + (tny * st)) / (2 * r)) + cx;
              my = (((tny * frmst) - (st * tnx)) / (2 * r)) + cy;
            }

            out_file << "    " << mx << " " << my << std::endl;
          }
          else {
            out_file << "   0" << std::endl;
            out_file << "    " << CGAL::to_double(pt.x()) << " "
                     << CGAL::to_double(pt.y()) << std::endl;
          }
          pt = cit->target();
          out_file << "    " << CGAL::to_double(pt.x()) << " " << CGAL::to_double(pt.y()) << std::endl;
        }
      }
      rOK = true;
    }
  }
  return rOK;
}

bool save_polyline(QString aFileName, Polyline_polygon_set &rSet) {
  bool rOK = false;

  Polyline_polygon_set pps;
  typedef typename Kernel::Point_2 Point;
  std::list <Polygon_2> pgns;

  std::ofstream out_file(qPrintable(aFileName));
  if (out_file) {
    out_file << "DOUBLE" << std::endl;
    std::vector <Polyline_polygon_with_holes> lpwh_container;
    rSet.polygons_with_holes(std::back_inserter(lpwh_container));

    out_file << lpwh_container.size() << std::endl;

    for (auto rit = lpwh_container.begin(); rit != lpwh_container.end(); ++rit)
    {
      Polyline_polygon_with_holes lpwh = *rit;

      out_file << " " << (1 + lpwh.number_of_holes()) << std::endl;

      int cc = lpwh.outer_boundary().size();
      int lc = cc - 1;
      out_file << "  " << cc << std::endl;
      int i = 0;

      for (auto cit = lpwh.outer_boundary().curves_begin();
           cit != lpwh.outer_boundary().curves_end(); ++cit, ++i) {
        auto it = cit->points_begin();
        out_file << " " << CGAL::to_double((*it).x()) << " "
                 << CGAL::to_double((*it).y()) << std::endl;
        /*while (it != cit->points_end()) {
          auto pt = *it;
          cout<<pt<<endl;
          (*it)++;
          }
          cout<<endl;*/
      }

      for (auto hit = lpwh.holes_begin(); hit != lpwh.holes_end(); ++hit) {
        int cc = hit->size();
        int lc = cc - 1;
        out_file << "  " << cc << std::endl;
        int i = 0;

        for (auto cit = hit->curves_begin(); cit != hit->curves_end();
             ++cit, ++i) {
          auto it = cit->points_begin();
          out_file << " " << CGAL::to_double((*it).x()) << " "
                   << CGAL::to_double((*it).y()) << std::endl;
          /*while (it != cit->points_end()) {
            auto pt = *it;
            out_file << "   " << CGAL::to_double(pt.x()) << " " << CGAL::to_double(pt.y()) << std::endl;
            *it++;
            }*/
        }
      }
      rOK = true;
    }
  }
  return rOK;
}

//bezier file in saved in two parts
void save_bezier_polygon(std::ostream &out_file, Bezier_polygon const &aBP) {
  typedef std::vector<Bezier_rat_point> Bezier_rat_point_vector;

  int cc = aBP.size();
  int lc = cc - 1;
  out_file << "  " << cc << std::endl;
  Bezier_rat_point lFirstP, lPrevP;
  int i = 0;

  for (auto cit = aBP.curves_begin(); cit != aBP.curves_end(); ++cit, ++i) {
    Bezier_rat_point_vector lQ;
    CGAL::Qt::Bezier_helper::clip(*cit, lQ);
    out_file << "   " << lQ.size() << std::endl;

    if (i == 0)
      lFirstP = lQ.front();

    if (i == lc)
      lQ.back() = lFirstP;

    for (auto pit = lQ.begin(); pit != lQ.end(); ++pit) {
      Bezier_rat_point lP = pit == lQ.begin() && i > 0 ? lPrevP : *pit;
      out_file << "    " << CGAL::to_double(lP.x()) << " " << CGAL::to_double(lP.y()) << std::endl;
      lPrevP = lP;
    }
  }
}

bool save_bezier_result(QString aFileName, Bezier_polygon_set const& aSet) {
  bool rOK = false;
  std::ofstream out_file(qPrintable(aFileName));
  if (out_file) {
    out_file << "DOUBLE" << std::endl;
    std::vector <Bezier_polygon_with_holes> bpwh_container;
    aSet.polygons_with_holes(std::back_inserter(bpwh_container));
    out_file << bpwh_container.size() << std::endl;

    for (auto rit = bpwh_container.begin();
         rit != bpwh_container.end(); ++rit) {
      Bezier_polygon_with_holes bpwh = *rit;
      out_file << " " << (1 + bpwh.number_of_holes()) << std::endl;
      save_bezier_polygon(out_file, bpwh.outer_boundary());
      for (auto hit = bpwh.holes_begin(); hit != bpwh.holes_end(); ++hit)
        save_bezier_polygon(out_file, *hit);

      rOK = true;
    }
  }
  return rOK;
}

//"file:save active bucket" button pressed
void MainWindow::on_actionSaveCurrentBucket_triggered() {
  if (m_circular_active) {
    if (!save_circular(QFileDialog::getSaveFileName(this,
                                                    tr("Save Result Circular Polygon Set"), "/data/index.cps",
                                                    tr("Circular Curve files (*.cps)")),
                       states_stack.back().active_set(m_color_active).circular()
                       )
        ) {
      show_error("Cannot save circular polygon set.");
    }
  }
  else if (m_bezier_active) {
    if (!save_bezier_result(QFileDialog::getSaveFileName(this,
                                                         tr("Save Result Bezier Polygon Set"), "/data/index.bps",
                                                         tr("Bezier Curve files (*.bps)")),
                            states_stack.back().active_set(m_color_active).bezier()
                            )
        ) {
      show_error("Cannot save bezier polygon set.");
    }
  }
  else if (m_polyline_active) {
    if (!save_polyline(QFileDialog::getSaveFileName(this,
                                                    tr("Save Result Polyline Polygon Set"), "/data/index.mps",
                                                    tr("Polyline Curve files (*.mps)")),
                       states_stack.back().active_set(m_color_active).polyline()
                       )
        ) {
      show_error("Cannot save polyline polygon set.");
    }
  }
  else {
    if (!save_linear(QFileDialog::getSaveFileName(this,
                                                  tr("Save Result Linear Polygon Set"), "/data/index.lps",
                                                  tr("Linear Curve files (*.lps)")),
                     states_stack.back().active_set(m_color_active).linear()
                     )
        ) {
      show_error("Cannot save linear polygon set.");
    }
  }
}

//check out
void MainWindow::switch_set_type(Curve_set &aSet, int aType)
{
  unlink_GI(aSet.gi());
  aSet.reset_type(aType);
  link_GI(aSet.gi());
  modelChanged();
}

void MainWindow::switch_sets_type(int aType) {
  switch_set_type(states_stack.back().blue_set(), aType);
  switch_set_type(states_stack.back().red_set(), aType);
  switch_set_type(states_stack.back().black_set(), aType);
  switch_set_type(states_stack.back().brown_set(), aType);
  switch_set_type(states_stack.back().yellow_set(), aType);
  switch_set_type(states_stack.back().magenta_set(), aType);
  switch_set_type(states_stack.back().aqua_set(), aType);
  switch_set_type(states_stack.back().result_set(), aType);
}

//switching from one polygon to another
bool MainWindow::ensure_circular_mode() {
  if (! m_circular_active) {
    bool lProceed = states_stack.back().blue_set().is_empty() && states_stack.back().red_set().is_empty() &&
      states_stack.back().black_set().is_empty() && states_stack.back().brown_set().is_empty() &&
      states_stack.back().yellow_set().is_empty() && states_stack.back().magenta_set().is_empty() &&
      states_stack.back().aqua_set().is_empty() && states_stack.back().result_set().is_empty();

    if (! lProceed)
      lProceed = ask_user_yesno("Circular mode switch",
                                "You are about to load a circular poygon, but there are linear/bezier/polyline curves already loaded.\n" \
                                "Both types are not interoperable. In order to proceed, the polygons must be removed first.\n" \
                                "Once deleted you cannot undo the action.\n" \
                                "Yes to remove and proceed?\n");

    if (lProceed) {
      switch_sets_type(2);
      m_linear_input->Reset();
      m_bezier_input->Reset();
      m_polyline_input->Reset();
      m_circular_active = true;
      m_bezier_active = false;
      m_polyline_active = false;
    }
  }
  return m_circular_active;
}

bool MainWindow::ensure_bezier_mode() {
  if (! m_bezier_active) {
    bool lProceed = states_stack.back().blue_set().is_empty() && states_stack.back().red_set().is_empty() &&
      states_stack.back().black_set().is_empty() && states_stack.back().brown_set().is_empty() &&
      states_stack.back().yellow_set().is_empty() && states_stack.back().magenta_set().is_empty() &&
      states_stack.back().aqua_set().is_empty() && states_stack.back().result_set().is_empty();

    if (!lProceed)
      lProceed = ask_user_yesno("Bezier mode switch",
                                "You are about to load a Bezier curve, but there are linear/bezier/polyline polygons already loaded.\n" \
                                "Both types are not interoperable. In order to proceed, the polygons must be removed first.\n" \
                                "Once deleted you cannot undo the action.\n" \
                                "Yes to remove and proceed?\n");

    if (lProceed) {
      switch_sets_type(3);
      m_circular_input->Reset();
      m_linear_input->Reset();
      m_polyline_input->Reset();
      m_bezier_active = true;
      m_circular_active = false;
      m_polyline_active = false;
    }
  }
  return m_bezier_active;
}

bool MainWindow::ensure_linear_mode() {
  if (m_circular_active || m_bezier_active || m_polyline_active) {
    bool lProceed = states_stack.back().blue_set().is_empty() && states_stack.back().red_set().is_empty() &&
      states_stack.back().black_set().is_empty() && states_stack.back().brown_set().is_empty() &&
      states_stack.back().yellow_set().is_empty() && states_stack.back().magenta_set().is_empty() &&
      states_stack.back().aqua_set().is_empty() && states_stack.back().result_set().is_empty();

    if (! lProceed)
      lProceed = ask_user_yesno("Linear mode switch",
                                "You are about to load a linear poygon, but there are circular/bezier/polyline polygons already loaded.\n" \
                                "Both types are not interoperable. In order to proceed, the polygons must be removed first.\n" \
                                "Once deleted you cannot undo the action.\n"
                                "Yes to remove and proceed?\n");

    // cout<<"Linear Mode"<<endl;
    if (lProceed) {
      //cout<<"fault 1"<<endl;
      switch_sets_type(1);
      m_circular_input->Reset();
      m_bezier_input->Reset();
      m_polyline_input->Reset();
      m_circular_active = false;
      m_bezier_active = false;
      //cout<<"fault 3"<<endl;
      m_polyline_active = false;
      //cout<<"End of ensure_linear_polygon()"<<endl;
    }
    //cout<<"end of this function"<<endl;
  }
  return ! m_circular_active;
}

bool MainWindow::ensure_polyline_mode() {
  if (!m_polyline_active) {
    bool lProceed = states_stack.back().blue_set().is_empty() && states_stack.back().red_set().is_empty() &&
      states_stack.back().black_set().is_empty() && states_stack.back().brown_set().is_empty() &&
      states_stack.back().yellow_set().is_empty() && states_stack.back().magenta_set().is_empty() &&
      states_stack.back().aqua_set().is_empty() && states_stack.back().result_set().is_empty();

    if (!lProceed)
      lProceed = ask_user_yesno("Polyline mode switch",
                                "You are about to load a polyline polygon, but there are circular/bezier/linear polygons already loaded.\n" \
                                "Both types are not interoperable. In order to proceed, the polygons must be removed first.\n" \
                                "Once deleted you cannot undo the action.\n" \
                                "Yes to remove and proceed?\n");

    if (lProceed) {
      switch_sets_type(4);
      m_circular_input->Reset();
      m_bezier_input->Reset();
      m_linear_input->Reset();
      m_circular_active = false;
      m_bezier_active = false;
      m_polyline_active = true;
    }
  }
  return m_polyline_active;
}

//load polygons from file that is "file:open" clicked
void MainWindow::open(QString fileName) {
  if (!fileName.isEmpty()) {
    bool lRead = false;
    if (fileName.endsWith(".lps")) {
      if (ensure_linear_mode())
        lRead = read_linear(fileName, states_stack.back().active_set(m_color_active).linear(),
                            states_stack.back().active_linear_sources(m_color_active));
    }
    else if (fileName.endsWith(".cps")) {
      if (ensure_circular_mode())
        lRead = read_circular(fileName, states_stack.back().active_set(m_color_active).circular(),
                              states_stack.back().active_circular_sources(m_color_active));
    }
    else if (fileName.endsWith(".bps")) {
      if (ensure_bezier_mode())
        lRead = read_bezier(fileName);
    }
    else if (fileName.endsWith(".mps")) {
      if (ensure_polyline_mode())
        lRead = read_polyline(fileName, states_stack.back().active_set(m_color_active).polyline(),
                              states_stack.back().active_polyline_sources(m_color_active));
    }
    if (lRead) {
      modelChanged();
      // zoomToFit();
      this->addToRecentFiles(fileName);
    }
  }
}

void MainWindow::on_actionInsertCircular_toggled(bool aChecked) {
  if (aChecked) {
    this->graphicsView->setDragMode(QGraphicsView::NoDrag);
    if (ensure_circular_mode()) {
      if (!m_pan && !m_reset) {
        while (m_state_num > 1) {
          on_actionUndo_triggered();
        }
      }

      actionPAN->setChecked(false);
      m_pan = false;
      actionOpenLinear->setEnabled(false);
      actionOpenDXF->setEnabled(true);
      actionOpenBezier->setEnabled(false);
      actionOpenPolyline->setEnabled(false);
      actionInsertLinear->setChecked(false);
      actionInsertBezier->setChecked(false);
      actionInsertPolyline->setChecked(false);
      m_scene.installEventFilter(m_circular_input);
      on_actionDeleteResult();
      actionMinkowski_SumH->setEnabled(false);

      m_circular_input->mOngoingPieceGI->setPen(sPens[m_color_active]);
      m_circular_input->mCircularGI->setPen(sPens[m_color_active]);
      m_circular_input->mHandleGI->setPen(sPens[m_color_active]);
    } else {
      actionInsertCircular->setChecked(false);
    }
  }
}

void MainWindow::on_actionInsertBezier_toggled(bool aChecked) {
  if (aChecked) {
    this->graphicsView->setDragMode(QGraphicsView::NoDrag);
    if (ensure_bezier_mode()) {
      if (!m_pan && !m_reset) {
        while (m_state_num > 1) {
          on_actionUndo_triggered();
        }
      }

      actionPAN->setChecked(false);
      m_pan = false;
      actionOpenLinear->setEnabled(false);
      actionOpenDXF->setEnabled(false);
      actionOpenBezier->setEnabled(true);
      actionOpenPolyline->setEnabled(false);
      actionInsertLinear->setChecked(false);
      actionInsertCircular->setChecked(false);
      actionInsertPolyline->setChecked(false);
      m_scene.installEventFilter(m_bezier_input);
      on_actionDeleteResult();

      actionMinkowski_SumH->setEnabled(false);

      m_bezier_input->mOngoingPieceGI->setPen(sPens[m_color_active]);
      m_bezier_input->mBezierGI->setPen(sPens[m_color_active]);
      m_bezier_input->mHandle0GI->setPen(sPens[m_color_active]);
      m_bezier_input->mHandle1GI->setPen(sPens[m_color_active]);
    } else {
      actionInsertBezier->setChecked(false);
    }
  }
}

void MainWindow::on_actionInsertPolyline_toggled(bool aChecked) {
  if (aChecked) {
    this->graphicsView->setDragMode(QGraphicsView::NoDrag);
    if (ensure_polyline_mode()) {
      if (!m_pan && !m_reset) {
        while (m_state_num > 1) {
          on_actionUndo_triggered();
        }
      }

      actionPAN->setChecked(false);
      m_pan = false;
      actionOpenLinear->setEnabled(false);
      actionOpenDXF->setEnabled(false);
      actionOpenBezier->setEnabled(false);
      actionOpenPolyline->setEnabled(true);
      actionInsertLinear->setChecked(false);
      actionInsertBezier->setChecked(false);

      actionInsertCircular->setChecked(false);
      m_scene.installEventFilter(m_polyline_input);
      on_actionDeleteResult();
      actionMinkowski_SumH->setEnabled(false);

      m_polyline_input->mOngoingPieceGI->setPen(sPens[m_color_active]);
      m_polyline_input->mPolylineGI->setPen(sPens[m_color_active]);
      m_polyline_input->mHandleGI->setPen(sPens[m_color_active]);
    } else {
      actionInsertPolyline->setChecked(false);
    }
  }
}

void MainWindow::on_actionInsertLinear_toggled(bool aChecked) {
  if (aChecked) {
    this->graphicsView->setDragMode(QGraphicsView::NoDrag);
    if (ensure_linear_mode()) {
      if (!m_pan && !m_reset) {
        while (m_state_num > 1) {
          on_actionUndo_triggered();
        }
      }
      actionPAN->setChecked(false);
      m_pan = false;
      actionOpenLinear->setEnabled(true);
      actionOpenDXF->setEnabled(false);
      actionOpenBezier->setEnabled(false);
      actionOpenPolyline->setEnabled(false);
      actionInsertCircular->setChecked(false);
      actionInsertBezier->setChecked(false);
      actionInsertPolyline->setChecked(false);
      m_scene.installEventFilter(m_linear_input);
      on_actionDeleteResult();

      actionMinkowski_SumH->setEnabled(true);
    }
    else {
      actionInsertLinear->setChecked(false);
    }
  }
}

//complement
void MainWindow::on_actionComplementH_toggled(bool aChecked) {
  if (actionComplementH->isChecked()) {
    bool lDone = false;
    QCursor old = this->cursor();
    this->setCursor(Qt::WaitCursor);
    actionUnionH->setChecked(false);
    actionIntersectionH->setChecked(false);
    actionDifferenceH->setChecked(false);
    actionSymmetric_DifferenceH->setChecked(false);
    actionMinkowski_SumH->setChecked(false);
    actionCopyH->setChecked(false);
    actionMoveH->setChecked(false);

    get_new_state(0);
    /*operation_name
    //operation_name
    // COMPLEMENT_OP = 0
    // INTERSECTION_OP = 1
    // UNION_OP = 2
    // DIFFERENCE_OP = 3
    // SYMMETRIC_DIFFERENCE_OP = 4
    // MINKOWSKI_SUM_OP = 5
    // RESET_OP = 6
    // COPY_OP = 7
    // MOVE_OP = 8
    // CLEAR_OP = 9
    // DELETEALL_OP = 10
    // START_OP = 11*/

    states_stack.back().result_set().clear();
    states_stack.back().result_linear_sources().clear();
    states_stack.back().result_circular_sources().clear();
    states_stack.back().result_bezier_sources().clear();

    if (! m_circular_active && ! m_bezier_active && ! m_polyline_active) {
      m_linear_input->get_BoundingRect();
    }
    else if (! m_bezier_active && !m_polyline_active) {
      m_circular_input->get_BoundingRect();
    }
    else if (! m_polyline_active) {
      m_bezier_input->get_BoundingRect();
    }
    else {
      m_polyline_input->get_BoundingRect();
    }
    //m_circular_input->get_BoundingRect();

    if (!states_stack.back().active_set(m_color_active).is_empty()) {
      if (empty_warn) {
        show_not_empty_warning();
      }
    }

    switch (m_color_active) {
     case 0:
      switch (m_color_complement) {
       case 0:
        if (!states_stack.back().blue_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().blue_set());
        }
        break;
       case 1:
        if (!states_stack.back().red_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().red_set());
        }
        break;
       case 2:
        if (!states_stack.back().black_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().black_set());
        }
        break;
       case 3:
        if (!states_stack.back().brown_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().brown_set());
        }
        break;
       case 4:
        if (!states_stack.back().yellow_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().yellow_set());
        }
        break;
       case 5:
        if (!states_stack.back().magenta_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        }
        break;
       case 6:
        if (!states_stack.back().aqua_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().aqua_set());
        }
        break;
      }
      states_stack.back().blue_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      break;

     case 1:
      switch (m_color_complement) {
       case 0:
        if (!states_stack.back().blue_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().blue_set());
        }
        break;
       case 1:
        if (!states_stack.back().red_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().red_set());
        }
        break;
       case 2:
        if (!states_stack.back().black_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().black_set());
        }
        break;
       case 3:
        if (!states_stack.back().brown_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().brown_set());
        }
        break;
       case 4:
        if (!states_stack.back().yellow_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().yellow_set());
        }
        break;
       case 5:
        if (!states_stack.back().magenta_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        }
        break;
       case 6:
        if (!states_stack.back().aqua_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().aqua_set());
        }
        break;
      }
      states_stack.back().red_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      break;

     case 2:
      switch (m_color_complement) {
       case 0:
        if (!states_stack.back().blue_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().blue_set());
        }
        break;
       case 1:
        if (!states_stack.back().red_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().red_set());
        }
        break;
       case 2:
        if (!states_stack.back().black_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().black_set());
        }
        break;
       case 3:
        if (!states_stack.back().brown_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().brown_set());
        }
        break;
       case 4:
        if (!states_stack.back().yellow_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().yellow_set());
        }
        break;
       case 5:
        if (!states_stack.back().magenta_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        }
        break;
       case 6:
        if (!states_stack.back().aqua_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().aqua_set());
        }
        break;
      }
      states_stack.back().black_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      break;

     case 3:
      switch (m_color_complement) {
       case 0:
        if (!states_stack.back().blue_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().blue_set());
        }
        break;
       case 1:
        if (!states_stack.back().red_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().red_set());
        }
        break;
       case 2:
        if (!states_stack.back().black_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().black_set());
        }
        break;
       case 3:
        if (!states_stack.back().brown_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().brown_set());
        }
        break;
       case 4:
        if (!states_stack.back().yellow_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().yellow_set());
        }
        break;
       case 5:
        if (!states_stack.back().magenta_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        }
        break;
       case 6:
        if (!states_stack.back().aqua_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().aqua_set());
        }
        break;
      }
      states_stack.back().brown_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      break;

     case 4:
      switch (m_color_complement) {
       case 0:
        if (!states_stack.back().blue_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().blue_set());
        }
        break;
       case 1:
        if (!states_stack.back().red_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().red_set());
        }
        break;
       case 2:
        if (!states_stack.back().black_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().black_set());
        }
        break;
       case 3:
        if (!states_stack.back().brown_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().brown_set());
        }
        break;
       case 4:
        if (!states_stack.back().yellow_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().yellow_set());
        }
        break;
       case 5:
        if (!states_stack.back().magenta_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        }
        break;
       case 6:
        if (!states_stack.back().aqua_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().aqua_set());
        }
        break;
      }
      states_stack.back().yellow_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      break;

     case 5:
      switch (m_color_complement) {
       case 0:
        if (!states_stack.back().blue_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().blue_set());
        }
        break;
       case 1:
        if (!states_stack.back().red_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().red_set());
        }
        break;
       case 2:
        if (!states_stack.back().black_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().black_set());
        }
        break;
       case 3:
        if (!states_stack.back().brown_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().brown_set());
        }
        break;
       case 4:
        if (!states_stack.back().yellow_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().yellow_set());
        }
        break;
       case 5:
        if (!states_stack.back().magenta_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        }
        break;
       case 6:
        if (!states_stack.back().aqua_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().aqua_set());
        }
        break;
      }
      states_stack.back().magenta_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      break;

     case 6:
      switch (m_color_complement) {
       case 0:
        if (!states_stack.back().blue_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().blue_set());
        }
        break;
       case 1:
        if (!states_stack.back().red_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().red_set());
        }
        break;
       case 2:
        if (!states_stack.back().black_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().black_set());
        }
        break;
       case 3:
        if (!states_stack.back().brown_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().brown_set());
        }
        break;
       case 4:
        if (!states_stack.back().yellow_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().yellow_set());
        }
        break;
       case 5:
        if (!states_stack.back().magenta_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        }
        break;
       case 6:
        if (!states_stack.back().aqua_set().is_empty()) {
          states_stack.back().result_set().difference(states_stack.back().aqua_set());
        }
        break;
      }
      states_stack.back().aqua_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      break;

     default :break;
    }
    actionComplementH->setChecked(false);
    lDone = true;
    this->setCursor(old);
    if (lDone) modelChanged();
  }
}

//intersection
void MainWindow::on_actionIntersectionH_toggled(bool aChecked) {
  if (actionIntersectionH->isChecked()) {
    //std::cout<<"Intersection"<<endl;
    bool lDone = false;
    QCursor old = this->cursor();
    this->setCursor(Qt::WaitCursor);

    actionComplementH->setChecked(false);
    actionUnionH->setChecked(false);
    actionDifferenceH->setChecked(false);
    actionSymmetric_DifferenceH->setChecked(false);
    actionMinkowski_SumH->setChecked(false);
    actionCopyH->setChecked(false);
    actionMoveH->setChecked(false);
    get_new_state(1);

    /*operation_name
    // COMPLEMENT_OP = 0
    // INTERSECTION_OP = 1
    // UNION_OP = 2
    // DIFFERENCE_OP = 3
    // SYMMETRIC_DIFFERENCE_OP = 4
    // MINKOWSKI_SUM_OP = 5
    // RESET_OP = 6
    // COPY_OP = 7
    // MOVE_OP = 8
    // CLEAR_OP = 9
    // DELETEALL_OP = 10
    // START_OP = 11*/

    states_stack.back().result_set().clear();
    states_stack.back().result_linear_sources().clear();
    states_stack.back().result_circular_sources().clear();
    states_stack.back().result_bezier_sources().clear();
    states_stack.back().result_polyline_sources().clear();

    if (!states_stack.back().active_set(m_color_active).is_empty()) {
      if (empty_warn) {
        show_not_empty_warning();
      }
    }
    switch (m_color_active) {
     case 0:
      if (!states_stack.back().blue_set().is_empty() && m_blue_int)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_int)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_int)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_int)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_int)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_int)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_int)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (m_blue_int) states_stack.back().result_set().intersect(states_stack.back().blue_set());
      if (m_red_int) states_stack.back().blue_set().intersect(states_stack.back().red_set());
      if (m_black_int) states_stack.back().blue_set().intersect(states_stack.back().black_set());
      if (m_brown_int) states_stack.back().result_set().intersect(states_stack.back().brown_set());
      if (m_yellow_int) states_stack.back().result_set().intersect(states_stack.back().yellow_set());
      if (m_magenta_int) states_stack.back().result_set().intersect(states_stack.back().magenta_set());
      if (m_aqua_int) states_stack.back().result_set().intersect(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().blue_set());
      states_stack.back().blue_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      //std::cout<<"Blue"<<endl;
      break;

     case 1:
      if (!states_stack.back().blue_set().is_empty() && m_blue_int)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_int)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_int)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_int)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_int)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_int)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_int)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (m_blue_int) states_stack.back().result_set().intersect(states_stack.back().blue_set());
      if (m_red_int) states_stack.back().result_set().intersect(states_stack.back().red_set());
      if (m_black_int) states_stack.back().result_set().intersect(states_stack.back().black_set());
      if (m_brown_int) states_stack.back().result_set().intersect(states_stack.back().brown_set());
      if (m_yellow_int) states_stack.back().result_set().intersect(states_stack.back().yellow_set());
      if (m_magenta_int) states_stack.back().result_set().intersect(states_stack.back().magenta_set());
      if (m_aqua_int) states_stack.back().result_set().intersect(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().red_set());
      states_stack.back().red_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      //std::cout<<"Red"<<endl;
      break;

     case 2:
      if (!states_stack.back().blue_set().is_empty() && m_blue_int)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_int)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_int)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_int)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_int)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_int)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_int)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (m_blue_int) states_stack.back().result_set().intersect(states_stack.back().blue_set());
      if (m_red_int) states_stack.back().result_set().intersect(states_stack.back().red_set());
      if (m_black_int) states_stack.back().result_set().intersect(states_stack.back().black_set());
      if (m_brown_int) states_stack.back().result_set().intersect(states_stack.back().brown_set());
      if (m_yellow_int) states_stack.back().result_set().intersect(states_stack.back().yellow_set());
      if (m_magenta_int) states_stack.back().result_set().intersect(states_stack.back().magenta_set());
      if (m_aqua_int) states_stack.back().result_set().intersect(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().black_set());
      states_stack.back().black_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      //std::cout<<"Black"<<endl;
      break;

     case 3:
      if (!states_stack.back().blue_set().is_empty() && m_blue_int)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_int)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_int)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_int)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_int)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_int)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_int)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (m_blue_int) states_stack.back().result_set().intersect(states_stack.back().blue_set());
      if (m_red_int) states_stack.back().result_set().intersect(states_stack.back().red_set());
      if (m_black_int) states_stack.back().result_set().intersect(states_stack.back().black_set());
      if (m_brown_int) states_stack.back().result_set().intersect(states_stack.back().brown_set());
      if (m_yellow_int) states_stack.back().result_set().intersect(states_stack.back().yellow_set());
      if (m_magenta_int) states_stack.back().result_set().intersect(states_stack.back().magenta_set());
      if (m_aqua_int) states_stack.back().result_set().intersect(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().brown_set());
      states_stack.back().brown_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      //std::cout<<"Brown"<<endl;
      break;

     case 4:
      if (!states_stack.back().blue_set().is_empty() && m_blue_int)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_int)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_int)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_int)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_int)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_int)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_int)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (m_blue_int) states_stack.back().result_set().intersect(states_stack.back().blue_set());
      if (m_red_int) states_stack.back().result_set().intersect(states_stack.back().red_set());
      if (m_black_int) states_stack.back().result_set().intersect(states_stack.back().black_set());
      if (m_brown_int) states_stack.back().result_set().intersect(states_stack.back().brown_set());
      if (m_yellow_int) states_stack.back().result_set().intersect(states_stack.back().yellow_set());
      if (m_magenta_int) states_stack.back().result_set().intersect(states_stack.back().magenta_set());
      if (m_aqua_int) states_stack.back().result_set().intersect(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().yellow_set());
      states_stack.back().yellow_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      //std::cout<<"Yellow"<<endl;
      break;

     case 5:
      if (!states_stack.back().blue_set().is_empty() && m_blue_int)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_int)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_int)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_int)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_int)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_int)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_int)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (m_blue_int) states_stack.back().result_set().intersect(states_stack.back().blue_set());
      if (m_red_int) states_stack.back().result_set().intersect(states_stack.back().red_set());
      if (m_black_int) states_stack.back().result_set().intersect(states_stack.back().black_set());
      if (m_brown_int) states_stack.back().result_set().intersect(states_stack.back().brown_set());
      if (m_yellow_int) states_stack.back().result_set().intersect(states_stack.back().yellow_set());
      if (m_magenta_int) states_stack.back().result_set().intersect(states_stack.back().magenta_set());
      if (m_aqua_int) states_stack.back().result_set().intersect(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().magenta_set());
      states_stack.back().magenta_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      //std::cout<<"Magenta"<<endl;
      break;

     case 6:
      if (!states_stack.back().blue_set().is_empty() && m_blue_int)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_int)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_int)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_int)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_int)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_int)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_int)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (m_blue_int) states_stack.back().result_set().intersect(states_stack.back().blue_set());
      if (m_red_int) states_stack.back().result_set().intersect(states_stack.back().red_set());
      if (m_black_int) states_stack.back().result_set().intersect(states_stack.back().black_set());
      if (m_brown_int) states_stack.back().result_set().intersect(states_stack.back().brown_set());
      if (m_yellow_int) states_stack.back().result_set().intersect(states_stack.back().yellow_set());
      if (m_magenta_int) states_stack.back().result_set().intersect(states_stack.back().magenta_set());
      if (m_aqua_int) states_stack.back().result_set().intersect(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().aqua_set());
      states_stack.back().aqua_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      //std::cout<<"Aqua"<<endl;
      break;

     default:break;
    }

    actionIntersectionH->setChecked(false);
    //cout<<"End of Intersection"<<endl;
    lDone = true;
    this->setCursor(old);
    if (lDone) modelChanged();
  }
}

//blue - red
//difference
void MainWindow::on_actionDifferenceH_toggled(bool aChecked) {
  if (actionDifferenceH->isChecked()) {
    bool lDone = false;
    QCursor old = this->cursor();
    this->setCursor(Qt::WaitCursor);

    actionComplementH->setChecked(false);
    actionUnionH->setChecked(false);
    actionIntersectionH->setChecked(false);
    actionSymmetric_DifferenceH->setChecked(false);
    actionMinkowski_SumH->setChecked(false);
    actionCopyH->setChecked(false);
    actionMoveH->setChecked(false);

    size_t count = 0;
    if (showBlueDiff->isChecked()) count++;
    if (showRedDiff->isChecked()) count++;
    if (showBlackDiff->isChecked()) count++;
    if (showBrownDiff->isChecked()) count++;
    if (showYellowDiff->isChecked()) count++;
    if (showMagentaDiff->isChecked()) count++;
    if (showAquaDiff->isChecked()) count++;

    if (count == 2) {
      size_t color1 = 111;
      size_t color2 = 1111;
      if (showBlueDiff->isChecked()) color1 = 0;
      if (showRedDiff->isChecked()) {
        if (color1 < 1) color2 = 1;
        else color1 = 1;
      }

      if (showBlackDiff->isChecked()) {
        if (color1 < 2) color2 = 2;
        else color1 = 2;
      }

      if (showBrownDiff->isChecked()) {
        if (color1 < 3) color2 = 3;
        else color1 = 3;
      }

      if (showYellowDiff->isChecked()) {
        if (color1 < 4) color2 = 4;
        else color1 = 4;
      }

      if (showMagentaDiff->isChecked()) {
        if (color1 < 5) color2 = 5;
        else color1 = 5;
      }

      if (showAquaDiff->isChecked()) {
        if (color1 < 6) color2 = 6;
        else color1 = 6;
      }

      get_new_state(3);

      /*peration_name
      // COMPLEMENT_OP = 0
      // INTERSECTION_OP = 1
      // UNION_OP = 2
      // DIFFERENCE_OP = 3
      // SYMMETRIC_dIFFERENCE_OP = 4
      // MINKOWSKI_SUM_OP = 5
      // RESET_OP = 6
      // COPY_OP = 7
      // MOVE_OP = 8
      // CLEAR_OP = 9
      // DELETEALL_OP = 10
      // START_OP = 11*/

      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      states_stack.back().result_polyline_sources().clear();

      if (!states_stack.back().active_set(m_color_active).is_empty()) {
        if (empty_warn) {
          show_not_empty_warning();
        }
      }

      switch (m_color_active) {
       case 0:
        if (color1 == 0) states_stack.back().result_set().assign(states_stack.back().blue_set());
        else if (color1 == 1) states_stack.back().result_set().assign(states_stack.back().red_set());
        else if (color1 == 2) states_stack.back().result_set().assign(states_stack.back().black_set());
        else if (color1 == 3) states_stack.back().result_set().assign(states_stack.back().brown_set());
        else if (color1 == 4) states_stack.back().result_set().assign(states_stack.back().yellow_set());
        else if (color1 == 5) states_stack.back().result_set().assign(states_stack.back().magenta_set());
        else if (color1 == 6) states_stack.back().result_set().assign(states_stack.back().aqua_set());

        if (color2 == 1) states_stack.back().result_set().difference(states_stack.back().red_set());
        else if (color2 == 2) states_stack.back().result_set().difference(states_stack.back().black_set());
        else if (color2 == 3) states_stack.back().result_set().difference(states_stack.back().brown_set());
        else if (color2 == 4) states_stack.back().result_set().difference(states_stack.back().yellow_set());
        else if (color2 == 5)
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        else if (color2 == 6) states_stack.back().result_set().difference(states_stack.back().aqua_set());

        states_stack.back().result_set().difference(states_stack.back().blue_set());
        states_stack.back().blue_set().join(states_stack.back().result_set());
        states_stack.back().result_set().clear();
        states_stack.back().result_linear_sources().clear();
        states_stack.back().result_circular_sources().clear();
        states_stack.back().result_polyline_sources().clear();
        states_stack.back().result_bezier_sources().clear();
        break;

       case 1:
        if (color1 == 0) states_stack.back().result_set().assign(states_stack.back().blue_set());
        else if (color1 == 1) states_stack.back().result_set().assign(states_stack.back().red_set());
        else if (color1 == 2) states_stack.back().result_set().assign(states_stack.back().black_set());
        else if (color1 == 3) states_stack.back().result_set().assign(states_stack.back().brown_set());
        else if (color1 == 4) states_stack.back().result_set().assign(states_stack.back().yellow_set());
        else if (color1 == 5) states_stack.back().result_set().assign(states_stack.back().magenta_set());

        if (color2 == 0) states_stack.back().result_set().difference(states_stack.back().blue_set());
        else if (color2 == 1) states_stack.back().result_set().difference(states_stack.back().red_set());
        else if (color2 == 2) states_stack.back().result_set().difference(states_stack.back().black_set());
        else if (color2 == 3) states_stack.back().result_set().difference(states_stack.back().brown_set());
        else if (color2 == 4) states_stack.back().result_set().difference(states_stack.back().yellow_set());
        else if (color2 == 5)
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        else if (color2 == 6) states_stack.back().result_set().difference(states_stack.back().aqua_set());
        states_stack.back().result_set().difference(states_stack.back().red_set());
        states_stack.back().red_set().join(states_stack.back().result_set());
        states_stack.back().result_set().clear();
        states_stack.back().result_linear_sources().clear();
        states_stack.back().result_circular_sources().clear();
        states_stack.back().result_polyline_sources().clear();
        states_stack.back().result_bezier_sources().clear();
        break;

       case 2:
        if (color1 == 0) states_stack.back().result_set().assign(states_stack.back().blue_set());
        else if (color1 == 1) states_stack.back().result_set().assign(states_stack.back().red_set());
        else if (color1 == 2) states_stack.back().result_set().assign(states_stack.back().black_set());
        else if (color1 == 3) states_stack.back().result_set().assign(states_stack.back().brown_set());
        else if (color1 == 4) states_stack.back().result_set().assign(states_stack.back().yellow_set());
        else if (color1 == 5) states_stack.back().result_set().assign(states_stack.back().magenta_set());

        if (color2 == 0) states_stack.back().result_set().difference(states_stack.back().blue_set());
        else if (color2 == 1) states_stack.back().result_set().difference(states_stack.back().red_set());
        else if (color2 == 2) states_stack.back().result_set().difference(states_stack.back().black_set());
        else if (color2 == 3) states_stack.back().result_set().difference(states_stack.back().brown_set());
        else if (color2 == 4) states_stack.back().result_set().difference(states_stack.back().yellow_set());
        else if (color2 == 5)
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        else if (color2 == 6) states_stack.back().result_set().difference(states_stack.back().aqua_set());
        states_stack.back().result_set().difference(states_stack.back().black_set());
        states_stack.back().black_set().join(states_stack.back().result_set());
        states_stack.back().result_set().clear();
        states_stack.back().result_linear_sources().clear();
        states_stack.back().result_circular_sources().clear();
        states_stack.back().result_polyline_sources().clear();
        states_stack.back().result_bezier_sources().clear();
        break;

       case 3:
        if (color1 == 0) states_stack.back().result_set().assign(states_stack.back().blue_set());
        else if (color1 == 1) states_stack.back().result_set().assign(states_stack.back().red_set());
        else if (color1 == 2) states_stack.back().result_set().assign(states_stack.back().black_set());
        else if (color1 == 3) states_stack.back().result_set().assign(states_stack.back().brown_set());
        else if (color1 == 4) states_stack.back().result_set().assign(states_stack.back().yellow_set());
        else if (color1 == 5) states_stack.back().result_set().assign(states_stack.back().magenta_set());

        if (color2 == 0) states_stack.back().result_set().difference(states_stack.back().blue_set());
        else if (color2 == 1) states_stack.back().result_set().difference(states_stack.back().red_set());
        else if (color2 == 2) states_stack.back().result_set().difference(states_stack.back().black_set());
        else if (color2 == 3) states_stack.back().result_set().difference(states_stack.back().brown_set());
        else if (color2 == 4) states_stack.back().result_set().difference(states_stack.back().yellow_set());
        else if (color2 == 5)
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        else if (color2 == 6) states_stack.back().result_set().difference(states_stack.back().aqua_set());
        states_stack.back().result_set().difference(states_stack.back().brown_set());
        states_stack.back().brown_set().join(states_stack.back().result_set());
        states_stack.back().result_set().clear();
        states_stack.back().result_linear_sources().clear();
        states_stack.back().result_circular_sources().clear();
        states_stack.back().result_polyline_sources().clear();
        states_stack.back().result_bezier_sources().clear();
        break;

       case 4:
        if (color1 == 0) states_stack.back().result_set().assign(states_stack.back().blue_set());
        else if (color1 == 1) states_stack.back().result_set().assign(states_stack.back().red_set());
        else if (color1 == 2) states_stack.back().result_set().assign(states_stack.back().black_set());
        else if (color1 == 3) states_stack.back().result_set().assign(states_stack.back().brown_set());
        else if (color1 == 4) states_stack.back().result_set().assign(states_stack.back().yellow_set());
        else if (color1 == 5) states_stack.back().result_set().assign(states_stack.back().magenta_set());

        if (color2 == 0) states_stack.back().result_set().difference(states_stack.back().blue_set());
        else if (color2 == 1) states_stack.back().result_set().difference(states_stack.back().red_set());
        else if (color2 == 2) states_stack.back().result_set().difference(states_stack.back().black_set());
        else if (color2 == 3) states_stack.back().result_set().difference(states_stack.back().brown_set());
        else if (color2 == 4) states_stack.back().result_set().difference(states_stack.back().yellow_set());
        else if (color2 == 5)
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        else if (color2 == 6) states_stack.back().result_set().difference(states_stack.back().aqua_set());
        states_stack.back().result_set().difference(states_stack.back().yellow_set());
        states_stack.back().yellow_set().join(states_stack.back().result_set());
        states_stack.back().result_set().clear();
        states_stack.back().result_linear_sources().clear();
        states_stack.back().result_circular_sources().clear();
        states_stack.back().result_polyline_sources().clear();
        states_stack.back().result_bezier_sources().clear();
        break;

       case 5:
        if (color1 == 0) states_stack.back().result_set().assign(states_stack.back().blue_set());
        else if (color1 == 1) states_stack.back().result_set().assign(states_stack.back().red_set());
        else if (color1 == 2) states_stack.back().result_set().assign(states_stack.back().black_set());
        else if (color1 == 3) states_stack.back().result_set().assign(states_stack.back().brown_set());
        else if (color1 == 4) states_stack.back().result_set().assign(states_stack.back().yellow_set());
        else if (color1 == 5) states_stack.back().result_set().assign(states_stack.back().magenta_set());

        if (color2 == 0) states_stack.back().result_set().difference(states_stack.back().blue_set());
        else if (color2 == 1) states_stack.back().result_set().difference(states_stack.back().red_set());
        else if (color2 == 2) states_stack.back().result_set().difference(states_stack.back().black_set());
        else if (color2 == 3) states_stack.back().result_set().difference(states_stack.back().brown_set());
        else if (color2 == 4) states_stack.back().result_set().difference(states_stack.back().yellow_set());
        else if (color2 == 5)
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        else if (color2 == 6) states_stack.back().result_set().difference(states_stack.back().aqua_set());
        states_stack.back().result_set().difference(states_stack.back().magenta_set());
        states_stack.back().magenta_set().join(states_stack.back().result_set());
        states_stack.back().result_set().clear();
        states_stack.back().result_linear_sources().clear();
        states_stack.back().result_circular_sources().clear();
        states_stack.back().result_polyline_sources().clear();
        states_stack.back().result_bezier_sources().clear();
        break;

       case 6:
        if (color1 == 0) states_stack.back().result_set().assign(states_stack.back().blue_set());
        else if (color1 == 1) states_stack.back().result_set().assign(states_stack.back().red_set());
        else if (color1 == 2) states_stack.back().result_set().assign(states_stack.back().black_set());
        else if (color1 == 3) states_stack.back().result_set().assign(states_stack.back().brown_set());
        else if (color1 == 4) states_stack.back().result_set().assign(states_stack.back().yellow_set());
        else if (color1 == 5) states_stack.back().result_set().assign(states_stack.back().magenta_set());

        if (color2 == 0) states_stack.back().result_set().difference(states_stack.back().blue_set());
        else if (color2 == 1) states_stack.back().result_set().difference(states_stack.back().red_set());
        else if (color2 == 2) states_stack.back().result_set().difference(states_stack.back().black_set());
        else if (color2 == 3) states_stack.back().result_set().difference(states_stack.back().brown_set());
        else if (color2 == 4) states_stack.back().result_set().difference(states_stack.back().yellow_set());
        else if (color2 == 5)
          states_stack.back().result_set().difference(states_stack.back().magenta_set());
        else if (color2 == 6) states_stack.back().result_set().difference(states_stack.back().aqua_set());
        states_stack.back().result_set().difference(states_stack.back().aqua_set());
        states_stack.back().aqua_set().join(states_stack.back().result_set());
        states_stack.back().result_set().clear();
        states_stack.back().result_linear_sources().clear();
        states_stack.back().result_circular_sources().clear();
        states_stack.back().result_polyline_sources().clear();
        states_stack.back().result_bezier_sources().clear();
        break;
      }
      lDone = true;
    } else {
      ask_user_ok("Difference Operation Error", "Operation valid for 2 colored polygon set\n");
    }
    actionDifferenceH->setChecked(false);
    this->setCursor(old);
    if (lDone) modelChanged();
  }
}

//symmetric difference
void MainWindow::on_actionSymmetric_DifferenceH_toggled(bool aChecked) {
  if (actionSymmetric_DifferenceH->isChecked()) {
    bool lDone = false;
    QCursor old = this->cursor();
    this->setCursor(Qt::WaitCursor);

    actionComplementH->setChecked(false);
    actionUnionH->setChecked(false);
    actionIntersectionH->setChecked(false);
    actionDifferenceH->setChecked(false);
    actionMinkowski_SumH->setChecked(false);
    actionCopyH->setChecked(false);
    actionMoveH->setChecked(false);

    get_new_state(4);

    /*operation_name
    // COMPLEMENT_OP = 0
    // INTERSECTION_OP = 1
    // UNION_OP = 2
    // DIFFERENCE_OP = 3
    // SYMMETRIC_dIFFERENCE_OP = 4
    // MINKOWSKI_SUM_OP = 5
    // RESET_OP = 6
    // COPY_OP = 7
    // MOVE_OP = 8
    // CLEAR_OP = 9
    // DELETEALL_OP = 10
    // START_OP = 11*/

    states_stack.back().result_set().clear();
    states_stack.back().result_linear_sources().clear();
    states_stack.back().result_circular_sources().clear();
    states_stack.back().result_bezier_sources().clear();
    states_stack.back().result_polyline_sources().clear();

    if (!states_stack.back().active_set(m_color_active).is_empty()) {
      if (empty_warn) {
        show_not_empty_warning();
      }
    }
    switch (m_color_active) {
     case 0:
      if (!states_stack.back().blue_set().is_empty() && m_blue_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().red_set());
      if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().black_set());
      if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().brown_set());
      if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().yellow_set());
      if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().magenta_set());
      if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().blue_set());
      states_stack.back().blue_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;

     case 1:
      if (!states_stack.back().blue_set().is_empty() && m_blue_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().red_set());
      if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().black_set());
      if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().brown_set());
      if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().yellow_set());
      if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().magenta_set());
      if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().red_set());
      states_stack.back().red_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;

     case 2:
      if (!states_stack.back().blue_set().is_empty() && m_blue_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().red_set());
      if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().black_set());
      if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().brown_set());
      if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().yellow_set());
      if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().magenta_set());
      if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().black_set());
      states_stack.back().black_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;

     case 3:
      if (!states_stack.back().blue_set().is_empty() && m_blue_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().red_set());
      if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().black_set());
      if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().brown_set());
      if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().yellow_set());
      if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().magenta_set());
      if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().aqua_set());\

      states_stack.back().result_set().difference(states_stack.back().brown_set());
      states_stack.back().brown_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;

     case 4:
      if (!states_stack.back().blue_set().is_empty() && m_blue_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().red_set());
      if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().black_set());
      if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().brown_set());
      if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().yellow_set());
      if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().magenta_set());
      if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().yellow_set());
      states_stack.back().yellow_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;

     case 5:
      if (!states_stack.back().blue_set().is_empty() && m_blue_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().red_set());
      if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().black_set());
      if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().brown_set());
      if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().yellow_set());
      if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().magenta_set());
      if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().magenta_set());
      states_stack.back().magenta_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;

     case 6:
      if (!states_stack.back().blue_set().is_empty() && m_blue_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      else if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().red_set());
      else if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().black_set());
      else if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      else if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      else if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      else if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().assign(states_stack.back().aqua_set());

      if (!states_stack.back().red_set().is_empty() && m_red_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().red_set());
      if (!states_stack.back().black_set().is_empty() && m_black_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().black_set());
      if (!states_stack.back().brown_set().is_empty() && m_brown_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().brown_set());
      if (!states_stack.back().yellow_set().is_empty() && m_yellow_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().yellow_set());
      if (!states_stack.back().magenta_set().is_empty() && m_magenta_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().magenta_set());
      if (!states_stack.back().aqua_set().is_empty() && m_aqua_sym_diff)
        states_stack.back().result_set().symmetric_difference(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().aqua_set());
      states_stack.back().aqua_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;
    }

    actionSymmetric_DifferenceH->setChecked(false);
    lDone = true;
    this->setCursor(old);
    if (lDone) modelChanged();
  }
}

//union
void MainWindow::on_actionUnionH_toggled(bool aChecked) {
  if (actionUnionH->isChecked()) {
    bool lDone = false;
    QCursor old = this->cursor();
    this->setCursor(Qt::WaitCursor);

    actionComplementH->setChecked(false);
    //actionUnion->setChecked(false);
    actionIntersectionH->setChecked(false);
    actionDifferenceH->setChecked(false);
    actionSymmetric_DifferenceH->setChecked(false);
    actionMinkowski_SumH->setChecked(false);
    actionCopyH->setChecked(false);
    actionMoveH->setChecked(false);

    get_new_state(2);

    /*operation_name
    // COMPLEMENT_OP = 0
    // INTERSECTION_OP = 1
    // UNION_OP = 2
    // DIFFERENCE_OP = 3
    // SYMMETRIC_dIFFERENCE_OP = 4
    // MINKOWSKI_SUM_OP = 5
    // RESET_OP = 6
    // COPY_OP = 7
    // MOVE_OP = 8
    // CLEAR_OP = 9
    // DELETEALL_OP = 10
    // START_OP = 11*/

    states_stack.back().result_set().clear();
    states_stack.back().result_linear_sources().clear();
    states_stack.back().result_circular_sources().clear();
    states_stack.back().result_bezier_sources().clear();
    states_stack.back().result_polyline_sources().clear();

    if (!states_stack.back().active_set(m_color_active).is_empty()) {
      if (empty_warn) {
        show_not_empty_warning();
      }
    }

    switch (m_color_active) {
     case 0:if (m_red_union) states_stack.back().result_set().assign(states_stack.back().red_set());
      if (m_blue_union) states_stack.back().result_set().join(states_stack.back().blue_set());
      if (m_black_union) states_stack.back().result_set().join(states_stack.back().black_set());
      if (m_brown_union) states_stack.back().result_set().join(states_stack.back().brown_set());
      if (m_magenta_union) states_stack.back().result_set().join(states_stack.back().magenta_set());
      if (m_yellow_union) states_stack.back().result_set().join(states_stack.back().yellow_set());
      if (m_aqua_union) states_stack.back().result_set().join(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().blue_set());
      states_stack.back().blue_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;

     case 1:if (m_red_union) states_stack.back().result_set().assign(states_stack.back().red_set());
      if (m_blue_union) states_stack.back().result_set().join(states_stack.back().blue_set());
      if (m_black_union) states_stack.back().result_set().join(states_stack.back().black_set());
      if (m_brown_union) states_stack.back().result_set().join(states_stack.back().brown_set());
      if (m_magenta_union) states_stack.back().result_set().join(states_stack.back().magenta_set());
      if (m_yellow_union) states_stack.back().result_set().join(states_stack.back().yellow_set());
      if (m_aqua_union) states_stack.back().result_set().join(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().red_set());
      states_stack.back().red_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;

     case 2:if (m_red_union) states_stack.back().result_set().assign(states_stack.back().red_set());
      if (m_blue_union) states_stack.back().result_set().join(states_stack.back().blue_set());
      if (m_black_union) states_stack.back().result_set().join(states_stack.back().black_set());
      if (m_brown_union) states_stack.back().result_set().join(states_stack.back().brown_set());
      if (m_magenta_union) states_stack.back().result_set().join(states_stack.back().magenta_set());
      if (m_yellow_union) states_stack.back().result_set().join(states_stack.back().yellow_set());
      if (m_aqua_union) states_stack.back().result_set().join(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().black_set());
      states_stack.back().black_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;

     case 3:if (m_red_union) states_stack.back().result_set().assign(states_stack.back().red_set());
      if (m_blue_union) states_stack.back().result_set().join(states_stack.back().blue_set());
      if (m_black_union) states_stack.back().result_set().join(states_stack.back().black_set());
      if (m_brown_union) states_stack.back().result_set().join(states_stack.back().brown_set());
      if (m_magenta_union) states_stack.back().result_set().join(states_stack.back().magenta_set());
      if (m_yellow_union) states_stack.back().result_set().join(states_stack.back().yellow_set());
      if (m_aqua_union) states_stack.back().result_set().join(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().brown_set());
      states_stack.back().brown_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;

     case 4:if (m_red_union) states_stack.back().result_set().assign(states_stack.back().red_set());
      if (m_blue_union) states_stack.back().result_set().join(states_stack.back().blue_set());
      if (m_black_union) states_stack.back().result_set().join(states_stack.back().black_set());
      if (m_brown_union) states_stack.back().result_set().join(states_stack.back().brown_set());
      if (m_magenta_union) states_stack.back().result_set().join(states_stack.back().magenta_set());
      if (m_yellow_union) states_stack.back().result_set().join(states_stack.back().yellow_set());
      if (m_aqua_union) states_stack.back().result_set().join(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().yellow_set());
      states_stack.back().yellow_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;

     case 5:if (m_red_union) states_stack.back().result_set().assign(states_stack.back().red_set());
      if (m_blue_union) states_stack.back().result_set().join(states_stack.back().blue_set());
      if (m_black_union) states_stack.back().result_set().join(states_stack.back().black_set());
      if (m_brown_union) states_stack.back().result_set().join(states_stack.back().brown_set());
      if (m_magenta_union) states_stack.back().result_set().join(states_stack.back().magenta_set());
      if (m_yellow_union) states_stack.back().result_set().join(states_stack.back().yellow_set());
      if (m_aqua_union) states_stack.back().result_set().join(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().magenta_set());
      states_stack.back().magenta_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;

     case 6:if (m_red_union) states_stack.back().result_set().assign(states_stack.back().red_set());
      if (m_blue_union) states_stack.back().result_set().join(states_stack.back().blue_set());
      if (m_black_union) states_stack.back().result_set().join(states_stack.back().black_set());
      if (m_brown_union) states_stack.back().result_set().join(states_stack.back().brown_set());
      if (m_magenta_union) states_stack.back().result_set().join(states_stack.back().magenta_set());
      if (m_yellow_union) states_stack.back().result_set().join(states_stack.back().yellow_set());
      if (m_aqua_union) states_stack.back().result_set().join(states_stack.back().aqua_set());

      states_stack.back().result_set().difference(states_stack.back().aqua_set());
      states_stack.back().aqua_set().join(states_stack.back().result_set());
      states_stack.back().result_set().clear();
      states_stack.back().result_linear_sources().clear();
      states_stack.back().result_circular_sources().clear();
      states_stack.back().result_polyline_sources().clear();
      states_stack.back().result_bezier_sources().clear();
      break;
    }

    actionUnionH->setChecked(false);
    lDone = true;
    this->setCursor(old);
    if (lDone) modelChanged();
  }
}

//copy
void MainWindow::on_actionCopyH_toggled(bool aChecked) {
  if (actionCopyH->isChecked()) {
    bool lDone = false;
    QCursor old = this->cursor();
    this->setCursor(Qt::WaitCursor);

    actionUnionH->setChecked(false);
    actionIntersectionH->setChecked(false);
    actionDifferenceH->setChecked(false);
    actionSymmetric_DifferenceH->setChecked(false);
    actionMinkowski_SumH->setChecked(false);
    actionMoveH->setChecked(false);

    get_new_state(7);

    /*operation_name
    // COMPLEMENT_OP = 0
    // INTERSECTION_OP = 1
    // UNION_OP = 2
    // DIFFERENCE_OP = 3
    // SYMMETRIC_dIFFERENCE_OP = 4
    // MINKOWSKI_SUM_OP = 5
    // RESET_OP = 6
    // COPY_OP = 7
    // MOVE_OP = 8
    // CLEAR_OP = 9
    // DELETEALL_OP = 10
    // START_OP = 11*/

    states_stack.back().result_set().clear();
    states_stack.back().result_linear_sources().clear();
    states_stack.back().result_circular_sources().clear();
    states_stack.back().result_bezier_sources().clear();
    states_stack.back().result_polyline_sources().clear();

    switch (m_color_copy) {
     case 0:
      if (!states_stack.back().blue_set().is_empty())
        states_stack.back().result_set().assign(states_stack.back().blue_set());
      break;
     case 1:
      if (!states_stack.back().red_set().is_empty())
        states_stack.back().result_set().assign(states_stack.back().red_set());
      break;
     case 2:
      if (!states_stack.back().black_set().is_empty())
        states_stack.back().result_set().assign(states_stack.back().black_set());
      break;
     case 3:
      if (!states_stack.back().brown_set().is_empty())
        states_stack.back().result_set().assign(states_stack.back().brown_set());
      break;
     case 4:
      if (!states_stack.back().yellow_set().is_empty())
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
      break;
     case 5:
      if (!states_stack.back().magenta_set().is_empty())
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
      break;
     case 6:
      if (!states_stack.back().aqua_set().is_empty())
        states_stack.back().result_set().assign(states_stack.back().aqua_set());
      break;
     default:show_warning("Please select a Bucket!!!");
      break;
    }

    actionCopyH->setChecked(false);
    lDone = true;
    m_color_cm = 0; //copy
    if (!states_stack.back().result_set().is_empty())
      states_stack.back().active_set(m_color_active).join(states_stack.back().result_set());
    states_stack.back().result_set().clear();
    states_stack.back().result_linear_sources().clear();
    states_stack.back().result_circular_sources().clear();
    states_stack.back().result_bezier_sources().clear();
    states_stack.back().result_polyline_sources().clear();
    this->setCursor(old);
    if (lDone) modelChanged();
  }
}

//move
void MainWindow::on_actionMoveH_toggled(bool aChecked) {
  if (actionMoveH->isChecked()) {
    bool lDone = false;
    QCursor old = this->cursor();
    this->setCursor(Qt::WaitCursor);

    actionUnionH->setChecked(false);
    actionIntersectionH->setChecked(false);
    actionDifferenceH->setChecked(false);
    actionSymmetric_DifferenceH->setChecked(false);
    actionMinkowski_SumH->setChecked(false);
    actionCopyH->setChecked(false);

    get_new_state(8);

    /*operation_name
    // COMPLEMENT_OP = 0
    // INTERSECTION_OP = 1
    // UNION_OP = 2
    // DIFFERENCE_OP = 3
    // SYMMETRIC_dIFFERENCE_OP = 4
    // MINKOWSKI_SUM_OP = 5
    // RESET_OP = 6
    // COPY_OP = 7
    // MOVE_OP = 8
    // CLEAR_OP = 9
    // DELETEALL_OP = 10
    // START_OP = 11*/

    states_stack.back().result_set().clear();
    states_stack.back().result_linear_sources().clear();
    states_stack.back().result_circular_sources().clear();
    states_stack.back().result_bezier_sources().clear();
    states_stack.back().result_polyline_sources().clear();

    switch (m_color_move_ext) {
     case 0:
      if (!states_stack.back().blue_set().is_empty()) {
        states_stack.back().result_set().assign(states_stack.back().blue_set());
        states_stack.back().active_set(m_color_active).join(states_stack.back().result_set());
        states_stack.back().blue_set().clear();
      }
      break;
     case 1:
      if (!states_stack.back().red_set().is_empty()) {
        states_stack.back().result_set().assign(states_stack.back().red_set());
        states_stack.back().active_set(m_color_active).join(states_stack.back().result_set());
        states_stack.back().red_set().clear();
      }
      break;
     case 2:
      if (!states_stack.back().black_set().is_empty()) {
        states_stack.back().result_set().assign(states_stack.back().black_set());
        states_stack.back().active_set(m_color_active).join(states_stack.back().result_set());
        states_stack.back().black_set().clear();
      }
      break;
     case 3:
      if (!states_stack.back().brown_set().is_empty()) {
        states_stack.back().result_set().assign(states_stack.back().brown_set());
        states_stack.back().active_set(m_color_active).join(states_stack.back().result_set());
        states_stack.back().brown_set().clear();
      }
      break;
     case 4:
      if (!states_stack.back().yellow_set().is_empty()) {
        states_stack.back().result_set().assign(states_stack.back().yellow_set());
        states_stack.back().active_set(m_color_active).join(states_stack.back().result_set());
        states_stack.back().yellow_set().clear();
      }
      break;
     case 5:
      if (!states_stack.back().magenta_set().is_empty()) {
        states_stack.back().result_set().assign(states_stack.back().magenta_set());
        states_stack.back().active_set(m_color_active).join(states_stack.back().result_set());
        states_stack.back().magenta_set().clear();
      }
      break;
     case 6:
      if (!states_stack.back().aqua_set().is_empty()) {
        states_stack.back().result_set().assign(states_stack.back().aqua_set());
        states_stack.back().active_set(m_color_active).join(states_stack.back().result_set());
        states_stack.back().aqua_set().clear();
      }
      break;
     default:show_warning("Please select a Bucket!!!");
      break;
    }

    actionMoveH->setChecked(false);
    states_stack.back().result_set().clear();
    states_stack.back().result_linear_sources().clear();
    states_stack.back().result_circular_sources().clear();
    states_stack.back().result_bezier_sources().clear();
    states_stack.back().result_polyline_sources().clear();

    lDone = true;
    m_color_cm = 1; //copy
    this->setCursor(old);
    if (lDone) modelChanged();
  }
}

//minkowski polygon(only for linear polygon)
Polygon_with_holes_2 MainWindow::getMinkInputPolygon(size_t color) {
  Linear_polygon_set lps;
  typedef std::vector <Linear_polygon_with_holes> Pgn_with_holes_container;
  typedef typename Kernel::Point_2 Point;
  std::list <Polygon_2> pgns;

  switch (color) {
   case 0:lps = states_stack.back().blue_set().linear();
    break;
   case 1:lps = states_stack.back().red_set().linear();
    break;
   case 2:lps = states_stack.back().black_set().linear();
    break;
   case 3:lps = states_stack.back().brown_set().linear();
    break;
   case 4:lps = states_stack.back().yellow_set().linear();
    break;
   case 5:lps = states_stack.back().magenta_set().linear();
    break;
   case 6:lps = states_stack.back().aqua_set().linear();
    break;
  }

  Pgn_with_holes_container res(lps.number_of_polygons_with_holes());
  res.begin() = lps.polygons_with_holes(res.begin());
  typename Pgn_with_holes_container::const_iterator it;

  Linear_polygon_with_holes temp;
  Polygon_with_holes_2 pwh;
  Polygon_2 p;
  //show_warning("Number of polygons with holes: "+to_string(lps.number_of_polygons_with_holes()));

  if (lps.number_of_polygons_with_holes() > 1) {
    m_disjoint = true;
  }

  if (!m_disjoint) {
    for (it = res.begin(); it != res.end(); ++it) {
      temp = *it;
      for (auto vit = temp.outer_boundary().curves_begin(); vit != temp.outer_boundary().curves_end(); vit++) {
        auto pt = vit->source();
        //show_warning(to_string(CGAL::to_double(pt.x()))+" , "+to_string(CGAL::to_double(pt.y())));
        p.push_back(Point(CGAL::to_double(pt.x()), CGAL::to_double(pt.y())));
      }
      if (p.orientation() == CGAL::CLOCKWISE) {
        p.reverse_orientation();
      }

      pgns.push_back(p);

      for (auto hi = temp.holes_begin(); hi != temp.holes_end(); ++hi) {
        Polygon_2 hole;
        auto vti = *hi;
        for (auto hit = vti.curves_begin(); hit != vti.curves_end(); hit++) {
          auto pti = hit->source();
          hole.push_back(Point(CGAL::to_double(pti.x()),
                               CGAL::to_double(pti.y())));
        }
        if (hole.orientation() == CGAL::COUNTERCLOCKWISE) {
          hole.reverse_orientation();
        }
        pgns.push_back(hole);
      }
    }
    std::list<Polygon_2>::iterator pit = pgns.begin();
    return Polygon_with_holes_2(pgns.front(), ++pit, pgns.end());
  }
  return pwh;
}

//minkowski button clicked
void MainWindow::on_actionMinkowski_SumH_toggled(bool aChecked) {
  if (actionMinkowski_SumH->isChecked()) {
    bool lDone = false;
    QCursor old = this->cursor();
    this->setCursor(Qt::WaitCursor);

    if (!m_circular_active && !m_bezier_active && !m_polyline_active) {
      actionComplementH->setChecked(false);
      actionUnionH->setChecked(false);
      actionIntersectionH->setChecked(false);
      actionDifferenceH->setChecked(false);
      actionSymmetric_DifferenceH->setChecked(false);

      get_new_state(5);

      /*operation_name
      // COMPLEMENT_OP = 0
      // INTERSECTION_OP = 1
      // UNION_OP = 2
      // DIFFERENCE_OP = 3
      // SYMMETRIC_dIFFERENCE_OP = 4
      // MINKOWSKI_SUM_OP = 5
      // RESET_OP = 6
      // COPY_OP = 7
      // MOVE_OP = 8
      // CLEAR_OP = 9
      // DELETEALL_OP = 10
      // START_OP = 11*/

      if (!states_stack.back().active_set(m_color_active).is_empty()) {
        if (empty_warn) {
          show_not_empty_warning();
        }
      }

      size_t count = 0;
      if (showBlueMink_Sum->isChecked()) count++;
      if (showRedMink_Sum->isChecked()) count++;
      if (showBlackMink_Sum->isChecked()) count++;
      if (showBrownMink_Sum->isChecked()) count++;
      if (showYellowMink_Sum->isChecked()) count++;
      if (showMagentaMink_Sum->isChecked()) count++;
      if (showAquaMink_Sum->isChecked()) count++;

      if (count == 2) {
        size_t color1 = 111;
        size_t color2 = 1111;

        if (showBlueMink_Sum->isChecked()) color1 = 0;
        if (showRedMink_Sum->isChecked()) {
          if (color1 < 1) color2 = 1;
          else color1 = 1;
        }

        if (showBlackMink_Sum->isChecked()) {
          if (color1 < 2) color2 = 2;
          else color1 = 2;
        }

        if (showBrownMink_Sum->isChecked()) {
          if (color1 < 3) color2 = 3;
          else color1 = 3;
        }

        if (showYellowMink_Sum->isChecked()) {
          if (color1 < 4) color2 = 4;
          else color1 = 4;
        }

        if (showMagentaMink_Sum->isChecked()) {
          if (color1 < 5) color2 = 5;
          else color1 = 5;
        }

        if (showAquaMink_Sum->isChecked()) {
          color2 = 6;
        }

        typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
        typedef Kernel::Point_2 Point_2;

        Polygon_with_holes_2 p0, p1;

        // Compute the Minkowski sum using the decomposition approach.
        CGAL::Small_side_angle_bisector_decomposition_2 <Kernel> ssab_decomp;

        if (color1 == 0 && !states_stack.back().blue_set().is_empty()) {
          p0 = getMinkInputPolygon(color1);
          if (color2 == 1 && !states_stack.back().red_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
          else if (color2 == 2 && !states_stack.back().black_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
          else if (color2 == 3 && !states_stack.back().brown_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
          else if (color2 == 4 && !states_stack.back().yellow_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
          else if (color2 == 5 && !states_stack.back().magenta_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
          else if (color2 == 6 && !states_stack.back().aqua_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
        }
        else if (color1 == 1 && !states_stack.back().red_set().is_empty()) {
          p0 = getMinkInputPolygon(color1);
          if (color2 == 2 && !states_stack.back().black_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
          else if (color2 == 3 && !states_stack.back().brown_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
          else if (color2 == 4 && !states_stack.back().yellow_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
          else if (color2 == 5 && !states_stack.back().magenta_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
          else if (color2 == 6 && !states_stack.back().aqua_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
        }
        else if (color1 == 2 && !states_stack.back().black_set().is_empty()) {
          p0 = getMinkInputPolygon(color1);
          if (color2 == 3 && !states_stack.back().brown_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          } else if (color2 == 4 && !states_stack.back().yellow_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          } else if (color2 == 5 && !states_stack.back().magenta_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          } else if (color2 == 6 && !states_stack.back().aqua_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
        } else if (color1 == 3 && !states_stack.back().brown_set().is_empty()) {
          p0 = getMinkInputPolygon(color1);
          if (color2 == 4 && !states_stack.back().yellow_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          } else if (color2 == 5 && !states_stack.back().magenta_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          } else if (color2 == 6 && !states_stack.back().aqua_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
        } else if (color1 == 4 && !states_stack.back().yellow_set().is_empty()) {
          p0 = getMinkInputPolygon(color1);
          if (color2 == 5 && !states_stack.back().magenta_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          } else if (color2 == 6 && !states_stack.back().aqua_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
        } else if (color1 == 5 && !states_stack.back().magenta_set().is_empty()) {
          p0 = getMinkInputPolygon(color1);
          if (color2 == 6 && !states_stack.back().aqua_set().is_empty()) {
            p1 = getMinkInputPolygon(color2);
          }
        }

        if (! m_disjoint) {
          Polygon_with_holes_2 mink_sum_res = CGAL::minkowski_sum_2(p0, p1);
          if (!mink_sum_res.is_unbounded()) {
            m_linear_input->m_is_mink = true;
            m_linear_input->get_Minkowski_result(mink_sum_res, p0);
            m_linear_input->m_hole = true;
            m_linear_input->get_Minkowski_holes(mink_sum_res, p0);
            states_stack.back().result_set().clear();
            states_stack.back().result_linear_sources().clear();
            m_linear_input->m_hole = false;
            m_linear_input->m_is_mink = false;
          } else ask_user_ok("Minkowski Sum Operation Error", "resultant polygon is unbounded\n");
          lDone = true;
          minkowski_sum_operated = true;
        }
        else {
          on_actionUndo_triggered();
          show_error(
                     "DISJOINT POLYGON SET ERROR\n\nCannot perform Minkowski Sum operation since Input Bucket contains more than disjoint polygons!!!");
          m_disjoint = false;
        }
      }
      else {
        ask_user_ok("Minkowski Sum Operation Error", "Function supports 2 polygon as input\n");
      }
    }

    actionMinkowski_SumH->setChecked(false);
    this->setCursor(old);
    if (lDone) {
      //zoomToFit();
      modelChanged();
    }
  }
}

//axis button clicked
void MainWindow::on_actionAxis_triggered() {
  if (! m_grid) {
    m_scene.addItem(xAxis);
    m_scene.addItem(yAxis);
    m_grid = true;
  }
  else {
    m_scene.removeItem(xAxis);
    m_scene.removeItem(yAxis);
    m_grid = false;
  }
  actionAxis->setChecked(false);
  m_scene.update();
  modelChanged();
}

//pan button
void MainWindow::on_actionPAN_triggered() {
  if (! m_pan) {
    //get_new_state(7);
    if (!m_circular_active && !m_bezier_active && !m_polyline_active) {
      m_scene.removeEventFilter(m_linear_input);
      m_scene.removeEventFilter(m_bezier_input);
      m_scene.removeEventFilter(m_circular_input);
      m_scene.removeEventFilter(m_polyline_input);
      m_linear_input->Reset();
      m_circular_input->Reset();
      m_bezier_input->Reset();
      m_polyline_input->Reset();
      actionInsertLinear->setChecked(false);
      this->graphicsView->setDragMode(QGraphicsView::ScrollHandDrag);
    }
    else if (!m_bezier_active && !m_polyline_active) {
      m_scene.removeEventFilter(m_linear_input);
      m_scene.removeEventFilter(m_bezier_input);
      m_scene.removeEventFilter(m_circular_input);
      m_scene.removeEventFilter(m_polyline_input);
      m_linear_input->Reset();
      m_circular_input->Reset();
      m_bezier_input->Reset();
      m_polyline_input->Reset();
      actionInsertCircular->setChecked(false);
      this->graphicsView->setDragMode(QGraphicsView::ScrollHandDrag);
    }
    else if (!m_polyline_active) {
      m_scene.removeEventFilter(m_linear_input);
      m_scene.removeEventFilter(m_bezier_input);
      m_scene.removeEventFilter(m_circular_input);
      m_scene.removeEventFilter(m_polyline_input);
      m_linear_input->Reset();
      m_circular_input->Reset();
      m_bezier_input->Reset();
      m_polyline_input->Reset();
      actionInsertBezier->setChecked(false);
      this->graphicsView->setDragMode(QGraphicsView::ScrollHandDrag);
    }
    else {
      m_scene.removeEventFilter(m_linear_input);
      m_scene.removeEventFilter(m_bezier_input);
      m_scene.removeEventFilter(m_circular_input);
      m_scene.removeEventFilter(m_polyline_input);
      m_linear_input->Reset();
      m_circular_input->Reset();
      m_bezier_input->Reset();
      m_polyline_input->Reset();
      actionInsertPolyline->setChecked(false);
      this->graphicsView->setDragMode(QGraphicsView::ScrollHandDrag);
    }

    this->graphicsView->setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
    m_scene.update();
    m_pan = true;
    modelChanged();
  }
  else {
    if (! (m_circular_active || m_bezier_active) || m_polyline_active) {
      actionInsertLinear->setChecked(true);
    }
    else if (m_circular_active) {
      actionInsertCircular->setChecked(true);
    }
    else if (m_bezier_active) {
      actionInsertBezier->setChecked(true);
    }
    else if (m_polyline_active) {
      actionInsertPolyline->setChecked(true);
    }
    m_pan = false;
    modelChanged();
  }
}

//add colors(represented by wheel)
void MainWindow::wheelEvent(QWheelEvent *event) {

  this->graphicsView->setTransformationAnchor(QGraphicsView::AnchorViewCenter);
  static const double scaleFactor = 1.15;
  static double currentScale = 1.0;  // stores the current scale value.
  static const double scaleMin = .1; // defines the min scale limit.

  if (m_scene.width() >= 1500000 && event->angleDelta().y() > 0) {
    // currentScale = 0.0;
    this->graphicsView->setDragMode(QGraphicsView::NoDrag);
    this->graphicsView->setSceneRect(-1500000, -1000000, 3100000, 2000000);
    this->graphicsView->fitInView(m_scene.sceneRect());
  }
  if (m_scene.width() <= 1500000) {
    if (event->angleDelta().y() > 0) {
      m_scene.setSceneRect(-m_scene.sceneRect().width() * scaleFactor / 2,
                           -m_scene.sceneRect().height() * scaleFactor / 2,
                           m_scene.sceneRect().width() * scaleFactor,
                           m_scene.sceneRect().height() * scaleFactor);
      this->graphicsView->fitInView(m_scene.sceneRect());
      this->graphicsView->scale(2.5, 2.5);
      // show_warning(to_string(m_scene.sceneRect().height())+" "+to_string(m_scene.sceneRect().width()));
    } else if (currentScale > scaleMin) {
      m_scene.setSceneRect(-m_scene.sceneRect().width() / (scaleFactor * 2),
                           -m_scene.sceneRect().height() / (scaleFactor * 2),
                           m_scene.sceneRect().width() / scaleFactor,
                           m_scene.sceneRect().height() / scaleFactor);
      this->graphicsView->fitInView(m_scene.sceneRect());
      this->graphicsView->scale(2.5, 2.5);
    }
  }
  modelChanged();
}

//zoom to fit button
void MainWindow::zoomToFit() {
    std::optional <QRectF> lTotalRect;
    for (auto si = states_stack.back().m_curve_sets.begin(); si != states_stack.back().m_curve_sets.end(); ++si) {
        if (!si->is_empty()) {
            QRectF lRect = si->bounding_rect();
            if (lTotalRect) lTotalRect = *lTotalRect | lRect;
            else lTotalRect = lRect;
        }
    }
    if (lTotalRect) {
        this->graphicsView->setSceneRect(*lTotalRect);
        this->graphicsView->fitInView(*lTotalRect, Qt::KeepAspectRatio);
    }
    modelChanged();
}

//function to detect whether the bucket is occupied or not and throw warning
void MainWindow::show_not_empty_warning() {
  QMessageBox msgBox;
  msgBox.setText("Selected Bucket is Not Empty!!!");
  msgBox.setIcon(QMessageBox::Warning);
  msgBox.addButton(QMessageBox::Ok);
  QCheckBox dontShowCheckBox("don't show this message again");
  dontShowCheckBox.blockSignals(true);
  msgBox.addButton(&dontShowCheckBox, QMessageBox::ResetRole);
  int32_t userReply = msgBox.exec();
  if (dontShowCheckBox.checkState() == Qt::Checked) {
    empty_warn = false;
  } else {
    empty_warn = true;
  }
}

//resize button
void MainWindow::resizeEvent(QResizeEvent *e) {
  this->graphicsView->fitInView(m_scene.sceneRect());
  this->graphicsView->scale(2.5, 2.5);
  modelChanged();
}

//all the information in the help/about section
void MainWindow::on_aboutCGAL_triggered() {
  QString title = "About CGAL";
  QString html_resource_name = "help/about_CGAL.html";

  QFile about_CGAL(html_resource_name);
  about_CGAL.open(QIODevice::ReadOnly);
  QString about_CGAL_txt = QTextStream(&about_CGAL).readAll();
#ifdef CGAL_VERSION_STR
  QString cgal_version(CGAL_VERSION_STR);
#  ifdef CGAL_FAKE_PUBLIC_RELEASE
  cgal_version.replace(QRegExp("-Ic?.*"), "");
#  endif
  about_CGAL_txt.replace("<!--CGAL_VERSION-->",
                         QString(" (version %1)")
                         .arg(cgal_version));
#endif
  QMessageBox mb;
  mb.setStyleSheet("QLabel{min-width: 550px; min-height: 250px;}");
  mb.setWindowTitle(title);
  mb.setText(about_CGAL_txt);
  mb.addButton(QMessageBox::Ok);

  QLabel *mb_label = mb.findChild<QLabel *>("qt_msgbox_label");
  if (mb_label) {
    mb_label->setTextInteractionFlags(mb_label->textInteractionFlags() |
                                      ::Qt::LinksAccessibleByMouse |
                                      ::Qt::LinksAccessibleByKeyboard);
  } else {
    std::cerr << "Cannot find child \"qt_msgbox_label\" in QMessageBox\n"
              << "  with Qt version " << QT_VERSION_STR << "!\n";
  }
  mb.exec();
}

void MainWindow::on_aboutMinkop_triggered() {
  QString title = "About 2D Minkowski sums";
  QString html_resource_name = "help/About2DMinkowskisums.html";
  QFile about_CGAL(html_resource_name);
  about_CGAL.open(QIODevice::ReadOnly);
  QString about_CGAL_txt = QTextStream(&about_CGAL).readAll();
#ifdef CGAL_VERSION_STR
  QString cgal_version(CGAL_VERSION_STR);
#  ifdef CGAL_FAKE_PUBLIC_RELEASE
  cgal_version.replace(QRegExp("-Ic?.*"), "");
#  endif
  about_CGAL_txt.replace("<!--CGAL_VERSION-->",
                         QString(" (version %1)")
                         .arg(cgal_version));
#endif
  QMessageBox mb;
  mb.setStyleSheet("QLabel{min-width: 700px; min-height: 300px;}");
  mb.setWindowTitle(title);
  mb.setText(about_CGAL_txt);
  mb.addButton(QMessageBox::Ok);

  QLabel *mb_label = mb.findChild<QLabel *>("qt_msgbox_label");
  if (mb_label) {
    mb_label->setTextInteractionFlags(mb_label->textInteractionFlags() |
                                      ::Qt::LinksAccessibleByMouse |
                                      ::Qt::LinksAccessibleByKeyboard);
  }
  else {
    std::cerr << "Cannot find child \"qt_msgbox_label\" in QMessageBox\n"
              << "  with Qt version " << QT_VERSION_STR << "!\n";
  }
  mb.exec();
}

void MainWindow::on_aboutBolop_triggered() {

  QString title = "About 2D Minkowski sums";
  QString html_resource_name = "help/About2DRegularizedBooleanOperations.html";
  QFile about_CGAL(html_resource_name);
  about_CGAL.open(QIODevice::ReadOnly);
  QString about_CGAL_txt = QTextStream(&about_CGAL).readAll();
#ifdef CGAL_VERSION_STR
  QString cgal_version(CGAL_VERSION_STR);
#  ifdef CGAL_FAKE_PUBLIC_RELEASE
  cgal_version.replace(QRegExp("-Ic?.*"), "");
#  endif
  about_CGAL_txt.replace("<!--CGAL_VERSION-->",
                         QString(" (version %1)")
                         .arg(cgal_version));
#endif
  QMessageBox mb;
  mb.setStyleSheet("QLabel{min-width: 700px; min-height: 300px;}");
  mb.setWindowTitle(title);
  mb.setText(about_CGAL_txt);
  mb.addButton(QMessageBox::Ok);

  QLabel *mb_label = mb.findChild<QLabel *>("qt_msgbox_label");
  if (mb_label) {
    mb_label->setTextInteractionFlags(mb_label->textInteractionFlags() |
                                      ::Qt::LinksAccessibleByMouse |
                                      ::Qt::LinksAccessibleByKeyboard);
  } else {
    std::cerr << "Cannot find child \"qt_msgbox_label\" in QMessageBox\n"
              << "  with Qt version " << QT_VERSION_STR << "!\n";
  }
  mb.exec();
}

void MainWindow::on_aboutDemo_triggered() {
  QString title = "About Demo";
  QString html_resource_name = "help/aboutThisDemo.html";
  QFile about_CGAL(html_resource_name);
  about_CGAL.open(QIODevice::ReadOnly);
  QString about_CGAL_txt = QTextStream(&about_CGAL).readAll();
#ifdef CGAL_VERSION_STR
  QString cgal_version(CGAL_VERSION_STR);
#  ifdef CGAL_FAKE_PUBLIC_RELEASE
  cgal_version.replace(QRegExp("-Ic?.*"), "");
#  endif
  about_CGAL_txt.replace("<!--CGAL_VERSION-->",
                         QString(" (version %1)")
                         .arg(cgal_version));
#endif
  QMessageBox mb;
  mb.setStyleSheet("QLabel{min-width: 700px; min-height: 300px;}");
  mb.setWindowTitle(title);
  mb.setText(about_CGAL_txt);
  mb.addButton(QMessageBox::Ok);

  QLabel *mb_label = mb.findChild<QLabel *>("qt_msgbox_label");
  if (mb_label) {
    mb_label->setTextInteractionFlags(mb_label->textInteractionFlags() |
                                      ::Qt::LinksAccessibleByMouse |
                                      ::Qt::LinksAccessibleByKeyboard);
  } else {
    std::cerr << "Cannot find child \"qt_msgbox_label\" in QMessageBox\n"
              << "  with Qt version " << QT_VERSION_STR << "!\n";
  }
  mb.exec();
}

//exception handling while switching modes
void MainWindow::exception_handler() {
  try {
    if (ensure_linear_mode()) {
      m_linear_input->Reset();
      m_linear_input->mState = m_linear_input->Start;
    } else if (ensure_circular_mode()) {
      m_circular_input->Reset();
      m_circular_input->mState = m_circular_input->Start;
    } else if (ensure_bezier_mode()) {
      m_bezier_input->Reset();
      m_bezier_input->mState = m_bezier_input->Start;
    } else if (ensure_polyline_mode()) {
      m_polyline_input->Reset();
      m_polyline_input->mState = m_polyline_input->Start;
    }
  }
  catch (const std::exception &e) {
    //todo::add catch statement
  }
}

//process all the inputs
void MainWindow::processInput(CGAL::Object o) {
  std::pair <Bezier_polygon, Bezier_boundary_source> lBI;
  Linear_polygon lLI;
  Circular_polygon lCI;
  Polyline_polygon lPI;
  //cout<<"in process input"<<endl;
  try {
    if (CGAL::assign(lBI, o)) {
      if (ensure_bezier_mode()) {
        CGAL::Orientation o = lBI.first.orientation();
        if (o == CGAL::CLOCKWISE)
          lBI.first.reverse_orientation();

        if (!m_bezier_input->isboundingRect()) {
          get_new_state(11);
          states_stack.back().active_set(m_color_active).bezier().join(Bezier_polygon_with_holes(lBI.first));
          Bezier_region_source br;
          br.push_back(lBI.second);
          states_stack.back().active_bezier_sources(m_color_active).push_back(br);
        } else {
          states_stack.back().result_set().bezier().join(Bezier_polygon_with_holes(lBI.first));
          Bezier_region_source br;
          br.push_back(lBI.second);
          states_stack.back().result_bezier_sources().push_back(br);
        }
      }
    } else if (CGAL::assign(lLI, o)) {
      if (ensure_linear_mode()) {
        CGAL::Orientation orient = lLI.orientation();
        if (orient == CGAL::CLOCKWISE) {
          lLI.reverse_orientation();
        }
        Linear_polygon_with_holes lCPWH(lLI);

        if (!m_linear_input->isboundingRect() && !m_linear_input->ishole()) {
          if (!m_linear_input->is_mink()) {
            get_new_state(11);
          }
          states_stack.back().active_set(m_color_active).linear().join(lCPWH);
          states_stack.back().active_linear_sources(m_color_active).push_back(lCPWH);
        } else {
          states_stack.back().result_set().linear().join(lCPWH);
          states_stack.back().result_linear_sources().push_back(lCPWH);
          if (m_linear_input->ishole()) {
            states_stack.back().active_set(m_color_active).difference(states_stack.back().result_set());
            states_stack.back().result_set().clear();
            states_stack.back().result_linear_sources().clear();
          }
        }
      }
    } else if (CGAL::assign(lCI, o)) {
      if (ensure_circular_mode()) {
        CGAL::Orientation o = lCI.orientation();
        if (o == CGAL::CLOCKWISE)
          lCI.reverse_orientation();

        Circular_polygon_with_holes lCPWH(lCI);

        if (!m_circular_input->isboundingRect()) {
          get_new_state(11);
          states_stack.back().active_set(m_color_active).circular().join(lCPWH);
          states_stack.back().active_circular_sources(m_color_active).push_back(lCPWH);
        } else {
          states_stack.back().result_set().circular().join(lCPWH);
          states_stack.back().result_circular_sources().push_back(lCPWH);
        }
      }
    } else if (CGAL::assign(lPI, o)) {
      if (ensure_polyline_mode()) {
        CGAL::Orientation o = lPI.orientation();
        if (o == CGAL::CLOCKWISE)
          lPI.reverse_orientation();

        Polyline_polygon_with_holes lPPWH(lPI);

        if (!m_polyline_input->isboundingRect()) {
          get_new_state(11);
          states_stack.back().active_set(m_color_active).polyline().join(lPPWH);
          states_stack.back().active_polyline_sources(m_color_active).push_back(lPPWH);
        } else {
          states_stack.back().result_set().polyline().join(lPPWH);
          states_stack.back().result_polyline_sources().push_back(lPPWH);
        }
      }
    }
  }
  catch (const std::exception &e) {
    on_actionUndo_triggered();
  }
  modelChanged();
}

#include "Boolean_set_operations_2.moc"
#include <CGAL/Qt/resources.h>

//the main function
int main(int argc, char *argv[]) {
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Boolean_operations_2 demo");
  CGAL_QT_INIT_RESOURCES;
  try {
    MainWindow w;
    w.show();
    return app.exec();
  }
  catch (const std::exception &e) {
    std::string s = e.what();
    show_error("Exception thrown during run of the program:\n" + s);
  }
}
