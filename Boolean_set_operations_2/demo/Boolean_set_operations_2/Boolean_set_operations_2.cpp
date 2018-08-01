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
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Apurva Bhatt <response2apurva@gmail.com>
//The demo contains no error handling

#include <QApplication>
#include <qmessagebox.h>

#include <QMainWindow>
#include <QGraphicsScene>
#include <QActionGroup>
#include <QtGui>
#include <QString>
#include <QFileDialog>
#include <QInputDialog>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QSlider>
#include <QProgressBar>
#include <QMessageBox>
#include <fstream>   
#include <string>
#include <sstream>
#include <iomanip>
#include <list>

#include <boost/shared_ptr.hpp>   

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
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/approximated_offset_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#ifdef CGAL_USE_GMP
  #include <CGAL/Gmpq.h>
#else
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>
#endif

#include <QT5/Circular_polygons.h>
#include <QT5/Linear_polygons.h>
#include <QT5/Graphics_view_circular_polygon_input.h>
#include <QT5/Graphics_view_linear_polygon_input.h>

#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/Qt/utility.h>
#include <CGAL/IO/Dxf_bsop_reader.h>

#include <CGAL/Qt/GraphicsViewNavigation.h>

#include "ui_Boolean_set_operations_2.h"

#include "Typedefs.h"

using namespace std;

typedef CGAL::Qt::Circular_set_graphics_item<Circular_polygon_set,Circular_traits> Circular_GI;
typedef CGAL::Qt::Linear_set_graphics_item<Linear_polygon_set,Linear_traits>     Linear_GI;

//Functions to show errors

void show_warning(std::string aS)
{
  QMessageBox::warning(NULL, "Warning", QString(aS.c_str()));
}

void show_error(std::string aS)
{
  QMessageBox::critical(NULL, "Critical Error", QString(aS.c_str()));
}

void error(std::string aS)
{
  show_error(aS);

  throw std::runtime_error(aS);
}

void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  std::ostringstream ss ;
  
  ss << "CGAL error: " << what << " violation!" << std::endl
     << "Expr: " << expr << std::endl
     << "File: " << file << std::endl
     << "Line: " << line << std::endl;
  if ( msg != 0)
    ss << "Explanation:" << msg << std::endl;
    
  error(ss.str());
    
}

//****************************************************

//A way to maintain 3 set of polygons namely red,blue and result for all boolean operations
enum { BLUE_GROUP, RED_GROUP, RESULT_GROUP } ;

//A way to maintain 2 category of polygons namely linear,circular
//enum genrates errors so, we wil use LINEAR_TYPE=1, CIRCULAR_TYPE=2
//enum { LINEAR_TYPE, CIRCULAR_TYPE } ;

//dawing tools
QPen   sPens   [] = { QPen(QColor(0,0,255),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin)
                    , QPen(QColor(255,0,0),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin)
                    , QPen(QColor(0,255,0),0,Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin) 
                    } ;
                    
QBrush sBrushes[] = { QBrush(QColor(0,0,255,32 ))
                    , QBrush(QColor(255,0,0,32 ))
                    , QBrush(QColor(0,255,0,220))
                    } ;
//**************************************

//A base call for rep class
struct Rep_base
{
  virtual ~Rep_base() {}
  
  virtual int type () const = 0 ;
   
  virtual CGAL::Qt::GraphicsItem* gi() const = 0 ;
  virtual CGAL::Qt::GraphicsItem* gi()       = 0 ;
  
  virtual void set_pen  ( QPen   const& aPen   ) = 0 ;
  virtual void set_brush( QBrush const& aBrush ) = 0 ;
  
  virtual QRectF bounding_rect() const { return gi()->boundingRect() ; }
  
  virtual bool is_empty() const = 0 ;
  
  virtual void clear               ()                         = 0 ;
  virtual void complement          ()                         = 0 ;
  //virtual void self_minkowski_sum  ()                         = 0 ;
  virtual void assign              ( Rep_base const& aOther ) = 0 ;
  virtual void intersect           ( Rep_base const& aOther ) = 0 ;
  virtual void join                ( Rep_base const& aOther ) = 0 ;
  virtual void difference          ( Rep_base const& aOther ) = 0 ;
  virtual void symmetric_difference( Rep_base const& aOther ) = 0 ;
  
} ;


//Class for initializing 
template<class GI_, class Set_,class Gps_traits>
class Rep : public Rep_base
{
public:

  typedef GI_  GI  ;
  typedef Set_ Set ;
  typedef Rep<GI,Set,Gps_traits> Self ;
  typedef Gps_traits m_tratis;
  
  Rep() { m_GI = new GI(&m_set,m_traits) ; }
  
  Set const& set() const { return m_set ; }
  Set      & set()       { return m_set ; }
 
  virtual CGAL::Qt::GraphicsItem* gi() const { return m_GI; }
  virtual CGAL::Qt::GraphicsItem* gi()       { return m_GI; }
  
  virtual void set_pen  ( QPen   const& aPen   ) { m_GI->setPen  (aPen);   } 
  virtual void set_brush( QBrush const& aBrush ) { m_GI->setBrush(aBrush); }
      
  virtual bool is_empty() const { return m_set.is_empty() ; }
  
  virtual void clear()                         
  { 
    try
    {
      m_set.clear() ; 
    }
    catch(...)
    {
      show_error("Exception thrown during boolean operation clear");
    }
  }
  
  virtual void complement()                         
  { 
    try
    {
      m_set.complement(); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation complement");
    } 
  }
/*
  virtual void self_minkowski_sum()
  {
    std::cout<<"dummy function"<<std::endl;
  }
  */
  virtual void assign( Rep_base const& aOther ) 
  { 
    try
    {
      m_set = cast(aOther).m_set; 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation assign");
    } 
  }
  
  virtual void intersect( Rep_base const& aOther ) 
  { 
    try
    {
      m_set.intersection( cast(aOther).m_set); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation intersect");
    } 
  }
  
  virtual void join( Rep_base const& aOther ) 
  { 
    try
    {
      m_set.join( cast(aOther).m_set); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation union");
    } 
  }
  
  virtual void difference( Rep_base const& aOther ) 
  { 
    try
    {
      m_set.difference( cast(aOther).m_set); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation difference");
    } 
  }
  
  virtual void symmetric_difference( Rep_base const& aOther ) 
  { 
    try
    {
      m_set.symmetric_difference( cast(aOther).m_set); 
    } 
    catch(...)
    {
      show_error("Exception thrown during boolean operation symmetric difference");
    } 
  }
        
  static Self const& cast( Rep_base const& aOther ) { return dynamic_cast<Self const&>(aOther); }
  static Self      & cast( Rep_base      & aOther ) { return dynamic_cast<Self      &>(aOther); }
  
private:

  //For maintaining all drawing operations 
  GI* m_GI;
  //Storage for all polygons of one type. It is used as a base to perform all boolean operations
  Set m_set;
protected:  
  //pass it
  Gps_traits m_traits;  
  
} ;

//A class for connecting GUI and this file
class Circular_rep : public Rep<Circular_GI, Circular_polygon_set,Circular_traits>
{
  typedef Rep<Circular_GI, Circular_polygon_set,Circular_traits> Base ;
  
public:
  
  Circular_rep () : Base() {} 
  
  virtual int type() const { return 2 ; }
  
} ;

//A class for connecting GUI and this file
class Linear_rep : public Rep<Linear_GI, Linear_polygon_set,Linear_traits>
{

typedef Rep<Linear_GI, Linear_polygon_set,Linear_traits> Base ;
public:
  
  Linear_rep () : Base() {}
  
  virtual int type() const { return 1 ; }
  
} ;

class Curve_set
{
  //a conatiner which deletes an object when last shared_ptr gets deleted or re-initiated
  typedef boost::shared_ptr<Rep_base> Rep_ptr ;
  
public:
  //constructor
  Curve_set( int aType, QPen aPen, QBrush aBrush ) : m_pen(aPen), m_brush(aBrush) 
  {
    reset_type(aType);
  }
  void reset_type( int aType ) 
  {
    //cout<<aType<<endl;
    //setting shared_ptr for repective polygon
    if(aType==1)
    {
      m_rep = Rep_ptr(new Linear_rep());
      //m_traits=Linear_traits;
    }
    else
    {
      m_rep=Rep_ptr(new Circular_rep());
      //m_traits=Circular_traits;
    }
    //setting pen and brush
    m_rep->set_pen  (m_pen);
    m_rep->set_brush(m_brush);
  }
  
  CGAL::Qt::GraphicsItem const* gi() const { return m_rep->gi() ; }
  CGAL::Qt::GraphicsItem*       gi()       { return m_rep->gi() ; }
  
  QRectF bounding_rect() const { return m_rep->bounding_rect() ; }
  
  bool is_empty() const { return !m_rep || m_rep->is_empty(); }
  
  void clear      () { m_rep->clear() ; }
  //boolean operations
  void complement () { m_rep->complement() ; }
  /*
  void self_minkowski_sum()
  {
    if ( is_circular() )
    {
      /*
        std::vector<Circular_polygon_with_holes> circular_vector;
      Circular_polygon_set circular_set;
      circular_set.polygons_with_holes( std::back_inserter(circular_vector));
      Circular_polygon_with_holes circular_holes;
      */
      /*
        CGAL::Polygon_with_holes<Kernel> pw1=new Polygon_with_holes();
    //circular_holes=CGAL::minkowski_sum_2(circular_vector.at(0),circular_vector.at(0));
    }
  }
  */
  void assign ( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->assign( *aOther.get_circular_rep() ) ;
    }
    else
    {
      get_linear_rep()->assign( *aOther.get_linear_rep() ) ;
    }
  }
  
  void intersect( Curve_set const& aOther ) 
  {    
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->intersect( *aOther.get_circular_rep() ) ;
    }
    else
    {
      get_linear_rep()->intersect( *aOther.get_linear_rep() ) ;
    }
  }
  
  void join( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->join( *aOther.get_circular_rep() ) ;
    }
    else
    {
      get_linear_rep()->join( *aOther.get_linear_rep() ) ;
    }
  }
  
  void difference( Curve_set const& aOther ) 
  {      
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->difference( *aOther.get_circular_rep() ) ;
    }
    else
    {
      get_linear_rep()->difference( *aOther.get_linear_rep() ) ;
    }
  }
  
  void symmetric_difference( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      get_circular_rep()->symmetric_difference( *aOther.get_circular_rep() ) ;
    }
    else
    {
      get_linear_rep()->symmetric_difference( *aOther.get_linear_rep() ) ;
    } 
  }
  /*
  void minkowski_sum( Curve_set const& aOther ) 
  {
    if ( is_circular() && aOther.is_circular() )
    {
      Rep_base pr;
      pr=CGAL::minkowski_sum_2(get_circular_rep(), *aOther.get_circular_rep() ) ;
    }
    
    else
    {
      get_linear_rep()->CGAL::minkowski_sum_2( *aOther.get_linear_rep() ) ;
    } 
  } 
  */
  
  //see its need keep it for now
  Rep_base const& rep() const { return *m_rep ; }
  Rep_base&       rep()       { return *m_rep ; }
  
  bool is_circular() const { return m_rep->type() == 2 ; }  
  bool is_linear  () const { return m_rep->type() == 1 ; } 
  
  //to get rep for circualr polygons
  Circular_rep const* get_circular_rep() const { //cout<<"get const circular_rep"<<endl;
return dynamic_cast<Circular_rep const*>( boost::get_pointer(m_rep) ); }
  Circular_rep      * get_circular_rep()       { //cout<<"get normal circular_rep"<<endl;
return dynamic_cast<Circular_rep*  >( boost::get_pointer(m_rep) ); }
  
  //to get Circular_polygon_set
  Circular_polygon_set const& circular() const { return get_circular_rep()->set(); }
  Circular_polygon_set      & circular()       { //cout<<"get normal circular_polygon_set"<<endl;
return get_circular_rep()->set(); }
  
   //to get rep for linear polygons
  Linear_rep const* get_linear_rep() const { //cout<<"get const linear_rep"<<endl;return 
    dynamic_cast<Linear_rep const*>( boost::get_pointer(m_rep) ); }
  Linear_rep      * get_linear_rep()       { //cout<<"get normal linear_rep"<<endl;
 return dynamic_cast<Linear_rep*  >( boost::get_pointer(m_rep) ); }
  
  //to get Linear_polygon_set
  Linear_polygon_set const& linear() const { return get_linear_rep()->set(); }
  Linear_polygon_set      & linear()       { //cout<<"get normal linear_polygon_set"<<endl;
return get_linear_rep()->set(); }
  
private:

  //drawing tools
  QPen                        m_pen ;
  QBrush                      m_brush ;
  //a conatiner which deletes an object when last shared_ptr gets deleted or re-initiated
  boost::shared_ptr<Rep_base> m_rep ;
  
  //Linear_polygon_with_holes lpw;
  //Circular_polygon_with_holes cpw;
  
} ;

typedef std::vector<Curve_set> Curve_set_container ;

typedef Curve_set_container::const_iterator Curve_set_const_iterator ;
typedef Curve_set_container::iterator       Curve_set_iterator ;


class MainWindow : public CGAL::Qt::DemosMainWindow ,
                   public Ui::Boolean_set_operations_2
{
  Q_OBJECT// removing it gives error ui not declared
  
private:  

  QGraphicsScene        m_scene;
  //keep it intact for now check it out
  bool                  m_circular_active;
  //which type is currently active now
  bool                  m_blue_active ;
  Curve_set_container   m_curve_sets ;
  //container for curves
  Circular_region_source_container   m_blue_circular_sources ;
  Circular_region_source_container   m_red_circular_sources ;
  Linear_region_source_container     m_blue_linear_sources ; 
  Linear_region_source_container     m_red_linear_sources ; 
  
  //typedefs of classes used to draw circular and linear polygon
  CGAL::Qt::Graphics_view_linear_polygon_input<Kernel>*     m_linear_input ;
  CGAL::Qt::Graphics_view_circular_polygon_input<Kernel>*   m_circular_input ;
   
public:

  MainWindow();

private:
  
  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);
  void zoomToFit();
  
protected slots:
  
  void open( QString filename ) ;//for file handling

public slots:
  
  void processInput(CGAL::Object o);
  //make it better
	void on_actionNew_triggered();
  //make it work
	void on_actionRecenter_triggered();
  
  //boolean operations
  void on_actionComplement_triggered();
  void on_actionUnion_triggered();
  void on_actionIntersection_triggered();
  void on_actionDifference_triggered();
  void on_actionSymmetric_Difference_triggered();
  void on_actionMinkowski_Sum_triggered();

  //insert polygon operations
  void on_actionInsertLinear_triggered();
  void on_actionInsertCircular_triggered();
    
  //file handling operations
  void on_actionOpenLinear_triggered();
  void on_actionOpenCircular_triggered();
  
  void on_showBlue_toggled  (bool a_check);
  void on_showRed_toggled   (bool a_check);
  void on_showResult_toggled(bool a_check);
  
  void on_drawBlue_toggled(bool a_check);
  void on_drawRed_toggled (bool a_check);
	
	void on_actionAdd_new_polygon_triggered();
	void on_actionPAN_triggered();
  
signals:

	void changed();
  
private:
  
  void modelChanged()
  {
    emit(changed());
  }
  
  
  //warning message for user
  bool ask_user_yesno( const char* aTitle, const char* aQuestion )
  {
    return QMessageBox::warning(this
                               ,aTitle
                               ,QString(aQuestion)
                               ,"&Yes"
                               ,"&No"
                               ,QString::null
                               , 1
                               , 1 
                               ) == 0 ;
  }
  
  //for setting Curve_set of aGroup type an int representing a set of polygon of a specific type
  Curve_set& set( int aGroup ) { //cout<<"set function"<<endl;
  return m_curve_sets[aGroup] ; }
  
  //setting curve
  Curve_set& blue_set  () { return set(BLUE_GROUP)  ; }
  Curve_set& red_set   () { return set(RED_GROUP)   ; }
  Curve_set& result_set() { return set(RESULT_GROUP); }

  //gets which group is currently active now
  int active_group() const { return m_blue_active ? BLUE_GROUP : RED_GROUP ; }
  
  //sets the current active group
  Curve_set& active_set()   { return set(active_group()) ; }

  //returns circular containers
  Circular_region_source_container const& blue_circular_sources() const { return m_blue_circular_sources ; }
  Circular_region_source_container      & blue_circular_sources()       { return m_blue_circular_sources ; }

  Circular_region_source_container const& red_circular_sources () const { return m_red_circular_sources ; }
  Circular_region_source_container      & red_circular_sources ()       { return m_red_circular_sources ; }
  
  //returns linear containers
  Linear_region_source_container const& blue_linear_sources() const { return m_blue_linear_sources ; }
  Linear_region_source_container      & blue_linear_sources()       { return m_blue_linear_sources ; }

  Linear_region_source_container const& red_linear_sources () const { return m_red_linear_sources ; }
  Linear_region_source_container      & red_linear_sources ()       { return m_red_linear_sources ; }

  //returns active blue container
  Circular_region_source_container const& active_circular_sources() const { return m_blue_active ? m_blue_circular_sources : m_red_circular_sources ; }
  Circular_region_source_container      & active_circular_sources()       { return m_blue_active ? m_blue_circular_sources : m_red_circular_sources ; }

  //returns active linear container
  Linear_region_source_container const& active_linear_sources() const { return m_blue_active ? m_blue_linear_sources : m_red_linear_sources ; }
  Linear_region_source_container      & active_linear_sources()       { return m_blue_active ? m_blue_linear_sources : m_red_linear_sources ; }  
  
	void SetViewBlue  ( bool a_check ) { showBlue  ->setChecked(a_check); }  
  void SetViewRed   ( bool a_check ) { showRed   ->setChecked(a_check); }  
  void SetViewResult( bool a_check ) { showResult->setChecked(a_check); }  
	
  //changes the set of polygons of a specific type
  void ToogleView( int aGROUP, bool a_check );
  
  void link_GI( CGAL::Qt::GraphicsItem* aGI )
  {
    QObject::connect(this, SIGNAL(changed()), aGI, SLOT(modelChanged()));
    m_scene.addItem( aGI );
  }
  
  void unlink_GI( CGAL::Qt::GraphicsItem* aGI )
  {
    m_scene.removeItem( aGI );
    QObject::disconnect(this, SIGNAL(changed()), aGI, SLOT(modelChanged()));
  }
  
  void switch_set_type( Curve_set& aSet, int aType );
  
  void switch_sets_type( int aType );
  
  bool ensure_circular_mode();
  
  bool ensure_linear_mode();//see if it is need
};


MainWindow::MainWindow()
  : DemosMainWindow()
  , m_circular_active(false)//default
  , m_blue_active(true)    //default
{
  CGAL::set_error_handler  (error_handler);
  CGAL::set_warning_handler(error_handler);
		
	setupUi(this);

  setAcceptDrops(true);
  //cout<<"elementry setups"<<endl;
  //default setups
  m_curve_sets.push_back( Curve_set(1, sPens[BLUE_GROUP]  , sBrushes[BLUE_GROUP]  ) ) ;
  m_curve_sets.push_back( Curve_set(1, sPens[RED_GROUP]   , sBrushes[RED_GROUP]   ) ) ;
  m_curve_sets.push_back( Curve_set(1, sPens[RESULT_GROUP], sBrushes[RESULT_GROUP]) ) ;
  
  //cout<<"curve setups"<<endl;
  for( Curve_set_iterator si = m_curve_sets.begin(); si != m_curve_sets.end() ; ++ si )
  { //cout<<"setting curves"<<endl;
    link_GI(si->gi()) ;
  }
  //
  // Setup the m_scene and the view
  //
  m_scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  m_scene.setSceneRect(-100, -100, 100, 100);
  this->graphicsView->setScene(&m_scene);
  this->graphicsView->setMouseTracking(true);
	//this->on_actionInsertLinear_triggered();

  // Turn the vertical axis upside down 
  this->graphicsView->scale(1, -1);
    //cout<<"UI setup"<<endl;
    
  //adding basic setups
    
  //need to be given finishing touch
	// The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);
  //setting the menus
  this->setupStatusBar();
  this->setupOptionsMenu();
	//link to a page describing
  this->addAboutDemo(":/cgal/help/index.html");
  //link for about page of CGAL
  this->addAboutCGAL();

  this->addRecentFiles(this->menuFile, this->actionQuit);
  //cout<<"extra setup"<<endl;
  
  //initializing classes to draw respective polygons using mouse
  m_linear_input  =new CGAL::Qt::Graphics_view_linear_polygon_input<Kernel>(this, &m_scene);
  m_circular_input=new CGAL::Qt::Graphics_view_circular_polygon_input<Kernel>(this, &m_scene);
    
  //connecting GUI and the code base
  QObject::connect(m_linear_input  , SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));
  QObject::connect(m_circular_input, SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));
  //QObject::connect(this->actionAdd_new_polygon, SIGNAL(triggered()), this, SLOT(on_actionAdd_new_polygon_triggered()));
  m_scene.installEventFilter(m_linear_input);  
    
  QObject::connect(this->actionQuit, SIGNAL(triggered()), this, SLOT(close()));
  //QObject::connect(this, SIGNAL(openRecentFile(QString)), this, SLOT(open(QString)));//for file handling
  QObject::connect(drawBlue, SIGNAL(toggled(bool)), this, SLOT(on_drawBlue_toggled (bool)));
  QObject::connect(drawRed , SIGNAL(toggled(bool)), this, SLOT(on_drawRed_toggled(bool)));
  //QObject::connect(actionPAN,SIGNAL(triggered()), this ,SLOT(on_actionPAN_triggered()));
		
  QObject::connect(showBlue  , SIGNAL(toggled(bool)), this, SLOT(on_showBlue_toggled (bool)));
  QObject::connect(showRed   , SIGNAL(toggled(bool)), this, SLOT(on_showRed_toggled  (bool)));
  QObject::connect(showResult, SIGNAL(toggled(bool)), this, SLOT(on_showResult_toggled (bool)));
		
	//cout<<"connecting stuff"<<endl;
}

void MainWindow::on_showBlue_toggled  (bool a_check) { ToogleView(BLUE_GROUP  ,a_check); }
void MainWindow::on_showRed_toggled   (bool a_check) { ToogleView(RED_GROUP   ,a_check); }
void MainWindow::on_showResult_toggled(bool a_check) { ToogleView(RESULT_GROUP,a_check); }
  
void MainWindow::on_drawBlue_toggled(bool a_check) { m_blue_active =  a_check ; }
void MainWindow::on_drawRed_toggled (bool a_check) { m_blue_active = !a_check ; }
//keep it no use for now
void MainWindow::on_actionNew_triggered() 
{
  for( Curve_set_iterator si = m_curve_sets.begin(); si != m_curve_sets.end() ; ++ si )
    si->clear();
 //cout<<"In new Polygon"<<endl;
 blue_circular_sources().clear();
 red_circular_sources().clear();
 blue_linear_sources().clear();
 red_circular_sources().clear();
	
 SetViewBlue  (true);
 SetViewRed   (true);
 SetViewResult(true);
 
 m_circular_active = false ;
		
 drawBlue->setChecked(true);
  
	modelChanged();
}

//extra utilities
void MainWindow::on_actionRecenter_triggered()
{
  zoomToFit();
}

void MainWindow::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("text/uri-list"))
    event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent *event)
{
  QString filename = event->mimeData()->urls().at(0).path();
  open(filename);
  event->acceptProposedAction();
}

void MainWindow::on_actionAdd_new_polygon_triggered()
{
  //cout<<"added new polygon"<<endl;
  //ToogleView(BLUE_GROUP, false);
  m_blue_active=false;
  ToogleView(RED_GROUP, true);
}

void MainWindow::on_actionOpenLinear_triggered()
{
  open(QFileDialog::getOpenFileName(this, tr("Open Linear Polygon"), "../data", tr("Linear Curve files (*.lps)") ));
}

void MainWindow::on_actionOpenCircular_triggered()
{
  open(QFileDialog::getOpenFileName(this, tr("Open Circular Polygon"), "../data", tr("Circular curve files (*.dxf)") ));
}

//for converting linear part of circular polygon to circular part
Circular_polygon linear_2_circ( Circular_Linear_polygon const& pgn )
{
  CGAL::Cartesian_converter<Kernel,Kernel> convert ;
  
  Circular_polygon rCP;
  
  for( Circular_Linear_polygon::Edge_const_iterator ei = pgn.edges_begin(); ei != pgn.edges_end(); ++ei )
  {
    if  ( ei->source() != ei->target() )
      rCP.push_back( Circular_X_monotone_curve( convert(ei->source()), convert(ei->target())) );
  }  

  return rCP;
}

//for converting linear part of circular polygon with holes to circular part
Circular_polygon_with_holes linear_2_circ( Circular_Linear_polygon_with_holes const& pwh )
{
  Circular_polygon_with_holes rCP( linear_2_circ(pwh.outer_boundary()) ) ;
  
  for( Circular_Linear_polygon_with_holes::Hole_const_iterator hi = pwh.holes_begin(); hi != pwh.holes_end(); ++ hi )
    rCP.add_hole( linear_2_circ(*hi)  );

  return rCP;
}

//check out
void MainWindow::switch_set_type( Curve_set& aSet, int aType )
{
  unlink_GI( aSet.gi() ) ;
  
  aSet.reset_type(aType);
  
  link_GI( aSet.gi() ) ;
  
  modelChanged();
}

void MainWindow::switch_sets_type( int aType )
{
  switch_set_type( blue_set  (), aType ) ; 
  switch_set_type( red_set   (), aType ) ; 
  switch_set_type( result_set(), aType ) ; 
}

bool MainWindow::ensure_circular_mode()
{
  
  if ( ! m_circular_active )
  {
    bool lProceed = blue_set().is_empty() && red_set().is_empty() ;
    
    if ( ! lProceed )
      lProceed = ask_user_yesno("Linear/Circular mode switch"
                               ,"You are about to load a circular poygon, but there are linear curves already loaded.\n" \
                                "Both types are not interoperable. In order to proceed, the circular curves must be removed first.\n" \
                                "Yes to remove and proceed?\n"
                               ) ;
      
    if ( lProceed )
    {
      switch_sets_type(2);
      m_circular_active = true ;
    }
  }
  return m_circular_active ;
}

bool MainWindow::ensure_linear_mode()
{
 
  if ( m_circular_active )
  {
    bool lProceed = blue_set().is_empty() && red_set().is_empty() ;
    
    if ( ! lProceed )
      lProceed = ask_user_yesno("Linear/Circular mode switch"
                               ,"You are about to load a linear poygon, but there are circular curves already loaded.\n" \
                                "Both types are not interoperable. In order to proceed, the linear curves must be removed first.\n" \
                                "Yes to remove and proceed?\n"
                               ) ;
      
    if ( lProceed )
    {
      switch_sets_type(1);
      m_circular_active = false ;
    }
  }
  return !m_circular_active ;
}
//check out
bool read_linear ( QString aFileName, Linear_polygon_set& rSet, Linear_region_source_container& rSources )
{
  bool rOK = false ;
  return rOK ;
}

bool read_circular ( QString aFileName, Circular_polygon_set& rSet, Circular_region_source_container& rSources )
{
  bool rOK = false ;
  return rOK ;
}


void MainWindow::open( QString fileName )
{
  if(! fileName.isEmpty())
  {
    bool lRead = false ;
    
    if(fileName.endsWith(".lps"))
    {
      if ( ensure_linear_mode() )
        lRead = read_linear(fileName,active_set().linear(), active_linear_sources() ) ;
    }
    else if (fileName.endsWith(".dxf"))
    {
      if ( ensure_circular_mode() )
        lRead = read_circular(fileName,active_set().circular(), active_circular_sources() ) ;
    }
     
    if ( lRead )
    {
      modelChanged();
      zoomToFit();
      this->addToRecentFiles(fileName);
      
    }
  }  
}


void MainWindow::on_actionInsertCircular_triggered()
{
  //cout<<"signal circular triggered"<<endl;
  /*
	if(!m_circular_active)
	{   
		m_scene.removeEventFilter(m_linear_input);
		m_circular_active=true;
	}
	*/
	this->graphicsView->setDragMode(QGraphicsView::NoDrag);
	if(ensure_circular_mode())	
		m_scene.installEventFilter(m_circular_input);    
  //else
		
}

void MainWindow::on_actionInsertLinear_triggered()
{
  //cout<<"signal linear triggered"<<endl;
  /*
	if(m_circular_active)
	{
		m_scene.removeEventFilter(m_circular_input);
		m_circular_active=false;
	}
	*/
	this->graphicsView->setDragMode(QGraphicsView::NoDrag);
	if(ensure_linear_mode())
		m_scene.installEventFilter(m_linear_input);
	//else
		
}

void MainWindow::processInput(CGAL::Object o )
{
  //m_blue_active =  true ;
  
  Linear_polygon   lLI ;
  Circular_polygon lCI ;
     
  //cout<<"process input"<<endl;  
  if(CGAL::assign(lLI, o))
  {
    //cout<<"came to linear"<<endl;
    if ( ensure_linear_mode() )
    {
      //cout<<"inside linear"<<endl;
      CGAL::Orientation o = lLI.orientation();
      //return;
      //cout<<"set linear's orientation"<<endl;
      if( o == CGAL::CLOCKWISE )
      {
        ////cout<<"passed if"<<endl;
        lLI.reverse_orientation();
      }
      //cout<<"oriented"<<endl;
      Linear_polygon_with_holes lCPWH(lLI);
      //cout<<"l l l l"<<endl;
      active_set().linear().join(lCPWH) ;  
      //cout<<"hi linear"<<endl;
      active_linear_sources().push_back(lCPWH);
      //cout<<"processed linear"<<endl;
    }
  }

  else if ( CGAL::assign(lCI, o) )
  {
    //cout<<"came to circular"<<endl;
    if ( ensure_circular_mode() )
    {
      //cout<<"inside circular"<<endl;
      CGAL::Orientation o = lCI.orientation();
      //cout<<"set circular's orientation"<<endl;
      if ( o == CGAL::CLOCKWISE )
        lCI.reverse_orientation();

      //cout<<"oriented"<<endl;
      Circular_polygon_with_holes lCPWH(lCI);
      //cout<<"c c c c"<<endl;
      active_set().circular().join(lCPWH) ;  
      //cout<<"hi circular"<<endl;
      active_circular_sources().push_back(lCPWH);
      //cout<<"processed circualar"<<endl;
    }
  }
  modelChanged();  
    
}

void MainWindow::on_actionMinkowski_Sum_triggered()
{
  bool lDone = false ;
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  if ( !blue_set().is_empty() )
  {
    /*
		result_set().assign( blue_set() ) ;
    result_set().complement();
    */
		lDone = true ;
  }
  this->setCursor(old);
  if ( lDone )
  {
    modelChanged();
  }
}

//only blue complement
void MainWindow::on_actionComplement_triggered()
{
  bool lDone = false ;
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  if ( !blue_set().is_empty() )
  {
    result_set().assign( blue_set() ) ;
    result_set().complement();
    lDone = true ;
  }
  this->setCursor(old);
  if ( lDone )
  {
    modelChanged();
  }
}

void MainWindow::on_actionIntersection_triggered()
{
  bool lDone = false ;
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  if ( !blue_set().is_empty() && !red_set().is_empty() )
  {
    result_set().assign( red_set() ) ;
    result_set().intersect(blue_set());
    lDone = true ;
  }
  this->setCursor(old);
  if ( lDone )
  {
    modelChanged();
  }
}  

//blue - red
void MainWindow::on_actionDifference_triggered()
{
  bool lDone = false ;  
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  if ( !blue_set().is_empty() && !red_set().is_empty() )
  {
    result_set().assign( blue_set() ) ;
    result_set().difference(red_set());
    lDone = true ;
  }
  this->setCursor(old);
    
  if ( lDone )
  {
    modelChanged();
  }
}

void MainWindow::on_actionSymmetric_Difference_triggered()
{
  bool lDone = false ;
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  if ( !blue_set().is_empty() && !red_set().is_empty() )
  {
    result_set().assign( red_set() ) ;
    result_set().symmetric_difference(blue_set());
    lDone = true ;
  }
  this->setCursor(old);
  if ( lDone )
  {
    modelChanged();
  }
}


void MainWindow::on_actionUnion_triggered()
{
  bool lDone = false ;
  //cout<<"came to union"<<endl;
  QCursor old = this->cursor();
  this->setCursor(Qt::WaitCursor);
  if ( !blue_set().is_empty() && !red_set().is_empty() )
  {
    //cout<<"ready for union"<<endl;
    result_set().clear();
    result_set().assign( red_set() ) ;
    result_set().join(blue_set());
    lDone = true ;
  }
  this->setCursor(old);
  
  if ( lDone )
  {    
    //cout<<"changed model"<<endl;
    modelChanged();
  }
  //cout<<"union is done"<<endl;
}

//to change which polygons to see on the screen
void MainWindow::ToogleView( int aGROUP, bool a_check )
{
  if ( a_check )
  {
    set(aGROUP).gi()->show();
    //cout<<"if triggered"<<endl;
  }
  else 
  {
    set(aGROUP).gi()->hide();
    //cout<<"else triggered"<<endl;
  }
}

void MainWindow::on_actionPAN_triggered()
{
	//cout<<"PAN triggered"<<endl;
	if(!m_circular_active)
		m_scene.removeEventFilter(m_linear_input);
		//QObject::disconnect(m_linear_input  , SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));
	else
		m_scene.removeEventFilter(m_circular_input);
		//QObject::disconnect(m_circular_input  , SIGNAL(generate(CGAL::Object)), this, SLOT(processInput(CGAL::Object)));
	this->graphicsView->setDragMode(QGraphicsView::ScrollHandDrag);
	
}
void MainWindow::zoomToFit()
{
  boost::optional<QRectF> lTotalRect ;
  
  for ( Curve_set_const_iterator si = m_curve_sets.begin() ; si != m_curve_sets.end() ; ++ si )
  {
    if ( !si->is_empty() ) 
    {
      QRectF lRect = si->bounding_rect();
      if ( lTotalRect )
           lTotalRect = *lTotalRect | lRect ;
      else lTotalRect = lRect ;  
    }
  }
                   
  if ( lTotalRect )
  {
    this->graphicsView->setSceneRect(*lTotalRect);
    this->graphicsView->fitInView(*lTotalRect, Qt::KeepAspectRatio);  
  }                 
}

//Main part
#include "Boolean_set_operations_2.moc"
#include <CGAL/Qt/resources.h>
int main(int argc, char *argv[])
{
  //QApplication a(argc, argv);
  QApplication app(argc, argv);

  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Boolean_operations_2 demo");
  CGAL_QT_INIT_RESOURCES;
  try
  {
//std::cout<<"hello";    
    MainWindow w;
    w.show();

    return app.exec();
  }
  catch (const std::exception e)
  {
    std::string s = e.what();
    show_error("Exception throne during run of the program:\n" + s);
  }
}
