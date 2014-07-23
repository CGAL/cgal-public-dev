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
// $URL: $
// $Id: $
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#ifndef ARRANGEMENT_DEMO_WINDOW_H
#define ARRANGEMENT_DEMO_WINDOW_H

#include "ArrangementGraphicsItem.h"
#include "ArrangementTypes.h"
#include "DeleteCurveCallback.h"
#include "PointLocationCallback.h"
#include "VerticalRayShootCallback.h"
#include "MergeEdgeCallback.h"
#include "SplitEdgeCallback.h"
#include "EnvelopeCallback.h"
#include "ArrangementDemoTab.h"

#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/Qt/DemosMainWindow.h>
#include <CGAL/IO/pixmaps/hand.xpm>

#include <Qt>

#include <boost/variant.hpp>
#include <boost/variant/static_visitor.hpp>

#include "ui_ArrangementDemoWindow.h"

//#include <QFileDialog>
//#include <QInputDialog>
//#include <QMessageBox>
//#include <QtGui>

namespace Ui { class ArrangementDemoWindow; }

class QActionGroup;
class ArrangementDemoTabBase;

class ArrangementDemoWindow : public CGAL::Qt::DemosMainWindow
{
  Q_OBJECT
  public:
#if 0
  typedef Seg_traits::Point_2 Point;
  typedef Seg_traits::Segment_2 Segment;
#endif
  typedef enum TraitsType {
    SEGMENT_TRAITS,
    POLYLINE_TRAITS,
    CONIC_TRAITS,
    LINEAR_TRAITS,
    CIRCULAR_ARC_TRAITS,
    BEZIER_TRAITS,
    ALGEBRAIC_TRAITS
  } TraitsType;

  private:
  typedef boost::variant< Seg_arr*,
    Pol_arr*,
    Conic_arr*,
    Lin_arr*,
    Arc_arr*,
    Bezier_arr*,
    Alg_seg_arr*
  > SomeArrPtrType;

  typedef boost::variant< std::pair< Seg_arr*, Seg_arr* >,
    std::pair< Pol_arr*, Pol_arr* >,
    std::pair< Conic_arr*, Conic_arr* >,
    std::pair< Lin_arr*, Lin_arr* >,
    std::pair< Arc_arr*, Arc_arr* >,
    std::pair< Bezier_arr*, Bezier_arr* >,
    std::pair< Alg_seg_arr*, Alg_seg_arr* >
  > SomePairOfArrPtrType;

  /**
  Used to make a new overlay tab from a pair of arrangements of matching type.

  Usage
  -----
  1. Construct with a reference to ArrangementDemoWindow.
  2. Apply it to SomePairOfArrPtrType with boost::apply_visitor.
  */
  class MakeOverlayVisitor :
    public boost::static_visitor< >
  {
  public:
    MakeOverlayVisitor( ArrangementDemoWindow& parent );

    void operator()( std::pair< Seg_arr*, Seg_arr* >& pa );
    void operator()( std::pair< Pol_arr*, Pol_arr* >& pa );
    void operator()( std::pair< Conic_arr*, Conic_arr* >& pa );
    void operator()( std::pair< Lin_arr*, Lin_arr* >& pa );
    void operator()( std::pair< Arc_arr*, Arc_arr* >& pa );
    void operator()( std::pair< Bezier_arr*, Bezier_arr* >& pa );
    void operator()( std::pair< Alg_seg_arr*, Alg_seg_arr* >& pa );

    ArrangementDemoWindow& m_parent;
  };
  friend class MakeOverlayVisitor;

  public:
  ArrangementDemoWindow(QWidget* parent = 0);
  ~ArrangementDemoWindow();

  ArrangementDemoTabBase* makeTab( TraitsType tt );
  ArrangementDemoTabBase* getTab( unsigned int tabIndex ) const;
  ArrangementDemoTabBase* getCurrentTab( ) const;

  std::vector< QString > getTabLabels( ) const;
  std::vector< CGAL::Object > getArrangements( ) const;

  void makeOverlayTab( SomePairOfArrPtrType arrs );

  /**
  \param[out] arrPair - the selected pair of arrangement pointers
  \return whether the conversion was successful
  */
  static bool ToPairOfArr( const std::vector< CGAL::Object >& arrs,
    SomePairOfArrPtrType* arrPair );

public slots:
  void updateMode( QAction* a );
  void updateEnvelope( QAction* a );
  void updateSnapping( QAction* a );
  void updateConicType( QAction* a );
  void on_actionNewTab_triggered( );
  void on_actionSaveAs_triggered( );
  void on_actionOpen_triggered( );
  void on_actionQuit_triggered( );
  void on_tabWidget_currentChanged( );
  void on_actionOverlay_triggered( );
  void on_actionCloseTab_triggered( );
  void on_actionPrintConicCurves_triggered( );
  void on_actionZoomIn_triggered( );
  void on_actionZoomOut_triggered( );
  void on_actionPreferences_triggered( );
  void on_actionFillColor_triggered( );


signals:
  void modelChanged( );

protected: // fields
  void setupUi( );
  void resetCallbackState( unsigned int tabIndex );
  void removeCallback( unsigned int tabIndex );
  void updateFillColorSwatch( );

  void openArrFile( QString filename );
  void openDatFile( QString filename );

  std::vector< ArrangementDemoTabBase* > tabs;
  std::vector< CGAL::Object > arrangements;
  std::vector< QAction* > activeModes; // for the current tab; always size 1
  unsigned int lastTabIndex;

public: // expose the UI for convenience
  Ui::ArrangementDemoWindow* ui;
  QActionGroup* modeGroup;
  QActionGroup* envelopeGroup;
  QActionGroup* snapGroup;
  QActionGroup* conicTypeGroup;
};

#endif // ARRANGEMENT_DEMO_WINDOW_H
