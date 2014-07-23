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

#ifndef ARRANGEMENT_DEMO_TAB_H
#define ARRANGEMENT_DEMO_TAB_H

#include <QWidget>

#include "ArrangementGraphicsItem.h"
#include "ArrangementDemoGraphicsView.h"
#include "ArrangementCurveInputCallback.h"
#include "DeleteCurveCallback.h"
#include "PointLocationCallback.h"
#include "VerticalRayShootCallback.h"
#include "MergeEdgeCallback.h"
#include "SplitEdgeCallback.h"
#include "EnvelopeCallback.h"
#include "FillFaceCallback.h"
#include "ArrangementDemoWindow.h"

class QGridLayout;
class ArrangementDemoWindow;

class ArrangementDemoTabBase : public QWidget
{
  Q_OBJECT

  signals:
  void modelChanged( );

public:
  ArrangementDemoTabBase( QWidget* parent );
  virtual ~ArrangementDemoTabBase( );

  virtual QGraphicsScene* getScene( ) const;
  virtual ArrangementDemoGraphicsView* getView( ) const;

  virtual CGAL::Qt::ArrangementGraphicsItemBase* getArrangementGraphicsItem( )
    const;
  virtual CGAL::Qt::GraphicsViewCurveInputBase* getCurveInputCallback( ) const;
  virtual CGAL::Qt::Callback* getDeleteCurveCallback( ) const;
  virtual CGAL::Qt::Callback* getPointLocationCallback( ) const;
  virtual VerticalRayShootCallbackBase* getVerticalRayShootCallback( ) const;
  virtual CGAL::Qt::Callback* getMergeEdgeCallback( ) const;
  virtual SplitEdgeCallbackBase* getSplitEdgeCallback( ) const;
  virtual EnvelopeCallbackBase* getEnvelopeCallback( ) const;
  virtual FillFaceCallbackBase* getFillFaceCallback( ) const;

  /**
  Set up the toolbar in an arrangement-specific way when this tab becomes
  active.
  */
  // TODO: Use tag dispatch to handle template-based setup
  virtual void setupToolbar( ArrangementDemoWindow* parent );

protected:
  virtual void setupUi( );

  ArrangementDemoGraphicsView* graphicsView;
  QGraphicsScene* scene;
  QGridLayout* layout;

  CGAL::Qt::ArrangementGraphicsItemBase* arrangementGraphicsItem;
  CGAL::Qt::GraphicsViewCurveInputBase* curveInputCallback;
  CGAL::Qt::Callback* deleteCurveCallback;
  CGAL::Qt::Callback* pointLocationCallback;
  VerticalRayShootCallbackBase* verticalRayShootCallback;
  CGAL::Qt::Callback* mergeEdgeCallback;
  SplitEdgeCallbackBase* splitEdgeCallback;
  EnvelopeCallbackBase* envelopeCallback;
  FillFaceCallbackBase* fillFaceCallback;

}; // class ArrangementDemoTabBase

template < class Arr_ >
class ArrangementDemoTab : public ArrangementDemoTabBase
{
public:
  typedef ArrangementDemoTabBase Superclass;
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Geometry_traits_2 TraitsType;

  ArrangementDemoTab( Arrangement* arrangement_, QWidget* parent = 0 ):
    Superclass( parent ),
    arrangement( arrangement_ )
  {
    // std::cout << this->scene->views( ).size( ) << std::endl;
    // set up demo components
    this->arrangementGraphicsItem =
      new CGAL::Qt::ArrangementGraphicsItem<Arrangement>(this->arrangement);
    this->curveInputCallback =
      new ArrangementCurveInputCallback<Arrangement>(this->arrangement, this);
    this->deleteCurveCallback =
      new DeleteCurveCallback<Arrangement>( this->arrangement, this );
    this->pointLocationCallback =
      new PointLocationCallback<Arrangement>( this->arrangement, this );
    this->verticalRayShootCallback =
      new VerticalRayShootCallback<Arrangement>( this->arrangement, this );
    this->mergeEdgeCallback =
      new MergeEdgeCallback<Arrangement>( this->arrangement, this );
    this->splitEdgeCallback =
      new SplitEdgeCallback<Arrangement>( this->arrangement, this );
    this->envelopeCallback =
      new EnvelopeCallback<Arrangement>( this->arrangement, this );
    this->fillFaceCallback =
      new FillFaceCallback<Arrangement>( this->arrangement, this );

    this->scene->addItem( this->arrangementGraphicsItem );
    this->curveInputCallback->setScene( this->scene );
    this->deleteCurveCallback->setScene( this->scene );
    this->pointLocationCallback->setScene( this->scene );
    this->verticalRayShootCallback->setScene( this->scene );
    this->mergeEdgeCallback->setScene( this->scene );
    this->splitEdgeCallback->setScene( this->scene );
    this->envelopeCallback->setScene( this->scene );
    this->fillFaceCallback->setScene( this->scene );

    // set up callbacks
    this->scene->installEventFilter( this->curveInputCallback );
    QObject::connect(this->curveInputCallback, SIGNAL(modelChanged()), this,
                     SIGNAL(modelChanged()));
    QObject::connect(this->deleteCurveCallback, SIGNAL(modelChanged()), this,
                     SIGNAL(modelChanged()));
    QObject::connect(this->fillFaceCallback, SIGNAL(modelChanged()), this,
                     SIGNAL(modelChanged()));
    QObject::connect(this, SIGNAL(modelChanged()),
                     this->arrangementGraphicsItem, SLOT(modelChanged()));
    QObject::connect(this, SIGNAL(modelChanged()), this->envelopeCallback,
                     SLOT(slotModelChanged()));
    // TODO: Add a connection to update the demo window when the fill color
    //       changes
  }

  void setArrangement( Arrangement* newArr )
  {
    this->scene->removeItem( this->arrangementGraphicsItem );
    delete this->arrangementGraphicsItem;
    delete this->curveInputCallback;
    delete this->deleteCurveCallback;
    delete this->pointLocationCallback;
    delete this->verticalRayShootCallback;
    delete this->mergeEdgeCallback;
    delete this->splitEdgeCallback;
    delete this->envelopeCallback;
    delete this->fillFaceCallback;

    this->arrangement = newArr;

    this->arrangementGraphicsItem =
      new CGAL::Qt::ArrangementGraphicsItem<Arrangement>( this->arrangement );

    this->curveInputCallback =
      new ArrangementCurveInputCallback<Arrangement>(this->arrangement, this);
    this->deleteCurveCallback =
      new DeleteCurveCallback<Arrangement>( this->arrangement, this );
    this->pointLocationCallback =
      new PointLocationCallback<Arrangement>( this->arrangement, this );
    this->verticalRayShootCallback =
      new VerticalRayShootCallback<Arrangement>( this->arrangement, this );
    this->mergeEdgeCallback =
      new MergeEdgeCallback<Arrangement>( this->arrangement, this );
    this->splitEdgeCallback =
      new SplitEdgeCallback<Arrangement>( this->arrangement, this );
    this->envelopeCallback =
      new EnvelopeCallback<Arrangement>( this->arrangement, this );
    this->fillFaceCallback =
      new FillFaceCallback<Arrangement>( this->arrangement, this );

    this->scene->addItem( this->arrangementGraphicsItem );
    this->curveInputCallback->setScene( this->scene );
    this->deleteCurveCallback->setScene( this->scene );
    this->pointLocationCallback->setScene( this->scene );
    this->verticalRayShootCallback->setScene( this->scene );
    this->mergeEdgeCallback->setScene( this->scene );
    this->splitEdgeCallback->setScene( this->scene );
    this->envelopeCallback->setScene( this->scene );
    this->fillFaceCallback->setScene( this->scene );

    this->scene->installEventFilter(this->curveInputCallback);
    QObject::connect(this->curveInputCallback, SIGNAL(modelChanged()), this,
                     SIGNAL(modelChanged()));
    QObject::connect(this->deleteCurveCallback, SIGNAL(modelChanged()), this,
                     SIGNAL(modelChanged()));
    QObject::connect(this->fillFaceCallback, SIGNAL(modelChanged()), this,
                     SIGNAL(modelChanged()));
    QObject::connect(this, SIGNAL(modelChanged()),
                     this->arrangementGraphicsItem, SLOT(modelChanged()));
    QObject::connect(this, SIGNAL(modelChanged()), this->envelopeCallback,
                     SLOT(slotModelChanged()));
    // TODO: Add a connection to update the demo window when the fill color
    //       changes

    emit modelChanged( );
  }

  virtual void setupToolbar( ArrangementDemoWindow* parent )
  {
    this->setupToolbar( parent, (TraitsType*) NULL );
  }

protected: // methods
  template < class TGeomTraits >
  void setupToolbar( ArrangementDemoWindow* parent,
    TGeomTraits* /*unused*/ )
  {
    parent->ui->actionConicSegment->setChecked( true );

    parent->ui->actionCurveRay->setVisible( false );
    parent->ui->actionCurveLine->setVisible( false );

    parent->ui->actionConicCircle->setVisible( false );
    parent->ui->actionConicEllipse->setVisible( false );
    parent->ui->actionConicThreePoint->setVisible( false );
    parent->ui->actionConicFivePoint->setVisible( false );

    parent->conicTypeGroup->setEnabled( true );
  }

  void setupToolbar( ArrangementDemoWindow* parent,
    Conic_arr::Geometry_traits_2* /*unused*/ )
  {
    parent->ui->actionConicSegment->setChecked( true );

    parent->ui->actionCurveRay->setVisible( false );
    parent->ui->actionCurveLine->setVisible( false );

    parent->ui->actionConicCircle->setVisible( true );
    parent->ui->actionConicEllipse->setVisible( true );
    parent->ui->actionConicThreePoint->setVisible( true );
    parent->ui->actionConicFivePoint->setVisible( true );

    parent->conicTypeGroup->setEnabled( true );
  }

  void setupToolbar( ArrangementDemoWindow* parent,
    Lin_arr::Geometry_traits_2* /*unused*/ )
  {
    parent->ui->actionConicSegment->setChecked( true );

    parent->ui->actionCurveRay->setVisible( true );
    parent->ui->actionCurveLine->setVisible( true );

    parent->ui->actionConicCircle->setVisible( false );
    parent->ui->actionConicEllipse->setVisible( false );
    parent->ui->actionConicThreePoint->setVisible( false );
    parent->ui->actionConicFivePoint->setVisible( false );

    parent->conicTypeGroup->setEnabled( true );
  }

protected:
  Arrangement* arrangement;

}; // class ArrangementDemoTab

#endif // ARRANGEMENT_DEMO_TAB_H
