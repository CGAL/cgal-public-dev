// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#include "ArrangementDemoTab.h"
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

#include <QGridLayout>

ArrangementDemoTabBase::ArrangementDemoTabBase( QWidget* parent ) :
  QWidget( parent ),
  graphicsView( new ArrangementDemoGraphicsView( this ) ),
  scene( new QGraphicsScene( ) ),
  layout( new QGridLayout( this ) ),
  arrangementGraphicsItem( nullptr ),
  curveInputCallback( nullptr ),
  deleteCurveCallback( nullptr ),
  pointLocationCallback( nullptr ),
  verticalRayShootCallback( nullptr ),
  mergeEdgeCallback( nullptr ),
  splitEdgeCallback( nullptr ),
  envelopeCallback( nullptr ),
  fillFaceCallback( nullptr ),
  activeCallback( nullptr )
{
  this->setupUi( );
}

ArrangementDemoTabBase::~ArrangementDemoTabBase( )
{
  this->unhookCallbacks();
}

void ArrangementDemoTabBase::setupUi( )
{
  this->layout->addWidget( this->graphicsView, 0, 0 );
  this->graphicsView->setScene( this->scene );
  // TODO: Find suitable values
  double xymin = -std::numeric_limits<double>::max() / 1048576;
  double wh = std::numeric_limits<double>::max() / 524288;
  this->scene->setSceneRect(xymin, xymin, wh, wh);
}

QGraphicsScene* ArrangementDemoTabBase::getScene( ) const
{
  return this->scene;
}

ArrangementDemoGraphicsView* ArrangementDemoTabBase::getView() const
{
  return this->graphicsView;
}

CGAL::Qt::ArrangementGraphicsItemBase*
ArrangementDemoTabBase::getArrangementGraphicsItem( ) const
{
  return this->arrangementGraphicsItem.get();
}

//! Determining the points of the arrangement
/*!
	\return call to the ArrangementCurveInputCallback
*/
CGAL::Qt::GraphicsViewCurveInputBase*
ArrangementDemoTabBase::getCurveInputCallback( ) const
{
  return this->curveInputCallback.get();
}

//! eraser option i.e. to delete the selected curve.
/*!
	\return the drawing after the selected curve has been removed
*/
CGAL::Qt::Callback* ArrangementDemoTabBase::getDeleteCurveCallback( ) const
{
  return this->deleteCurveCallback.get();
}

//! Returns the point where the mouse is selected.
/*!
	\return the point from the curve originates
*/
CGAL::Qt::Callback* ArrangementDemoTabBase::getPointLocationCallback( ) const
{
  return this->pointLocationCallback.get();
}

//! Vertical ray offshoot feedback
/*!
	\return the ray in the direction closest to the edge of the screen
*/
VerticalRayShootCallbackBase*
ArrangementDemoTabBase::getVerticalRayShootCallback( ) const
{
  return this->verticalRayShootCallback.get();
}

//! Merging the segments
/*!
	\return the curves after merging them back together
*/
CGAL::Qt::Callback* ArrangementDemoTabBase::getMergeEdgeCallback( ) const
{
  return this->mergeEdgeCallback.get();
}

//! Splitting the curves drawn in the screen with points.
/*!
	\return the points of splitting
*/
SplitEdgeCallbackBase* ArrangementDemoTabBase::getSplitEdgeCallback( ) const
{
  return this->splitEdgeCallback.get();
}

//! feedback after the envelope call.
/*!
	\return result of the envelope call
*/
EnvelopeCallbackBase* ArrangementDemoTabBase::getEnvelopeCallback( ) const
{
  return this->envelopeCallback.get();
}

//! member function to fill the viewport
/*!
	\return result after calling the fill color option
*/
FillFaceCallbackBase* ArrangementDemoTabBase::getFillFaceCallback( ) const
{
  return this->fillFaceCallback.get();
}

void ArrangementDemoTabBase::activateCurveInputCallback(CGAL::Qt::CurveType type)
{
  this->unhookCallbacks();

  this->curveInputCallback->setCurveType(type);
  this->getScene()->installEventFilter(this->curveInputCallback.get());
  this->activeCallback = this->curveInputCallback.get();
}

void ArrangementDemoTabBase::unhookCallbacks()
{
  if (this->activeCallback)
  {
    this->getScene()->removeEventFilter(this->activeCallback);

    // TODO(Ahmed Essam): This is ugly. Fix it.
    // GraphicsViewCurveInputBase should inherit from Callback
    auto callback = dynamic_cast<CGAL::Qt::Callback*>(this->activeCallback);
    if (callback) { callback->reset(); }
    else
    {
      auto curveInput = dynamic_cast<CGAL::Qt::GraphicsViewCurveInputBase*>(
        this->activeCallback);
      if (curveInput) curveInput->reset();
    }

    this->activeCallback = nullptr;
  }
}

void ArrangementDemoTabBase::unhookAndInstallEventFilter(QObject* obj)
{
  this->unhookCallbacks();
  this->getScene()->installEventFilter(obj);
  this->activeCallback = obj;
}

void ArrangementDemoTabBase::activateDeleteCurveCallback()
{
  // TODO(Ahmed Essam): Create different button for modes of delete
  if (
    this->activeCallback ==
    static_cast<QObject*>(this->deleteCurveCallback.get()))
  {
    auto deleteMode = this->deleteCurveCallback->getDeleteMode();
    if (deleteMode == DeleteMode::DeleteOriginatingCuve)
      this->deleteCurveCallback->setDeleteMode(DeleteMode::DeleteEdge);
    else
      this->deleteCurveCallback->setDeleteMode(
        DeleteMode::DeleteOriginatingCuve);
  }
  else
  {
    this->unhookAndInstallEventFilter(this->deleteCurveCallback.get());
  }
}

void ArrangementDemoTabBase::activatePointLocationCallback()
{
  this->unhookAndInstallEventFilter(this->pointLocationCallback.get());
}

void ArrangementDemoTabBase::activateVerticalRayShootCallback(bool shootingUp)
{
  this->verticalRayShootCallback->setShootingUp(shootingUp);
  this->unhookAndInstallEventFilter(this->verticalRayShootCallback.get());
}

void ArrangementDemoTabBase::activateMergeEdgeCallback()
{
  this->unhookAndInstallEventFilter(this->mergeEdgeCallback.get());
}

void ArrangementDemoTabBase::activateSplitEdgeCallback()
{
  this->unhookAndInstallEventFilter(this->splitEdgeCallback.get());
}

void ArrangementDemoTabBase::activateFillFaceCallback()
{
  this->unhookAndInstallEventFilter(this->fillFaceCallback.get());
}

template <class Arr_>
ArrangementDemoTab<Arr_>::ArrangementDemoTab(
  QWidget* parent, std::unique_ptr<Arrangement> arrangement_) :
    Superclass(parent),
    arrangement(std::move(arrangement_))
{
  if (!this->arrangement) this->initArrangement();
  this->initComponents();
  this->setupCallbacks();
}

template <class Arr_>
ArrangementDemoTab<Arr_>::~ArrangementDemoTab()
{
  // make unique_ptr handle deletion instead of Qt
  this->scene->removeItem(this->arrangementGraphicsItem.get());
}

template <class Arr_>
void ArrangementDemoTab<Arr_>::initArrangement()
{
  this->arrangement = std::make_unique<Arrangement>();
}

template <class Arr_>
void ArrangementDemoTab<Arr_>::initComponents()
{
  this->arrangementGraphicsItem =
    std::make_unique<CGAL::Qt::ArrangementGraphicsItem<Arrangement>>(
      this->arrangement.get());
  this->curveInputCallback =
    std::make_unique<ArrangementCurveInputCallback<Arrangement>>(
      this->arrangement.get(), this, this->scene);
  this->deleteCurveCallback =
    std::make_unique<DeleteCurveCallback<Arrangement>>(
      this->arrangement.get(), this);
  this->pointLocationCallback =
    std::make_unique<PointLocationCallback<Arrangement>>(
      this->arrangement.get(), this);
  this->verticalRayShootCallback =
    std::make_unique<VerticalRayShootCallback<Arrangement>>(
      this->arrangement.get(), this);
  this->mergeEdgeCallback = std::make_unique<MergeEdgeCallback<Arrangement>>(
    this->arrangement.get(), this);
  this->splitEdgeCallback = std::make_unique<SplitEdgeCallback<Arrangement>>(
    this->arrangement.get(), this);
  this->envelopeCallback = std::make_unique<EnvelopeCallback<Arrangement>>(
    this->arrangement.get(), this);
  this->fillFaceCallback = std::make_unique<FillFaceCallback<Arrangement>>(
    this->arrangement.get(), this);

  this->scene->addItem(this->arrangementGraphicsItem.get());
  this->arrangementGraphicsItem->setScene(this->scene);
  this->curveInputCallback->setScene(this->scene);
  this->deleteCurveCallback->setScene(this->scene);
  this->pointLocationCallback->setScene(this->scene);
  this->verticalRayShootCallback->setScene(this->scene);
  this->mergeEdgeCallback->setScene(this->scene);
  this->splitEdgeCallback->setScene(this->scene);
  this->envelopeCallback->setScene(this->scene);
  this->fillFaceCallback->setScene(this->scene);
}

template <class Arr_>
void ArrangementDemoTab<Arr_>::setupCallbacks()
{
  // set up callbacks
  QObject::connect(
    this->curveInputCallback.get(), SIGNAL(modelChanged()), this,
    SIGNAL(modelChanged()));
  QObject::connect(
    this->deleteCurveCallback.get(), SIGNAL(modelChanged()), this,
    SIGNAL(modelChanged()));
  QObject::connect(
    this->fillFaceCallback.get(), SIGNAL(modelChanged()), this,
    SIGNAL(modelChanged()));
  QObject::connect(
    this, SIGNAL(modelChanged()), this->arrangementGraphicsItem.get(),
    SLOT(modelChanged()));
  QObject::connect(
    this, SIGNAL(modelChanged()), this->envelopeCallback.get(),
    SLOT(slotModelChanged()));
  QObject::connect(
    this->splitEdgeCallback.get(), SIGNAL(modelChanged()), this,
    SIGNAL(modelChanged()));
  QObject::connect(
    this, SIGNAL(modelChanged()), this, SLOT(slotModelChanged()));
}

template <class Arr_>
CGAL::Object ArrangementDemoTab<Arr_>::getArrangement() const
{
  return CGAL::make_object(this->arrangement.get());
}

template <class Arr_>
void ArrangementDemoTab<Arr_>::slotModelChanged()
{
}

static CGAL::Bbox_2 makeFinite(const CGAL::Bbox_2& box)
{
  double xmin = std::numeric_limits<double>::infinity();
  double ymin = std::numeric_limits<double>::infinity();
  double xmax = -std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();
  if (!std::isinf(box.xmin())) xmin = box.xmin();
  if (!std::isinf(box.ymin())) ymin = box.ymin();
  if (!std::isinf(box.xmax())) xmax = box.xmax();
  if (!std::isinf(box.ymax())) ymax = box.ymax();
  return {xmin, ymin, xmax, ymax};
}

static bool isFinite(const CGAL::Bbox_2& box)
{
  return !std::isinf(box.xmin()) && !std::isinf(box.xmax()) &&
         !std::isinf(box.ymin()) && !std::isinf(box.ymax());
}

static CGAL::Bbox_2 addMargins(const CGAL::Bbox_2& box)
{
  // add margin to bounding box
  float x_margin;
  float y_margin;
  if (box.xmin() == box.xmax() || box.ymin() == box.ymax())
  {
    static constexpr float const_margin = 50;
    x_margin = const_margin;
    y_margin = const_margin;
  }
  else
  {
    static constexpr float prop_margin = 0.10;
    x_margin = (box.xmax() - box.xmin()) * prop_margin;
    y_margin = (box.ymax() - box.ymin()) * prop_margin;
  }
  return {
    box.xmin() - x_margin, box.ymin() - y_margin,
    box.xmax() + x_margin, box.ymax() + y_margin};
}

static const auto& getXyCurves()
{
  using Traits = Alg_seg_traits;
  static std::vector<typename Traits::X_monotone_curve_2> xy_curves;
  if (xy_curves.empty())
  {
    Traits traits{};
    typedef typename Traits::Polynomial_2 Polynomial_2;
    auto construct_curve = traits.construct_curve_2_object();
    auto make_x_monotone = traits.make_x_monotone_2_object();

    Polynomial_2 x = CGAL::shift(Polynomial_2(1), 1, 0);
    Polynomial_2 y = CGAL::shift(Polynomial_2(1), 1, 1);
    auto x_cv = construct_curve(x);
    auto y_cv = construct_curve(y);

    std::vector<CGAL::Object> arcs;
    make_x_monotone(x_cv, std::back_inserter(arcs));
    make_x_monotone(y_cv, std::back_inserter(arcs));
    for (auto& arc_obj : arcs)
    {
      typename Traits::X_monotone_curve_2 arc;
      CGAL::assign(arc, arc_obj);
      xy_curves.push_back(arc);
    }
  }
  return xy_curves;
}

template <class Arr_>
CGAL::Bbox_2
findOtherInterestingPoints(const std::unique_ptr<Arr_>& arr)
{
  return {};
}

template <>
CGAL::Bbox_2 findOtherInterestingPoints<Alg_seg_arr>(
  const std::unique_ptr<Alg_seg_arr>& arr)
{
  using Traits = Alg_seg_traits;
  CGAL::Bbox_2 bb = {};
  std::vector<CGAL::Object> intersections;
  for (auto it = arr->edges_begin(); it != arr->edges_end(); ++it)
  {
    for (auto& arc : getXyCurves())
      if (arc.is_vertical() != it->curve().is_vertical())
        it->curve().intersections(arc, std::back_inserter(intersections));
  }
  for (auto it = intersections.begin(); it != intersections.end(); it++)
  {
    std::pair<typename Traits::Point_2, unsigned int> point_multiplicity;
    CGAL::assign(point_multiplicity, *it);
    auto& point = point_multiplicity.first;
    if (point.location() == CGAL::ARR_INTERIOR)
    {
      auto xy = point.to_double();
      bb += makeFinite({xy.first, xy.second, xy.first, xy.second});
    }
  }
  return bb;
}

template <class Arr_>
void ArrangementDemoTab<Arr_>::adjustViewport()
{
  CGAL::Bbox_2 bb = {};
  for (auto it = arrangement->vertices_begin();
       it != arrangement->vertices_end(); it++)
  {
    double x = CGAL::to_double(it->point().x());
    double y = CGAL::to_double(it->point().y());
    bb += makeFinite({x, y, x, y});
  }

  for (auto it = arrangement->edges_begin(); it != arrangement->edges_end();
       ++it)
  {
    // can throws "CGAL::internal::Zero_resultant_exception"
    try { bb += makeFinite(it->curve().bbox()); }
    catch (...) { }
  }

  if (!isFinite(bb))
    bb += findOtherInterestingPoints(this->arrangement);

  // this should happen only if the arrangement is empty
  if (!isFinite(bb))
    bb += {0, 0, 0, 0};

  bb = addMargins(bb);

  auto viewportRect =
    QRectF(bb.xmin(), bb.ymin(), bb.xmax() - bb.xmin(), bb.ymax() - bb.ymin());

  this->graphicsView->resetTransform();
  // TODO: Find suitable values
  double xmin = viewportRect.x() - std::numeric_limits<double>::max() / 1048576;
  double ymin = viewportRect.y() - std::numeric_limits<double>::max() / 1048576;
  double wh = std::numeric_limits<double>::max() / 524288;
  this->scene->setSceneRect(xmin, ymin, wh, wh);
  this->graphicsView->fitInView(viewportRect, Qt::KeepAspectRatio);
}

template class ArrangementDemoTab<Seg_arr>;
template class ArrangementDemoTab<Pol_arr>;
template class ArrangementDemoTab<Conic_arr>;
template class ArrangementDemoTab<Lin_arr>;
template class ArrangementDemoTab<Alg_seg_arr>;
