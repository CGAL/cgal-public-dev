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
// Author(s) : Kabir Kedia   <kabirkedia0111@gmail.com>
//             Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_GRAPHICS_VIEW_POLYLINE_INPUT_H
#define CGAL_GRAPHICS_VIEW_POLYLINE_INPUT_H

#include <iostream>
#include <CGAL/auto_link/Qt.h>
#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Arr_polyline_traits_2.h>

#include <QPolygonF>
#include <QPointF>
#include <QGraphicsLineItem>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

#include "QT5/Polyline_curves.h"

namespace CGAL {
namespace Qt {

  template <typename Kernel_>
  class Graphics_view_polyline_input : public GraphicsViewInput {
  public:
    int print_err(int x) {
      std::cout << "here " << x << std::endl;
      return x;
    }

    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel1;
    typedef Kernel1::Point_2                                  Point_2;
    typedef CGAL::Polygon_2<Kernel1>                          Polygon_2;
    typedef CGAL::Polygon_with_holes_2<Kernel1>               Polygon_with_holes_2;
    typedef std::list <Polygon_with_holes_2>                  Pgn_with_holes_2_container;

    typedef Kernel_ Kernel;
    typedef Polyline_traits Gps_traits;
    typedef typename Gps_traits::Curve_2                        Polyline_curve;
    typedef typename Gps_traits::X_monotone_curve_2             Polyline_X_monotone_curve;
    typedef typename Gps_traits::Polygon_2                      Polyline_polygon;
    typedef typename Gps_traits::Point_2                        Polyline_point;
    typedef typename Gps_traits::General_polygon_with_holes_2   Polyline_polygon_with_holes;
    typedef typename Kernel::FT                                 FT;
    typedef typename Kernel::Vector_2                           Vector;
    typedef typename Kernel::Point_2                            Point;
    typedef std::vector <Polyline_curve>                        Polyline_curve_vector;

    typedef typename Polyline_curve_vector::const_iterator
    const_polyline_curve_iterator;


    typedef Polyline_boundary_pieces_graphics_item <Polyline_curve_vector> GI;

    typedef CGAL::Arr_segment_traits_2<Kernel>                                          Segment_traits_2;
    typedef CGAL::Gps_traits_2<CGAL::Arr_polycurve_traits_2<Segment_traits_2> >         Base;
    typedef typename Base::X_monotone_subcurve_2                                        X_monotone_subcurve_2;


    //constructor
    Graphics_view_polyline_input(QObject *aParent, QGraphicsScene *aScene)
      :
      GraphicsViewInput(aParent), mScene(aScene)
      , mOngoingPieceGI(new GI(&mOngoingPieceCtr)), mHandleGI(new QGraphicsLineItem()),
      mPolylinePolygonPen(QColor(0, 255, 0)), mOngoingCurvePen(QColor(255, 0, 0)),
      mHandlePen(QColor(0, 0, 255)), mState(Start), m_bound_rect(true), m_last_polyline(false),
      m_last(false) {
      mOngoingPieceGI->setPen(mOngoingCurvePen);
      mHandleGI->setPen(mHandlePen);

      mHandleGI->setLine(0, 0, 1, 1);
      mHandleGI->hide();

      mPolylineGI = new GI(&mPolylinePolygonPieces);

      mPolylineGI->setPen(mPolylinePolygonPen);

      mScene->addItem(mOngoingPieceGI);
      mScene->addItem(mHandleGI);
      mScene->addItem(mPolylineGI);
    }

    ~Graphics_view_polyline_input() {}

    bool eventFilter(QObject *obj, QEvent *aEvent) {
      bool rHandled = false;
      if (aEvent->type() == QEvent::GraphicsSceneMousePress) {
        rHandled = mousePressEvent(static_cast<QGraphicsSceneMouseEvent *>(aEvent));
      } else if (aEvent->type() == QEvent::GraphicsSceneMouseRelease) {
        rHandled = mouseReleaseEvent(static_cast<QGraphicsSceneMouseEvent *>(aEvent));
      } else if (aEvent->type() == QEvent::GraphicsSceneMouseMove) {
        rHandled = mouseMoveEvent(static_cast<QGraphicsSceneMouseEvent *>(aEvent));
      } else if (aEvent->type() == QEvent::KeyPress) {
        rHandled = keyPressEvent(static_cast<QKeyEvent *>(aEvent));
      }

      if (!rHandled)
        rHandled = QObject::eventFilter(obj, aEvent);

      return rHandled;
    }

  public:
    enum State {
      Start, PieceStarted, PieceOngoing, HandleOngoing, PieceEnded, CurveEnded
    };

    Point cvt(QPointF const &aP) const { return Point(aP.x(), aP.y()); }

    //All functions related to mouse activity
    bool mousePressEvent(QGraphicsSceneMouseEvent *aEvent) {
      bool rHandled = false;
      m_bound_rect = false;
      Point lP = cvt(aEvent->scenePos());
      //left click to complete polyline
      if (aEvent->button() == ::Qt::LeftButton) {
        switch (mState) {
         case Start:
          mP0 = lP;
          mState = PieceStarted;
          rHandled = true;
          break;

         case PieceOngoing:
          mP1 = lP;
          mState = HandleOngoing;
          rHandled = true;
          break;

         default:
          break; //! \todo handle default case
        }
      }
      //right click to complete the polygon
      else if (aEvent->button() == ::Qt::RightButton) {
        switch (mState) {
         case PieceOngoing:
          // allowing user to complete polygon
          m_last = true;
          mState = HandleOngoing;
          rHandled = true;
          break;

         default:
          break; //! \todo handle default case
        }
      }
      return rHandled;
    }

    bool mouseMoveEvent(QGraphicsSceneMouseEvent *aEvent) {
      bool rHandled = false;

      Point lP = cvt(aEvent->scenePos());

      switch (mState) {
       case PieceOngoing:
        mP1 = lP;
        UpdateOngoingPiece();
        rHandled = true;
        break;

       case HandleOngoing:
        if (m_last) {
          Point ps(mPolylinePolygonPieces.front()[0].source().x(),
                   mPolylinePolygonPieces.front()[0].source().y());
          mP1 = ps;
          m_last_polyline = true;
        }
        UpdateHandle(lP);
        UpdateOngoingPiece();
        rHandled = true;
        break;

       case PieceEnded:
        mState = PieceOngoing;
        rHandled = true;
        break;

       default:
        break; //! \todo handle default case
      }

      return rHandled;
    }

    bool mouseReleaseEvent(QGraphicsSceneMouseEvent *aEvent) {
      bool rHandled = false;
      Point lP = cvt(aEvent->scenePos());
      if (aEvent->button() == ::Qt::LeftButton) {
        switch (mState) {
         case PieceStarted:
          mState = PieceOngoing;
          rHandled = true;
          break;

         case HandleOngoing:
          HideHandle();
          CommitOngoingPiece(lP);
          mState = PieceEnded;
          rHandled = true;
          break;

         default:
          break; //! \todo handle default case
        }
      } else if (aEvent->button() == ::Qt::RightButton) {
        switch (mState) {
         case HandleOngoing:
          if (m_last_polyline) {
            HideHandle();
            CommitOngoingPiece(lP);
          }
          m_bound_rect = false;
          CommitCurrPolylinePolygon();
          ReStart();
          rHandled = true;
          //  cout << "right click over" << endl;
          break;

         default:
          break; //! \todo handle default case
        }
      }
      return rHandled;
    }

    bool keyPressEvent(QKeyEvent *aEvent) {
      bool rHandled = false;

      if ((aEvent->key() == ::Qt::Key_Delete) ||
          (aEvent->key() == ::Qt::Key_Backspace)) {
        RemoveLastPiece();
        mState = (mPolylinePolygonPieces.size() > 0) ? PieceEnded : Start;
        rHandled = true;
      } else if (aEvent->key() == ::Qt::Key_Escape) {
        Reset();
        mState = Start;
        rHandled = true;
      }
      return rHandled;
    }

  public:

    Polyline_curve const *ongoing_piece() const {
      return (mOngoingPieceCtr.size() == 1) ? &mOngoingPieceCtr[0] : NULL;
    }

    void ReStart() {
      mH = std::optional<Point>();
      mState = Start;
    }

    void Reset() {
      mPolylinePolygonPieces.clear();
      mOngoingPieceCtr.clear();
      mPolylineGI->modelChanged();
      mOngoingPieceGI->modelChanged();
      ReStart();
    }

    void HideHandle() {
      mHandleGI->hide();
    }

    //change this
    Polyline_curve CreatePiece() {
      Gps_traits poly_tr;
      auto construct_poly = poly_tr.construct_curve_2_object();
      return Polyline_curve(construct_poly(mP0, mP1));
    }

    void RemoveLastPiece() {
      mPolylinePolygonPieces.pop_back();
      mOngoingPieceCtr.clear();
      mPolylineGI->modelChanged();
      mOngoingPieceGI->modelChanged();
      if (mPolylinePolygonPieces.size() > 0) {
        auto xx = mPolylinePolygonPieces.back().number_of_subcurves();
        Point pt(mPolylinePolygonPieces.back()[xx - 1].target().x(),
                 mPolylinePolygonPieces.back()[xx - 1].target().y());
        mP0 = pt;
        UpdateOngoingPiece();
      }
      mH = std::optional<Point>();
    }

    void UpdateOngoingPiece() {
      if (mOngoingPieceCtr.size() > 0) mOngoingPieceCtr.clear();
      mOngoingPieceCtr.push_back(CreatePiece());
      //cout<<"hi"<<endl;
      mOngoingPieceGI->modelChanged();
    }

    void CommitOngoingPiece(Point const &aP) {
      if (ongoing_piece()) {
        mPolylinePolygonPieces.push_back(*ongoing_piece());
        mPolylineGI->modelChanged();
        mOngoingPieceCtr.clear();
        mOngoingPieceGI->modelChanged();
        mP0 = mP1;
        mP1 = aP;
        mH = std::optional<Point>();
      }
    }

    void UpdateHandle(Point const &aP) {
      if (squared_distance(mP1, aP) >= 4) {
        mH = aP;
        mHandleGI->setLine(to_double(mP1.x()), to_double(mP1.y()),
                           to_double(mH->x()), to_double(mH->y()));
        mHandleGI->show();
      } else {
        HideHandle();
        mH = std::optional<Point>();
      }
    }

    Point cvt(QPointF const &aP) { return Point(aP.x(), aP.y()); }

    void CommitCurrPolylinePolygon() {
      GeneratePolylinePolygon();
      mOngoingPieceCtr.clear();

      mOngoingPieceGI->modelChanged();

      mPolylinePolygonPieces.clear();
      mPolylineGI->modelChanged();

      mH = std::optional<Point>();

      HideHandle();
    }

    //main function to generate polyline polygon
    void GeneratePolylinePolygon() {
      if (mPolylinePolygonPieces.size() > 0) {
        Gps_traits traits;
        auto make_x_monotone = traits.make_x_monotone_2_object();
        using Pnt = typename Gps_traits::Point_2;
        using Xcv = typename Gps_traits::X_monotone_curve_2;
        using Make_x_monotone_result = std::variant<Pnt, Xcv>;

        std::vector <Polyline_X_monotone_curve> xcvs;
        for (auto it = mPolylinePolygonPieces.begin();
             it != mPolylinePolygonPieces.end(); ++it) {
          std::vector<Make_x_monotone_result> x_objs;
          make_x_monotone(*it, std::back_inserter(x_objs));

          for (auto i = 0; i < x_objs.size(); ++i) {
            auto* xcv = std::get_if<Polyline_X_monotone_curve > (&x_objs[i]);
            CGAL_assertion(xcv != nullptr);
            xcvs.push_back(*xcv);
          }
        }
        if (xcvs.size() > 0) {
          if (! m_last_polyline) {
            Polyline_point const &first_point = xcvs.front()[0].source();
            Polyline_point const &last_point = xcvs.back()[xcvs.back().number_of_subcurves() -1].target();
            FT fxs = first_point.x();
            FT fys = first_point.y();
            FT lxs = last_point.x();
            FT lys = last_point.y();
            //cout<<"point"<<endl;
            Gps_traits x;
            X_monotone_subcurve_2 seg = x.subcurve_traits_2()->
              construct_x_monotone_curve_2_object()(Point(lxs,lys),Point(fxs,fys));//

            xcvs.push_back(Polyline_X_monotone_curve(seg));
          }
          m_last = false;
          m_last_polyline = false;
          /*for(auto &i:xcvs)
            cout<<i<<endl;
            cout<<endl;*/
          Polyline_polygon pp(xcvs.begin(), xcvs.end());
          emit(generate(std::variant<Polyline_polygon>(pp)));
        }
      }
    }

    //for complement operations
    void CommitComplementCurrPolylinePolygon()
    {
      GenerateComplementPolylinePolygon();
      mOngoingPieceCtr.clear();
      mOngoingPieceGI->modelChanged();

      mPolylinePolygonPieces.clear();
      mPolylineGI->modelChanged();

      mH = std::optional<Point>();

      HideHandle();
    }

    //for complement operations
    void GenerateComplementPolylinePolygon() {
      if (mPolylinePolygonPieces.size() > 0) {
        Gps_traits traits;
        auto make_x_monotone = traits.make_x_monotone_2_object();
        using Pnt = typename Gps_traits::Point_2;
        using Xcv = typename Gps_traits::X_monotone_curve_2;
        using Make_x_monotone_result = std::variant<Pnt, Xcv>;

        std::vector <Polyline_X_monotone_curve> xcvs;
        for (auto it = mPolylinePolygonPieces.begin();
             it != mPolylinePolygonPieces.end(); ++it) {
          std::vector<Make_x_monotone_result> x_objs;
          make_x_monotone(*it, std::back_inserter(x_objs));

          for(auto i = 0; i < x_objs.size(); ++i) {
            auto* xcv = std::get_if<Polyline_X_monotone_curve>(&x_objs[i]);
            CGAL_assertion(xcv != nullptr);
            xcvs.push_back(*xcv);
          }
        }
        if (xcvs.size() > 0) {
          if (!m_last_polyline) {
            Polyline_point const &first_point = xcvs.front()[0].source();
            Polyline_point const &last_point = xcvs.back()[xcvs.back().number_of_subcurves() -1].target();
            FT fxs = first_point.x();
            FT fys = first_point.y();
            FT lxs = last_point.x();
            FT lys = last_point.y();
            //cout<<fxs<<" "<<fys<<endl<<lxs<<" "<<lys<<"\n\n";
            //comment out the below two lines to see why they are important
            fxs=fxs*-1;
            lxs=lxs*-1;
            Gps_traits x;
            X_monotone_subcurve_2 seg = x.subcurve_traits_2()->
              construct_x_monotone_curve_2_object()(Point(fxs,fys),Point(lxs,lys));

            xcvs.push_back(Polyline_X_monotone_curve(seg));
          }

          m_last = false;
          m_last_polyline = false;

          Polyline_polygon pp(xcvs.begin(), xcvs.end());
          emit(generate(std::variant<Polyline_polygon>(pp)));
        }
      }
    }

    void get_BoundingRect() {
      m_bound_rect = true;

      mP0 = Point(-15500000, -10000000);
      mState = PieceStarted;

      mState = PieceOngoing;
      mP1 = Point(-15500000, 10000000);
      UpdateOngoingPiece();

      mP1 = Point(-15500000, 10000000);
      mState = HandleOngoing;
      HideHandle();
      CommitOngoingPiece(Point(-15500000, 10000000));
      mState = PieceEnded;

      mState = PieceOngoing;
      mP1 = Point(15500000, 10000000);
      UpdateOngoingPiece();

      mP1 = Point(15500000, 10000000);
      mState = HandleOngoing;
      HideHandle();
      CommitOngoingPiece(Point(15500000, 10000000));
      mState = PieceEnded;

      mState = PieceOngoing;
      mP1 = Point(15500000, -10000000);
      UpdateOngoingPiece();

      mP1 = Point(15500000, -10000000);
      mState = HandleOngoing;
      HideHandle();
      CommitOngoingPiece(Point(15500000, -10000000));
      mState = PieceEnded;

      mState = PieceOngoing;
      mP1 = Point(-9000000, -9000000);
      UpdateOngoingPiece();
      CommitComplementCurrPolylinePolygon();
      ReStart();
    }

    bool isboundingRect() {
      return m_bound_rect;
    }

  public:
    QGraphicsScene *mScene;
    GI *mPolylineGI;
    GI *mOngoingPieceGI;
    QGraphicsLineItem *mHandleGI;

    QPen mPolylinePolygonPen;
    QPen mOngoingCurvePen;
    QPen mHandlePen;

    bool m_bound_rect;
    bool m_last_polyline;
    bool m_last;

    Polyline_curve_vector mPolylinePolygonPieces;
    Polyline_curve_vector mOngoingPieceCtr;

    int mState;

    Point mP0;
    Point mP1;

    std::optional <Point> mH;
  };

} // namespace Qt
} // namespace CGAL

#endif //CGAL_GRAPHICS_VIEW_POLYLINE_INPUT_H
