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
// Author(s) : Apurva Bhatt <response2apurva@gmail.com>
//             Ronnie Gandhi <ronniegandhi19999@gmail.com>
//             Kabir Kedia   <kabirkedia0111@gmail.com>
//             Efi Fogel <efifogel@gmain.com>

#ifndef CGAL_GRAPHICS_VIEW_POLYLINE_INPUT_H
#define CGAL_GRAPHICS_VIEW_POLYLINE_INPUT_H

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
        template<typename Kernel_>
        class Graphics_view_polyline_input : public GraphicsViewInput {
        public:

            typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel1;
            typedef Kernel1::Point_2 Point_2;
            typedef CGAL::Polygon_2<Kernel1> Polygon_2;
            typedef CGAL::Polygon_with_holes_2<Kernel1> Polygon_with_holes_2;
            typedef std::list <Polygon_with_holes_2> Pgn_with_holes_2_container;

            typedef Kernel_ Kernel;

            typedef Polyline_traits Gps_traits;
            typedef typename Gps_traits::Curve_2 Polyline_curve;
            typedef typename Gps_traits::X_monotone_curve_2 Polyline_X_monotone_curve;
            typedef typename Gps_traits::Polygon_2 Polyline_polygon;
            typedef typename Gps_traits::Point_2 Polyline_point;
            typedef typename Kernel::FT FT;
            typedef typename Kernel::Vector_2 Vector;
            typedef typename Kernel::Point_2 Point;
            typedef std::vector <Polyline_curve> Polyline_curve_vector;

            typedef typename Polyline_curve_vector::const_iterator
                    const_polyline_curve_iterator;

            typedef Polyline_boundary_pieces_graphics_item <Polyline_curve_vector> GI;

            //constructor
            Graphics_view_polyline_input(QObject* aParent, QGraphicsScene* aScene)
            :
            GraphicsViewInput( aParent         )
            , mScene             ( aScene          )
            , mState             ( Start           )
            , mPolylinePolygonPen( QColor(0,255,0) )
            , mOngoingCurvePen   ( QColor(255,0,0) )
            , mHandlePen         ( QColor(0,0,255) )
            , mPolylineGI        ( 0               )
            , m_bound_rect       ( true            )
            , m_last_polyline    ( false           )
            , m_last             ( false           )
            {
                mOngoingPieceGI = new GI(&mOngoingPieceCtr) ;
                mHandle0GI      = new QGraphicsLineItem();
                mHandle1GI      = new QGraphicsLineItem();

                mOngoingPieceGI->setPen(mOngoingCurvePen);
                mHandle0GI     ->setPen(mHandlePen);
                mHandle1GI     ->setPen(mHandlePen);

                mHandle0GI->setLine(0,0,1,1);
                mHandle1GI->setLine(0,0,1,1);
                mHandle0GI->hide();
                mHandle1GI->hide();

                mPolylineGI = new GI(&mPolylinePolygonPieces) ;

                mPolylineGI->setPen(mPolylinePolygonPen);

                mScene->addItem(mOngoingPieceGI);
                mScene->addItem(mHandle0GI);
                mScene->addItem(mHandle1GI);
                mScene->addItem(mPolylineGI);
            }

            ~Graphics_view_polyline_input(){}

            bool eventFilter(QObject *obj, QEvent *aEvent)
            {
                bool rHandled = false ;

                if (aEvent->type() == QEvent::GraphicsSceneMousePress)
                {
                    rHandled = mousePressEvent( static_cast<QGraphicsSceneMouseEvent *>(aEvent) ) ;
                }
                else if (aEvent->type() == QEvent::GraphicsSceneMouseRelease)
                {
                    rHandled = mouseReleaseEvent( static_cast<QGraphicsSceneMouseEvent *>(aEvent) ) ;
                }
                else if (aEvent->type() == QEvent::GraphicsSceneMouseMove)
                {
                    rHandled = mouseMoveEvent( static_cast<QGraphicsSceneMouseEvent *>(aEvent) ) ;
                }
                else if (aEvent->type() == QEvent::KeyPress)
                {
                    rHandled = keyPressEvent( static_cast<QKeyEvent *>(aEvent) ) ;
                }

                if ( !rHandled )
                    rHandled = QObject::eventFilter(obj, aEvent);

                return rHandled ;
            }

        public:

            enum State { Start, PieceOrFirstHandleStarted, PieceOngoing, FirstHandleOngoing, HandleOngoing, PieceEnded, CurveEnded } ;

            Point cvt ( QPointF const& aP ) const { return Point(aP.x(),aP.y()) ; }

            bool mousePressEvent(QGraphicsSceneMouseEvent *aEvent)
            {

                bool rHandled = false;
                m_bound_rect = false;

                Point lP = cvt(aEvent->scenePos());

                if ( aEvent->button() == ::Qt::LeftButton )
                {
                    switch (mState)
                    {
                        case Start:
                            mP0      = lP;
                            mState   = PieceOrFirstHandleStarted;
                            rHandled = true;
                            break;

                        case PieceOngoing:
                            mP1      = lP;
                            mState   = HandleOngoing;
                            rHandled = true;
                            break;

                        default: break; //!todo handle default case
                    }
                }

                else  if (aEvent->button() == ::Qt::RightButton) {
                    switch (mState) {
                        case PieceOngoing:
                            // allowing user to curve last piece as well
                            m_last = true;
                            mState = HandleOngoing;
                            rHandled = true;
                            break;

                        default: break; //! \todo handle default case
                    }
                }

                return rHandled ;
            }


            bool mouseMoveEvent(QGraphicsSceneMouseEvent *aEvent)
            {
                bool rHandled = false ;
                Point lP = cvt(aEvent->scenePos());

                switch (mState)
                {
                    case PieceOrFirstHandleStarted:
                        mState   = FirstHandleOngoing;
                        rHandled = true;
                        break;

                    case PieceOngoing:
                        mP1 = lP;
                        UpdateOngoingPiece();
                        rHandled = true ;
                        break;

                    case FirstHandleOngoing:
                        UpdateVeryFirstHandle(lP);
                        rHandled = true ;
                        break;

                    case HandleOngoing:
                        if(m_last)
                        {
                            mP1 = mPolylinePolygonPieces.front().points(0);
                            m_last_polyline = true;
                        }
                        UpdateHandles(lP);
                        UpdateOngoingPiece();
                        rHandled = true ;
                        break;

                    case PieceEnded:
                        mState   = PieceOngoing;
                        rHandled = true;
                        break;
                }

                return rHandled ;
            }

            bool mouseReleaseEvent(QGraphicsSceneMouseEvent *aEvent)
            {
                bool rHandled = false ;

                Point lP = cvt(aEvent->scenePos());

                if ( aEvent->button() == ::Qt::LeftButton )
                {
                    switch (mState)
                    {
                        case PieceOrFirstHandleStarted:
                            mState   = PieceOngoing;
                            rHandled = true;
                            break;

                        case FirstHandleOngoing:
                            UpdateVeryFirstHandle(lP);
                            mPrevH0  = mH1 ;
                            mH1      = boost::optional<Point>();
                            mState   = PieceOngoing;
                            rHandled = true;
                            break;

                        case HandleOngoing:
                            UpdateHandles(lP);
                            CommitOngoingPiece(lP);
                            mState   = PieceEnded;
                            rHandled = true;
                            break;

                        default: break; //!todo add default case handling
                    }
                }
                else if ( aEvent->button() == ::Qt::RightButton )
                {
                    switch (mState)
                    {
                        case HandleOngoing:
                            m_bound_rect = false;
                            if(m_last_polyline)
                            {
                                HideHandles();
                                CommitOngoingPiece(lP);
                            }
                            CloseCurrBoundary();
                            CommitCurrPolylinePolygon();
                            ReStart();
                            rHandled = true;
                            break;

                        default: break; //!todo handle default case
                    }
                }
                return rHandled ;
            }

            bool keyPressEvent(QKeyEvent *aEvent)
            {
                bool rHandled = false ;

                if( aEvent->key() == ::Qt::Key_Delete || aEvent->key() == ::Qt::Key_Backspace )
                {
                    RemoveLastPiece();
                    mState   = mPolylinePolygonPieces.size() > 0 ? PieceEnded : Start ;
                    rHandled = true;
                }
                else if( aEvent->key() == ::Qt::Key_Escape)
                {
                    Reset();
                    mState   = Start;
                    rHandled = true;
                }
                return rHandled ;
            }


        public:

            Polyline_curve const* ongoing_piece() const { return mOngoingPieceCtr.size() == 1 ? &mOngoingPieceCtr[0] : NULL ; }

            void ReStart()
            {
                mPrevH0 = mH0 = mH1 = boost::optional<Point>();
                mState = Start ;
            }

            void Reset()
            {
                mPolylinePolygonPieces.clear();
                mOngoingPieceCtr      .clear();
                mPolylineGI    ->modelChanged();
                mOngoingPieceGI->modelChanged();
                ReStart();
            }

            void RemoveLastPiece()
            {
                mPolylinePolygonPieces.pop_back();
                mOngoingPieceCtr      .clear();
                mPolylineGI      ->modelChanged();
                mOngoingPieceGI->modelChanged();
                if ( mPolylinePolygonPieces.size() > 0 )
                {
                    mP0 = mPolylinePolygonPieces.back().control_point(mPolylinePolygonPieces.back().number_of_control_points()-1);
                    UpdateOngoingPiece();
                }
                mPrevH0 = mH0 = mH1 = boost::optional<Point>();
            }

            void HideHandles()
            {
                mHandle0GI->hide();
                mHandle1GI->hide();
            }

            Polyline_curve CreatePiece()
            {
                if ( mPrevH0 && mH1 && *mPrevH0 != *mH1 && *mPrevH0 != mP0 && *mH1 != mP1 )
                {
                    Point lControlPoints[4] = { mP0
                            , *mPrevH0
                            , *mH1
                            , mP1
                    } ;
                    return Polyline_curve( lControlPoints, lControlPoints + 4 ) ;
                }
                else if ( mPrevH0 && !mH1 && *mPrevH0 != mP0 && *mPrevH0 != mP1 )
                {
                    Point lControlPoints[3] = { mP0
                            , *mPrevH0
                            , mP1
                    } ;
                    return Polyline_curve ( lControlPoints, lControlPoints + 3 ) ;
                }
                else if ( !mPrevH0 && mH1 && *mH1 != mP0 && *mH1 != mP1 )
                {
                    Point lControlPoints[3] = { mP0
                            , *mH1
                            , mP1
                    } ;
                    return Polyline_curve ( lControlPoints, lControlPoints + 3 ) ;
                }
                else
                {
                    Point lControlPoints[2] = { mP0
                            , mP1
                    } ;
                    return Polyline_curve( lControlPoints, lControlPoints + 2 ) ;
                }
            }

            void UpdateOngoingPiece()
            {
                if ( mOngoingPieceCtr.size() > 0 )
                    mOngoingPieceCtr.clear();
                mOngoingPieceCtr.push_back(CreatePiece());
                mOngoingPieceGI->modelChanged();
            }

            void CommitOngoingPiece( Point const& aP )
            {
                if ( ongoing_piece() )
                {
                    mPolylinePolygonPieces.push_back( *ongoing_piece() );
                    mPolylineGI->modelChanged();
                    mOngoingPieceCtr.clear();
                    mOngoingPieceGI->modelChanged();
                    mP0 = mP1 ;
                    mP1 = aP ;
                    mPrevH0 = mH0 ;
                    mH0 = mH1 = boost::optional<Point>();
                }
            }

            void UpdateVeryFirstHandle(Point const& aP)
            {
                if ( squared_distance(mP0,aP) >= 9 )
                {
                    mH1 = aP ;
                    mHandle1GI->setLine( to_double(mP0.x()), to_double(mP0.y()),
                                         to_double(mH1->x()), to_double(mH1->y()));
                    mHandle1GI->show();

                    mH0 = boost::optional<Point>();
                    mHandle0GI->hide();
                }
                else
                {
                    HideHandles();
                    mH0 = mH1 = boost::optional<Point>();
                }
            }

            void UpdateHandles(Point const& aP)
            {
                if ( squared_distance(mP1,aP) >= 9 )
                {
                    mH0 = aP ;
                    mH1 = mP1 - (aP - mP1);

                    mHandle0GI->setLine( to_double(mP1.x()), to_double(mP1.y()), to_double(mH0->x()), to_double(mH0->y()));
                    mHandle1GI->setLine( to_double(mP1.x()), to_double(mP1.y()), to_double(mH1->x()), to_double(mH1->y()));
                    mHandle0GI->show();
                    mHandle1GI->show();
                }
                else
                {
                    HideHandles();
                    mH0 = mH1 = boost::optional<Point>();
                }
            }

            void CloseCurrBoundary()
            {
                if ( mPolylinePolygonPieces.size() > 0 && ongoing_piece()!= NULL  && !m_last_polyline)
                {
                    std::vector<Point> lControlPoints(ongoing_piece()->control_points_begin(),ongoing_piece()->control_points_end());

                    lControlPoints.back() = mPolylinePolygonPieces.front().control_points(0);

                    mPolylinePolygonPieces.push_back( Polyline_curve( lControlPoints.begin(), lControlPoints.end() ) ) ;

                    mPolylineGI->modelChanged() ;
                }
            }

            void CommitCurrPolylinePolygon()
            {
                GeneratePolylinePolygon();

                mOngoingPieceCtr.clear();
                mOngoingPieceGI->modelChanged();

                mPolylinePolygonPieces.clear();
                mPolylineGI->modelChanged() ;

                mPrevH0 = mH0 = mH1 = boost::optional<Point>();

                HideHandles();
            }

            void GeneratePolylinePolygon()
            {
                Gps_traits traits ;
                typename Gps_traits::Make_x_monotone_2 make_x_monotone = traits.make_x_monotone_2_object();

                std::vector<Polyline_X_monotone_curve> xcvs;

                for ( const_polyline_curve_iterator it = mPolylinePolygonPieces.begin() ; it != mPolylinePolygonPieces.end() ; ++ it )
                {
                    std::vector<CGAL::Object>                 x_objs;
                    std::vector<CGAL::Object>::const_iterator xoit;

                    make_x_monotone ( *it, std::back_inserter (x_objs));

                    Polyline_X_monotone_curve xcv;
                    xoit = x_objs.begin();
                    CGAL::assign(xcv,*xoit);

                    for (xoit = x_objs.begin(); xoit != x_objs.end(); ++xoit)
                    {
                        if (CGAL::assign (xcv, *xoit))
                            xcvs.push_back (xcv);
                    }
                }

                if ( xcvs.size() > 0 )
                {
                    m_last = false;
                    m_last_polyline = false;
                    Polyline_polygon bp(xcvs.begin(), xcvs.end());
                    emit(generate(CGAL::make_object( std::make_pair(bp,mPolylinePolygonPieces))));
                }
            }

            void get_BoundingRect()
            {

                m_bound_rect = true;

                mP0 = Point(-15500000,-10000000);
                mState  = PieceOrFirstHandleStarted;

                mState   = PieceOngoing;
                mP1 = Point(-15500000,10000000);
                UpdateOngoingPiece();

                mP1 = Point(-15500000,10000000);
                mState = HandleOngoing;
                UpdateHandles(Point(-15500000,10000000));
                CommitOngoingPiece(Point(-15500000,10000000));
                mState   = PieceEnded;

                mState   = PieceOngoing;
                mP1 = Point(15500000,10000000);
                UpdateOngoingPiece();

                mP1 = Point(15500000,10000000);
                mState = HandleOngoing;
                UpdateHandles(Point(15500000,10000000));
                CommitOngoingPiece(Point(15500000,10000000));
                mState   = PieceEnded;

                mState   = PieceOngoing;
                mP1 = Point(15500000,-10000000);
                UpdateOngoingPiece();

                mP1 = Point(15500000,-10000000);
                mState = HandleOngoing;
                UpdateHandles(Point(15500000,-10000000));
                CommitOngoingPiece(Point(15500000,-10000000));
                mState   = PieceEnded;

                mState   = PieceOngoing;
                mP1 = Point(-9000000,-9000000);
                UpdateOngoingPiece();

                CloseCurrBoundary();
                CommitCurrPolylinePolygon();
                ReStart();
            }

            bool isboundingRect()
            {
                return m_bound_rect;
            }

        public:
            QGraphicsScene*    mScene ;
            GI*                mPolylineGI ;
            GI*                mOngoingPieceGI ;
            QGraphicsLineItem* mHandle0GI ;
            QGraphicsLineItem* mHandle1GI ;

            QPen mPolylinePolygonPen ;
            QPen mOngoingCurvePen ;
            QPen mHandlePen ;

            bool m_bound_rect;
            bool m_last_polyline;
            bool m_last;

            Polyline_curve_vector  mPolylinePolygonPieces ;
            Polyline_curve_vector  mOngoingPieceCtr ;

            int mState;

            Point mP0;
            Point mP1;

            boost::optional<Point> mPrevH0;
            boost::optional<Point> mH0;
            boost::optional<Point> mH1;

        }; // end class Graphics_view_polyline_input

    }//namespace Qt
}//namespace CGAL


#endif //CGAL_GRAPHICS_VIEW_POLYLINE_INPUT_H
