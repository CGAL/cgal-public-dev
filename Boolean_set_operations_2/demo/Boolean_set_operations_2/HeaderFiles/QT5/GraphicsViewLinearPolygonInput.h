//keep it as it for now
//it is still not used,so it isn't giving any errors

#ifndef CGAL_QT_GRAPHICS_VIEW_LINEAR_POLYGON_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_LINEAR_POLYGON_INPUT_H

#include <CGAL/auto_link/Qt.h>

#include <QPolygonF>
#include <QPointF>
#include <QGraphicsLineItem> 
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <QT5/LinearPolygons.h>
#include "Typedefs.h"

namespace CGAL {

namespace Qt {

  template <class K>
  class GraphicsViewLinearPolygonInput : public GraphicsViewInput
  {
  public:

    typedef K Kernel ;
    
    typedef CGAL::Gps_segement_traits_2_apurva<K> Gps_traits;
    
    typedef typename Gps_traits::Curve_2            Linear_curve;
    typedef typename Gps_traits::X_monotone_curve_2 Linear_X_monotone_curve;
    typedef typename Gps_traits::Polygon_2          Linear_polygon;
    typedef typename Kernel::Vector_2               Vector ;
    typedef typename Kernel::Point_2            Point;
    //typedef CGAL::Point_2<Linear_kernel>                  Point;
    typedef std::vector<Linear_curve> Linear_curve_vector ;
    
    typedef typename Linear_curve_vector::const_iterator const_linear_curve_iterator ;
    
    typedef Linear_boundary_pieces_graphics_item<Linear_curve_vector> GI ;

    GraphicsViewLinearPolygonInput(QObject* aParent, QGraphicsScene* aScene)
      :
        GraphicsViewInput  ( aParent         )
      , mScene             ( aScene          )
      , mState             ( Start           )
      , mLinearPolygonPen  ( QColor(0,255,0) )
      , mOngoingCurvePen   ( QColor(255,0,0) )
      , mHandlePen         ( QColor(0,0,255) )
      , mLinearGI          ( 0               )
    {
      mOngoingPieceGI = new GI(&mOngoingPieceCtr) ;
      mHandleGI       = new QGraphicsLineItem();
      
      mOngoingPieceGI->setPen(mOngoingCurvePen);
      mHandleGI      ->setPen(mHandlePen);
      
      mHandleGI->setLine(0,0,1,1);
      mHandleGI->hide();
      
      mLinearGI = new GI(&mLinearPolygonPieces) ;
      
      mLinearGI->setPen(mLinearPolygonPen);
      
      mScene->addItem(mOngoingPieceGI);
      mScene->addItem(mHandleGI);
      mScene->addItem(mLinearGI);
    }
    
    ~GraphicsViewLinearPolygonInput()
    {
    }
    
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
    
  protected:

    enum State { Start, PieceStarted, PieceOngoing, HandleOngoing, PieceEnded, CurveEnded } ;
    
    Point cvt ( QPointF const& aP ) const { return Point(aP.x(),aP.y()) ; }

    bool mousePressEvent(QGraphicsSceneMouseEvent *aEvent)
    {
      bool rHandled = false ;
      
      Point lP = cvt(aEvent->scenePos());
      
      if ( aEvent->button() == ::Qt::LeftButton )
      {
        switch (mState)
        {
          case Start: 
            mP0      = lP;
            mState   = PieceStarted;
            rHandled = true;
            break;

          case PieceOngoing: 
            mP1      = lP;
            mState   = HandleOngoing;
            rHandled = true;
            break;
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
        case PieceOngoing: 
          mP1 = lP;
          UpdateOngoingPiece();
          rHandled = true ;
          break;

          
        case HandleOngoing:
        
          UpdateHandle(lP);
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
          case PieceStarted:
            mState   = PieceOngoing;
            rHandled = true;
            break;
            
          case HandleOngoing: 
            HideHandle();
            CommitOngoingPiece(lP);
            mState   = PieceEnded;
            rHandled = true;
            break;
        }
      }
      else if ( aEvent->button() == ::Qt::RightButton )
      {
        switch (mState)
        {
          case PieceOngoing: 
            CommitCurrLinearPolygon();
            ReStart();
            rHandled = true;
            break;
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
        mState   = mLinearPolygonPieces.size() > 0 ? PieceEnded : Start ;
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

    
    
  private:

    Linear_curve const* ongoing_piece() const { return mOngoingPieceCtr.size() == 1 ? &mOngoingPieceCtr[0] : NULL ; }

    void ReStart()
    {
      mH = boost::optional<Point>();
      mState = Start ;     
    }
    
    void Reset()
    {
      mLinearPolygonPieces.clear();
      mOngoingPieceCtr      .clear();
      mLinearGI->modelChanged();
      mOngoingPieceGI->modelChanged();
      ReStart();
    }
    
    void HideHandle()
    {
      mHandleGI->hide();
    }  

    Linear_curve CreatePiece()
    {
      if ( mH )
      {
        Vector lD = *mH - mP1 ;
        Vector lU = lD * 1.5 ;
        Point  lH = mP1 - lU ;
        return Linear_curve(mP0,lH,mP1); 
      }
      else
      {
        return Linear_curve(mP0,mP1); 
      }
    }
    
    
    void RemoveLastPiece()
    {
      mLinearPolygonPieces.pop_back();
      mOngoingPieceCtr.clear();
      mLinearGI->modelChanged();
      mOngoingPieceGI->modelChanged();
      if ( mLinearPolygonPieces.size() > 0 )
      {
        mP0 = cvt(mLinearPolygonPieces.back().target());
        UpdateOngoingPiece();
      }
      mH = boost::optional<Point>();
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
        mLinearPolygonPieces.push_back( *ongoing_piece() ) ;
        mLinearGI->modelChanged();
        mOngoingPieceCtr.clear();
        mOngoingPieceGI->modelChanged();
        mP0 = mP1 ;
        mP1 = aP ;
        mH = boost::optional<Point>();
      }
    }      

    void UpdateHandle(Point const& aP)
    {
      if ( squared_distance(mP1,aP) >= 4 )
      {
        mH = aP ;
        
        mHandleGI->setLine( to_double(mP1.x()), to_double(mP1.y()), to_double(mH->x()), to_double(mH->y()));
        mHandleGI->show();
      }
      else
      {
        HideHandle();
        mH = boost::optional<Point>();
      }
    }
          
    Point cvt ( typename Linear_curve::Point_2 const& aP ) { return Point( to_double(aP.x()), to_double(aP.y()) ) ; } 
        
    void CommitCurrLinearPolygon()
    {
      GenerateLinearPolygon();

      mOngoingPieceCtr.clear();
      mOngoingPieceGI->modelChanged();
      
      mLinearPolygonPieces.clear();
      mLinearGI->modelChanged() ;
      
      mH = boost::optional<Point>();
      
      HideHandle();
    }
    
    void GenerateLinearPolygon() 
    {
      if ( mLinearPolygonPieces.size() >  0 )
      {
        Gps_traits traits ;
        typename Gps_traits::Make_x_monotone_2 make_x_monotone = traits.make_x_monotone_2_object();
        
        std::vector<Linear_X_monotone_curve> xcvs;

        for ( const_linear_curve_iterator it = mLinearPolygonPieces.begin() ; it != mLinearPolygonPieces.end() ; ++ it )
        {       
          std::vector<CGAL::Object>                 x_objs;
          std::vector<CGAL::Object>::const_iterator xoit;
          
          make_x_monotone ( *it, std::back_inserter (x_objs));
          
          for (xoit = x_objs.begin(); xoit != x_objs.end(); ++xoit) 
          {
            Linear_X_monotone_curve xcv;
            if (CGAL::assign (xcv, *xoit))
              xcvs.push_back (xcv);
          }    
        }
        
      }
    }
    
  private:
  
    QGraphicsScene*    mScene ;
    GI*                mLinearGI ; 
    GI*                mOngoingPieceGI ; 
    QGraphicsLineItem* mHandleGI ;          

    QPen mLinearPolygonPen ;
    QPen mOngoingCurvePen ;
    QPen mHandlePen ;    
    
    Linear_curve_vector mLinearPolygonPieces ;
    Linear_curve_vector mOngoingPieceCtr ;  

    int mState;
    
    Point mP0;
    Point mP1;
    
    boost::optional<Point> mH;
  
  };

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_LINEAR_POLYGON_INPUT_H
