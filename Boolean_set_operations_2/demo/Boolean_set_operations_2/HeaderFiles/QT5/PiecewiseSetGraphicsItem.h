#ifndef CGAL_QT_PIECEWISE_SET_GRAPHICS_ITEM_H
#define CGAL_QT_PIECEWISE_SET_GRAPHICS_ITEM_H

#include <CGAL/Gps_circle_segment_traits_2.h>
#include <QT5/PiecewiseRegionGraphicsItem.h>

namespace CGAL {

namespace Qt {

  namespace internal
  {
    template<class Piecewise_set> 
    struct Piecewise_set_traits 
    {
      typedef typename Piecewise_set::Base base ;
      
      typedef typename base::Polygon_with_holes_2 Region ;
    } ;
    
    template<class K, class C, class D>
    struct Piecewise_set_traits< Polygon_set_2<K,C,D> >
    {
      typedef Polygon_set_2<K,C,D> PS ;
      
      typedef typename PS::Polygon_with_holes_2 Region ;
    } ;
  }

template <class Piecewise_set_, class Draw_piece_, class Piece_bbox_>
class Piecewise_set_graphics_item : public Piecewise_region_graphics_item< typename internal::Piecewise_set_traits<Piecewise_set_>::Region, Draw_piece_, Piece_bbox_ > 
{
  typedef Piecewise_set_ Piecewise_set ;
  typedef Draw_piece_    Draw_piece ;
  typedef Piece_bbox_    Piece_bbox ;
  
  typedef typename internal::Piecewise_set_traits<Piecewise_set_>::Region Region ;
  
  typedef Piecewise_region_graphics_item<Region, Draw_piece, Piece_bbox> Base ;
 
  typedef std::vector<Region> Region_vector ;

  typedef typename Region_vector::const_iterator Region_const_iterator ;
  
public:

  Piecewise_set_graphics_item( Piecewise_set* aSet, Draw_piece const& aPieceDrawer = Draw_piece(), Piece_bbox const& aPieceBBox = Piece_bbox() )
    :
     Base(aPieceDrawer,aPieceBBox)
    ,mSet(aSet)
  {}  

public:

  virtual bool isModelEmpty() const { return !mSet || mSet->is_empty() ; }
  
protected:
  
  virtual void update_bbox( Piecewise_graphics_item_base::Bbox_builder& aBboxBuilder)
  {
    if ( mSet ) 
      update_set_bbox(*mSet, aBboxBuilder ) ;
  }    

  virtual void draw_model ( QPainterPath& aPath ) 
  {
    if ( mSet )
      draw_set(*mSet,aPath);  
  }

  void update_set_bbox( Piecewise_set const& aSet, Piecewise_graphics_item_base::Bbox_builder& aBboxBuilder ) ;
  void draw_set       ( Piecewise_set const& aSet, QPainterPath& aPath ) ;
  
protected:

  Piecewise_set* mSet;
};

template <class S, class D, class P>
void Piecewise_set_graphics_item<S,D,P>::update_set_bbox( Piecewise_set const& aSet, Piecewise_graphics_item_base::Bbox_builder& aBboxBuilder )
{
  Region_vector vec ;
  
  aSet.polygons_with_holes( std::back_inserter(vec) ) ;
  
  for( Region_const_iterator rit = vec.begin(); rit != vec.end() ; ++ rit )
    this->update_region_bbox(*rit,aBboxBuilder);
}

template <class S, class D, class P>
void Piecewise_set_graphics_item<S,D,P>::draw_set( Piecewise_set const& aSet, QPainterPath& aPath )
{
  Region_vector vec ;
  
  aSet.polygons_with_holes( std::back_inserter(vec) ) ;
  
  for( Region_const_iterator rit = vec.begin(); rit != vec.end() ; ++ rit )
    this->draw_region(*rit,aPath);
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PIECEWISE_SET_GRAPHICS_ITEM_H
