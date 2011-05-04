// Self_intersection_handler_2.h
#ifndef Self_intersection_handler_2_h
#define Self_intersection_handler_2_h

#include <CGAL/Offset_statistics_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Circle_approximation_2.h>
#include <CGAL/Minkowski_sum_construction_2.h>

namespace CGAL {

/*! \classf
 * A class implementing detection and handling of self-intersections
 * that might be caused by replacement of the kgon sloped edges with
 * corresponding arcs.
 * Implemented via arrangement observer, used to receive notifications of
 * arc/segment or arc/arc intersections and handle them
 */

// declare template class existence
template <typename CPolygon_2>
class Minkowski_sum_construction_2;

template <typename CPolygon_2>
class Self_intersection_handler_2 : public Offset_statistics_2<CPolygon_2>
{
public:

  typedef Offset_statistics_2<CPolygon_2> Base;
  typedef Minkowski_sum_construction_2<CPolygon_2> Min_sum_construction_2;
  typedef typename Min_sum_construction_2::Circle_app_2 Circle_app_2;
  
  typedef typename Base::Types Types;
  typedef typename Types::Polygon_2 Polygon_2;
  typedef typename Types::Polygon_traits_2 Polygon_traits_2;
  typedef typename Types::Input_rational Input_rational;

  typedef typename Polygon_traits_2::Circle_2 Circle_2;
  typedef typename std::list<Circle_2> Circle_list_2;
  
  typedef typename Types::Approximate_offset_polygon_2 Approximate_offset_polygon_2;
  typedef typename Types::Approximate_polygon_2 Approximate_polygon_2;
  typedef typename Types::Approximate_polygon_list_2 Approximate_polygon_list_2;

  typedef typename CGAL::Polygon_with_holes_2<Polygon_traits_2>  Kgon_sum_polygon_2;
  
  typedef typename Min_sum_construction_2::Arc_center_multimap_2 Arc_center_multimap_2;
  typedef typename Arc_center_multimap_2::iterator Ac_multimap_2_iterator;

  typedef typename Types::Rat_traits_2 Rat_traits_2;
  typedef typename Rat_traits_2::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Rat_traits_2::Point_2  Arc_point_2;
  typedef typename Rat_traits_2::CoordNT One_root_nt;
  typedef typename CGAL::Cartesian<One_root_nt> One_root_kernel_2;
  typedef typename One_root_kernel_2::Point_2 One_root_point_2;
  
  typedef typename CGAL::Arrangement_2<Rat_traits_2>  Arrangement_2;
  typedef typename CGAL::Arr_observer<Arrangement_2>  Arr_observer_2;
  
  typedef typename Arrangement_2::Point_2 Point_2;  // the point type.
  //typedef typename Arrangement_2::X_monotone_curve_2  X_monotone_curve_2;   // the x-monotone curve type.

  typedef typename Arrangement_2::Vertex_handle Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle Halfedge_handle;
  typedef typename Arrangement_2::Face_handle Face_handle;
  
  typedef typename std::vector<typename Arrangement_2::Face_const_iterator> Face_vector;
 // typedef typename std::vector<typename Arrangement_2::Face*> Face_vector;
    
  static const Input_rational* smp_offset;
  static Circle_list_2* smp_circles;
  static Arc_center_multimap_2* smp_ac_multimap;
  
  Self_intersection_handler_2()
  {
  }

  ~Self_intersection_handler_2(void)
  {
  }

class Intersection_observer_2 : public Arr_observer_2
{
public:

  Intersection_observer_2 (Arrangement_2& i_arr) :
    Arr_observer_2 (i_arr),
    m_arr(i_arr),
      m_sqr_offset((*smp_offset)*(*smp_offset))
  {
    
    LOG_DEBUG << "observer" << std::endl;
  }

  void sweep(const Approximate_polygon_2& polygon)
  {
    CGAL::insert(m_arr, polygon.curves_begin(), polygon.curves_end());
  }

  bool no_intersections()
  {
    LOG_DEBUG << "no of faces " << m_arr.number_of_faces() << std::endl;
    return (m_arr.number_of_faces() <= 2); // supposed to be one closed contour
  }

  void detect_intersections()
  {
    // For each original arc we need to find a list of sub-arcs
    // that got "inside" the contour, causing self-intersections
    // For this we:
    // 1) Find face(s) corresponding to the contour with wiped(merged)
    // self-intersections:
    // a) in case of the outer contour - the face(s) that is(are) touching
    // the unbounded face
    // b) in case of the inner contour - the unbounded face itself
    // 2) Go over all edges:
    // circular arcs that are "whole" as before the sweep are not intersecting
    // the concave sub-curves, or missing parts of convex curves, or missing
    // altogether curves are those that need to be handled
    // (replaced by corresponding preimage cover chain)
    // * concave = clockwise both for an outer contour and for a hole

    // TODO:
    // we use the global map of circular arcs to check for each arc
    // if it has changed
    typename Arrangement_2::Face_const_iterator            fit;
    typename Arrangement_2::Ccb_halfedge_const_circulator  curr;

    LOG_DEBUG << m_arr.number_of_faces() << " faces:" << std::endl;

    Face_vector outer_faces;
    int i = 0;
    for (fit = m_arr.faces_begin(); fit != m_arr.faces_end(); ++fit)
    {
      ++i;
      if (fit->is_unbounded())
      {
        LOG_DEBUG << "Unbounded Face no. " << i << " with " << fit->number_of_holes() <<" holes: ";
        typename Arrangement_2::Hole_const_iterator            hit;
        int                                 index = 1;
        for (hit = fit->holes_begin(); hit != fit->holes_end(); ++hit)
        {
          LOG_DEBUG << "\tHole #" << index << ": ";
          traverse_ccb_neighbours(*hit, outer_faces);
        }
      }
    }

    i = 0;
    typename Face_vector::iterator fpit;
    for (fpit = outer_faces.begin(); fpit != outer_faces.end(); ++fpit)
    {
      typename Arrangement_2::Face_const_iterator face = (*fpit);
      assert(!face->is_unbounded());

      LOG_DEBUG << "Outer Face #" << ++i << ": ";
      curr = face->outer_ccb();
      traverse_ccb(curr);
      
/* 
      typename Arrangement_2::Hole_const_iterator            hit;
      int                                 index = 1;
      for (hit = face->holes_begin(); hit != face->holes_end(); ++hit)
      {
        LOG_DEBUG << "\tHole #" << index << ": ";
        traverse_ccb(*hit);
      }
 */
    }
    
  }

  static bool is_convex_arc(const Point_2& start_point,
                            const Point_2& end_point,
                            const typename Polygon_traits_2::Point_2& arc_center)
  {
    // since it is at most half-arc it is convex if the points
    // source-target-center make a left turn
    // or, equivalently, if the points center-source-target make a left turn

    // the orientation predicate is_left_turn is true if the sign of the determinant
    // for 3 points a, b, c is positive, when the determinant is
    // bx - ax   by - ay
    // cx - ax   cy - ay
    //
    // if we take a=center all numbers in the determinant are one-root numbers

    CGAL::Comparison_result start_center = CGAL::compare(start_point.y(), arc_center.y());
    CGAL::Comparison_result end_center = CGAL::compare(end_point.y(), arc_center.y());
    CGAL::Comparison_result start_end = CGAL::compare(start_point.x(), end_point.x());
    assert(start_end != EQUAL); //  SMALLER, EQUAL, LARGER
    bool is_start_before_end = (start_end == SMALLER);
    
    assert((start_center == SMALLER && end_center != LARGER) || // can't be not x-monotone
          (start_center == LARGER && end_center != SMALLER) || // can't be not x-monotone
          (start_center == EQUAL && end_center != EQUAL) ); // can't be full half-arc
    bool is_below_center = (start_center == SMALLER) || (end_center == SMALLER);
    bool is_above_center = (start_center == LARGER) || (end_center == LARGER);

    return (is_start_before_end && is_below_center) ||
             (!is_start_before_end && is_above_center);
  }


  // find all faces that are neighbours of given hole ccb
  static void traverse_ccb_neighbours(typename Arrangement_2::Ccb_halfedge_const_circulator curr,
                                      Face_vector& o_faces)
  {
    typename Arrangement_2::Ccb_halfedge_const_circulator end = curr;
    do
    {
      typename Arrangement_2::Ccb_halfedge_const_circulator opp = curr->twin();
      assert(!opp->is_on_hole());
      typename Arrangement_2::Face_const_iterator pFace = opp->face();
      typename Face_vector::const_iterator fit = std::find(o_faces.begin(), o_faces.end(), pFace);
      if(fit == o_faces.end())
      {
        o_faces.push_back(pFace);
      }
      ++curr;
    } while (curr != end);
  }

  static void traverse_ccb(typename Arrangement_2::Ccb_halfedge_const_circulator  curr)
  {
    typename Arrangement_2::Ccb_halfedge_const_circulator end = curr;
        LOG_DEBUG << curr->source()->point();
        do
        {
          if(curr->curve().is_circular())
          {
            typename Polygon_traits_2::Point_2 arc_center = curr->curve().supporting_circle().center();
            //One_root_point_2 center(arc_center.x(), arc_center.y());
            //One_root_point_2 source(curr->source()->point().x(), curr->source()->point().y());
            //One_root_point_2 target(curr->target()->point().x(), curr->target()->point().y());
            // since curves are x-monotone there can't be an arc longer then half-circle,
            // so checking whether arc is above/below center and is start before end or vice versa
            // is equivalent to checking whether and arc is convex/concave
            //CGAL::Orientation orientation = curr->curve().orientation();
            //bool is_right = curr->curve().is_directed_right();
            bool is_convex = is_convex_arc(curr->source()->point(), curr->target()->point(), arc_center);
            LOG_DEBUG << " ~"<< (is_convex? "<ccw":"cw>>") <<"~> " << std::endl;

            if(!is_convex)
            {
              Input_rational sqr_offset = (*smp_offset) * (*smp_offset);
              Circle_2 circle(arc_center, sqr_offset);
              typename Circle_list_2::const_iterator cit = std::find(smp_circles->begin(), smp_circles->end(), circle);
              if(cit == smp_circles->end())
              {
                LOG_DEBUG << "add new concave arc" << std::endl;
                smp_circles->push_back(circle);
              }
              else
              {
                LOG_DEBUG << "existing concave arc" << std::endl;
              }

              typename std::pair<Ac_multimap_2_iterator, Ac_multimap_2_iterator> range =
                smp_ac_multimap->equal_range(arc_center);
              for(Ac_multimap_2_iterator it = range.first; it != range.second; ++it)
              {
                LOG_DEBUG << "\t c(" << (*it).first <<") s("<<(*it).second.source<<") - t("<<(*it).second.target<<")" << std::endl;
              }             
            }
            else
            {
              // add arc kgon end-points to "allowed" arc list
             typename std::pair<Ac_multimap_2_iterator, Ac_multimap_2_iterator> range =
                smp_ac_multimap->equal_range(arc_center);
              for(Ac_multimap_2_iterator it = range.first; it != range.second; ++it)
              {
                LOG_DEBUG << "\t c(" << (*it).first <<") s("<<(*it).second.source<<") - t("<<(*it).second.target<<")" << std::endl;
              }                           
            }
          }
          else
          {
            LOG_DEBUG << " --> ";
          }
          LOG_DEBUG << curr->target()->point();
          ++curr;
        } while (curr != end);
        LOG_DEBUG << std::endl;
  }

  virtual void before_split_face (Face_handle,
                                  Halfedge_handle e)
  {
    LOG_DEBUG << "-> The insertion of :  [ " << e->curve()
              << " ]  causes a face to split." << std::endl;

    if(e->curve().is_circular())
    {
      typename Polygon_traits_2::Point_2 arc_center = e->curve().supporting_circle().center();
//      Circle_2 circle(arc_center, m_sqr_offset);
//      smp_circles->push_back(circle);

      typename std::pair<Ac_multimap_2_iterator, Ac_multimap_2_iterator> range =
        smp_ac_multimap->equal_range(arc_center);
      for(Ac_multimap_2_iterator it = range.first; it != range.second; ++it)
      {
        LOG_DEBUG << "\t c(" << (*it).first <<") s("<<(*it).second.source<<") - t("<<(*it).second.target<<")" << std::endl;
      }
    }
  }

  // issued immediately after the existing face f1 has been split,
  // such that a portion of it now forms a new face f2. The flag
  // is_hole designates whether f2 forms a hole inside f1.

  virtual void after_split_face ( Face_handle f1,
                                  Face_handle f2,
                                  bool is_hole)
  {
    LOG_DEBUG << "-> Face was split, f2 is hole inside f1 :  [ " << f2->is_unbounded()
              << " ]  is hole: " << is_hole << " f1 is unbounded: " << f1->is_unbounded() << std::endl;
    
  }


/*
  virtual void before_merge_face (Face_handle,
                                  Face_handle,
                                  Halfedge_handle e)
  {
    LOG_DEBUG << "-> The removal of :  [ " << e->curve()
              << " ]  causes two faces to merge." << std::endl;
  }
*/
  // issued just before an edge e is split into two edges that should be 
  // associated with the x-monotone curves c1 and c2. The vertex v 
  // corresponds to the split point, and will be used to separate 
  // the two resulting edges.
  virtual void before_split_edge ( 	Halfedge_handle e,
    Vertex_handle v,
    X_monotone_curve_2 c1,
    X_monotone_curve_2 c2)
  {
    LOG_DEBUG << "-> Before the split of :  [ " << e->curve()
              << " ]  by [ " << v->point() << "]." << std::endl;
  }

  // issued immediately after an existing edge has been split into 
  // the two given edges e1 and e2. 
  virtual void 	after_split_edge ( Halfedge_handle e1, 
    Halfedge_handle e2)
  {
    LOG_DEBUG << "-> After the split to :  [ " << e1->curve()
              << " ]  and [ " << e2->curve() << "]." << std::endl;

  }
  private:

  Arrangement_2& m_arr;
  const Input_rational m_sqr_offset;
};

class Contour_handler_2
{
  public:
    Contour_handler_2(const Approximate_polygon_2& i_contour,
                      const Polygon_2& i_segment_contour,
                      Approximate_polygon_2& o_contour):
                      m_circle_contour(i_contour),
                      m_segment_contour(i_segment_contour),
                      m_good_contour(o_contour)
             
    {
    }

    void handle_self_intersections()
    {
      Arrangement_2 arr;
      Intersection_observer_2 observer(arr);
      observer.sweep(m_circle_contour);
      if(!observer.no_intersections())
      {
        LOG_DEBUG << "self-intersecting contour" << std::endl;
        observer.detect_intersections();  // in the meanwhile only prints detected data
        overwrite(m_segment_contour, m_good_contour);
      }
      else
      {
        LOG_DEBUG << "no self-intersection in contour" << std::endl;
      }

    }
    
  private:

     // TODO: fix "bad" arcs instead whole contour replacement
    void overwrite(const Polygon_2& p1, Approximate_polygon_2& p2)
    {
      LOG_DEBUG << "OVERWRITE" << std::endl;
      p2.clear();
      if(p1.is_empty()) return;

      typename Polygon_2::Edge_const_circulator ecirc = p1.edges_circulator();
      typename Polygon_2::Edge_const_circulator ecirc_end = ecirc;
      /*
      typename Approximate_polygon_2::Curve_const_iterator circ = p2.curves_begin();
      typename Approximate_polygon_2::Curve_const_iterator circ_end = p2.curves_end();

      Point_2 start_point = circ->source()->point();
      Point_2 orig_point = ecirc.source();
      while (start_point != orig_point)
      {
        ++ecirc;
        assert(ecirc != ecirc_end);
        orig_point = ecirc.source();
      }

      LOG_DEBUG << "App_point" << start_point <<std:endl;
      LOG_DEBUG << "Seg_point" << orig_point <<std:endl;
      */
      // TODO:
      // find the same point in two contours
      // iterate over both contours and whenever there is an arc that should be frozen
      // replace the arc with corresponding segments
      
      do
      {
        p2.push_back(X_monotone_curve_2 ((*ecirc).source(), (*ecirc).target()));
        ++ecirc;
      }while(ecirc != ecirc_end);
    }
    
    const Approximate_polygon_2& m_circle_contour;
    const Polygon_2& m_segment_contour;
    Approximate_polygon_2& m_good_contour;

};

  /// Detect self-intersections and replace self-intersecting contours with the minkowski sum results
   void handle_self_intersections(Approximate_polygon_list_2& kgon_sum_with_arcs,
                                 const Kgon_sum_polygon_2& kgon_sum,
                                 Approximate_offset_polygon_2& kgon_offset_polygon)
  {
    LOG_DEBUG << "handle_self_intersections()" << std::endl;
    assert(!kgon_sum_with_arcs.empty());
    Approximate_polygon_2 outer_contour = *kgon_sum_with_arcs.begin();
    kgon_sum_with_arcs.pop_front(); // pop outer contour

    LOG_DEBUG << "outer_contour"<< std::endl;
    Contour_handler_2 handler(outer_contour, kgon_sum.outer_boundary(), outer_contour);
    handler.handle_self_intersections();

    int i = 1;
    typename Approximate_polygon_list_2::iterator cit = kgon_sum_with_arcs.begin();
    for(typename Kgon_sum_polygon_2::Hole_const_iterator hit = kgon_sum.holes_begin();
        hit != kgon_sum.holes_end();
        ++hit, ++cit, ++i)
    {
      LOG_DEBUG << "hole " << i << std::endl;
      Contour_handler_2 handler(*cit, *hit, *cit);
      handler.handle_self_intersections();
    }
    
    kgon_offset_polygon = Approximate_offset_polygon_2(outer_contour,
      kgon_sum_with_arcs.begin(), kgon_sum_with_arcs.end());
  }


private:

};

template <typename CPolygon_2>
const typename Self_intersection_handler_2<CPolygon_2>::
Input_rational* Self_intersection_handler_2<CPolygon_2>::smp_offset;

template <typename CPolygon_2>
typename Self_intersection_handler_2<CPolygon_2>::
Circle_list_2* Self_intersection_handler_2<CPolygon_2>::smp_circles;

template <typename CPolygon_2>
typename Self_intersection_handler_2<CPolygon_2>::
Arc_center_multimap_2* Self_intersection_handler_2<CPolygon_2>::smp_ac_multimap;

} // namespace CGAL


#endif // Self_intersection_handler_2_h
