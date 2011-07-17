// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_on_surface_2/include/CGAL/Arr_new_landmarks_point_location.h $
// $Id: Arr_new_landmarks_point_location.h 56667 2010-06-09 07:37:13Z sloriot $
// 
//
// Author(s)     : Michael Hemmer <hemmer@googlemail.com> 
//                 Ophir Setter <ophir.setter@cs.tau.ac.il>

#ifndef CGAL_ARR_LANDMARKS_POINT_LOCATION_NO_CONSTRUCTION_H
#define CGAL_ARR_LANDMARKS_POINT_LOCATION_NO_CONSTRUCTION_H

/*! \file
 * Definition of the Arr_new_landmarks_point_location<Arrangement> template.
 */

//#define CGAL_DEBUG_LM

#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <set>

#include <boost/foreach.hpp>

namespace CGAL {

  /*! \class Arr_landmarks_point_location_no_construction
   * A class that answers point-location queries on an arrangement using the
   * landmarks algorithm, namely by locating the (approximately) nearest
   * landmark point to the qury point and walking from it toward the query
   * point.
   * This class-template has two parameters:
   * Arrangement corresponds to an arrangement-on-surface instantiation.
   * Generator is a class that generates the set of landmarks.
   */

  template <class Arrangement_, 
    class Generator_ = Arr_landmarks_vertices_generator<Arrangement_> >
    class Arr_landmarks_point_location_no_construction
    {
    public:
 
    typedef Arrangement_                                Arrangement_2;
    typedef typename Arrangement_2::Geometry_traits_2   Geometry_traits_2;
    typedef Generator_                                  Generator;
 
    protected:
 
    typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
    typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
    
    typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;
    typedef typename Arrangement_2::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
    typedef typename Arrangement_2::Outer_ccb_const_iterator
    Outer_ccb_const_iterator;
    typedef typename Arrangement_2::Inner_ccb_const_iterator
    Inner_ccb_const_iterator;
    typedef typename Arrangement_2::Isolated_vertex_const_iterator
    Isolated_vertex_const_iterator;
    
    typedef typename Arrangement_2::Point_2             Point_2;
    typedef typename Arrangement_2::X_monotone_curve_2  X_monotone_curve_2;


    typedef Arr_traits_basic_adaptor_2<Geometry_traits_2>    Traits_adaptor_2;

    // Data members:
    const Arrangement_2     *p_arr;     // The associated arrangement.
    const Traits_adaptor_2  *m_traits;  // Its associated traits object.
    Generator               *lm_gen;    // The associated landmark generator.
    bool                     own_gen;   // Indicates whether the generator
    // has been locally allocated.

    public:

    /*! Default constructor. */
    Arr_landmarks_point_location_no_construction (): 
    p_arr (NULL),
    m_traits (NULL),
    lm_gen(NULL),
    own_gen (false)
    {}

    /*! Constructor given an arrangement only. */
    Arr_landmarks_point_location_no_construction (const Arrangement_2& arr) :
    p_arr (&arr)
    {
      // Allocate the landmarks generator.
      m_traits = static_cast<const Traits_adaptor_2*> (p_arr->geometry_traits());
      lm_gen = new Generator(arr);
      own_gen = true;
    }

    /*! Constructor given an arrangement, and landmarks generator. */
    Arr_landmarks_point_location_no_construction (const Arrangement_2& arr, 
                                      Generator *gen) :
    p_arr (&arr),
    lm_gen (gen),
    own_gen (false)
    {
      m_traits = static_cast<const Traits_adaptor_2*> (p_arr->geometry_traits());
    }

    /*! Destructor. */
    ~Arr_landmarks_point_location_no_construction () 
    {
      if (own_gen) 
        delete lm_gen;
    }
   
    /*! Attach an arrangement object (and a generator, if supplied). */
    void attach (const Arrangement_2& arr, Generator *gen = NULL)
    {
      // Keep a pointer to the associated arrangement.
      p_arr = &arr;
      m_traits = static_cast<const Traits_adaptor_2*> (p_arr->geometry_traits());

      // Update the landmarks generator.
      if (gen != NULL)
        {
          // In case a generator is given, keep a pointer to it.
          CGAL_assertion (lm_gen == NULL);
          lm_gen = gen;
          own_gen = false;
        }
      else if (lm_gen != NULL)
        {
          // In case a generator exists internally, make sure it is attached to
          // the given arrangement.
          Arrangement_2 &non_const_arr = const_cast<Arrangement_2&>(*p_arr);
          lm_gen->attach(non_const_arr); 
        }
      else
        {
          // Allocate a new generator, attached to the given arrangement.
          lm_gen = new Generator(arr);
          own_gen = true;
        }
    }

    /*! Detach the instance from the arrangement object. */
    void detach () 
    {
      p_arr = NULL;
      m_traits = NULL;

      CGAL_assertion(lm_gen != NULL);
      if (lm_gen)
        lm_gen->detach();
    }
  
    /*!
     * Locate the arrangement feature containing the given point.
     * \param p The query point.
     * \return An object representing the arrangement feature containing the
     *         query point. This object is either a Face_const_handle or a
     *         Halfedge_const_handle or a Vertex_const_handle.
     */
    Object locate (const Point_2& p) const {
      
      // generator is empty, start blind walk
      if(lm_gen->is_empty())
        return _walk(p_arr->faces_begin(),p); 
      
      // Use the generator and to find the closest landmark to the query point.
      Object         lm_location_obj; 
      const Point_2& landmark_point = lm_gen->closest_landmark (p, lm_location_obj);
      
      // easy exit, querry is a landmark point
      if (landmark_point == p)
        return lm_location_obj;

  
      // Walk from the nearest_vertex to the point p, using walk algorithm,
      // and find the location of the query point p.

      // Locate the arrangement feature that contains the landmark.

      Vertex_const_handle     vh;
      Halfedge_const_handle   eh;
      Face_const_handle       fh;

      CGAL_assertion_msg (! lm_location_obj.is_empty(),
                          "lm_location_obj is empty.");
  
      if (assign(vh, lm_location_obj)) {
        
        CGAL_assertion_msg(vh->is_at_open_boundary() == false, "Generator is not allowed to generate a fictitious point");
         
        if(vh->is_isolated())
          return _walk(vh->face(), p);
        else
          return _walk(vh->incident_halfedges()->face(),p);
      }
      else if (assign(eh, lm_location_obj))
        return _walk(eh->face(),p);
      else if (assign(fh, lm_location_obj))
        return _walk(fh, p);
      
      // never reached 
      CGAL_assertion_msg (! lm_location_obj.is_empty(),
                          "lm_location_obj of an unknown type.");
      assert(false);
      return CGAL::Object(); 
    }
    

    protected:
    
    /** 
     * Returns a face to a specific direction of vh.
     * 
     * @param vh The vertex.
     * @param direction SMALLER - left, LARGER - right, EQUAL is not legal.
     * 
     * @return Returns a face to a specific direction of vh.
     */
    Face_const_handle _get_face_at_direction(Vertex_const_handle vh, CGAL::Comparison_result direction) const {
      CGAL_precondition (direction != EQUAL);

      // if we want to move right, we first have to look at edges to our right.
      Arr_halfedge_direction h_dir = (direction == LARGER) ? ARR_RIGHT_TO_LEFT : ARR_LEFT_TO_RIGHT;

      Halfedge_around_vertex_const_circulator hec = vh->incident_halfedges();
      do {
        CGAL_precondition(vh == hec->target());
        if (hec->direction() == h_dir) {
          if (!hec->face()->is_fictitious())
            return hec->face();
          else
            return hec->twin()->face();
        }
        hec++;
      } while(hec != vh->incident_halfedges());
      
      // There are no edges going to the proper direction, so we take the upper/lower edge of the other direction.
      Halfedge_around_vertex_const_circulator left_ext_edge = hec; hec++;
      while(hec != vh->incident_halfedges()) {
        if(m_traits->compare_y_at_x_left_2_object()
           (hec->curve(), left_ext_edge->curve(), vh->point()) == direction) {
          left_ext_edge = hec;
        }
        hec++;
      }
      
      if (!left_ext_edge->face()->is_fictitious())
        return left_ext_edge->face();
      else
        return left_ext_edge->twin()->face();
    }

    CGAL::Object _walk(Face_const_handle fh, const Point_2& p) const {

      CGAL_postcondition(!fh->is_fictitious());
    
      // if querry is in xrange start vertical walk 
      // otherwise continue walk in x direction by selecting 
      // the x-extremal vertex 
    
      std::vector<Ccb_halfedge_const_circulator> ccbs; 
      std::copy(fh->outer_ccbs_begin(), fh->outer_ccbs_end(), std::back_inserter(ccbs));
      std::copy(fh->inner_ccbs_begin(), fh->inner_ccbs_end(), std::back_inserter(ccbs));
    
      // If there are no CCBs, we still need to check isolated points.
      if(ccbs.size() == 0) {
        return _vertical_walk(fh, p); 
      }
    
      Vertex_const_handle vminh = (*ccbs.begin())->source();
      Vertex_const_handle vmaxh = (*ccbs.begin())->target();

      if(compare_x(vminh,vmaxh) == LARGER)
        std::swap(vminh, vmaxh);   
      CGAL_postcondition(compare_x(vminh, vmaxh) != LARGER);
    
    
      BOOST_FOREACH(Ccb_halfedge_const_circulator start, ccbs) {
        Ccb_halfedge_const_circulator hec = start;
        do {
          if (compare_x(hec->target(), vminh) == SMALLER)
            vminh = hec->target();
          if (compare_x(hec->target(), vmaxh) == LARGER)
            vmaxh = hec->target();
          hec++;
        } while(hec != start);
      }  
    
      // a valid  face must have a non singular x-range 
      CGAL_postcondition(compare_x(vminh, vmaxh) == SMALLER);
 
      // There are cases where there is no face to the right, and there is no boundary.
      // For example, in the regular (bounded) segment traits there is no outer ccb which is 
      // unbounded. If _get_face_at_direction returns to us the same face we are currently using
      // it means that there is no other face in that direction, so we go vertical.
      // 
      if (vmaxh->parameter_space_in_x() == ARR_INTERIOR && 
          compare_x(vmaxh, p) == SMALLER) {
        Face_const_handle next_f = _get_face_at_direction(vmaxh, LARGER);
        if (fh != next_f)
          return _walk(next_f, p);
      }
      
      if (vminh->parameter_space_in_x() == ARR_INTERIOR && 
          compare_x(p, vminh) == SMALLER) {
        Face_const_handle next_f = _get_face_at_direction(vminh, SMALLER);
        if (fh != next_f)
          return _walk(next_f, p);
      }
      
      return _vertical_walk(fh,p);
    }
  
    template<class OutputIterator>
    OutputIterator _get_edges_in_x_range(Ccb_halfedge_const_circulator begin, 
                                         const Point_2& p,OutputIterator oit) const {
      Ccb_halfedge_const_circulator hec = begin;
      do {
        if (hec->is_fictitious() == false &&
            m_traits->is_in_x_range_2_object() (hec->curve(), p)) {
          *oit++=hec;
        }
        hec++;
      } while(hec != begin);
    
      return oit;
    }

  
    CGAL::Object _vertical_walk(Face_const_handle fh, const Point_2& p) const {
    
      CGAL_postcondition(!fh->is_fictitious());
    
      std::vector<Ccb_halfedge_const_circulator> edges;
    
      for(Outer_ccb_const_iterator oit = fh->outer_ccbs_begin(); 
          oit != fh->outer_ccbs_end(); oit++){
        _get_edges_in_x_range(*oit, p, std::back_inserter(edges));    
      }

      for(Inner_ccb_const_iterator hit = fh->inner_ccbs_begin(); 
          hit != fh->inner_ccbs_end(); hit++){
        _get_edges_in_x_range(*hit, p, std::back_inserter(edges));
      }
      
      Ccb_halfedge_const_circulator ubounde; // the smallest edge larger  than p 
      Ccb_halfedge_const_circulator lbounde; // the largest  edge smaller than p 
    
      BOOST_FOREACH(Ccb_halfedge_const_circulator hec, edges) {
        CGAL_precondition(!hec->is_fictitious()); 
        CGAL_precondition(m_traits->is_in_x_range_2_object()(hec->curve(),p));

        CGAL::Comparison_result comp = CGAL::opposite(m_traits->compare_y_at_x_2_object()(p, hec->curve()));
        
        switch(comp) {

        case EQUAL: {
          if(hec->target()->parameter_space_in_x() == ARR_INTERIOR && 
             hec->target()->parameter_space_in_y() == ARR_INTERIOR) {
            
            if(m_traits->compare_xy_2_object() (p, hec->target()->point()) == EQUAL)
              return CGAL::make_object(hec->target());
          }
          
          if (hec->source()->parameter_space_in_x() == ARR_INTERIOR && 
              hec->source()->parameter_space_in_y() == ARR_INTERIOR) {
            if(m_traits->compare_xy_2_object()(p, hec->source()->point()) == EQUAL)
              return  CGAL::make_object(hec->source());
          }
          return CGAL::make_object(Halfedge_const_handle(hec));
        }

          // Curve is below the point
        case SMALLER: { 
          if( lbounde == Ccb_halfedge_const_circulator() ||
              m_traits->compare_y_position_2_object() (hec->curve(), lbounde->curve()) == LARGER)
            lbounde = hec;
          break;
        }
          
          // Cureve is above the point
        case LARGER: {
          if( ubounde == Ccb_halfedge_const_circulator()||
              m_traits->compare_y_position_2_object()(hec->curve(),ubounde->curve()) == SMALLER)
            ubounde = hec;  
        
          break;
        }   
        }
      }

      if (lbounde != Ccb_halfedge_const_circulator()) {
        CGAL_postcondition(m_traits->compare_y_at_x_2_object()(p, lbounde->curve()) == LARGER);
        
        // If we are not the lowest face, we move to the next face.
        if(lbounde->direction() == ARR_LEFT_TO_RIGHT || lbounde->twin()->face() == fh)
          return _check_isolated_vertices(fh,p);
        else{       
          return _vertical_walk(lbounde->twin()->face(),p); 
        }
      }
      if (ubounde != Ccb_halfedge_const_circulator()) {
        CGAL_postcondition(m_traits->compare_y_at_x_2_object()(p,ubounde->curve()) == SMALLER);

        // If we are not the top most face we move to the next face.
        if (ubounde->direction() == ARR_RIGHT_TO_LEFT || ubounde->twin()->face() == fh)
          return _check_isolated_vertices(fh,p);
        else     
          return _vertical_walk(ubounde->twin()->face(),p); 
      }
      
      // The only way to arrive here is if the face has no edges in the x-range of the point
      // (empty vertical area).
      CGAL_assertion(edges.size() == 0);
      return _check_isolated_vertices(fh, p);
    }
  
    CGAL::Object _check_isolated_vertices(Face_const_handle fh, const Point_2& p) const{    
      for(Isolated_vertex_const_iterator vit = fh->isolated_vertices_begin();
          vit != fh->isolated_vertices_end();vit++){
        if(m_traits->compare_xy_2_object()(p,vit->point()) == EQUAL){
          return CGAL::make_object(Vertex_const_handle(vit));
        }
      }
      return make_object(fh);
    };

  
    CGAL::Comparison_result compare_x(Point_2 p , Vertex_const_handle vh) const { 
      return  CGAL::opposite(compare_x(vh,p));
    }
  
    CGAL::Comparison_result compare_x(Vertex_const_handle vh, Point_2 p) const {    
      if(vh->parameter_space_in_x()==ARR_INTERIOR){
        if(vh->parameter_space_in_y() == ARR_INTERIOR){
          return m_traits->compare_x_2_object()(vh->point(),p);
        }else{
          Halfedge_around_vertex_const_circulator hec = vh->incident_halfedges();
          while(hec->is_fictitious()) hec++; 
          const X_monotone_curve_2& xcurve = hec->curve();
          Arr_curve_end ce = (hec->direction()==ARR_LEFT_TO_RIGHT)?ARR_MAX_END:ARR_MIN_END; 
          return CGAL::opposite((m_traits->compare_x_point_curve_end_2_object()(p, xcurve, ce)));
        }
      }
      if(vh->parameter_space_in_x() == ARR_LEFT_BOUNDARY) return SMALLER; 
      return LARGER; 
    }
  
    
    CGAL::Comparison_result compare_x(Vertex_const_handle v1, Vertex_const_handle v2) const {

      if(v1 == v2) return EQUAL;
    
      if( v1->parameter_space_in_x()==ARR_INTERIOR && 
          v2->parameter_space_in_x()==ARR_INTERIOR){
        if(v1->parameter_space_in_y() == ARR_INTERIOR){
          if(v2->parameter_space_in_y() == ARR_INTERIOR){
            return m_traits->compare_x_2_object()(v1->point(),v2->point()); 
          }else{
            Halfedge_around_vertex_const_circulator hec = v2->incident_halfedges();
            while(hec->is_fictitious()) hec++; 
            const X_monotone_curve_2& xcurve2 = hec->curve();
            Arr_curve_end ce2 = (hec->direction()==ARR_LEFT_TO_RIGHT)?ARR_MAX_END:ARR_MIN_END; 
            return m_traits->compare_x_point_curve_end_2_object() (v1->point(), xcurve2, ce2);
          }
        }else{
          if(v2->parameter_space_in_y() == ARR_INTERIOR){
            Halfedge_around_vertex_const_circulator hec = v1->incident_halfedges();
            while(hec->is_fictitious()) hec++; 
            const X_monotone_curve_2& xcurve1 = hec->curve();
            Arr_curve_end ce1 = (hec->direction()==ARR_LEFT_TO_RIGHT)?ARR_MAX_END:ARR_MIN_END; 
            return CGAL::opposite(m_traits->compare_x_point_curve_end_2_object()(v2->point(),xcurve1,ce1));
          }else{
            Halfedge_around_vertex_const_circulator hec;
            hec = v1->incident_halfedges();
            while(hec->is_fictitious()) hec++; 
            const X_monotone_curve_2& xcurve1 = hec->curve();
            Arr_curve_end ce1 = (hec->direction()==ARR_LEFT_TO_RIGHT)?ARR_MAX_END:ARR_MIN_END; 
            hec = v2->incident_halfedges();
            while(hec->is_fictitious()) hec++; 
            const X_monotone_curve_2& xcurve2 = hec->curve();
            Arr_curve_end ce2 = (hec->direction()==ARR_LEFT_TO_RIGHT)?ARR_MAX_END:ARR_MIN_END; 
            return m_traits->compare_x_curve_ends_2_object()(xcurve1,ce1,xcurve2,ce2);
          }
        }
      }
      if(v1->parameter_space_in_x() == v2->parameter_space_in_x()) return EQUAL; 
      if(v1->parameter_space_in_x() == ARR_LEFT_BOUNDARY ) return SMALLER;
      if(v2->parameter_space_in_x() == ARR_RIGHT_BOUNDARY) return SMALLER;
      return LARGER;   
    }
    };
} //namespace CGAL


#endif // CGAL_ARR_NEW_LANDMARKS_POINT_LOCATION_H
