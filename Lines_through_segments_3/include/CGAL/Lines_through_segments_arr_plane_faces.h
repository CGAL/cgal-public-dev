// Copyright (c) 2010  Tel-Aviv University (Israel).
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
// $URL: $
// $Id: $
// 
//
// Author(s)     : Asaf Porat          <asafpor1@post.tau.ac.il>

#ifndef LINE_THROUGH_SEGMENTS_ARR_PLANE_FACES_H
#define LINE_THROUGH_SEGMENTS_ARR_PLANE_FACES_H

/*************************************************************
 * This class represent a vector of arrangements.
 * Each entry at the vector holds an arrangement that was created from 3 lines
 * on the same plane that are contained in a single plane.
 * After all of the lines were inserted, an overlay of all of the arrangements
 * is done.
 * If two faces overlap, the created face will represent a plane that passes 
 * through 4 lines.
 *
 * After this process is completed the overlaid arrangement is overlaid again
 * with arr_on_plane,
 * Each edge that falls inside the faces of the overlaid arrangement represents
 * a plane that passes through 4 lines.
 *************************************************************/

#include <CGAL/basic.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

namespace CGAL {

template <typename Arrangement_2,typename Observer>
class Lines_through_segments_arr_plane_faces
{
private:
  class Overlay_traits_private : 
    public CGAL::Arr_default_overlay_traits<Arrangement_2>
  {
    typedef typename Arrangement_2::Face_const_handle     Face_handle_A;
    typedef typename Arrangement_2::Face_const_handle     Face_handle_B;
    typedef typename Arrangement_2::Face_handle           Face_handle_R;
    typedef typename Arrangement_2::Halfedge_const_handle Halfedge_handle_A;
    typedef typename Arrangement_2::Halfedge_const_handle Halfedge_handle_B;
    typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle_R;
    typedef typename Arrangement_2::Dcel::const_iterator  segs_const_it;
   
  public:
    /* constructs the face f induced by the an overlap between the faces
       f1 and f2. */
    virtual void create_face (Face_handle_A f1,
                              Face_handle_B f2,
                              Face_handle_R f) const
    {
      segs_const_it f1_it,f2_it;
           
      for (f1_it = f1->segs_begin();
           f1_it != f1->segs_end();
           ++f1_it)
      {
        f->add_segment(*f1_it);
      }
           
      for (f2_it = f2->segs_begin();
           f2_it != f2->segs_end();
           ++f2_it)
      {
        f->add_segment(*f2_it);
      }
    }

    /*!
     * Create an edge e that matches the overlap between e1 and e2.
     */
    virtual void create_edge (Halfedge_handle_A e1,
                              Halfedge_handle_B e2,
                              Halfedge_handle_R e) const
    {
      /* The segments lists of e1 and e2 must be strictly distinct.
       */
      CGAL_assertion_code (bool are_distinct = true;
                           segs_const_it assert_it1;
                           for (assert_it1 = e1->segs_begin(); 
                                assert_it1 != e1->segs_end();
                                ++assert_it1)
                           {
                             segs_const_it assert_it2;
                             for (assert_it2 = e2->segs_begin(); 
                                  assert_it2 != e2->segs_end();
                                  ++assert_it2)
                             {
                               if (*assert_it2 == *assert_it1)
                                 are_distinct = false;
                             }
                           }
                           );
      CGAL_assertion(are_distinct == true);
           
      segs_const_it it;
      for (it = e1->segs_begin(); it != e1->segs_end(); ++it)
      {
        e->add_segment(*it);
        e->twin()->add_segment(*it);
      }
           
      for (it = e2->segs_begin(); it != e2->segs_end(); ++it)
      {
        e->add_segment(*it);
        e->twin()->add_segment(*it);
      }
    }
            
    /*!
     * Create an edge e that matches the edge e1, contained in the face f2.
     */
    virtual void create_edge (Halfedge_handle_A e1,
                              Face_handle_B f2,
                              Halfedge_handle_R e) const
    {
      segs_const_it it;
      for (it = e1->segs_begin(); it != e1->segs_end(); ++it)
      {
        e->add_segment(*it);
        e->twin()->add_segment(*it);
      }
    }
            
    /*!
     * Create an edge e that matches the edge e2, contained in the face f1.
     */
    virtual void create_edge (Face_handle_A f1,
                              Halfedge_handle_B e2,
                              Halfedge_handle_R e) const
    {
      segs_const_it it;
      for (it = e2->segs_begin(); it != e2->segs_end(); ++it)
      {
        e->add_segment(*it);
        e->twin()->add_segment(*it);
      }
    }
  };
      
private:
  std::vector<Arrangement_2*> arr_vector; 
  std::vector<Observer*> obs_arr_vector;
      
public:
  Lines_through_segments_arr_plane_faces()
  {
         
  }
      
  ~Lines_through_segments_arr_plane_faces()
  {
    typename std::vector<Observer*>::iterator oit;
    for (oit = obs_arr_vector.begin(); oit != obs_arr_vector.end(); ++oit)
    {
      delete (*oit);
    }
  }

  void add_element(Arrangement_2* arr,Observer* obs)
  {
    arr_vector.push_back(arr);
    obs_arr_vector.push_back(obs);
  }
      
  unsigned int size()
  {
    return arr_vector.size();
  }
      
  Arrangement_2* get_overlaid_arr(Arrangement_2* arr_on_plane)
  {
    Arrangement_2 *ov_arr = merge_arrangements(0,this->size()-1);
         
    Arrangement_2* res_arr = new Arrangement_2();
    Overlay_traits_private OT;
            
    CGAL::overlay(*ov_arr,
                  *arr_on_plane,
                  *res_arr,
                  OT);

    delete ov_arr;
    return res_arr;
  }
      
private:
  /*************************************************************
   * The following function merges all of the arrangements of arr_vector
   * to one arrangement.
   **************************************************************/
  Arrangement_2* merge_arrangements(int start,int end)
  {
    if (start == end)
    {
      return arr_vector[start];
    }
    /* Only 2 cells at the vector */
    else if (start + 1 == end)
    {
      Arrangement_2* res_arr = new Arrangement_2();
      Overlay_traits_private OT;
            
      CGAL::overlay(*arr_vector[start], *arr_vector[end], *res_arr, OT);
            
      delete arr_vector[start];
      delete arr_vector[end];
            
      return res_arr;
    }
    else
    {
      Arrangement_2* res_arr = new Arrangement_2();
      Arrangement_2 *arr1,*arr2;
      Overlay_traits_private OT;
            
      arr1 = merge_arrangements(start, start + (end - start)/2);
      arr2 = merge_arrangements(start + (end - start)/2 + 1, end);
            
      CGAL::overlay(*arr1, *arr2, *res_arr, OT);
            
      delete arr1;
      delete arr2;
            
      return res_arr;
    }
  }
};

} //namespace CGAL

#endif
