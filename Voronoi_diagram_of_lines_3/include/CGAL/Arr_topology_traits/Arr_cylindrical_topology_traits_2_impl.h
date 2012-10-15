// Copyright (c) 2006, Tel-Aviv University (Israel).
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
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Arrangement_on_surface_2/include/CGAL/Arr_topology_traits/Arr_cylindrical_topology_traits_2_impl.h $
// $Id: Arr_cylindrical_topology_traits_2_impl.h 50366 2009-07-05 12:56:48Z efif $
// 
// Author(s)     : Michael Hemmer       <hemmer@mpi-inf.mpg.de>

#ifndef CGAL_ARR_CYLINDRICAL_TOPOLOGY_TRAITS_2_IMPL_H
#define CGAL_ARR_CYLINDRICAL_TOPOLOGY_TRAITS_2_IMPL_H

/*! \file
 * Member-function definitions for the
 * Arr_cylindrical_topology_traits_2<GeomTraits> class.
 */

namespace CGAL {

/*! \brief no default constructor */
/*! \brief constructs with a geometry-traits class */
template <class GeomTraits, class Dcel>
Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::
Arr_cylindrical_topology_traits_2(const Geometry_traits_2 * traits) :
  m_cylindrical_face(NULL),
  m_l_fictitious_face(NULL),
  m_r_fictitious_face(NULL),
  m_l_fictitious_vertex(NULL),
  m_r_fictitious_vertex(NULL),
  n_inf_verts(0)
{
  m_boundary_vertices = Vertex_map(Vertex_key_comparer(m_traits));
  m_traits = static_cast<const Traits_adaptor_2*>(traits);
}

/*! \brief assigns the contents of another topology-traits class */
template <class GeomTraits, class Dcel>
void Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::
assign(const Self & other)
{
  std::cout << "Assign Arr_cylindrical_topology_traits_2 " << std::endl;
  // Clear the current DCEL and duplicate the other DCEL.
  m_dcel.delete_all();
  m_dcel.assign(other.m_dcel);
  m_traits = other.m_traits;
  // Update this after change of 
  dcel_updated();
  this->check();
  return;
}

template <class GeomTraits, class Dcel>
void Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::check() const
{
  // std::cout << "Arr_cylindrical_topology_traits_2: check " << std::endl;
  CGAL_assertion(this->m_cylindrical_face != NULL);
  CGAL_assertion(this->m_l_fictitious_face != NULL); 
  CGAL_assertion(this->m_r_fictitious_face != NULL); 
  CGAL_assertion(this->m_l_fictitious_vertex != NULL); 
  CGAL_assertion(this->m_r_fictitious_vertex != NULL);
  CGAL_assertion(this->m_traits != NULL);
  return; 
}


/*! \brief updates *this after change of dcel */
template <class GeomTraits, class Dcel>
void Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::dcel_updated()
{
  std::cout << "Arr_cylindrical_topology_traits_2: dcel_updated"<< std::endl;
  assert(false);
  typename Dcel::Vertex_iterator       vit;
  Arr_parameter_space                  psx, psy;
  
  // Go over the DCEL vertices and locate the special vertices and 
  // the vertices on the identification 
  m_r_fictitious_vertex = NULL;
  m_l_fictitious_vertex = NULL;
  n_inf_verts = 0 ;
  m_boundary_vertices.clear();

  for (vit = this->m_dcel.vertices_begin();
       vit != this->m_dcel.vertices_end(); ++vit)
    {
      psx = vit->parameter_space_in_x();
      psy = vit->parameter_space_in_x();
    
      // count number of infinit vertices 
      if(psy == ARR_LEFT_BOUNDARY || psy == ARR_RIGHT_BOUNDARY){
        n_inf_verts++;
      }

      if(psx == ARR_TOP_BOUNDARY || psx == ARR_BOTTOM_BOUNDARY){
        // get vertices on identification 
        if(psy == ARR_INTERIOR){
          const Point_2& key = vit->point();
          m_boundary_vertices.insert(Vertex_value(key, &(*vit)));
        }
        // get the LEFT fictitious vertex 
        if(psy == ARR_LEFT_BOUNDARY){
          CGAL_assertion(m_l_fictitious_vertex == NULL);
          m_l_fictitious_vertex = &(*vit);
        }
        // get the RIGHT fictitious vertex 
        if(psy == ARR_RIGHT_BOUNDARY){
          CGAL_assertion(m_r_fictitious_vertex == NULL);
          m_r_fictitious_vertex = &(*vit);
        }
      } 
    }
  // Since fictitious they have to exist 
  CGAL_assertion(m_r_fictitious_vertex != NULL);
  CGAL_assertion(m_l_fictitious_vertex != NULL);


  // Go over the DCEL faces and locate the special faces 
  m_cylindrical_face = NULL;
  typename Dcel::Face_iterator         fit;  
  for (fit = this->m_dcel.faces_begin(); fit != this->m_dcel.faces_end(); ++fit)
    {
      if (fit->number_of_outer_ccbs() == 0)
        {
          CGAL_assertion(m_cylindrical_face == NULL);
          m_cylindrical_face = &(*fit);
        }
    
    }
  CGAL_assertion(m_cylindrical_face != NULL);
  
  assert(false); // not sure what todo with these guys yet.. 
  m_l_fictitious_face = NULL;
  m_r_fictitious_face = NULL;
  CGAL_assertion(m_l_fictitious_face != NULL);
  CGAL_assertion(m_r_fictitious_face != NULL); 
  this->check();
  return;
}

/*! \brief initializes an empty DCEL structure. */
template <class GeomTraits, class Dcel>
void Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::init_dcel()
{
  std::cout << "Arr_cylindrical_topology_traits_2: init_dcel()" << std::endl;
  // Clear the current DCEL.
  m_dcel.delete_all();
  m_boundary_vertices.clear();
  m_cylindrical_face = NULL;
  m_l_fictitious_face = NULL; 
  m_r_fictitious_face = NULL; 
  m_l_fictitious_vertex = NULL; 
  m_r_fictitious_vertex = NULL;
  n_inf_verts = 0; 



  // Create the face.
  m_cylindrical_face = this->m_dcel.new_face();
  m_cylindrical_face->set_unbounded(true);
  m_cylindrical_face->set_fictitious(false);
 
  // init LEFT part 
  // init fictitious face 
  m_l_fictitious_face = this->m_dcel.new_face();
  m_l_fictitious_face->set_unbounded(false);
  m_l_fictitious_face->set_fictitious(true);
  // init fictitious vertex 
  m_l_fictitious_vertex = this->m_dcel.new_vertex();
  m_l_fictitious_vertex->set_boundary(ARR_LEFT_BOUNDARY,ARR_TOP_BOUNDARY);  
  // get Halfedges for loop 
  Halfedge  *he_l_i = this->m_dcel.new_edge(); // will form inner ccb of f
  Halfedge  *he_l_o = he_l_i->opposite(); // will form outer ccb of left_ff 
  // 'init' curve
  he_l_i->set_curve (NULL);
  he_l_o->set_curve (NULL);
  // set next edges (self)
  he_l_i->set_next (he_l_i);  
  he_l_o->set_next (he_l_o);
  // set target to fictitious vertex 
  he_l_i->set_vertex (m_l_fictitious_vertex);
  he_l_o->set_vertex (m_l_fictitious_vertex);
  // get new outer/inner ccbs 
  Outer_ccb   *oc_l = this->m_dcel.new_outer_ccb();
  Inner_ccb   *ic_l = this->m_dcel.new_inner_ccb();
  // init ccbs 
  he_l_i->set_inner_ccb (ic_l);       
  he_l_o->set_outer_ccb (oc_l);
  
  // Assign the incident halfedges of the two fictitious vertices.
  m_l_fictitious_vertex->set_halfedge (he_l_i);
  m_l_fictitious_vertex->set_halfedge (he_l_o);
  // Set the direction of the halfedges:
  he_l_i->set_direction (ARR_RIGHT_TO_LEFT);
  he_l_o->set_direction (ARR_LEFT_TO_RIGHT);
  // add inner ccb to m_cylindrical_face
  m_cylindrical_face ->add_inner_ccb (ic_l, he_l_i); 
  ic_l->set_face (m_cylindrical_face);
  // add outer ccb the m_l_fictitious_face.
  m_l_fictitious_face->add_outer_ccb (oc_l, he_l_o);
  oc_l->set_face (m_l_fictitious_face);

  CGAL_assertion(*(m_cylindrical_face->inner_ccbs_begin())  == he_l_i);
  CGAL_assertion(*(m_l_fictitious_face->outer_ccbs_begin()) == he_l_o);

  // init RIGHT part 
  // init fictitious face 
  m_r_fictitious_face = this->m_dcel.new_face();
  m_r_fictitious_face->set_unbounded(false);
  m_r_fictitious_face->set_fictitious(true);
  // init fictitious vertex 
  m_r_fictitious_vertex = this->m_dcel.new_vertex();
  m_r_fictitious_vertex->set_boundary(ARR_RIGHT_BOUNDARY,ARR_TOP_BOUNDARY);  
  // get Halfedges for loop 
  Halfedge  *he_r_i = this->m_dcel.new_edge(); // will form inner ccb of f
  Halfedge  *he_r_o = he_r_i->opposite(); // will form outer ccb of left_ff 
  // 'init' curve
  he_r_i->set_curve (NULL);
  he_r_o->set_curve (NULL);
  // set next edges (self)
  he_r_i->set_next (he_r_i);  
  he_r_o->set_next (he_r_o);
  // set target to fictitious vertex 
  he_r_i->set_vertex (m_r_fictitious_vertex);
  he_r_o->set_vertex (m_r_fictitious_vertex);
  // get new outer/inner ccbs 
  Outer_ccb   *oc_r = this->m_dcel.new_outer_ccb();
  Inner_ccb   *ic_r = this->m_dcel.new_inner_ccb();
  // init ccbs 
  he_r_i->set_inner_ccb (ic_r);       
  he_r_o->set_outer_ccb (oc_r);
  // Assign the incident halfedges of the two fictitious vertices.
  m_r_fictitious_vertex->set_halfedge (he_r_i);
  m_r_fictitious_vertex->set_halfedge (he_r_o);
  // Set the direction of the halfedges:
  he_r_i->set_direction (ARR_LEFT_TO_RIGHT);
  he_r_o->set_direction (ARR_RIGHT_TO_LEFT);
  // add inner ccb to m_cylindrical_face
  m_cylindrical_face ->add_inner_ccb (ic_r, he_r_i);
  ic_r->set_face (m_cylindrical_face);
  // add outer ccb the m_r_fictitious_face.
  m_r_fictitious_face->add_outer_ccb (oc_r, he_r_o);
  oc_r->set_face (m_r_fictitious_face);

  n_inf_verts = 2; 
  this->check();
  CGAL_assertion(this->is_empty_dcel());
  return;
}

/*! \brief determines whether a point lies in the interior of a given face. */
template <class GeomTraits, class Dcel>
bool Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::
is_in_face(const Face * f, const Point_2 & p, const Vertex * v) const
{ 
  this->check();
  // The principal idea is to shoot a ray upward an count the  number 
  // of intersections with the outer_ccb. However, since Top/Bottom 
  // are identified this is not enough. Thus we stop at the ARR_TOP_BOUNDARY
  // and shoot a ray to the left until we reach ARR_LEFT_BOUNDARY
  // If the total number of intersections is odd the point is inside the face. 
  
  // This is implemented in the following way:
  // The ray intersects the edge if the point is in x range and below the edge. 
  // The ray along the identification intersects an edge at each and endpoint 
  // that is at ARR_TOP_BOUNDARY and to the left of the point. 
  // Note that do not count vertices at ARR_BOTTOM_BOUNDARY since we actually 
  // want to count the total number of intersections with the outer_ccb
  // Note the code below must also take care about special cases such as vertical 
  // edges.
  
  // the followoing code is a modifaction of the code in 
  // Arr_topology_traits/Arr_planar_topology_traits_base_2.h 
  // That is, we also maintain a counter for every end point 
  // at the top and to the left of p. 

  CGAL_precondition (v == NULL || ! v->has_null_point());
  CGAL_precondition (v == NULL || 
                     m_traits->equal_2_object()(p, v->point()));

  // In case the face is unbounded and has no outer ccbs, this is the single
  // unbounded face of an arrangement of bounded curves. This face obviously
  // contains any point in its interior.
  if (f->number_of_outer_ccbs() == 0){
    CGAL_assertion(f->is_unbounded());
    return (true);
  }

  assert(false); // TODO 

  // Keep a counter of the number of x-monotone curves that intersect an upward
  // vertical emanating from p (except for some degenerate cases that are
  // explained below).
  unsigned int       n_ray_intersections = 0;

  // Get a halfedge along the outer CCB of the given face, go over all curves
  // of the boundary component, and count those which are above p.
  // We begin by comparing p to the source vertex of the first halfedge.
  // Note that if p coincides with this vertex, p is obviously not in the
  // interior of the component.
  const Halfedge    *first = *(f->outer_ccbs_begin());
  const Halfedge    *curr = first;
  Comparison_result  res_source;
  Comparison_result  res_target;
  Comparison_result  res_y_at_x;

  if (curr->opposite()->vertex() == v) return (false);
  
  res_source = 
    _compare_xy (p, curr->curve(), // curve end of the SOURCE !!  
        (curr->direction() == ARR_LEFT_TO_RIGHT)?ARR_MIN_END:ARR_MAX_END);

  if(res_source == EQUAL) return (false);
  
  do {
    if (curr->vertex() == v) return (false);

    res_target =  
      _compare_xy (p, curr->curve(), // curve end of the TARGET !!  
          (curr->direction() == ARR_LEFT_TO_RIGHT)?ARR_MAX_END:ARR_MIN_END);
    
    if(res_target == EQUAL) return (false);
    
//     // In case the current halfedge belongs to an "antenna", namely its
//     // incident face is the same as its twin's, we can simply skip it
//     // (in order not to count it twice).
//     // This is just an optimization 
//     if (! curr->opposite()->is_on_inner_ccb() &&
//         curr->outer_ccb()->face() == curr->opposite()->outer_ccb()->face())
//     {
//       curr = curr->next();
//       res_source = res_target;
//       continue;
//     }

    // Check that if we shoot a "tilted" vertical ray from p upward
    // (by "tilted" we mean the angle it forms with the x-axis is
    //  PI/2 + epsilon, where epsilon is arbitrarily small), then we hit
    // the x-monotone curve associated with curr once.
    if (res_source != res_target)
    {
      res_y_at_x = compare_y_at_x (p, curr);

      if (res_y_at_x == SMALLER)
      {
        n_ray_intersections++;
      }        
      else if (res_y_at_x == EQUAL)
      {
        // In this case p lies on the current edge, so it is obviously not
        // contained in the interior of the component.
        return (false);
      }
    }

    // Proceed to the next halfedge along the component boundary.
    // Note that the source vertex of this halfedge is the current target.
    curr = curr->next();
    res_source = res_target;

  } while (curr != first);

  // The query point lies inside the connected components if and only if the
  // ray we shoot from it intersects the boundary an odd number of time.
  // return ((n_ray_intersections % 2) != 0);
  assert(false);

  // The following counts the number of curve ends at ARR_TOP_BOUNDARY 
  // that are to the left of p. 
  // 
  unsigned int n_top_left_curve_ends = 0;
  {
    const Halfedge    *first = *(f->outer_ccbs_begin());
    const Halfedge    *curr = first;
    do{
      if( (this->m_traits->parameter_space_in_x_2_object()
              (curr->curve(),ARR_MIN_END) == ARR_TOP_BOUNDARY )
          &&
          (this->m_traits->compare_x_near_boundary_2_object()
              (p,curr->curve(),ARR_MIN_END) == LARGER)){
        n_top_left_curve_ends++;
      }
      if( (this->m_traits->parameter_space_in_x_2_object()
              (curr->curve(),ARR_MAX_END) == ARR_TOP_BOUNDARY )
          &&
          (this->m_traits->compare_x_near_boundary_2_object()
              (p,curr->curve(),ARR_MAX_END) == LARGER)){
        n_top_left_curve_ends++;
      }
      curr = curr ->next();  
    }while(curr != first);
  }

  std::cout << "n_top_left_curve_ends" << n_top_left_curve_ends << std::endl;
  
  
  std::cout <<" is_in_face: "  << std::endl;
  std::cout << p << std::endl;
  assert(false); // use is_in_face of ARR_unb_top_traits_base_2.. ?
}

/*! \brief compares the relative y-position of a point and a halfedge */
template <class GeomTraits, class Dcel>
Comparison_result
Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::
compare_y_at_x(const Point_2 & p, const Halfedge * he) const
{
  this->check();
  return m_traits->compare_y_at_x_2_object()(p, he->curve());
}

/*! \brief determine whether a vertex is associated with a curve end */
template <class GeomTraits, class Dcel>
bool Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::
are_equal(const Vertex * v,
    const X_monotone_curve_2 & xc, Arr_curve_end ind,
    Arr_parameter_space ps_x, Arr_parameter_space ps_y) const
{
  this->check();
  // In case the given boundary conditions do not match those of the given
  // vertex, v cannot represent the curve end.
  if (ps_x != v->parameter_space_in_x())  return false;   
  if(ps_y == ARR_INTERIOR){
    if(v->parameter_space_in_y() != ARR_INTERIOR) return false;
  }else{
    if(v->parameter_space_in_y() == ARR_INTERIOR) return false;
  }  

  // if v is on identification line it must have a point which we can compare to. 
  if(ps_y == ARR_TOP_BOUNDARY || ps_y == ARR_BOTTOM_BOUNDARY){
    CGAL_assertion(!v->has_null_point());
    return 
      EQUAL == m_traits->compare_xy_2_object()(v->point(),_get_end_point(xc,ind));
  }else{
    CGAL_assertion(ps_y == ARR_INTERIOR);
    CGAL_assertion(ps_x == ARR_LEFT_BOUNDARY || ps_x == ARR_RIGHT_BOUNDARY);
    // CGAL_assertion(v->degree()==3); // there must be exactly one incident curve 

    Arr_curve_end            v_ind;
    const X_monotone_curve_2& v_cv = _curve (v, v_ind);
    CGAL_assertion (v_ind == ind);
    return (EQUAL == 
        m_traits->compare_y_near_boundary_2_object() (xc, v_cv, v_ind));
  }
}

/*! \brief receives a notification on the creation of a new boundary vertex */
template <class GeomTraits, class Dcel>
void
Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::
notify_on_boundary_vertex_creation(Vertex * v,
    const X_monotone_curve_2 & xc,
    Arr_curve_end ind,
    Arr_parameter_space ps_x,
    Arr_parameter_space ps_y)
{
  this->check();
  CGAL_assertion(ps_x == ARR_INTERIOR || ps_y == ARR_INTERIOR);

  if(ps_y == ARR_TOP_BOUNDARY || ps_y == ARR_BOTTOM_BOUNDARY){
    const Point_2 & key = (ind == ARR_MIN_END)     ?
      m_traits->construct_min_vertex_2_object()(xc):
      m_traits->construct_max_vertex_2_object()(xc);
    m_boundary_vertices.insert(Vertex_value(key, v));
  }

  if(ps_x == ARR_LEFT_BOUNDARY){
    assert(false); // do I have/want to do something here ?
  } 
  if(ps_x == ARR_RIGHT_BOUNDARY){
    assert(false); // do I have/want to do something here ?
  }
}






/*! \brief given a curve end with boundary conditions and a face that contains
 * the interior of the curve, find a place for a boundary vertex that will
 * represent the curve end along the face boundary
 *
 * For ARR_LEFT_BOUNDARY / ARR_RIGHT_BOUNDARY we have to return the 
 * fictitious edge that contains the curve end (which will be split)
 * For ARR_TOP_BOUNDARY / ARR_BOTTOM_BOUNDARY we have to return a 
 * vertex on the indetification line, if it already exists
 * Otherwise we HAVE TO return the empty object, this is the main difference to the 
 * function locate_curve_end, which returns the face that contains the curve end. 
 */

template <class GeomTraits, class Dcel>
CGAL::Object
Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::
place_boundary_vertex(Face * CGAL_assertion_code(f),
    const X_monotone_curve_2 & xc, Arr_curve_end ind,
    Arr_parameter_space ps_x, Arr_parameter_space ps_y)
{ 
  this->check();
#if 1
  std::cout << " place_boundary_vertex " << std::endl;
  std::cout << xc << std::endl; 
  std::cout << ind << std::endl; 
  std::cout << ps_x << std::endl;
  std::cout << ps_y << std::endl;
#endif 

  if(ps_x == ARR_LEFT_BOUNDARY  || ps_x == ARR_RIGHT_BOUNDARY){
    Halfedge   *first = (ps_x == ARR_LEFT_BOUNDARY)? 
      *m_l_fictitious_face->outer_ccbs_begin(): // take the left face 
      *m_r_fictitious_face->outer_ccbs_begin(); // take the right face 
    Halfedge   *curr = first;
    bool        is_source, is_target;
    do {
      // Note we consider only fictitious halfedges and check whether they
      // contain the relevant curve end.
      if (_is_on_fictitious_edge (xc, ind, ps_x, ps_y, curr, is_source, is_target)){
        CGAL_assertion(! is_source);
        CGAL_assertion(! is_target);
        CGAL_assertion(f == _get_face(curr->opposite()));
        return (CGAL::make_object (curr->opposite()));
      }
      // Move to the next halfegde along the CCB.
      curr = curr->next();
    } while (curr != first);

    // If we reached here, we did not find a suitable halfegde, which should
    // never happen.
    CGAL_error();
    return CGAL::Object();
  }
  
  CGAL_assertion(ps_y == ARR_TOP_BOUNDARY || ps_y == ARR_BOTTOM_BOUNDARY);
  typename Vertex_map::iterator  it;
  Vertex                         *v = NULL;
  // return a feature on the identification 
  // Check if the given curve end is incident to a vertex on the line of
  // discontinuity. If so, return this vertex. Otherwise, locate the first
  // vertex above it.
  
  const Point_2 & key = _get_end_point(xc,ind);
  it = m_boundary_vertices.find(key);
  if (it != m_boundary_vertices.end()) {
    v = it->second;
    return CGAL::make_object(v);
  }
  // there is no vertex that equals the curve end;
  return CGAL::Object();
}

/*! \brief locates a DCEL feature that contains a given curve end. */
template <class GeomTraits, class Dcel>
CGAL::Object Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::
locate_curve_end(
    const X_monotone_curve_2 & xc, 
    Arr_curve_end ind,
    Arr_parameter_space ps_x, 
    Arr_parameter_space ps_y)
{ 
  this->check();
  std::cout << " located curve end: " << std::endl; 
  
  // std::cout << ind << std::endl; 
  std::cout << ps_x << std::endl;
  std::cout << ps_y << std::endl;
  std::cout << xc << std::endl; 

  if(ps_x == ARR_LEFT_BOUNDARY  || ps_x == ARR_RIGHT_BOUNDARY){
    Halfedge   *first = (ps_x == ARR_LEFT_BOUNDARY)? 
      *m_l_fictitious_face->outer_ccbs_begin(): // take the left face 
      *m_r_fictitious_face->outer_ccbs_begin(); // take the right face 
    Halfedge   *curr = first;
    bool        is_source, is_target;
    do{
      if (_is_on_fictitious_edge(xc,ind,ps_x,ps_y,curr,is_source,is_target)){
        if (is_source){
          std::cout <<" is_source "<< std::endl;
          // cv's end coincides with the source vertex of the current
          // fictitious halfedge. This means that cv overlaps the curve that
          // is associated with the only non-fictitious halfedge incident to
          // this vertex. We therefore return a pointer to this halfedge.
          Halfedge     *he = curr->opposite()->next();
          
          CGAL_assertion (! he->has_null_curve());
          return (CGAL::make_object (he));
        }
        if (is_target){
          std::cout <<" is_target "<< std::endl;
          // cv's end coincides with the target vertex of the current
          // fictitious halfedge. This means that cv overlaps the curve that
          // is associated with the only non-fictitious halfedge incident to
          // this vertex. We therefore return a pointer to this halfedge.
          Halfedge     *he = curr->opposite()->prev();

          CGAL_assertion (! he->has_null_curve());
          return (CGAL::make_object (he));
        }
        std::cout << " is_on fictitious edge " << std::endl;
        // The current ficitious edge contains cv's end in its interior.
        // Note we use curr's twin, whose incident face is a valid
        // face (whereas the incident face of curr is the fictitious face).
        // Note that we rely on the fact that this always belongs to an inner ccb
        Face      *uf = _get_face(curr->opposite());
        CGAL_assertion (!uf->is_fictitious());
        CGAL_assertion ( uf->is_unbounded());
        return (CGAL::make_object(uf));
      }
      curr = curr->next();
    } while (curr != first);
    // We should never reach here.
    CGAL_error();
    return Object();
  } 

  if(ps_y == ARR_TOP_BOUNDARY || ps_y == ARR_BOTTOM_BOUNDARY){

    typename Vertex_map::iterator  it;
    Vertex                         *v = NULL;

    // return a feature on the identification 
    // Check if the given curve end is incident to a vertex on the line of
    // discontinuity. If so, return this vertex. Otherwise, locate the first
    // vertex above it.

    const Point_2 & key = _get_end_point(xc,ind);

    it = m_boundary_vertices.find(key);
    if (it != m_boundary_vertices.end()) {
      v = it->second;
      return CGAL::make_object(v);
    }
    // there is no vertex that equals the curve end;
    CGAL_assertion(it == m_boundary_vertices.end());
    // find the vertex right below the curve end; 
    it = m_boundary_vertices.lower_bound(key);
    
    if(it == m_boundary_vertices.end()){
      std::cout << "m_l_fictitious_vertex " << std::endl; 
      v =  this->m_l_fictitious_vertex;
    }else{ 
      std::cout << "  it->second " << std::endl; 
      v  = it->second; 
    }
    CGAL_assertion(v != NULL);
    return make_object(_face_right_of_vertex_on_top_identification(v)); 
  }
  // never get here! 
  assert(false);
  return CGAL::Object();
  // return the face to the right of *it
}



/*! \brief determines whether a given boundary vertex is redundant */
template <class GeomTraits, class Dcel>
bool
Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::
is_redundant(const Vertex * v) const
{
  return (v->halfedge() == NULL);
}

/* \brief erases a given redundant vertex */
template <class GeomTraits, class Dcel>
typename Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::Halfedge *
Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::
erase_redundant_vertex(Vertex * v)
{
  this->check();
  assert(false);
  return NULL; 
/*
  const Arr_parameter_space ps_y = v->parameter_space_in_y();
  if (ps_y == ARR_BOTTOM_BOUNDARY) {
  m_south_pole = NULL;
  return NULL;
  }
  if (ps_y == ARR_TOP_BOUNDARY) {
  m_north_pole = NULL;
  return NULL;
  }
  CGAL_assertion_code(Arr_parameter_space ps_x = v->parameter_space_in_x());
  CGAL_assertion(ps_x != ARR_INTERIOR);
  m_boundary_vertices.erase(v->point());
  return NULL;
*/
}

/*! \brief locate the predecessor halfedge for the given curve around a given
 * vertex with boundary conditions. */
template <class GeomTraits, class Dcel>
typename Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::Halfedge * 
Arr_cylindrical_topology_traits_2<GeomTraits,Dcel>::
locate_around_boundary_vertex(Vertex * v,
    const X_monotone_curve_2 & xc,
    Arr_curve_end ind,
    Arr_parameter_space ps_x,
    Arr_parameter_space ps_y) const {

  this->check();
  
  CGAL_assertion(ps_y != ARR_INTERIOR);
  

  // If the vertex is isolated, there is no predecssor halfedge.
  if (v->is_isolated()){
    std::cout << "locate_around_boundary_vertex: "
              << "v->is_isolated()" << std::endl;
    return NULL;
  }

  // Get the first incident halfedge around v and the next halfedge.
  Halfedge * first = v->halfedge();
  Halfedge * curr = first;
  CGAL_assertion(curr != NULL);
  
  // If is only one halfedge incident to v, 
  // return this halfedge as xc's predecessor:
  
  if (curr == curr->next()->opposite()) {
    std::cout << "locate_around_boundary_vertex: "
              << "curr == curr->next()->opposite()" << std::endl;
    return curr;
  }
  

  assert(false);
  
  Halfedge* upper_bound = NULL;
  Halfedge* lower_bound = NULL;
  do{
    Arr_curve_end       ind_c  = _get_curve_end(curr);
    Arr_parameter_space ps_x_c = _parameter_space_in_x(curr);
    Arr_parameter_space ps_y_c = _parameter_space_in_y(curr);
//    if(_less_around_vertex_on_identification(curr))


    curr = curr->next()->opposite();
  }while(curr != first);
    
  assert(false);
  Halfedge * next = first;
  // Otherwise, we traverse the halfedges around v until we find the pair
  // of adjacent halfedges between which we should insert xc.
  typename Traits_adaptor_2::Is_between_cw_2 is_between_cw =
    m_traits->is_between_cw_2_object();
  bool eq_curr, eq_next;

  while (!is_between_cw(xc, (ind == ARR_MIN_END), curr->curve(), 
          (curr->direction() == ARR_RIGHT_TO_LEFT),
          next->curve(), 
          (next->direction() == ARR_RIGHT_TO_LEFT), v->point(),
          eq_curr, eq_next))
    {
// The curve must not be equal to one of the curves already incident to v.
      CGAL_assertion(!eq_curr && !eq_next);

// Move to the next pair of incident halfedges.
      curr = next;
      next = curr->next()->opposite();

// Make sure we have not completed a full traversal around v without
// locating a place for the new curve xc.
      CGAL_assertion(curr != first);
    }

// Return the halfedge we have located.

  assert(false);
  return curr;
}
  


/*! \brief Return the face that lies to the right of the given vertex
 *  \prec v is on the ARR_TOP_BOUNDARY/ARR_BOTTOM_BOUNDARY identification.
 */
template <class GeomTraits, class Dcel>
typename Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::Face *
Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::
_face_right_of_vertex_on_top_identification (Vertex* v) const {
  this->check();
  std::cout << "_face_right_of_vertex_on_top_identification begin " << std::endl;
  
  CGAL_assertion(v != m_r_fictitious_vertex); 
  CGAL_assertion(v->parameter_space_in_y() != ARR_INTERIOR);

  // If the vertex is isolated, just return the face that contains it.
  if (v->is_isolated())
    return (v->isolated_vertex()->face());
  
  // Get the first incident halfedge around v 
  Halfedge  *current = v->halfedge();  
  // If there is only one halfedge incident to v, return its incident face.
  if (current ==  current->next()->opposite()) return _get_face(current);  
  Halfedge  *first = current;
  Halfedge  *max_e = first;  
  current = current->next()->opposite();
  while(current != first){
    max_e = _max(max_e,current);
    current = current->next()->opposite();
  }
  CGAL_assertion(_get_face(max_e) != this->m_l_fictitious_face);
  CGAL_assertion(_get_face(max_e) != this->m_r_fictitious_face);
  return _get_face(max_e);
}

//-----------------------------------------------------------------------------
// Get the curve associated with a boundary vertex.
template <class GeomTraits, class Dcel_>
const typename 
Arr_cylindrical_topology_traits_2<GeomTraits, Dcel_>::X_monotone_curve_2& 
Arr_cylindrical_topology_traits_2<GeomTraits, Dcel_>::
_curve (const Vertex *v, Arr_curve_end& ind) const
{
  this->check();
  // Go over the incident halfedges of v until encountering the halfedge
  // associated with a valid curve (v should have three incident halfedges,
  // two of them are fictitious and one associated with a curve).
  const Halfedge         *he = v->halfedge();
  CGAL_assertion(he != NULL); // vertex under construction ? 

  while (he->has_null_curve())
    {
      he = he->next()->opposite();
      if (he == v->halfedge())
        CGAL_error_msg("boundary vertex without supporting curve");
    }

  // The halfedge he is directed toward v, so if it is directed from left to
  // right, v represents the maximal end of cv, otherwise it represents its
  // minimal end.
  ind = (he->direction() == ARR_LEFT_TO_RIGHT) ? ARR_MAX_END : ARR_MIN_END;
  
  // Return the x-monotone curve.
  CGAL_assertion(!he->has_null_curve());
  return he->curve();
}

/*! \brief determines whether prev1 will be incident to the newly created face
 * (which will become a hole in the other face), as the result of an insertion
 * of a new halfedge
 */
template <class GeomTraits, class Dcel>
bool
Arr_cylindrical_topology_traits_2<GeomTraits, Dcel>::
is_on_new_perimetric_face_boundary(
    const Halfedge * prev1,
    const Halfedge * prev2,
    const X_monotone_curve_2 & xc, 
    bool try_other_way) const
{
  this->check();
  /*! We need to maintain the variant that the face that contains everything,
   * and has no outer CCB's, also contains the north pole. In the degenerate
   * case, where the north pole coincides with a vertex, the face that
   * contains everythin is incident to the north pole.
   * We count the number of times that path from prev1 to prev2 crosses the
   * discontinuity arc from left and from right, that is from
   * ARR_LEFT_BOUNDARY to ARR_RIGHT_BOUNDARY, and the number of times it
   * crosses the other way around.
   */
  assert(false);
  return true; 
}


//-----------------------------------------------------------------------------
// Check whether the given infinite curve end lies on the given fictitious
// halfedge.
//
template <class GeomTraits, class Dcel_>
bool 
Arr_cylindrical_topology_traits_2<GeomTraits, Dcel_>::
_is_on_fictitious_edge (const X_monotone_curve_2& cv, Arr_curve_end ind,
                        Arr_parameter_space ps_x, Arr_parameter_space ps_y,
                        Halfedge *he,
                        bool& eq_source, bool& eq_target) const 
{
  this->check();
  CGAL_assertion(he != NULL);
  
  CGAL_assertion(
      _get_face(he) == m_l_fictitious_face || 
      _get_face(he) == m_r_fictitious_face );
 
  eq_source = false;
  eq_target = false;

  // Get the end-vertices of the edge.
  const Vertex      *source = he->opposite()->vertex();
  const Vertex      *target = he->vertex();
 
#if 0 
  Arr_curve_end dummy; 
  if(source == m_l_fictitious_vertex)
    std::cout << "source is m_l_fictitious_vertex" << std::endl;
  else
    std::cout << "source is on " << _curve(source,dummy) << std::endl; 

  if(target == m_l_fictitious_vertex)
    std::cout << "target is m_l_fictitious_vertex" << std::endl;
  else
    std::cout << "target is on " << _curve(target,dummy) << std::endl; 
#endif   

  Comparison_result  comp_source, comp_target;
  
  // This must be a vertical ficticious edge 
  CGAL_assertion(source != NULL);
  CGAL_assertion(target != NULL);
  CGAL_assertion(ps_y == ARR_INTERIOR);
  CGAL_assertion(ps_x == ARR_LEFT_BOUNDARY || ps_x == ARR_RIGHT_BOUNDARY);
  CGAL_assertion(source->parameter_space_in_x() == ps_x);
  CGAL_assertion(target->parameter_space_in_x() == ps_x);
  CGAL_assertion(ind == (ps_x == ARR_LEFT_BOUNDARY)?ARR_MIN_END:ARR_MAX_END);

  // Compare the y-position of the curve end to the source vertex.
  if (source == this->m_l_fictitious_vertex || source == this->m_r_fictitious_vertex){
    // if the curve is oriented from left to right the 
    // ficticious vertex as source vertex counts as the lowest vertex 
    comp_source = (he->direction()==ARR_LEFT_TO_RIGHT)?SMALLER:LARGER;
  }else{ 
    Arr_curve_end v_ind; 
    comp_source = this->m_traits->compare_y_near_boundary_2_object() 
      (_curve (source, v_ind), cv, ind);
    CGAL_assertion(v_ind == ind); 
    if (comp_source == EQUAL) {
      eq_source = true;
      return (true);
    }
  }
  
  // Compare the y-position of the curve end to the target vertex.
  if (target == this->m_l_fictitious_vertex || target == this->m_r_fictitious_vertex){
    // if the curve is oriented from left to right the 
    // ficticious vertex as target vertex counts as the highest vertex 
    comp_target = (he->direction()==ARR_LEFT_TO_RIGHT)?LARGER:SMALLER;
  }else{ 
    Arr_curve_end v_ind; 
    comp_target = this->m_traits->compare_y_near_boundary_2_object() 
      (_curve (target, v_ind), cv, ind);
    CGAL_assertion(v_ind == ind); 
    if (comp_target == EQUAL) {
      eq_target = true;
      return (true);
    }
  }
  if(source == target){ // the edge is a full loop. 
    CGAL_assertion(source == this->m_l_fictitious_vertex || source == this->m_r_fictitious_vertex);
    CGAL_assertion(comp_source != comp_target);
  }
//  std::cout << "comp_source: " << comp_source << std::endl;
//  std::cout << "comp_target: " << comp_target << std::endl;
  return (comp_source != comp_target);
}

 /*! Split a fictitious edge using the given vertex.
   * \param e The edge to split (one of the pair of halfedges).
   * \param v The split vertex.
   * \pre e is a fictitious halfedge.
   * \return A halfedge whose direction is the same as e's and whose target is
   *         the split vertex v.
   */
template <class GeomTraits, class Dcel_>
typename Arr_cylindrical_topology_traits_2<GeomTraits, Dcel_>::Halfedge * 
Arr_cylindrical_topology_traits_2<GeomTraits, Dcel_>::
split_fictitious_edge(Halfedge *e, Vertex *v)
{
  this->check();
  CGAL_precondition(v->parameter_space_in_x() != ARR_INTERIOR);
  CGAL_precondition(v->parameter_space_in_y() == ARR_INTERIOR);
  
  std::cout << "split_fictitious_edge at: " << std::endl; 
  std::cout << v->parameter_space_in_x() << std::endl; 
  std::cout << v->parameter_space_in_y() << std::endl; 

  // Increment the number of vertices at infinity.
  n_inf_verts++;

  // Get the split halfedge and its twin, and their incident faces.
  // Note that we initialized the fictitious edges as inner boundaries 
  // of the cylindrical face. However, it is likely that fictitious edges 
  // become part of an outer_ccb, as they belong to faces that are 
  // considered as holes in the cylindrical face. 
  // In any case, its opposite belongs to the outer_ccb of one of the 
  // fictitious faces. 

  // get edges to be split 
  Halfedge       *he1 = e;
  Halfedge       *he2 = he1->opposite();
  // Allocate a pair of new halfedges.
  Halfedge       *he3 = this->m_dcel.new_edge();
  Halfedge       *he4 = he3->opposite();

  // Connect the new halfedges:
  //
  //            he1      he3
  //         -------> ------->
  //       (.)      (.)v     (.)
  //         <------- <-------
  //            he2      he4
  // 
  //       m_[l,r]_fictitious_face
  //

  // Connect e3 between e1 and its successor.
  he3->set_next (he1->next());
  // Insert he4 between he2 and its predecessor.
  he2->prev()->set_next (he4);

  // Set the properties of the new halfedges.
  he3->set_vertex (he1->vertex());

  v->set_halfedge (he4);
  he4->set_vertex (v);
  he4->set_next (he2);

  if (he1->vertex()->halfedge() == he1)
    // If he1 is the incident halfedge to its target, he3 replaces it.
    he1->vertex()->set_halfedge (he3);

  // Update the properties of the twin halfedges we have just split.
  he1->set_next(he3);
  he1->set_vertex(v);

  // The direction of he3 is the same as he1's (and the direction of he4 is
  // the same as he2).
  he3->set_direction (he1->direction());


  // all inits related to inner/outer ccbs
  // he2 is on outer ccb of an fictitious face 
  CGAL_assertion (!he2->is_on_inner_ccb());
  Outer_ccb      *oc2 = he2->outer_ccb();
  CGAL_assertion (
      oc2->face() == m_l_fictitious_face || 
      oc2->face() == m_r_fictitious_face);
  he4->set_outer_ccb (oc2);

  if( he1->is_on_inner_ccb()){
    Inner_ccb      *ic1 = he1->inner_ccb();
    CGAL_assertion (ic1->face()->is_unbounded());
    CGAL_assertion (ic1->face() != m_l_fictitious_face);
    CGAL_assertion (ic1->face() != m_r_fictitious_face);
    he3->set_inner_ccb (ic1);
  }else{
    Outer_ccb      *oc1 = he1->outer_ccb();
    CGAL_assertion (oc1->face()->is_unbounded());
    CGAL_assertion (oc1->face() != m_l_fictitious_face);
    CGAL_assertion (oc1->face() != m_r_fictitious_face);
    he3->set_outer_ccb (oc1);
  }
  
  CGAL_assertion(he1->prev()->next() == he1);
  CGAL_assertion(he2->prev()->next() == he2);
  CGAL_assertion(he3->prev()->next() == he3);
  CGAL_assertion(he4->prev()->next() == he4);
  CGAL_assertion(he1->next()->prev() == he1);
  CGAL_assertion(he2->next()->prev() == he2);
  CGAL_assertion(he3->next()->prev() == he3);
  CGAL_assertion(he4->next()->prev() == he4);


  // Return a pointer to one of the existing halfedge that is incident to the
  // split vertex.
  return (he1);
}

// compare a point to curve end 
template <class GeomTraits, class Dcel_>
Comparison_result 
Arr_cylindrical_topology_traits_2<GeomTraits, Dcel_>::_compare_xy
( const Point_2& p, const X_monotone_curve_2& cv, Arr_curve_end ind ) const{
  
  this->check();
  Arr_parameter_space ps_x_p =  m_traits->parameter_space_in_x_2_object()(p);
  Arr_parameter_space ps_y_p =  m_traits->parameter_space_in_y_2_object()(p);
  Arr_parameter_space ps_x_ce = m_traits->parameter_space_in_x_2_object()(cv,ind);
  Arr_parameter_space ps_y_ce = m_traits->parameter_space_in_y_2_object()(cv,ind);
  
  CGAL_assertion(ps_x_p == ARR_INTERIOR);
  if(ps_x_ce == ARR_LEFT_BOUNDARY)  return LARGER; 
  if(ps_x_ce == ARR_RIGHT_BOUNDARY) return SMALLER; 
  CGAL_assertion(ps_x_ce == ARR_INTERIOR);
  
  Point_2 p2 = _get_end_point(cv,ind);
  if(ps_y_ce ==  ARR_BOTTOM_BOUNDARY){ 
    Comparison_result res = m_traits->compare_x_2_object()(p,p2);
    if(res == EQUAL){
      if(ps_y_p != ARR_INTERIOR) return EQUAL; 
      else return LARGER; 
    }
    return res; 
  }else{
    return  m_traits->compare_xy_2_object()(p,p2);
  }
}


} //namespace CGAL

#endif
