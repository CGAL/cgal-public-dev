// Copyright (c) 2005-2007 Tel-Aviv University (Israel).
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
// Author(s): Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//
//

// Helper function to construct a (ratinal) CK_2::Point_2 

#ifndef CGAL_CHECK_ENVELOP_DIAGRAM_H
#define CGAL_CHECK_ENVELOP_DIAGRAM_H

#include <CGAL/construct_point_2.h>
#include <CGAL/VDOL_3/vol_tools.h>

namespace CGAL {

template <class SVCET_3>
void check_opposite_diagrams(
    CGAL::Envelope_diagram_2<SVCET_3> & diag_pos,
    CGAL::Envelope_diagram_2<SVCET_3> & diag_neg){


  typedef CGAL::Envelope_diagram_2<SVCET_3> Diagram;
  typedef typename Diagram::Topology_traits Topology_traits;
  typedef typename Topology_traits::Vertex Vertex; 
  typedef typename Topology_traits::Halfedge Halfedge;
  typedef typename Topology_traits::Face Face;


  // The ARR_LEFT_BOUNDARY of the positive plane should match 
  // the ARR_RIGHT_BOUNDARY of the negative plane and vice versa. 

  // check positive left boundary with negative right boundary. 
  Vertex *v_bl_pos = diag_pos.topology_traits()->bottom_left_vertex();  
  Vertex *v_tl_neg = diag_neg.topology_traits()->top_left_vertex();  
  // get a halfedge 
  Halfedge *he_pos = v_bl_pos->halfedge();
  Halfedge *he_neg = v_tl_neg->halfedge();
  // get the proper one (i.e. non vertical)
  if(VDOL_3::is_vertical(he_pos)) he_pos = he_pos->next()->opposite();
  if(VDOL_3::is_vertical(he_neg)) he_neg = he_neg->next()->opposite();
  // both must be non vertical 
  CGAL_assertion(!VDOL_3::is_vertical(he_pos));
  CGAL_assertion(!VDOL_3::is_vertical(he_neg));
  CGAL_assertion(he_pos->vertex() == v_bl_pos );
  CGAL_assertion(he_neg->vertex() == v_tl_neg );
  he_pos = he_pos->opposite();
  he_neg = he_neg->opposite();
  
 
  while(he_pos->vertex() != diag_pos.topology_traits()->bottom_right_vertex()){
    // no events at discontinuity line 
    CGAL_assertion(VDOL_3::degree(he_pos->vertex()) == 3); 
    CGAL_assertion(VDOL_3::degree(he_pos->vertex()) == 3);
    // check for common curve (may cannonicalize / gcd )
    CGAL_assertion( he_pos->next()->curve().curve().polynomial_2() == 
        he_neg->opposite()->prev()->opposite()->curve().curve().polynomial_2());
      
    // walk along the boundary of the ficticious face 
    he_pos = he_pos->opposite()->prev()->opposite();   
    he_neg = he_neg->next(); 
  }
  CGAL_assertion(he_neg->vertex() 
      == diag_neg.topology_traits()->top_right_vertex()); 
}






template <class SVCET_3>
void check_envelop_diagram_by_random_points(
    const CGAL::Envelope_diagram_2<SVCET_3> & diag,
    const std::vector<typename SVCET_3::Line_3>& lines){   

  // The test generatos some randome rational points. 
  // For each point a point location determines the object it lies in. 
  // Thereafter, it is veryfied that this is the correct object using 
  // ray shooting. 

  typedef CGAL::Envelope_diagram_2<SVCET_3> Envelope_diagram_2;
  typedef typename SVCET_3::Point_2 Point_2;
  typedef typename SVCET_3::Line_3 Line_3;
  typedef typename SVCET_3::FT FT;

  typedef typename SVCET_3::Linear_kernel Linear_kernel; 
  
  CGAL::Arr_naive_point_location<Envelope_diagram_2> point_location(diag);

  const SVCET_3*  svcet_3 = diag.geometry_traits();
  
  //CGAL::Random random(153);

  for(int i = -200; i <=200; i+= 50){
    for(int j = -200; j <=200; j+= 50){
 
      std::pair<FT,FT> rat_point = std::make_pair(FT(i),FT(j));

      Point_2 point = 
        CGAL::construct_point_2(svcet_3 ,rat_point.first , rat_point.second);
      CGAL::Object obj(point_location.locate(point));
      
    
      // Get the line of the located entity 
      // IF THERE IS ONE !! 
      Line_3 line;
      bool found_line = false;
      // now a very nice example of code duplication 
      typename Envelope_diagram_2::Face_const_handle       f;
      typename Envelope_diagram_2::Halfedge_const_handle   e;
      typename Envelope_diagram_2::Vertex_const_handle     v;
      if (CGAL::assign (f, obj)) { 
        if(f->number_of_surfaces() >0){
          line = f->surface().line(); 
          found_line = true; 
        }
      }else if (CGAL::assign (e, obj)) { 
        if(e->number_of_surfaces() >0){
          line = e->surface().line(); 
          found_line = true; 
        }
      }
      else if (CGAL::assign (v, obj)){
        if(v->number_of_surfaces() >0){
          line = v->surface().line(); 
          found_line = true; 
        }
      }
      else { CGAL_assertion_msg (false, "Invalid object."); }
      
      typename std::vector<Line_3>::const_iterator it;
      for(it = lines.begin(); it != lines.end(); it++){
        if(found_line){
          // line must be less or equal to all other lines. 
          assert( compare_lines_at(svcet_3,line,*it,rat_point) != CGAL::LARGER );
        }else{
          // It seems that there is no bisector in the sample direction 
          // that is, all bisectors are "at infinity" here, lets check:
          assert(is_at_infinity(svcet_3,*it,rat_point));  
        }
      }
    }
  }
}



template <class SVCET_3>
void check_envelop_diagram_features(
    const CGAL::Envelope_diagram_2<SVCET_3> & diag,
    const std::vector<typename SVCET_3::Line_3>& lines){ 
  
  // test all edges:
  // take a point on the edge 
  // and compare with all other lines, result must be less or equal. 
  // number of equal must be the number of associated lines 
  
  typedef CGAL::Envelope_diagram_2<SVCET_3> Envelope_diagram_2;
  CGAL_SNAP_SVCET_3_TYPEDEFS;

  typedef typename Envelope_diagram_2::Topology_traits Topology_traits;
  typedef typename Topology_traits::Vertex Vertex; 
  typedef typename Topology_traits::Halfedge Halfedge;
  typedef typename Topology_traits::Face Face;
  
  typedef typename Envelope_diagram_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Envelope_diagram_2::Surface_const_iterator Surface_const_iterator;
  typedef typename Envelope_diagram_2::Edge_const_iterator Edge_const_iterator;
  
  const SVCET_3 *svcet_3 = diag.geometry_traits();

  std::vector<Xy_monotone_surface_3> xy_surfaces;
  typedef typename std::vector<typename SVCET_3::Line_3>::const_iterator LIT;
  for(LIT lit = lines.begin(); lit != lines.end();lit++){
    svcet_3->make_xy_monotone_3_object()(*lit,true, std::back_inserter(xy_surfaces));
  }
  

  for( Edge_const_iterator eit = diag.edges_begin(); eit != diag.edges_end(); eit++){
   
    
    Xy_monotone_surface_3  surface = eit->surface(); 
    Point_2 point = svcet_3->construct_interior_vertex_2_object()(eit->curve());
    
    for( Surface_const_iterator sit = eit->surfaces_begin(); sit != eit->surfaces_end();sit++){
      assert(svcet_3->compare_z_at_xy_3_object()(point,surface,*sit) == CGAL::EQUAL);
    }       
    
    std::set<Poly_int_3> bisectors_A;
    typedef typename std::vector<Xy_monotone_surface_3>::const_iterator SIT;
    for(SIT sit = xy_surfaces.begin(); sit != xy_surfaces.end();sit++){
      CGAL::Comparison_result comp = svcet_3->compare_z_at_xy_3_object()(point,surface,*sit);
      assert(comp != CGAL::LARGER);
      if(comp == CGAL::EQUAL)      
        bisectors_A.insert(sit->bisector());
    }
    std::set<Poly_int_3> bisectors_B;
    for(typename Envelope_diagram_2::Surface_const_iterator sit = eit->surfaces_begin(); 
        sit != eit->surfaces_end(); sit++){
      bisectors_B.insert(sit->bisector());
    }

    if(eit->number_of_surfaces() == 1){
      std::cerr << " suspicious edge number_of_surfaces() == 1" << std::endl;
    }
    if(bisectors_A.size()!=bisectors_B.size() || eit->number_of_surfaces() == 1 ){
      std::cerr << eit->curve() << std::endl; 
      std::cerr << eit->curve().curve().polynomial_2() << std::endl; 
      std::cerr << "number of xy surfaces: " << eit->number_of_surfaces()<< std::endl;
      std::cerr <<" A: (found in test) "<< std::endl;
      std::copy(bisectors_A.begin(), bisectors_A.end(), std::ostream_iterator<Poly_int_3>(std::cerr,"\n")); 
      std::cerr <<" B: (given in envelope)"<< std::endl;
      std::copy(bisectors_B.begin(), bisectors_B.end(), std::ostream_iterator<Poly_int_3>(std::cerr,"\n"));
    }
   
    
    assert(bisectors_A.size()==bisectors_B.size());    
  }
}

template <class SVCET_3>
void check_envelop_diagram(
    const CGAL::Envelope_diagram_2<SVCET_3> & diag,
    const std::vector<typename SVCET_3::Line_3>& lines){ 
  check_envelop_diagram_features(diag,lines);
  check_envelop_diagram_by_random_points(diag,lines);
}
} //namespace CGAL

#endif // CGAL_CHECK_ENVELOP_DIAGRAM_H
