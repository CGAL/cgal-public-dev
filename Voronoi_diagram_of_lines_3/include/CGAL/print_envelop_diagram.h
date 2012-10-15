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

#ifndef CGAL_PRINT_ENVELOP_DIAGRAM_H
#define CGAL_PRINT_ENVELOP_DIAGRAM_H


namespace CGAL {

template <class EnvelopeDiagram_2>
void print_envelop_diagram(const EnvelopeDiagram_2& diag){
  
  typedef EnvelopeDiagram_2 Envelope_diagram_2;
  
  //typedef typename Envelope_diagram_2::
  
  typename Envelope_diagram_2::Face_const_iterator            fit;
  typename Envelope_diagram_2::Edge_const_iterator           eit;
  typename Envelope_diagram_2::Vertex_const_iterator        vit;

  std::cerr << "Envelope_diagram_2 number of faces / edges / vertices :  " 
            << diag.number_of_faces ()    << " / " 
            << diag.number_of_edges ()    << " / " 
            << diag.number_of_vertices () << std::endl;
  
  std::cerr << "Number of xy_monotone_surfaces  for faces:     "  ;
    fit = diag.faces_begin();
  while(fit != diag.faces_end()){
    std::cerr << fit->number_of_surfaces();
    fit++;
  }
  std::cerr  << std::endl;
  
  std::cerr << "Number of xy_monotone_surfaces  for edges:     "  ;
  eit = diag.edges_begin();
  while(eit != diag.edges_end()){
    std::cerr << eit->number_of_surfaces();
    eit++;
  }
  std::cerr  << std::endl;
  
  std::cerr << "Number of xy_monotone_surfaces  for vertices:  "  ;
  vit = diag.vertices_begin();
  while(vit != diag.vertices_end()){
    std::cerr << vit->number_of_surfaces();
    vit++;
  }
  std::cerr  << std::endl;
}

} //namespace CGAL

#endif // CGAL_PRINT_ENVELOP_DIAGRAM_H
