// Copyright (c) 2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Number_types/include/CGAL/Sqrt_extension.h $
// $Id: Sqrt_extension.h 47264 2008-12-08 06:25:14Z hemmer $
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>

#include <CGAL/basic.h>
#include <CGAL/Voronoi_diagram_of_lines_3.h>


#ifndef CGAL_VORONOI_DIAGRAM_OF_LINES_HIERARCHY_3_H
#define CGAL_VORONOI_DIAGRAM_OF_LINES_HIERARCHY_3_H

namespace CGAL {

template <class VoronoiDiagramOfLines_3>
class Voronoi_diagram_of_lines_hierarchy_3{
public:
  typedef VoronoiDiagramOfLines_3 Voronoi_diagram_of_lines_3;
  
  typedef Voronoi_diagram_of_lines_hierarchy_3<Voronoi_diagram_of_lines_3> Self; 
  typedef Voronoi_diagram_of_lines_3 VDOL_3;
  
  typedef typename VDOL_3::Line_3  Line_3;
  typedef typename VDOL_3::Point_3 Point_3;
  typedef typename VDOL_3::FT      FT;
  typedef typename VDOL_3::RT      RT;
 
  
  std::vector<VDOL_3> m_hierarchy;
  const std::vector<VDOL_3>& hierarchy() const {return m_hierarchy;}
  std::vector<VDOL_3>& hierarchy()       {return m_hierarchy;}
  
  mutable boost::shared_ptr<VDOL_3> m_vdol_3; 
  const VDOL_3& vdol_3() const {return *m_vdol_3;}
  //VDOL_3& vdol_3()       {return m_vdol_3;}
  
 
  template<class InputIterator> 
  Voronoi_diagram_of_lines_hierarchy_3(InputIterator begin, InputIterator end, int k = 2){

    std::vector<Line_3> lines(begin,end);

   
    if(k>=2){
      timer_vdol_hierarchy_total.start();  

      long KK=k; 
      while((KK*=k)<lines.size()){
        std::cerr << KK << std::endl;
      }
      KK/=k;
      
      while(KK >= k){
        int i = lines.size()/KK;
        std::cerr <<i << " " <<  lines.size() << std::endl; 
        if(i>1){
          m_hierarchy.push_back(VDOL_3(lines.begin(),lines.begin()+i));
        }
        KK/=k; 
      }
      std::cerr << "SIZE" <<  m_hierarchy.size() << std::endl; 
      timer_vdol_hierarchy_total.stop();
    }
    
    timer_vdol_total.start();
    m_vdol_3 = boost::shared_ptr<VDOL_3>(
        new VDOL_3(lines.begin(),lines.end()));
    timer_vdol_total.stop();
  }
 
  template<class Point_3, class OutputIterator>
  OutputIterator locate(const Point_3& p, OutputIterator oit, int hint = 0) const{

    for(int i = 0; i < m_hierarchy.size(); i++){
      std::vector<int> indices;
      m_hierarchy[i].locate(p,std::back_inserter(indices),hint);
      CGAL_postcondition(indices.size()>=1);
      hint = indices[0];
    }
    return m_vdol_3->locate(p,oit, hint); 
  }
  
};

}


#endif // CGAL_VORONOI_DIAGRAM_OF_LINES_HIERARCHY_3_H
