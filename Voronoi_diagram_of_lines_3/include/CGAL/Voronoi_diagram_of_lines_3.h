// Copyright (c) 2006-2010 Max-Planck-Institute Saarbruecken (Germany),
 // INRIA Sophia-Antipolis (France). All rights reserved.
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
#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Get_arithmetic_kernel.h>
#include <CGAL/Voronoi_cell_of_line_3.h>
#include <boost/variant.hpp>
#include <boost/foreach.hpp>


#ifndef CGAL_VORONOI_DIAGRAM_OF_LINES_3_H
#define CGAL_VORONOI_DIAGRAM_OF_LINES_3_H

namespace CGAL {
// TODO add Algebraic_kernel as optional template argument. 
template <class LinearKernel, class ArithmeticKernel = CGAL::Arithmetic_kernel>
class Voronoi_diagram_of_lines_3{
public:
  typedef LinearKernel                                                Linear_kernel; 
  typedef ArithmeticKernel                                            Arithmetic_kernel; 
  typedef Voronoi_diagram_of_lines_3<Linear_kernel,Arithmetic_kernel> Self; 

  typedef typename Linear_kernel::FT       FT; 
  typedef typename Linear_kernel::RT       RT; 
  typedef typename Linear_kernel::Point_3  Point_3; 
  typedef typename Linear_kernel::Line_3   Line_3;
  
private:
  typedef typename Arithmetic_kernel::Integer           Integer; 
  typedef typename Arithmetic_kernel::Rational          Rational; 
  typedef typename Arithmetic_kernel::Bigfloat_interval BFI; 
  
public:
  typedef CGAL::Voronoi_cell_of_line_3<Linear_kernel,Arithmetic_kernel>  VCOL_3;   
  typedef typename VCOL_3::Poly_int_3   Poly_int_3;
  typedef typename VCOL_3::Poly_int_2   Poly_int_2;
  typedef typename VCOL_3::Poly_int_1   Poly_int_1;
  
private: 
// ##########################################

  std::vector<Line_3>       m_lines;
  std::vector<VCOL_3>       m_cells; 
  
  mutable std::map<int,int> m_line_index_map;
  boost::shared_ptr<CGAL::Random> mp_rand; 

public:
  typedef typename std::vector<Line_3>::const_iterator Lines_iterator;
  const std::vector<Line_3>& lines() const {return m_lines;}
        std::vector<Line_3>& lines()       {return m_lines;}
  

  typedef typename std::vector<VCOL_3>::const_iterator Cell_iterator;
  const std::vector<VCOL_3>& cells() const {return m_cells;}
        std::vector<VCOL_3>& cells()       {return m_cells;}
  

  template <class InputIterator>
  Voronoi_diagram_of_lines_3(InputIterator begin, InputIterator end, int vd_seed = 0)
    :m_lines(begin,end)
  {
    mp_rand = boost::shared_ptr<CGAL::Random>(new CGAL::Random(vd_seed));
    std::cerr << " number of lines " <<  m_lines.size() << std::endl;
    
    for(unsigned int i = 0; i < m_lines.size(); i++){
      m_line_index_map[compute_hash_value(m_lines[i])] = i;
    }

    CGAL_postcondition(m_line_index_map.size() == m_lines.size());

    for(unsigned int i = 0; i < m_lines.size(); i++){
      //for(int i = 0; i < 1 ; i++){
      std::cerr << " cell number " <<  i  << std::endl;
      std::swap(m_lines[0],m_lines[i]);

      long seed_range = 16; 
      bool done = false; 
      while(!done){
        std::cerr << "try with seed " << seed_range << std::flush; 
        try{
          m_cells.push_back(VCOL_3( mp_rand->get_int(0,seed_range)+1,m_lines.begin(),m_lines.end()));
          done = true ;
        }catch(CGAL::VDOL_3::Event_at_discontinuity_exception& e){
#ifdef CGAL_VDOL_USE_COUNTER
          count_catch_frame++;
#endif
          std::cerr << "catch:" << std::endl 
                    << e.what() << std::endl;
          seed_range*=2; 
          if(seed_range > 100) {
            std::cerr << "infinite try/catch unexpected special case ? " 
                      << std::endl;
            // assert(false);
          }
        } 
      }
      std::swap(m_lines[0],m_lines[i]);// swap back 
    }
 
#ifdef CGAL_VDOL_CELL_SAMPLE_CURVES
    init_points();
#endif
  }
  
  template<class OutputIterator>
  OutputIterator locate(const Point_3& p, OutputIterator oit, int hint = 0) const {
    CGAL_precondition(hint < m_cells.size());
    
    std::vector<Line_3> lines;
    bool found = false;
    while(!found){
      lines.clear();
      found = m_cells[hint].contains(p,std::back_inserter(lines));
      CGAL_postcondition(found || lines.size()>=0);
      hint = (lines.size())>=0?m_line_index_map[compute_hash_value(lines[0])]:0;
    }
    std::set<int> line_indices;
    for(int i = 0; i < lines.size();i++){
      // std::cerr << m_line_index_map[compute_hash_value(lines[i])] << " - "  ;
      line_indices.insert(m_line_index_map[compute_hash_value(lines[i])]);
    }
    return copy(line_indices.begin(),line_indices.end(),oit);
  } 

public: 
  
  template <class OutputIterator>
  OutputIterator generate_points ( 
      OutputIterator oit, const VDOL_3::Approximation_info& approximation_info) const {
    oit =  generate_vertex_points(oit,approximation_info); 
    oit =  generate_edge_points(oit,approximation_info); 
    oit =  generate_face_points(oit,approximation_info); 
   
    return oit;
  }
  
  template <class OutputIterator>
  OutputIterator 
  generate_vertex_points(OutputIterator oit, const VDOL_3::Approximation_info& approximation_info) const {
    for(Cell_iterator cit = cells().begin(); cit != cells().end(); cit++){
      oit = cit->generate_vertex_points(oit,approximation_info);
    }
    return oit; 
  }
  
  template <class OutputIterator>
  OutputIterator 
  generate_edge_points(OutputIterator oit, const VDOL_3::Approximation_info& approximation_info) const {
    for(Cell_iterator cit = cells().begin(); cit != cells().end(); cit++){
      oit = cit->generate_edge_points(oit,approximation_info);
    }
    return oit; 
  }


  
 template <class OutputIterator>
 OutputIterator 
 generate_face_points(OutputIterator oit, const VDOL_3::Approximation_info& approximation_info) const {  
   for(Cell_iterator cit = cells().begin(); cit != cells().end(); cit++){
     oit = cit->generate_face_points(oit,approximation_info);
   }
   return oit; 
 }
  // reports total number of vertices (first) and (second) number of vertices
  // that are really relevant for the voronoi diagram
  std::pair<int,int> number_of_vertices() const {
    std::pair<int,int> result(0,0); 
    std::pair<int,int> cresult(0,0);
    BOOST_FOREACH(const VCOL_3& cell, cells()){
      cresult = cell.number_of_vertices();
      result.first  += cresult.first; 
      result.second += cresult.second; 
    }
    return result;
  }
  
  template <class OutputIterator>
  OutputIterator report_bisectors(OutputIterator oit) const {
    BOOST_FOREACH(const VCOL_3& cell, cells()){
      oit = cell.report_bisectors(oit); 
    }
    return oit; 
  } 
  
  template <class OutputIterator>
  OutputIterator report_transformed_bisectors(OutputIterator oit) const {
    BOOST_FOREACH(const VCOL_3& cell, cells()){
      oit = cell.report_transformed_bisectors(oit); 
    }
    return oit; 
  }
  
  template <class OutputIterator>
  OutputIterator report_resultants(OutputIterator oit) const {
    BOOST_FOREACH(const VCOL_3& cell, cells()){
      oit = cell.report_resultants(oit); 
    }
    return oit; 
  }
};




} // namespace CGAL 
#endif // CGAL_VORONOI_DIAGRAM_OF_LINES_3_H 
