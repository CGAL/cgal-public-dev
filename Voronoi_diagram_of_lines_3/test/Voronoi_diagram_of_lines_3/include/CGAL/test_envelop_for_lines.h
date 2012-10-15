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

#ifndef CGAL_TEST_ENVELOP_FOR_LINES_H
#define CGAL_TEST_ENVELOP_FOR_LINES_H

#include <CGAL/check_envelop_diagram.h>
#include <CGAL/print_envelop_diagram.h>

namespace CGAL {

template<class SVCET_3> 
void test_envelop_for_lines( SVCET_3 svcet_3, std::vector<typename SVCET_3::Line_3> lines){
  typename std::vector<typename SVCET_3::Line_3>::iterator begin = lines.begin();
  typename std::vector<typename SVCET_3::Line_3>::iterator end = lines.end();
 
  std::cout<< "Start Envelope_3 for positive plane for "<< lines.size() << " lines... "<< std::endl;
  CGAL::Timer timer_1;
  timer_1.start();
  typedef typename CGAL::Envelope_diagram_2<SVCET_3> Envelope_diagram_2;
  Envelope_diagram_2 diag(&svcet_3); 
  CGAL::lower_envelope_3(begin, end, diag);
  timer_1.stop();
  print_envelop_diagram(diag);
  std::cout<< "Time: "<< timer_1.time() << std::endl;

  std::cout<< "Start Envelope_3 for negative plane for "<< lines.size()<< " lines... "<< std::endl;
  CGAL::Timer timer_2;
  timer_2.start();

  SVCET_3 svcet_3_mirror = svcet_3.mirror();

  Envelope_diagram_2 diag_mirror(&svcet_3_mirror); 
  CGAL::lower_envelope_3(begin, end, diag_mirror);
  timer_2.stop();
  print_envelop_diagram(diag_mirror);
  std::cout<< "Time: "<< timer_2.time() << std::endl;
  
  std::cout<< "Check positive plane ... ";
  std::cout<< std::flush;
  check_envelop_diagram(diag,lines);
  std::cout<< "-------------------- ok" << std::endl;
  
  std::cout<< "Check negative plane ... ";
  std::cout<< std::flush;
  check_envelop_diagram(diag_mirror,lines);
  std::cout<< "-------------------- ok" << std::endl;

  std::cout<< "Check consitency of both ... ";
  check_opposite_diagrams(diag,diag_mirror);
  std::cout<< "-------------------- ok" << std::endl;
  std::cout<<std::endl;
}

template<class SVCET_3, class InputIterator> 
void test_envelop_for_line_range( InputIterator begin, InputIterator end)
{
  typedef CGAL::VDOL_3::Event_at_discontinuity_exception
    Event_at_discontinuity_exception;
  CGAL_assertion(begin != end); 
  int seed = -1; 
  typename SVCET_3::Line_3 line_0(*begin++);
  std::vector<typename SVCET_3::Line_3> lines(begin,end);
  
  for(int i=0; i<3; i++){
    std::cout << ++seed << std::endl;
    bool done = false;
    while(!done){
      std::cout << "try ..." <<std::flush; 
      try{
        SVCET_3 svcet_3(line_0,seed,CGAL::POSITIVE);
        test_envelop_for_lines(svcet_3,lines);
        std::cout << "success" << std::endl;
        done = true ;
      }catch(Event_at_discontinuity_exception& e){
        std::cout << "catch:" << std::endl 
                  << e.what() << std::endl;
        seed++; 
        
        if(seed > 10) {
          std::cout << "infinite try/catch unexpected special case ? " 
                    << std::endl;
          assert(false);
        }
      } 
    }
  }
}


template<class SVCET_3> 
void test_envelop_for_lines(
    const typename SVCET_3::Line_3& line_0,
    const typename SVCET_3::Line_3& line_1
){
  std::vector<typename SVCET_3::Line_3> lines;
  lines.push_back(line_0);
  lines.push_back(line_1);
  typedef typename std::vector<typename SVCET_3::Line_3>::iterator IT;
  test_envelop_for_line_range<SVCET_3,IT>(lines.begin(),lines.end());
}
template<class SVCET_3> 
void test_envelop_for_lines(
    const typename SVCET_3::Line_3& line_0,
    const typename SVCET_3::Line_3& line_1,
    const typename SVCET_3::Line_3& line_2
){
  std::vector<typename SVCET_3::Line_3> lines;
  lines.push_back(line_0);
  lines.push_back(line_1);
  lines.push_back(line_2);
  typedef typename std::vector<typename SVCET_3::Line_3>::iterator IT;
  test_envelop_for_line_range<SVCET_3,IT>(lines.begin(),lines.end());
}
template<class SVCET_3> 
void test_envelop_for_lines(
    const typename SVCET_3::Line_3& line_0,
    const typename SVCET_3::Line_3& line_1,
    const typename SVCET_3::Line_3& line_2,
    const typename SVCET_3::Line_3& line_3
){
  std::vector<typename SVCET_3::Line_3> lines;
  lines.push_back(line_0);
  lines.push_back(line_1);
  lines.push_back(line_2);
  lines.push_back(line_3);
  typedef typename std::vector<typename SVCET_3::Line_3>::iterator IT;
  test_envelop_for_line_range<SVCET_3,IT>(lines.begin(),lines.end());
}
template<class SVCET_3> 
void test_envelop_for_lines(
    const typename SVCET_3::Line_3& line_0,
    const typename SVCET_3::Line_3& line_1,
    const typename SVCET_3::Line_3& line_2,
    const typename SVCET_3::Line_3& line_3,
    const typename SVCET_3::Line_3& line_4
){
  std::vector<typename SVCET_3::Line_3> lines;
  lines.push_back(line_0);
  lines.push_back(line_1);
  lines.push_back(line_2);
  lines.push_back(line_3);
  lines.push_back(line_4);
  typedef typename std::vector<typename SVCET_3::Line_3>::iterator IT;
  test_envelop_for_line_range<SVCET_3,IT>(lines.begin(),lines.end());
}





} //namespace CGAL

#endif // CGAL_TEST_ENVELOP_FOR_LINES_H
