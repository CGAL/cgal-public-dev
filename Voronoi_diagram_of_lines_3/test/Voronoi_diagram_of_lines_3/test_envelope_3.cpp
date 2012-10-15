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
// This file is provided AS IS with NO
// WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id:$
// 
//
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>
//            Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//
//

#define CGAL_USE_FICTITIOUS_MIRROR_SVCET 0
#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1 

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/macros.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Cached_algebraic_kernel_1.h>
#include <CGAL/Cached_algebraic_kernel_2.h>

#include <CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h> // template arg to CKvA_2
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h> // traits for Arr_2

#include <CGAL/VDOL_3/Single_voronoi_cell_envelope_traits_3.h>
// #include <CGAL/VDOL_3/Topology_traits_2.h>

#include <CGAL/envelope_3.h>
#include <CGAL/Timer.h>
#include <CGAL/Random.h>
 
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_unb_planar_topology_traits_2.h>

#include <CGAL/test_envelop_for_lines.h>


// get default arithmetic 
typedef CGAL::Arithmetic_kernel AK;
typedef AK::Integer Integer;
typedef AK::Rational Rational;

// define linear kernel 
typedef CGAL::Cartesian<Rational>                       Linear_kernel_3;
typedef Linear_kernel_3::Line_3                         Line_3;
typedef Linear_kernel_3::Point_3                        Point_3;

// define algebraic kerenel 
typedef CGAL::Algebraic_kernel_d_1<Integer>             AK_1;
typedef CGAL::Algebraic_curve_kernel_2<AK_1>            AK_2;

// define curved kernel 
typedef CGAL::Curved_kernel_via_analysis_2< AK_2 >      Curved_kernel_2;

// define traits for lower envelop diagram 
typedef CGAL::VDOL_3::Single_voronoi_cell_envelope_traits_3
<Linear_kernel_3, Curved_kernel_2>                 SVCET_3;

void copy_file_into(const char *infilename, const char *outfilename){
  std::ifstream     in_file (infilename);
  std::ofstream     out_file (outfilename);
  std::cout << "save " << infilename << " to " << outfilename << "... " << std::flush;
  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << infilename << "!" << std::endl;
    return; 
  }
  if (! out_file.is_open()) {
    std::cerr << "Failed to open " << outfilename << "!" << std::endl;
    return; 
  }
  out_file << in_file.rdbuf();
  std::cout << "  done "<< std::endl;
}

template <typename OutputIterator>
OutputIterator read_lines(const char *filename, OutputIterator oi){

  // Open the input file.
  std::ifstream     in_file (filename);
  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << filename << "!" << std::endl;
    return oi;
  }
  
  // Read the lines from the file
  // The input file format should be (all coordinate values are integers):
  // <n>                                       // number of lines.
  // <a1_x> <a1_y> <a1_z> <b1_x> <b1_y> <b1_z> // line #1.
  // <a2_x> <a2_y> <a2_z> <b2_x> <b2_y> <b2_z> // line #2.
  //   :      :       :      :
  // <an_x> <an_y> <an_z> <bn_x> <bn_y> <bn_z> // line #n.
  
  // read number of lines 
  unsigned int n; in_file >> n;

  for (int k = 0; k < n; ++k) {
    int a,b,c,d,e,f;
    in_file >> a >> b >> c>> d >> e >> f;
    typedef Point_3::FT RT; 
    *oi++ = Line_3(Point_3(RT(a),RT(b),RT(c)),Point_3(RT(d),RT(e),RT(f)));
  }
  in_file.close();
}

int main(int argc, char* argv[])
{
  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);
  
  // Get the name of the input file from the command line, or use the default
  // last.dat file if no command-line parameters are given.
  const char * filename = (argc > 1) ? argv[1] : "last.dat"; 
  if(argc > 1){
    copy_file_into(filename,"last.dat");
  }
  std::cout << "Run test with file: " << filename << std::endl;

  std::vector<SVCET_3::Surface_3> lines;
  read_lines(filename,std::back_inserter(lines));
  if(lines.size()==0) 
    std::cerr << "Failed to read file " << filename << std::endl;
  
  std::cout << "Number of lines (incl base line): "<< lines.size() << std::endl;
  
 
  CGAL::test_envelop_for_line_range<SVCET_3>(lines.begin(),lines.end());
  return 0; 
}
