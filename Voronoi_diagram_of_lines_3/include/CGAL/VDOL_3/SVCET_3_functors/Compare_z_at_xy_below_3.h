// Copyright (c) 2005-2009 Tel-Aviv University (Israel).
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
// $Id:
// 
//
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>
//            Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//

// SVCET_3 = Single_voronoi_cell_envelope_traits_3


#ifndef CGAL_VDL3_SVCET_33_COMPARE_Z_AT_XY_BELOW_3
#define CGAL_VDL3_SVCET_33_COMPARE_Z_AT_XY_BELOW_3

#include <CGAL/config.h>
#include <CGAL/VDOL_3/SVCET_3_functors/Compare_z_at_xy_near_arc_3.h>

namespace CGAL { 
namespace VDOL_3 {
namespace SVCET_3_functors {

template < class SVCET_3 >
class Compare_z_at_xy_below_3 
  : public Compare_z_at_xy_near_arc_3<SVCET_3>{
public:
  //! first template parameter
  typedef SVCET_3 Single_voronoi_cell_envelope_traits_3;
  //! the base type 
  typedef  Compare_z_at_xy_near_arc_3<SVCET_3> Base;
  
  CGAL_SNAP_SVCET_3_TYPEDEFS;
  Compare_z_at_xy_below_3(const SVCET_3 *traits): Base(traits,false) {};
};

} // namespace SVCET_3_functors
} // namespace VDOL_3 
} // namespace CGAL

#endif // CGAL_VDL3_SVCET_33_COMPARE_Z_AT_XY_BELOW_3
