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
// $Id:
// 
//
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>
//            Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//


#ifndef CGAL_VDL3_SVCET_33_FUNCTOR_BASE_H
#define CGAL_VDL3_SVCET_33_FUNCTOR_BASE_H

#include <CGAL/Sweep_line_2_algorithms.h>

namespace CGAL{ 
namespace VDOL_3 {
namespace SVCET_3_functors {

template < class SingleVoronoiCellEnvelopeTraits_3> 
class Single_voronoi_cell_envelope_traits_3_functor_base {
public: 
  //! this instance's template parameter
  typedef SingleVoronoiCellEnvelopeTraits_3
    Single_voronoi_cell_envelope_traits_3;
  
  // CGAL_SNAP_SVCET_3_TYPEDEFS;   
  
protected:
  const Single_voronoi_cell_envelope_traits_3* m_svcet_3;
  
protected:
  const Single_voronoi_cell_envelope_traits_3* svcet_3() const {
    return m_svcet_3;
  }
  
public:
  Single_voronoi_cell_envelope_traits_3_functor_base(
      const Single_voronoi_cell_envelope_traits_3 *traits) : m_svcet_3(traits) {
    CGAL_precondition(traits != NULL);
  }  
};

} // namespace SVCET_3_functors
} // namespace VDOL_3 
} // namespace CGAL


  

#endif // CGAL_VDL3_SVCET_33_FUNCTOR_BASE_H
