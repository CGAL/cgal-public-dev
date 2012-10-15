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
// Author(s): Ophir Setter          <ophirset@post.tau.ac.il>
//            Michael Hemmer        <hemmer@mpi-inf.mpg.de>
//
//
#ifndef CGAL_VDL3_EXCEPTIONS_H
#define CGAL_VDL3_EXCEPTIONS_H

#include <stdexcept>

namespace CGAL {

namespace VDOL_3
{

class Event_at_discontinuity_exception : public std::range_error
{
public:
  Event_at_discontinuity_exception(const std::string &s)
    : std::range_error(s) {}

  ~Event_at_discontinuity_exception() throw() {}
};

} // end of namespace VDOL_3

} //namespace CGAL

#endif // CGAL_VDL3_EXCEPTIONS_H
