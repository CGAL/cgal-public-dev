// Copyright (c) 2010  Tel-Aviv University (Israel).
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
// Author(s)     : Asaf Porat          <asafpor1@post.tau.ac.il>


#ifndef LINES_THROUGH_SEGMENTS_EXCEPTIONS_H
#define LINES_THROUGH_SEGMENTS_EXCEPTIONS_H

#include <exception>
namespace CGAL {

class Lines_through_segments_exp_2_lines_overlap: public std::exception
{
public:
  virtual const char* what() const throw()
  {
    return "Lines_through_segments_exp_2_lines_overlap";
  }
};

} //namespace CGAL

#endif
