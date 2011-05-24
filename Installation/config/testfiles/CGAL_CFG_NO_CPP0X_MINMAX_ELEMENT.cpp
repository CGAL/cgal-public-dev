// Copyright (c) 2011  GeometryFactory (France).  All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Phillip Möller

//| If a compiler does not support C++0x minmax_element
//| CGAL_CFG_NO_CPP0X_MINMAX_ELEMENT is set. 

#undef NDEBUG
#include <algorithm> 
#include <cassert>

int main()
{
  int f[5] = {1, 2, 3, 4, 5};
  std::pair<int*, int*> ret = std::minmax_element(f, f + 5);
  assert(*(ret.first) == 1);
  assert(*(ret.second) == 5);
}
