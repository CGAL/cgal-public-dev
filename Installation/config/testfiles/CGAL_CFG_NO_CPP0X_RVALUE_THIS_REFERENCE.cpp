// Copyright (c) 2011  INRIA Saclay (France).
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
// $URL$
// $Id$
//
// Author(s)     : Marc Glisse

//| If a compiler does not support rvalue references on this (from C++11)
//| CGAL_CFG_NO_RVALUE_THIS_REFERENCE is set. 

struct A {
	int i;
	A(int j):i(j){}
	A operator-()const&{return -i;}
	A operator-()&&{i=-i;return static_cast<A&&>(*this);}
};

int main()
{
  A a = 1;
  A const& b = a;
  -a;
  -b;
  - -b;
  -A(2);
  return 0;
}
