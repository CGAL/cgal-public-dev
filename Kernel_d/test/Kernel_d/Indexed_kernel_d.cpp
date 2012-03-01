// Copyright 2011-2012 National and Kapodistrian University of Athens,
// Greece.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://penarand@scm.gforge.inria.fr/svn/cgal/branches/features/Kernel_d-new_models-penarand_vfisikop/Kernel_d/test/Kernel_d/Indexed_kernel_d.cpp $
// $Id: Indexed_kernel_d.cpp 67920 2012-03-01 16:00:51Z penarand $
//
//
// Authors:     Vissarion Fisikopoulos <vissarion@di.uoa.gr>
//              Luis Peñaranda <luis.penaranda@gmx.com>

// This test program only creates some points and makes sure that:
// -their indices are generated sequentially when they are constructed from
// scratch, and
// -their indices are copied by the copy constructor.

#include <CGAL/Cartesian_d.h>
#include <CGAL/Indexed_kernel_d.h>
#include <CGAL/assertions.h>
#include <vector>

int main(){
        typedef CGAL::Cartesian_d<double>                       Base;
        typedef CGAL::Indexed_point_kernel_d<Base>              K;
        typedef K::Point_d                                      Point;
        Point p(1,2);
        std::vector<double> vq;
        vq.push_back(3);
        vq.push_back(4);
        vq.push_back(5);
        vq.push_back(6);
        vq.push_back(7);
        Point q(vq.size(),vq.begin(),vq.end());
        Point r;
        Point s(q);
        CGAL_assertion(p.index()+1==q.index());
        CGAL_assertion(q.index()+1==r.index());
        CGAL_assertion(s.index()==q.index());
        return 0;
}
