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
// $URL$
// $Id$
//
//
// Authors:     Vissarion Fisikopoulos <vissarion@di.uoa.gr>
//              Luis Peñaranda <luis.penaranda@gmx.com>

#include <CGAL/Cartesian_d.h>
#include <CGAL/Indexed_kernel_d.h>
#include <CGAL/Lifted_kernel_d.h>
#include <CGAL/assertions.h>

int main(){
        typedef CGAL::Cartesian_d<double>                       Base;
        typedef CGAL::Indexed_point_kernel_d<Base>              Indexed;
        typedef CGAL::Lifted_kernel_d<Indexed>                  K;
        typedef K::Point_d                                      Point;
        typedef CGAL::Convex_hull_d<K>                          CH;

        Point p(5,-1),q(3,-4),r(1,2);
        CH ch1(2);
        ch1.insert(p);
        ch1.insert(q);
        ch1.insert(r);

        K::set_lifting(q,4);
        K mykernel;
        std::vector<Point> points;
        points.push_back(p);
        points.push_back(q);
        points.push_back(r);
        CH ch2(2,mykernel);
        ch2.insert(points.begin(),points.end());

        return 0;
}
