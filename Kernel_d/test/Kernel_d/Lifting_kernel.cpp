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
#include <CGAL/Lifting_kernel_d.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/assertions.h>

int main(){
        typedef CGAL::Cartesian_d<double>                       Base;
        typedef CGAL::Lifting_kernel_d<Base>                    K;
        typedef K::Point_d                                      Point;
        typedef CGAL::Convex_hull_d<K>                          CH;
        typedef K::Orientation_d                                Ori;
        typedef K::Volume_d                                     Vol;

        Point p(5,-1),q(3,-4),r(1,2),s(4,2),t(3,0);
        CH ch1(2);
        ch1.insert(p);
        ch1.insert(q);
        ch1.insert(r);
        ch1.insert(s);
        ch1.insert(t);

        CGAL_assertion(ch1.is_valid());

        K::set_lifting(q,4);
        K mykernel;
        std::vector<Point> points;
        points.push_back(p);
        points.push_back(q);
        points.push_back(r);
        CH ch2(2,mykernel);
        ch2.insert(points.begin(),points.end());

        CGAL_assertion(ch2.is_valid());

        // Test predicates with d+1 points.
        CGAL_assertion(Ori()(points.begin(),points.end())==
                       CGAL::COUNTERCLOCKWISE);
        CGAL_assertion(Vol()(points.begin(),points.end())==14);

        // Test predicates with d points.
        points.pop_back();
        CGAL_assertion(Vol()(points.begin(),points.end())==2);
        CGAL_assertion(Ori()(points.begin(),points.end())==
                       CGAL::COUNTERCLOCKWISE);

        return 0;
}
