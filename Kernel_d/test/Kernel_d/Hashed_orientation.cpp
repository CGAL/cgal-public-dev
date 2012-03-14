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

#include <CGAL/Cartesian_d.h>
#include <CGAL/Lifting_kernel_d.h>
#include <CGAL/assertions.h>

int main(){
        typedef CGAL::Cartesian_d<double>                       Base;
        typedef CGAL::Lifting_kernel_d<Base>                    K;
        typedef K::Point_d                                      Point;
        Point p(5,-1),q(3,-4),r(1,2);
        Point s(q);
        CGAL_assertion(p.index()+1==q.index());
        CGAL_assertion(q.index()+1==r.index());
        CGAL_assertion(s.index()==q.index());
        std::vector<Point> points;
        points.push_back(r);
        points.push_back(s);
        points.push_back(p);
        CGAL::Orientation o1=
                K().orientation_d_object()(points.begin(),points.end());
        K::set_lifting(s,4);
        points.clear();
        points.push_back(r);
        points.push_back(s);
        points.push_back(p);
        CGAL::Orientation o2=
                K().orientation_d_object()(points.begin(),points.end());
        CGAL_assertion(o1!=o2);
        return 0;
}
