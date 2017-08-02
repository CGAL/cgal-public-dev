/*
 * chooseNewAxis.h
 *
 *  Created on: Jul 31, 2017
 *      Author: shahar
 */

#ifndef EXAMPLES_SET_MOVABLE_SEPARABILITY_3_CHOOSENEWAXIS_H_
#define EXAMPLES_SET_MOVABLE_SEPARABILITY_3_CHOOSENEWAXIS_H_
#include <math.h>       /* sqrt */
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Linear_algebraHd.h>
#include <CGAL/min_quadrilateral_2.h>
#include <CGAL/convex_hull_2.h>

template<typename Kernel>
class AxisConvertor
{

  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef CGAL::Aff_transformation_3<Kernel> transformation_3;
  transformation_3 m_toNew;
  transformation_3 m_toOriginal;
public:
  AxisConvertor( Vector_3 new_x,Vector_3 new_y, Vector_3 new_z)
  {

    m_toOriginal= CGAL::Aff_transformation_3<Kernel>(
	 new_x.x(),new_y.x(),new_z.x(),
	 new_x.y(),new_y.y(),new_z.y(),
	 new_x.z(),new_y.z(),new_z.z()
    );

    typename Kernel::FT determinant =    +new_x.x()*(new_y.y()*new_z.z()-new_z.y()*new_y.z())
                            -new_x.y()*(new_y.x()*new_z.z()-new_y.z()*new_z.x())
                            +new_x.z()*(new_y.x()*new_z.y()-new_y.y()*new_z.x());
    typename Kernel::FT invdet = 1/determinant;
    m_toNew = CGAL::Aff_transformation_3<Kernel>(
	 (new_y.y()*new_z.z()-new_z.y()*new_y.z())*invdet,-(new_y.x()*new_z.z()-new_y.z()*new_z.x())*invdet,(new_y.x()*new_z.y()-new_z.x()*new_y.y())*invdet,
	 -(new_x.y()*new_z.z()-new_x.z()*new_z.y())*invdet,(new_x.x()*new_z.z()-new_x.z()*new_z.x())*invdet,-(new_x.x()*new_z.y()-new_z.x()*new_x.y())*invdet,
	 (new_x.y()*new_y.z()-new_x.z()*new_y.y())*invdet,-(new_x.x()*new_y.z()-new_y.x()*new_x.z())*invdet,(new_x.x()*new_y.y()-new_y.x()*new_x.y())*invdet
    );

  }
  Point_3 toNewCords(Point_3 p)
  {
    return p.transform(m_toNew);

  }

  Point_3 toOriginalCords(Point_3 p)
  {
    return p.transform(m_toOriginal);
}
};

template<typename Polyhedron,typename Kernel>
AxisConvertor<Kernel> chooseNewAxis (const Polyhedron&p,const typename Polyhedron::Face_iterator &fit)
  {
  typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Point_3 Point_3;
    typedef typename Kernel::Vector_3 Vector_3;
    typedef CGAL::Polygon_2<Kernel>   Polygon;

    typename Polyhedron::Halfedge_const_handle h=fit->halfedge();
    Vector_3 new_x;
    Vector_3 new_y;
    Vector_3 new_z;
    Point_3 p1 = h->vertex()->point();
    h=h->next();
    Point_3 p2 = h->vertex()->point();
    h=h->next();
    Point_3 p3 = h->vertex()->point();
    new_z = CGAL::normal(p1,p2,p3);
    new_x = p1-p3;
     Point_3 a = p1+new_z;
    new_y = -CGAL::normal(p1,p3,a);
    std::cout<<new_x<<std::endl<<new_y<<std::endl<<new_z<<std::endl<<std::endl;
    double xlength =sqrt(CGAL::to_double(new_x.squared_length()));
    new_x =Vector_3(new_x.x() /xlength,new_x.y() /xlength,new_x.z() /xlength);

    double ylength =sqrt(CGAL::to_double(new_y.squared_length()));
    new_y =Vector_3(new_y.x() /ylength,new_y.y() /ylength,new_y.z() /ylength);

    double zlength =sqrt(CGAL::to_double(new_z.squared_length()));
    new_z =Vector_3(new_z.x() /zlength,new_z.y() /zlength,new_z.z() /zlength);
    AxisConvertor<Kernel> randomAxis(new_x,new_y,new_z);
    std::vector< Point_2>points;
    points.reserve(p.size_of_vertices());
   typename Polyhedron::Point_const_iterator pit = p.points_begin();

         for(;pit!= p.points_end();++pit)
   	{
   	  Point_3 randomAxisPit = randomAxis.toNewCords(*pit);
   	  points.push_back(Point_2(randomAxisPit.x(),randomAxisPit.y()));
   	    std::cout << "randomAxisPit="<<randomAxisPit << std::endl;

   	}

    Polygon convex_polygon;

    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(convex_polygon));
    std::cout << convex_polygon << std::endl;

    Polygon p_m;
    CGAL::min_rectangle_2(convex_polygon.vertices_begin(), convex_polygon.vertices_end() ,std::back_inserter(p_m));
    std::cout << p_m << std::endl;
    typename Polygon::Vertex_const_iterator vci =p_m.vertices_begin();
    Point_2 optimalR1= *(vci++);
    Point_2 optimalR2= *(vci++);
    Point_2 optimalR3= *(vci);
    Point_3 orgOptimalR1= randomAxis.toOriginalCords(Point_3(optimalR1.x(),optimalR1.y(),0));
    Point_3 orgOptimalR2= randomAxis.toOriginalCords(Point_3(optimalR2.x(),optimalR2.y(),0));
    Point_3 orgOptimalR3= randomAxis.toOriginalCords(Point_3(optimalR3.x(),optimalR3.y(),0));

    new_x  = orgOptimalR1-orgOptimalR2;
    xlength =sqrt(CGAL::to_double(new_x.squared_length()));
    new_x =Vector_3(new_x.x() /xlength,new_x.y() /xlength,new_x.z() /xlength);

    new_y  = orgOptimalR2-orgOptimalR3;
    ylength =sqrt(CGAL::to_double(new_y.squared_length()));
    new_y =Vector_3(new_y.x() /ylength,new_y.y() /ylength,new_y.z() /ylength);
    std::cout << "x="<<new_x<< " y="<<new_y<< " z="<<new_z << std::endl;

    return AxisConvertor<Kernel>(new_x,new_y,new_z);
}


#endif /* EXAMPLES_SET_MOVABLE_SEPARABILITY_3_CHOOSENEWAXIS_H_ */
