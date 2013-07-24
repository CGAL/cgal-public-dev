// Author(s) : 

#include <iostream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

typedef CGAL::Simple_cartesian<double> K;

//2D geometric data types 
typedef K::FT FT;
typedef K::Ray_2 Ray;
typedef K::Line_2 Line;
typedef K::Point_2 Point;
typedef K::Triangle_2 Triangle;

//AABB tree is instantiated with the 2D geometric data, notice that there is no need to explicitly define the dimension. 
typedef std::list<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

int main()
{
    Point a(0,0);
	Point b(1,-2);
	Point c(2,-1);
	Point d(2,1);
	Point e(1,2);
	Point f(-1,-2);
	Point g(-2,1);
	Point h(-2,-1);
	Point i(-1,-2);

    std::list<Triangle> triangles;
    list.push_back(Triangle(a,b,c));
	list.push_back(Triangle(a,d,e));
	list.push_back(Triangle(a,f,g));
	list.push_back(Triangle(a,h,i));

    // constructs 2D AABB tree
    Tree tree(triangles.begin(),triangles.end());

    // counts #intersections
    Ray ray_query(g,d);
    std::cout << tree.number_of_intersected_primitives(ray_query)
        << " intersections(s) with ray query" << std::endl;

    // compute closest point and squared distance
    Point point_query(1.0, 1.0);
    Point closest_point = tree.closest_point(point_query);
    std::cerr << "closest point is: " << closest_point << std::endl;
    FT sqd = tree.squared_distance(point_query);
    std::cout << "squared distance: " << sqd << std::endl;

    return EXIT_SUCCESS;
}
