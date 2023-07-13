// Author(s): Camille Wormser, Pierre Alliez
// An example of an AABB tree constructed with custom point and triangle types.

#include <iostream>
#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>


// custom triangle type with
// three pointers to points
template<class Point_3, class Vector_3, class FT>
struct triangle_trajectory {
    // Vertex positions of triangle
    const Point_3 & pa;
    const Point_3 & pb;
    const Point_3 & pc;

    // Vertex velocities of triangle
    const Vector_3 & va;
    const Vector_3 & vb;
    const Vector_3 & vc;

    // Timestep
    const FT & timestep;

    triangle_trajectory(){}
    triangle_trajectory(
        const Point_3& pa, 
        const Point_3& pb, 
        const Point_3& pc, 
        const Vector_3& va, 
        const Vector_3& vb, 
        const Vector_3& vc, 
        const FT& timestep
    ) : pa{pa}, pb{pb}, pc{pc}, va{va}, vb{vb}, vc{vc}, timestep{timestep} {}
};

// The following primitive provides the conversion facilities between
// the custom triangle and point types and the CGAL ones
template <class K>
struct AABB_tree_triangle_trajectory_primitive {
public:

    // CGAL types returned
    typedef typename K::FT                                                                  FT;
    typedef typename K::Point_3                                                             Point; // CGAL 3D point type
    typedef typename K::Bbox_3                                                              Datum; // CGAL 3D triangle type
    typedef typename K::Vector_3                                                            Vector;
    typedef typename K::Aff_transformation_3                                                Transform;      
    typedef typename std::vector<triangle_trajectory<Point, Vector, FT>>::const_iterator    Iterator;

    // this is the type of data that the queries returns. For this example
    // we imagine that, for some reasons, we do not want to store the iterators
    // of the vector, but raw pointers. This is to show that the Id type
    // does not have to be the same as the one of the input parameter of the
    // constructor.
    typedef const triangle_trajectory<Point, Vector, FT>* Id;


private:
    Id m_pt; // this is what the AABB tree stores internally

public:
    AABB_tree_triangle_trajectory_primitive() {} // default constructor needed

    // the following constructor is the one that receives the iterators from the
    // iterator range given as input to the AABB_tree
    AABB_tree_triangle_trajectory_primitive(Iterator it)
        : m_pt(&(*it)) {}

    const Id& id() const { return m_pt; }

    // on the fly conversion from the internal data to the CGAL types
    Datum datum() const
    {
        Transform transform_a(CGAL::TRANSLATION, va*timestep);
        Transform transform_b(CGAL::TRANSLATION, vb*timestep);
        Transform transform_c(CGAL::TRANSLATION, vc*timestep);

        Point next_a = a.transform(va);
        Point next_b = b.transform(vb);
        Point next_c = c.transform(vc);

        return Datum(
            (std::min)(a.x(), b.x(), c.x(), next_a.x(), next_b.x(), next_c.x()),
            (std::min)(a.y(), b.y(), c.y(), next_a.y(), next_b.y(), next_c.y()),
            (std::min)(a.z(), b.z(), c.z(), next_a.z(), next_b.z(), next_c.z()),
            (std::max)(a.x(), b.x(), c.x(), next_a.x(), next_b.x(), next_c.x()),
            (std::max)(a.y(), b.y(), c.y(), next_a.y(), next_b.y(), next_c.y()),
            (std::max)(a.z(), b.z(), c.z(), next_a.z(), next_b.z(), next_c.z())
        );
    }

    // returns a reference point which must be on the primitive
    Point reference_point() const
    { return convert(m_pt->pa); }
};



// typedef CGAL::AABB_traits<K, My_triangle_primitive> My_AABB_traits;
// typedef CGAL::AABB_tree<My_AABB_traits> Tree;

// int main()
// {
//     My_point a(1.0, 0.0, 0.0);
//     My_point b(0.0, 1.0, 0.0);
//     My_point c(0.0, 0.0, 1.0);
//     My_point d(0.0, 0.0, 0.0);

//     std::vector<My_triangle> triangles;
//     triangles.push_back(My_triangle(&a,&b,&c));
//     triangles.push_back(My_triangle(&a,&b,&d));
//     triangles.push_back(My_triangle(&a,&d,&c));

//     // constructs AABB tree
//     Tree tree(triangles.begin(),triangles.end());

//     // counts #intersections
//     K::Ray_3 ray_query(K::Point_3(1.0, 0.0, 0.0), K::Point_3(0.0, 1.0, 0.0));
//     std::cout << tree.number_of_intersected_primitives(ray_query)
//         << " intersections(s) with ray query" << std::endl;

//     // computes closest point
//     K::Point_3 point_query(2.0, 2.0, 2.0);
//     K::Point_3 closest_point = tree.closest_point(point_query);
//     std::cerr << "closest point is: " << closest_point << std::endl;

//     return EXIT_SUCCESS;
// }
