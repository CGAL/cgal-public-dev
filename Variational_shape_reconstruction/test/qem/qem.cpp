#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Variational_shape_reconstruction/internal/metrics.h>
#include <iostream>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Point_set_3/IO/XYZ.h>
// CGAL
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;


typedef Kernel::FT                  FT;
typedef Kernel::Point_3             Point;
typedef Kernel::Vector_3            Vector;
typedef Kernel::Triangle_3          Triangle;
typedef Kernel::Plane_3             Plane;


/// @brief test the distance between one point of one plane
/// @return tuple function name, resul
std::tuple<std::string,bool> test_qem_one_point_distance_one_plane()
{
    auto qem = QEM_metric();
    auto query = Point(0.,0.,0.);
    auto normal = Vector(0.,0.,1.);
    qem.init_qem_metrics_face(1., query,  normal);
    
    double distance = 2.;
    auto p = Eigen::Vector4d(0.,0.,distance,1.);
    double error_metric = p.transpose() * qem.get_4x4_matrix() * p;
      
    return {__func__ ,error_metric == (distance*distance)}; 
}
/// @brief test the distance between one point and one plane
/// @return tuple function name, resul
std::tuple<std::string,bool>  test_qem_one_point_on_one_plane()
{
    auto qem = QEM_metric();
    auto query = Point(0.,0.,0.);
    auto normal = Vector(0.,0.,1.);
    qem.init_qem_metrics_face(1., query,  normal);
    
    double distance = 0.;
    auto p = Eigen::Vector4d(0.,0.,distance,1.);
    double error_metric = p.transpose() * qem.get_4x4_matrix() * p;
    return {__func__ ,error_metric == 0}; 
}
/// @brief test the distance between one point at the intersection of 2 planes
/// @return tuple function name, resul
std::tuple<std::string,bool>  test_qem_one_point_on_two_plane()
{
    auto qem = QEM_metric();
    auto query = Point(0.,0.,0.);
    auto normal = Vector(0.,0.,1.);
    qem.init_qem_metrics_face(1., query,  normal);

    auto qem2 = QEM_metric();
    normal = Vector(1.,0.,0.);
    qem2.init_qem_metrics_face(1., query,  normal);
    
    qem= qem+qem2;
    auto p = Eigen::Vector4d(0.,0.,0.,1.);
    double error_metric = p.transpose() * qem.get_4x4_matrix() * p;
    return {__func__ ,error_metric == 0}; 
}
/// @brief test the distance between one point and 2 planes
/// @return tuple function name, resul
std::tuple<std::string,bool>  test_qem_one_point_distance_two_plane()
{
    auto qem = QEM_metric();
    auto query = Point(0.,0.,0.);
    auto normal = Vector(0.,0.,1.);
    qem.init_qem_metrics_face(1., query,  normal);

    auto qem2 = QEM_metric();
    normal = Vector(1.,0.,0.);
    qem2.init_qem_metrics_face(1., query,  normal);
    
    qem= qem+qem2;
    double distance_x = 2.;
    double distance_z = 3.;
    auto p = Eigen::Vector4d(distance_x,0.,distance_z,1.);
    double error_metric = p.transpose() * qem.get_4x4_matrix() * p;
    return {__func__ ,error_metric == (distance_x*distance_x + distance_z*distance_z)}; 
}
/// @brief test the distance between one point at the intersection of n planes
/// @return tuple function name, resul
std::tuple<std::string,bool>  test_qem_one_point_on_n_plane()
{
    std::vector<QEM_metric> qem_list;
    int n = 3;
    std::vector<Vector> normals
    {
        Vector(0.,0.,1.),
        Vector(0.,1.,0.),
        Vector(1.,0.,0.)
    };
    auto qem = QEM_metric();
    for(int i = 0 ; i< n ;i++)
    {
        auto qemv = QEM_metric();
        auto query = Point(0.,0.,0.);
        qemv.init_qem_metrics_face(1., query,  normals[i]);
        qem= qem+qemv;
    }

    auto p = Eigen::Vector4d(0.,0.,0.,1.);
    double error_metric = p.transpose() * qem.get_4x4_matrix() * p;
    return {__func__ ,error_metric == 0}; 
}
/// @brief test the distance between one point and n planes
/// @return tuple function name, resul
std::tuple<std::string,bool>  test_qem_one_point_distance_n_plane()
{
    std::vector<QEM_metric> qem_list;
    int n = 4;
    std::vector<Vector> normals
    {
        Vector(0.,0.,1.),
        Vector(0.,1.,0.),
        Vector(1.,0.,0.),
        Vector(-1.,0.,0.)
    };
    auto qem = QEM_metric();
    for(int i = 0 ; i< n ;i++)
    {
        auto qemv = QEM_metric();
        auto query = Point(0.,0.,0.);
        qemv.init_qem_metrics_face(1., query,  normals[i]);
        qem= qem+qemv;
    }
    double distance_y = 1.;
    double distance_x = 4.;
    double distance_z = 3.;
    auto p = Eigen::Vector4d(distance_x,distance_y,distance_z,1.);
    double error_metric = p.transpose() * qem.get_4x4_matrix() * p;
    // distance to p1 + distance to p2 +distance to p3+distance to p4 
    return {__func__ ,error_metric == (distance_x*distance_x+ distance_x*distance_x +distance_y*distance_y +distance_z*distance_z)}; 
}
/// @brief compute sign of double
/// @param x 
/// @return sign of x 1 positif, -1 négatif or 0
double sign(double x) {
    return x < 0 ? -1. : (x > 0 ? 1. : 0.);
}
/// @brief takes a function f and compute monotony over [a,b] using eps for numerical steps
/// @param f function
/// @param a lower bound
/// @param b upper bound
/// @param eps epsilon for each step
/// @return 1 increasing, -1 decreasing, 0 constant and 2 if not monotonic
int monotonic(const std::function< double(double) >& f,
                      double a, double b, double eps) {
    double x = a, y1 = f(x), y2 = f(x + eps);
    double d = y2 - y1;
    double s = sign(d);
    while( x < b ) {
        x += eps;
        y1 = y2;
        y2 = f(x + eps);
        d = y2 - y1;
        if( s == 0. ) {
            s = sign(d);
        }
        if( s * d < 0 ) {
            return 2;
        }
    }
    return s > 0 ? 1 :
          (s < 0 ? -1 :
                   0);
}

double test_qem_one_point_one_plane_with_distance(double distance)
{
    auto qem = QEM_metric();
    auto query = Point(0.,0.,0.);
    auto normal = Vector(0.,0.,1.);
    qem.init_qem_metrics_face(1., query,  normal);
    
    auto p = Eigen::Vector4d(0.,0.,distance,1.);
    double error_metric = p.transpose() * qem.get_4x4_matrix() * p;
    return error_metric;
}
double test_qem_one_point_two_planes_with_distance(double distance)
{
    auto query = Point(0.,0.,0.);

    auto qem = QEM_metric();
    auto normal = Vector(0.,0.,1.);
    qem.init_qem_metrics_face(1., query,  normal);

    auto qem2 = QEM_metric();
    normal = Vector(1.,0.,0.);
    qem2.init_qem_metrics_face(1., query,  normal);
    
    qem= qem+qem2;

    auto p = Eigen::Vector4d(distance,distance,distance,1.);
    double error_metric = p.transpose() * qem.get_4x4_matrix() * p;
    return error_metric;
}
double test_qem_one_point_n_planes_with_distance(double distance)
{
    std::vector<QEM_metric> qem_list;
    int n = 4;
    std::vector<Vector> normals
    {
        Vector(0.,0.,1.),
        Vector(0.,1.,0.),
        Vector(1.,0.,0.),
        Vector(-1.,0.,0.)
    };
    auto qem = QEM_metric();
    for(int i = 0 ; i< n ;i++)
    {
        auto qemv = QEM_metric();
        auto query = Point(0.,0.,0.);
        qemv.init_qem_metrics_face(1., query,  normals[i]);
        qem= qem+qemv;
    }
    auto p = Eigen::Vector4d(distance,distance,distance,1.);
    double error_metric = p.transpose() * qem.get_4x4_matrix() * p;
    return error_metric;
}
/// @brief test the qem when the intersection of the planes is not one point
/// @return tuple name, boolean value if the minimun qem is not 0
std::tuple<std::string,bool>  test_qem_one_point_n_planes_no_point_intersection_with_distance()
{
    std::vector<QEM_metric> qem_list;
    int n = 4;
    double distance = 2.;
    std::vector<Vector> normals
    {
        Vector(0.,0.,1.),
        Vector(0.,1.,0.),
        Vector(1.,0.,0.),
        Vector(-1.,0.,0.)
    };
    std::vector<Point> points
    {
        Point(0.,0.,-distance),
        Point(0.,distance,0.),
        Point(distance,0.,0.),
        Point(-distance,0.,0.)
    };
    // We set a box of 4 planes, the center of the box is at 0. 0. 0. 
    // we will have 4 (if distance set to 1)
    // as the sum of the minimun distance squared from all the planes 
    auto qem = QEM_metric();
    for(int i = 0 ; i< n ;i++)
    {
        auto qemv = QEM_metric();
        auto query = points[i];
        qemv.init_qem_metrics_face(1., query,  normals[i]);
        qem= qem+qemv;
    }
    auto p = Eigen::Vector4d(0.,0.,0.,1.);
    double error_metric = p.transpose() * qem.get_4x4_matrix() * p;
     return {__func__ ,error_metric == (n*(distance*distance))}; 
}
/// @brief Computes the qem for 6 triangles as a star
/// @param distance from the center vertex of the star
/// @return qem
double test_qem_one_point_star_triangles_with_distance(double distance)
{
    std::vector<QEM_metric> qem_list;
    int n = 6;
    std::vector<Point> points;
    std::vector<Vector> normals;
    std::vector<Triangle> triangles;
    double two_pi = 2*atan(1)*4;
    auto query = Point(0.,0.,0.5);
    for(int i = 0 ; i< n ;i++)
    {
        double a =i*(two_pi/n);
        auto x = cos(a);
        auto y = sin(a);
        points.push_back(Point(x,y,0.));
        if(i>0)
        {
            auto t = Triangle(points[i-1],points[i],query);
            Vector normal = CGAL::normal(points[i-1],
                                         points[i],
                                         query);
            normals.push_back(normal);
            triangles.push_back(t);
        }
    }
    auto t = Triangle(points.back(),points[0],query);
    Vector normal = CGAL::normal(points.back(),points[0],query);
    normals.push_back(normal);
    triangles.push_back(t);
    
    auto qem = QEM_metric();
    for(int i = 0 ; i< n ;i++)
    {
        auto qemv = QEM_metric();
        qemv.init_qem_metrics_face(sqrt(triangles[i].squared_area()), query,  normals[i]);
        qem= qem+qemv;
    }
    auto p = Eigen::Vector4d(0.,0.,distance+0.5,1.);
    double error_metric = p.transpose() * qem.get_4x4_matrix() * p;
    return error_metric;
}
/// @brief Test if the center vertex intersecting all triangles has a qem equals to 0
/// @return tuple function name, result
std::tuple<std::string,bool>  test_qem_one_point_star_triangles()
{
    return {__func__ ,test_qem_one_point_star_triangles_with_distance(0.)==0.}; 
}
/// @brief Test if the qem from a point at a distance from the center vertex in the star is monotonic
/// @return tuple function name, result
std::tuple<std::string,bool>  test_qem_one_point_star_triangles_monotonic()
{
    int type = monotonic(&test_qem_one_point_star_triangles_with_distance,0., 1., 10e-2);
    return {__func__ ,type==1}; 
}
/// @brief test the distance between one point and n planes
/// @return tuple function name, result
std::tuple<std::string,bool>  test_qem_one_point_on_one_plane_monotonic()
{
    int type = monotonic(&test_qem_one_point_one_plane_with_distance,0., 1., 10e-2);
    return {__func__ ,type==1}; 
}
/// @brief test the distance between one point and n planes
/// @return tuple function name, resul
std::tuple<std::string,bool>  test_qem_one_point_on_two_planes_monotonic()
{
    int type = monotonic(&test_qem_one_point_two_planes_with_distance,0., 1., 10e-2);
    return {__func__ ,type==1}; 
}
/// @brief test the distance between one point and n planes
/// @return tuple function name, resul
std::tuple<std::string,bool>  test_qem_one_point_on_n_planes_monotonic()
{
    int type = monotonic(&test_qem_one_point_n_planes_with_distance,0., 1., 10e-2);
    return {__func__ ,type==1}; 
}
/// @brief function that test if a function returns true, print its name in green otherwise in red
/// @param myfunction the function to test
void test(std::function<std::tuple<std::string,bool>(void)> myfunction)
{
    auto tuple = myfunction();
    if(std::get<1>(tuple))
    {
        std::cout<<"\x1B[32m"<<std::get<0>(tuple)<<" passed\033[0m\n";
    }
    else
    {
        std::cout<<"\x1B[31m"<<std::get<0>(tuple)<<" not passed\033[0m\n";

    }
}

int main(int argc, char** argv)
{	
    std::cout<<"Test egality of distance (0 or sum of squared distance) : \n";
    test(test_qem_one_point_on_one_plane);
    test(test_qem_one_point_distance_one_plane);
    std::cout<<"\n";
    test(test_qem_one_point_on_two_plane);
    test(test_qem_one_point_distance_two_plane);
    std::cout<<"\n";
    test(test_qem_one_point_on_n_plane);
    test(test_qem_one_point_distance_n_plane);
    std::cout<<"\n";
    test(test_qem_one_point_n_planes_no_point_intersection_with_distance);
    test(test_qem_one_point_star_triangles);
    std::cout<<"Test monotony of the QEM: \n";
    test(test_qem_one_point_on_one_plane_monotonic);
    test(test_qem_one_point_on_two_planes_monotonic);
    test(test_qem_one_point_on_n_planes_monotonic);

    test(test_qem_one_point_star_triangles_monotonic);
	return 0;
}


