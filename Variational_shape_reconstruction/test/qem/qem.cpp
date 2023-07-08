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

typedef Kernel::Point_2             Point_2;
typedef Kernel::Segment_2           Segment_2;

typedef CGAL::First_of_pair_property_map<std::pair<Point, Vector>>                     Point_map;
typedef CGAL::Second_of_pair_property_map<std::pair<Point, Vector>>                    Normal_map;
typedef CGAL::Point_set_3< Point, Vector > Pointset;

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
    //std::cout<<"metric : "<<error_metric<<"\n";
      
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
    test(test_qem_one_point_on_one_plane);
    test(test_qem_one_point_distance_one_plane);
    std::cout<<"\n";
    test(test_qem_one_point_on_two_plane);
    test(test_qem_one_point_distance_two_plane);
    std::cout<<"\n";
    test(test_qem_one_point_on_n_plane);
    test(test_qem_one_point_distance_n_plane);
    
	return 0;
}


