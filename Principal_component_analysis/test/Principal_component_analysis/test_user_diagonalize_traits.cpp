#include <CGAL/Simple_cartesian.h>

#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>

typedef CGAL::Simple_cartesian<double> Kernel;

// Dummy structure just to test compilation
// with something different than Eigen_diagonalize_traits
struct User_diagonalize_traits
{
  typedef std::vector<double> Vector;
  typedef std::vector<double> Matrix;
  typedef std::vector<double> Covariance_matrix;
  static bool diagonalize_selfadjoint_covariance_matrix(const Covariance_matrix&, Vector&)
  { return true; }
  static bool diagonalize_selfadjoint_covariance_matrix(const Covariance_matrix&,
                                                        Vector&, Matrix&)
  { return true; }
  static bool extract_largest_eigenvector_of_covariance_matrix(const Covariance_matrix&,
                                                               Vector&)
  { return true; }
};

// Generate default objects so that the test compiles *AND* runs fine

Kernel::Point_2 default_object(const Kernel::Point_2&)
{ return Kernel::Point_2(0,0); }
Kernel::Segment_2 default_object(const Kernel::Segment_2&)
{ return Kernel::Segment_2(Kernel::Point_2(0,0), Kernel::Point_2(0,1)); }
Kernel::Circle_2 default_object(const Kernel::Circle_2&)
{ return Kernel::Circle_2(Kernel::Point_2(0,0), Kernel::FT(1.0)); }
Kernel::Triangle_2 default_object(const Kernel::Triangle_2&)
{ return Kernel::Triangle_2(Kernel::Point_2(0,0), Kernel::Point_2(0,1), Kernel::Point_2(1,0)); }
Kernel::Iso_rectangle_2 default_object(const Kernel::Iso_rectangle_2&)
{ return Kernel::Iso_rectangle_2(Kernel::Point_2(0,0), Kernel::Point_2(1,1)); }
Kernel::Point_3 default_object(const Kernel::Point_3&)
{ return Kernel::Point_3(0,0,0); }
Kernel::Segment_3 default_object(const Kernel::Segment_3&)
{ return Kernel::Segment_3(Kernel::Point_3(0,0,0), Kernel::Point_3(0,0,1)); }
Kernel::Sphere_3 default_object(const Kernel::Sphere_3&)
{ return Kernel::Sphere_3(Kernel::Point_3(0,0,0), Kernel::FT(1.0)); }
Kernel::Triangle_3 default_object(const Kernel::Triangle_3&)
{ return Kernel::Triangle_3(Kernel::Point_3(0,0,0), Kernel::Point_3(0,0,1), Kernel::Point_3(0,1,0)); }
Kernel::Tetrahedron_3 default_object(const Kernel::Tetrahedron_3&)
{ return Kernel::Tetrahedron_3(Kernel::Point_3(0,0,0), Kernel::Point_3(0,0,1),
                                   Kernel::Point_3(0,1,0), Kernel::Point_3(1,0,0)); }
Kernel::Iso_cuboid_3 default_object(const Kernel::Iso_cuboid_3&)
{ return Kernel::Iso_cuboid_3(Kernel::Point_3(0,0,0), Kernel::Point_3(1,1,1)); }

template <typename Object, int dim>
void test_2d()
{
  std::array<Object, 1> dummy = { default_object(Object()) };
  Kernel::Line_2 line;
  Kernel::Point_2 centroid;
  CGAL::linear_least_squares_fitting_2 (dummy.begin(), dummy.end(), line, centroid,
                                        CGAL::Dimension_tag<dim>(), Kernel(),
                                        User_diagonalize_traits());
}

template <typename Object, int dim>
void test_3d()
{
  std::array<Object, 1> dummy = { default_object(Object()) };
  Kernel::Line_3 line;
  Kernel::Plane_3 plane;
  Kernel::Point_3 centroid;
  CGAL::linear_least_squares_fitting_3 (dummy.begin(), dummy.end(), line, centroid,
                                        CGAL::Dimension_tag<dim>(), Kernel(),
                                        User_diagonalize_traits());
  CGAL::linear_least_squares_fitting_3 (dummy.begin(), dummy.end(), plane, centroid,
                                        CGAL::Dimension_tag<dim>(), Kernel(),
                                        User_diagonalize_traits());
}


int main()
{
  test_2d<Kernel::Point_2, 0>();
  test_2d<Kernel::Segment_2, 1>();
  test_2d<Kernel::Segment_2, 0>();
  test_2d<Kernel::Circle_2, 2>();
  test_2d<Kernel::Circle_2, 1>();
  test_2d<Kernel::Triangle_2, 2>();
  test_2d<Kernel::Triangle_2, 1>();
  test_2d<Kernel::Triangle_2, 0>();
  test_2d<Kernel::Iso_rectangle_2, 2>();
  test_2d<Kernel::Iso_rectangle_2, 1>();
  test_2d<Kernel::Iso_rectangle_2, 0>();
  test_3d<Kernel::Point_3, 0>();
  test_3d<Kernel::Segment_3, 1>();
  test_3d<Kernel::Segment_3, 0>();
  test_3d<Kernel::Sphere_3, 3>();
  test_3d<Kernel::Sphere_3, 2>();
  test_3d<Kernel::Triangle_3, 2>();
  test_3d<Kernel::Triangle_3, 1>();
  test_3d<Kernel::Triangle_3, 0>();
  test_3d<Kernel::Tetrahedron_3, 3>();
  test_3d<Kernel::Tetrahedron_3, 2>();
  test_3d<Kernel::Tetrahedron_3, 1>();
  test_3d<Kernel::Tetrahedron_3, 0>();
  test_3d<Kernel::Iso_cuboid_3, 3>();
  test_3d<Kernel::Iso_cuboid_3, 2>();
  test_3d<Kernel::Iso_cuboid_3, 1>();
  test_3d<Kernel::Iso_cuboid_3, 0>();

  return EXIT_SUCCESS;
}
