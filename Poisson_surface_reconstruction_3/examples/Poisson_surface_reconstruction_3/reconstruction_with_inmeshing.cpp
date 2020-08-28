//includes in-meshing
#include <core/util.h>
#include <core/output.h>
#include <core/mesh_completion.h>
#include <core/adjust_geometry.h>
#include <core/dualcontouring/connectivity.h>
#include <core/dualcontouring/dual_contouring.h>
#include <poisson_recon/Bridge.h>

//includes CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/Polygon_mesh_processing/distance.h>

#include <utility>


// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Point_with_normal;
typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
typedef Kernel::Sphere_3 Sphere;
typedef std::vector<Point_with_normal> PointList;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

std::vector<oriented_point> read_xyz_points(const std::string& file)
{
  std::vector<oriented_point> points;

  std::ifstream is(file);
  float x, y, z, nx, ny, nz;

  while(!is.eof())
  {
    is >> x >> y >> z >> nx >> ny >> nz;
    points.emplace_back(Eigen::Vector3f(x, y, z), Eigen::Vector3f(nx, ny, nz));
  }

  return points;
}

/// Return an octree that correspond to the 'space nodes' of the reference Poisson reconstruction
template<typename T>
sorted_octree<T> space_nodes(const poisson_outputs& poisson, const T& def = T())
{
  std::vector<octree_id> snodes;
  poisson.tree.for_each([&](const octree_id& node, const poisson_outputs::node&)
                        {
                          if(node == poisson.FEM_root || poisson.FEM_root.is_ancestor_of(node))
                            snodes.emplace_back(node.relative_pos(poisson.FEM_root));
                        });
  return sorted_octree<T>(snodes, def);
}

struct Implicit_function
{
  poisson_outputs poisson;
  Eigen::Matrix4f transform;
  Eigen::Matrix4f back_transform;

  Implicit_function(poisson_outputs  poisson_, Eigen::Matrix4f  transform_, Eigen::Matrix4f  back_transform_):
      poisson(std::move(poisson_)),
      transform(std::move(transform_)),
      back_transform(std::move(back_transform_))
  {}

  FT operator()(const Point& p) const
  {
    Eigen::Vector3f v1 (p.x(), p.y(), p.z());
    Eigen::Vector3f v2 = from_hpos<float>(transform * to_hpos(v1));
    float x = poisson.implicit_value()(v2);
    if(x == -1) // exterior
    {
      return -50;
    }
    return -50 * x;
  }
};

// constructs a point inside the solid starting from the i-th point and following the opposite direction of the normal
Point incursion(const PointList& points, const Implicit_function& f, const FT& average_spacing, unsigned i)
{
  Point p = points[i].first;
  Vector n = points[i].second;
  p -= 0.5 * average_spacing * n;
  return p;
}

// constructs a point inside the solid
Point make_inner_point(const PointList& points, const Implicit_function& f, const FT& average_spacing)
{
  float threshold = 10;
  unsigned i = 0;
  Point p = incursion(points, f, average_spacing, i);

  // we make incursion until we find a point whose implicit fonction is greater than the threshold
  while(f(p) < threshold)
  {
    ++i;
    p = incursion(points, f, average_spacing, i);
  }

  return p;
}

FT make_bounding_sphere_radius(PointList points)
{
  FT squared_max_distance = 0;
  for(unsigned i = 0; i < points.size(); ++i)
  {
    for(unsigned j = i + 1; j < points.size(); ++j)
    {
      Vector u (points[i].first.x(), points[i].first.y(), points[i].first.z());
      Vector v (points[j].first.x(), points[j].first.y(), points[j].first.z());
      FT challenger = (u - v).squared_length();
      if(squared_max_distance < challenger)
      {
        squared_max_distance = challenger;
      }
    }
  }

  return 2 * squared_max_distance;
}

void in_meshing_reconstruction(const std::string& input_file, const std::string& output_file)
{
  // read and setup
  point_mesh mesh;
  mesh.vertices = read_xyz_points(input_file);
  // the poisson from in meshing will be performed in a unit box, we need to be able to switch between both
  Eigen::Matrix4f back_transform = mesh.transform_to_unit(1.25f);
  Eigen::Matrix4f transform = back_transform.inverse();

  // solve the poisson problem
  completed_mesh out_mesh = std::move(mesh);
  unsigned max_depth = 7;
  poisson_outputs poisson = run_poisson({ out_mesh.vertices, max_depth, poisson_boundary_type::dirichlet, poisson_output_tree });
  upscale_1(poisson.tree);

  // setup the inplicit function
  Implicit_function implicit_function (poisson, transform, back_transform);


  // Poisson options
  FT sm_angle = 20.0; // Min triangle angle in degrees.
  FT sm_radius = 30; // Max triangle size w.r.t. point set average spacing.
  FT sm_distance = 0.375; // Surface Approximation error w.r.t. point set average spacing.

  // Computes average spacing
  PointList points;
  for(auto& pwn : out_mesh.vertices)
  {
    pwn.position = from_hpos<float>(back_transform * to_hpos(pwn.position));
    pwn.normal = from_hpos<float>(back_transform * to_hpos(pwn.normal));
    Point p (pwn.position[0], pwn.position[1], pwn.position[2]);
    Vector n (pwn.normal[0], pwn.normal[1], pwn.normal[2]);
    points.emplace_back(p, n);
  }

  FT average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
      (points, 6 /* knn = 1 ring */,
       CGAL::parameters::point_map (Point_map()));

  // Gets one point inside the implicit surface
  // and computes implicit function bounding sphere radius.
  Point inner_point = make_inner_point(points, implicit_function, average_spacing);//= function.get_inner_point();
  Sphere bsphere (inner_point, 9);//= function.bounding_sphere();
  FT radius = std::sqrt(bsphere.squared_radius());

  // Defines the implicit surface: requires defining a
  // conservative bounding sphere centered at inner point.
  FT sm_sphere_radius = 5.0 * radius;
  FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
  Surface_3 surface(implicit_function,
                    Sphere(inner_point,sm_sphere_radius*sm_sphere_radius),
                    sm_dichotomy_error/sm_sphere_radius);

  // Defines surface mesh generation criteria
  CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
                                                      sm_radius*average_spacing,  // Max triangle size
                                                      sm_distance*average_spacing); // Approximation error

  // Generates surface mesh with manifold option
  STr tr; // 3D Delaunay triangulation for surface mesh generation
  C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
  CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
                          surface,                              // implicit surface
                          criteria,                             // meshing criteria
                          CGAL::Manifold_with_boundary_tag());  // require manifold mesh

  // saves reconstructed surface mesh
  std::ofstream out(output_file);
  Polyhedron output_mesh;
  CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, output_mesh);
  out << output_mesh;
  out.close();
}

int main(int argc, char** argv)
{
  in_meshing_reconstruction(argv[1], argv[2]);
}
