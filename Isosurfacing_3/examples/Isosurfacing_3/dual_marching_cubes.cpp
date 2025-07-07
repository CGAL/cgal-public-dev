#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>
#include <CGAL/Isosurfacing_3/Value_function_3.h>
#include <CGAL/Isosurfacing_3/Gradient_function_3.h>
#include <CGAL/Isosurfacing_3/internal/dual_marching_cube_functors.h>  

#include <CGAL/Bbox_3.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <CGAL/IO/write_off_points.h>


#include <cmath>
#include <iostream>
#include <vector>
#include <cmath>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = Kernel::FT;
using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
using Values = CGAL::Isosurfacing::Value_function_3<Grid>;
using Gradients = CGAL::Isosurfacing::Gradient_function_3<Grid>;
using Domain = CGAL::Isosurfacing::Dual_contouring_domain_3<Grid, Values, Gradients>;

using Pointrange = std::vector<Point>;
using Polygonrange = std::vector<std::vector<std::size_t> >;

using Mesh = CGAL::Surface_mesh<Point>;

// example function

// auto Sphere_value = [](const Point& p) {
//     const FT x = p.x(), y = p.y(), z = p.z();
//     return x*x + y*y + z*z - 1; 
// };

// auto Sphere_gradient = [](const Point& p) {
//     const FT x = p.x(), y = p.y(), z = p.z();
//     const FT gx = 2.0 * x;
//     const FT gy = 2.0 * y;
//     const FT gz = 2.0 * z;
//     Vector g(gx, gy, gz);
//     return g / std::sqrt(gx*gx + gy*gy + gz*gz);
// };

// auto Kusner_Schmitt_value = [](const Point& point)
// {
//   const FT x = point.x(), y = point.y(), z = point.z();
//   return x*x*y*y*z*z
//          + 3*x*x*y*y + 3*x*x*z*z + 9*x*x
//          + 3*y*y*z*z + 9*y*y + 9*z*z
//          - 32*x*y*z - 5.0;
// };

// auto Kusner_Schmitt_gradient = [](const Point& p) 
// {
//   const FT x = p.x(), y = p.y(), z = p.z();

//   const FT gx = 2*x*y*y*z*z + 6*x*y*y + 6*x*z*z + 18*x - 32*y*z;
//   const FT gy = 2*x*x*y*z*z + 6*x*x*y + 6*y*z*z + 18*y - 32*x*z;
//   const FT gz = 2*x*x*y*y*z + 6*x*x*z + 6*y*y*z + 18*z - 32*x*y;

//   Vector g(gx, gy, gz);
//   return g / std::sqrt(gx*gx + gy*gy + gz*gz);
// };

// auto S22_value = [](const Point& point)
// {
//   const FT x = point.x(), y = point.y(), z = point.z();
//   return 4*x*x*z + y*y*y + 4*y*x;
// };

// auto S22_gradient = [](const Point& point)
// {
//   const FT x = point.x(), y = point.y(), z = point.z();
  
//   const FT gx = 8*x*z + 4*y;
//   const FT gy = 3*y*y + 4*x;
//   const FT gz = 4*x*x;
//   Vector g(gx, gy, gz);
//   return g / std::sqrt(gx*gx + gy*gy + gz*gz);
// };

// auto Shallowtail_value = [](const Point& p) {
//     const FT x = p.x(), y = p.y(), z = p.z();
//     return -4 * z*z*z * y*y
//          - 27 * y*y*y*y
//          + 16 * x * z*z*z*z
//          - 128 * x*x * z*z
//          + 144 * x * y*y * z
//          + 256 * x*x*x;
// };

// auto Shallowtail_gradient = [](const Point& p) {
//     const FT x = p.x(), y = p.y(), z = p.z();
//     const FT gx = 16 * std::pow(z, 4)
//         - 256 * x * z * z
//         + 144 * y * y * z
//         + 768 * x * x;
//     const FT gy = -8 * z * z * z * y
//         - 108 * y * y * y
//         + 288 * x * y * z;
//     const FT gz = -12 * z * z * y * y
//         + 64 * x * z * z * z
//         - 256 * x * x * z
//         + 144 * x * y * y;
//     Vector g(gx, gy, gz);
//     return g / std::sqrt(gx*gx + gy*gy + gz*gz);
// };

auto rotated_wave_value = [](const Point& p) 
{
    const FT theta = 30; // rotation angle
    const FT x = p.x(), y = p.y(), z = p.z();
    const FT sinT = std::sin(theta);
    const FT cosT = std::cos(theta);

    FT xx = x * cosT + z * sinT;
    FT zz = -x * sinT + z * cosT;

    FT val1 = zz - std::sin(xx) - std::sin(y);
    FT val2 = val1 - 0.3;
    FT val3 = val1 + 0.3;
    return val1 * val2 * val3;
};

auto rotated_wave_gradient = [](const Point& p) 
{
    const FT theta = 30; // rotation angle
    const FT x = p.x(), y = p.y(), z = p.z();

    const FT sinT = std::sin(theta);
    const FT cosT = std::cos(theta);

    FT xp = x * cosT + z * sinT;
    FT zp = -x * sinT + z * cosT;

    FT sxp = std::sin(xp);
    FT cxp = std::cos(xp);
    FT sy = std::sin(y);
    FT cy = std::cos(y);

    FT u = zp - sxp - sy;

    FT u_m = u - 0.3;
    FT u_0 = u;
    FT u_p = u + 0.3;

    FT prod = u_m * u_0 * u_p;

    FT dudx = (-sinT) - cxp * cosT; 
    FT dudy = -cy; 
    FT dudz = (cosT) - cxp * sinT; 

    // product rule
    FT sum_factors = (u_0 * u_p) + (u_m * u_p) + (u_m * u_0);

    FT gx = dudx * sum_factors;
    FT gy = dudy * sum_factors;
    FT gz = dudz * sum_factors;

    Vector g(gx, gy, gz);
    return g / std::sqrt(gx*gx + gy*gy + gz*gz);
};


int main(int argc, char** argv)
{
    const FT isovalue = (argc > 1) ? std::stod(argv[1]) : 0;
    const FT box_c = (argc > 2) ? std::abs(std::stod(argv[2])) : 1.;
    const std::size_t grid_n = (argc > 3) ? std::stoi(argv[3]) : 50;

    // create bounding box and grid
    const CGAL::Bbox_3 bbox { -box_c, -box_c, -box_c, box_c, box_c, box_c };
    Grid grid { bbox, CGAL::make_array<std::size_t>(grid_n, grid_n, grid_n) };

    std::cout << "Span: " << grid.span() << std::endl;
    std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
    std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

    // Values values { Sphere_value, grid };
    // Values values { Kusner_Schmitt_value, grid };
    // Values values { S22_value, grid };
    // Values values { Shallowtail_value, grid };
    Values values { rotated_wave_value, grid};

    // Gradients gradients { Sphere_gradient, grid };
    // Gradients gradients { Kusner_Schmitt_gradient, grid };
    // Gradients gradients { S22_gradient, grid };
    // Gradients gradients { Shallowtail_gradient, grid };
    Gradients gradients { rotated_wave_gradient, grid};

    Domain domain = CGAL::Isosurfacing::create_dual_contouring_domain_3(grid, values, gradients);

    // Output containers
    Pointrange points;
    Polygonrange quads;

    std::cout << "Running Dual Marching Cubes with isovalue = " << isovalue << std::endl;

    CGAL::Isosurfacing::internal::dual_marching_cubes(domain, isovalue, points, quads, false, CGAL::Isosurfacing::Vertex_strategy::Centroid);
    // CGAL::Isosurfacing::internal::dual_marching_cubes_tmc(domain, isovalue, points, quads, true, CGAL::Isosurfacing::Vertex_strategy::QEM);

    // for (size_t i = 0; i < points.size(); ++i)
    // std::cout << "Vertex " << i << ": " << points[i] << std::endl;

    // for (const auto& poly : quads) 
    // {
    //     std::cout << "Face: ";
    //     for (auto idx : poly)
    //         std::cout << idx << " ";
    //     std::cout << std::endl;
    // }


    std::cout << "Soup #vertices: " << points.size() << std::endl;
    std::cout << "Soup #quads: " << quads.size() << std::endl;

    std::ofstream out("DMC-soup.off");

    CGAL::IO::write_OFF("DMC-soup.off", points, quads); // print soup to OFF file

    if(!CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(quads)) {
        std::cerr << "Warning: the soup is not a 2-manifold surface, non-manifoldness?..." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Create the mesh..." << std::endl;
    Mesh mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, quads, mesh);

    CGAL::IO::write_polygon_mesh("dual_marching_cubes.off", mesh, CGAL::parameters::stream_precision(17));
    std::cout << "Done" << std::endl;

    return EXIT_SUCCESS;
}
