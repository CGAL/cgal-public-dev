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
#include <CGAL/Image_3.h>

#include <CGAL/Isosurfacing_3/Finite_difference_gradient_3.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = Kernel::FT;
using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
using Values = CGAL::Isosurfacing::Value_function_3<Grid>;
using Gradient = CGAL::Isosurfacing::Finite_difference_gradient_3<Kernel>;
using Domain = CGAL::Isosurfacing::Dual_contouring_domain_3<Grid, Values, Gradient>;

using Pointrange = std::vector<Point>;
using Polygonrange = std::vector<std::vector<std::size_t>>;
using Mesh = CGAL::Surface_mesh<Point>;

int main(int argc, char** argv)
{
    const std::string raw_file = argv[1];
    const FT isovalue = std::stod(argv[2]);
    const std::size_t nx = static_cast<std::size_t>(std::stoul(argv[3]));
    const std::size_t ny = static_cast<std::size_t>(std::stoul(argv[4]));
    const std::size_t nz = static_cast<std::size_t>(std::stoul(argv[5]));
    const std::string out_file = (argc > 6) ? argv[6] : "dual_mc_raw.off";

    CGAL::Image_3 image;

    // const std::size_t nx = 256;
    // const std::size_t ny = 256;
    // const std::size_t nz = 98;
    const double vx = 1, vy = 1, vz = 1;
    // const std::size_t offset = 0;

    if (!image.read_raw(raw_file.c_str(), nx, ny, nz, vx, vy, vz)) 
    {
        std::cerr << "Error: failed to read raw file " << raw_file << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Read successful: "
            << image.xdim() << " x "
            << image.ydim() << " x "
            << image.zdim() << std::endl;

    // // Log raw image values
    // std::ofstream log("raw_image_values.log");
    // log << "# x y z value\n";
    // for (std::size_t x = 0; x < nx; ++x)
    //     for (std::size_t y = 0; y < ny; ++y)
    //         for (std::size_t z = 0; z < nz; ++z)
    //         {
    //             FT val = image.value(x, y, z);
    //             if (!CGAL::is_finite(val))
    //             {
    //                 std::cerr << "Non-finite value at (" << x << ", " << y << ", " << z
    //                         << "): " << val << "\n";
    //                 log.close();
    //                 return EXIT_FAILURE;
    //             }
    //             log << x << " " << y << " " << z << " " << val << "\n";
    //         }
    // log.close();
    // std::cout << "Raw voxel values logged to raw_image_values.log\n";

    const CGAL::Bbox_3 bbox(0.0, 0.0, 0.0, vx * nx, vy * ny, vz * nz);

    
    std::size_t factor = 3;
    // Grid grid(bbox, CGAL::make_array<std::size_t>(nx / factor, ny / factor, nz / factor));
    Grid grid(bbox, CGAL::make_array<std::size_t>(30, 30, 30));

    // make padded grid to resolve nonmanifold vertices (works but not for all and change original geometry (bad))
    // Grid grid1(bbox, CGAL::make_array<std::size_t>(30, 30, 30));
    // Grid grid = CGAL::Isosurfacing::internal::make_padded_grid(grid1);

    Values values([&](const Point& p) -> FT 
    {
        float x = std::clamp(static_cast<float>(p.x()), 0.0f, float(vx * (nx - 1)));
        float y = std::clamp(static_cast<float>(p.y()), 0.0f, float(vy * (ny - 1)));
        float z = std::clamp(static_cast<float>(p.z()), 0.0f, float(vz * (nz - 1)));

         float val = image.value(x, y, z);
        if (std::isnan(val)) {
            return std::numeric_limits<FT>::quiet_NaN(); // propagate NaN to domain
        }

        return static_cast<FT>(val);
    }, grid);

    Gradient gradients(values, FT(vx), FT(vy), FT(vz));

    Domain domain = CGAL::Isosurfacing::create_dual_contouring_domain_3(grid, values, gradients);

    Pointrange points;
    Polygonrange quads;

    std::cout << "Running Dual Marching Cubes on " << raw_file << " with isovalue = " << isovalue << std::endl;


    // CGAL::Isosurfacing::internal::dual_marching_cubes(domain, isovalue, points, quads, false, CGAL::Isosurfacing::Vertex_strategy::Centroid);
    // CGAL::Isosurfacing::internal::dual_marching_cubes_tmc(domain, isovalue, points, quads, false, CGAL::Isosurfacing::Vertex_strategy::Centroid, CGAL::Isosurfacing::internal::PostProcessOff());

    CGAL::Isosurfacing::internal::gridDim gridDims = {grid.xdim(), grid.ydim(), grid.zdim()};
    CGAL::Isosurfacing::internal::dual_marching_cubes_tmc(domain, isovalue, points, quads, false, CGAL::Isosurfacing::Vertex_strategy::Centroid, CGAL::Isosurfacing::internal::PostProcessOn(), gridDims);

    std::cout << "Soup #vertices: " << points.size() << std::endl;
    std::cout << "Soup #quads: " << quads.size() << std::endl;

    CGAL::IO::write_OFF(out_file, points, quads);

    if (!CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(quads)) {
        std::cerr << "Warning: the soup is not a 2-manifold surface, non-manifoldness?..." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Create the mesh..." << std::endl;
    Mesh mesh;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, quads, mesh);

    std::string cleaned_file = "cleaned_" + out_file;
    CGAL::IO::write_polygon_mesh(cleaned_file, mesh, CGAL::parameters::stream_precision(17));
    std::cout << "Done" << std::endl;

    return EXIT_SUCCESS;
}
