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
    const double vx = 1.0, vy = 1.0, vz = 1.0;
    // const std::size_t offset = 0;

    if (!image.read_raw(raw_file.c_str(), nx, ny, nz, vx, vy, vz)) 
    {
        std::cerr << "Error: failed to read raw fiel " << raw_file << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Read successful: "
            << image.xdim() << " x "
            << image.ydim() << " x "
            << image.zdim() << std::endl;

    const CGAL::Bbox_3 bbox(0.0, 0.0, 0.0, vx * nx, vy * ny, vz * nz);
    Grid grid(bbox, CGAL::make_array<std::size_t>(nx, ny, nz));

    // Values values([&](const Point& p) -> FT 
    // {
    //     float x = static_cast<float>(p.x());
    //     float y = static_cast<float>(p.y());
    //     float z = static_cast<float>(p.z());

    //     return static_cast<FT>(image.value(x, y, z));
    // }, grid);

    // -1 to get rid of the segfault
    Values values([&](const Point& p) -> FT 
    {
        float x = std::clamp(static_cast<float>(p.x()), 0.0f, float(vx * (nx - 1)));
        float y = std::clamp(static_cast<float>(p.y()), 0.0f, float(vy * (ny - 1)));
        float z = std::clamp(static_cast<float>(p.z()), 0.0f, float(vz * (nz - 1)));

        return static_cast<FT>(image.value(x, y, z));
    }, grid);

    Gradient gradients(values, FT(1), FT(1), FT(1)); 

    Domain domain = CGAL::Isosurfacing::create_dual_contouring_domain_3(grid, values, gradients);

    Pointrange points;
    Polygonrange quads;

    std::cout << "Running Dual Marching Cubes on " << raw_file << " with isovalue = " << isovalue << std::endl;

    // std::size_t non_zero_count = 0;
    // float min_val = 1e9, max_val = -1e9;

    // for (std::size_t k = 0; k < nz; ++k)
    // for (std::size_t j = 0; j < ny; ++j)
    //     for (std::size_t i = 0; i < nx; ++i) 
    // {
    //     float val = image.value(i, j, k);
    //     if (val != 0) {
    //         ++non_zero_count;
    //         min_val = std::min(min_val, val);
    //         max_val = std::max(max_val, val);
    //     }
    //     }

    // std::cout << "Non-zero voxels: " << non_zero_count << std::endl;
    // std::cout << "Actual non-zero min = " << min_val << ", max = " << max_val << std::endl;

    // exit(0);

    // CGAL::Isosurfacing::internal::dual_marching_cubes(domain, isovalue, points, quads, false, CGAL::Isosurfacing::Vertex_strategy::QEM);
    CGAL::Isosurfacing::internal::dual_marching_cubes_tmc(domain, isovalue, points, quads, true, CGAL::Isosurfacing::Vertex_strategy::Centroid);

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
