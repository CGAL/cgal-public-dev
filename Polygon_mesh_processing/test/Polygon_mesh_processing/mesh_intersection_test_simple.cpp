#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Real_timer.h> 
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <inputA.off> <inputB.off> <output.txt>" << std::endl;
        return 1;
    }

    std::string fileA = argv[1];
    std::string fileB = argv[2];
    std::string outFile = argv[3];

    CGAL::Real_timer total_timer;
    total_timer.start();

    // [Step 1&2] Load Meshes & Build Internal Data Structures
    CGAL::Real_timer t_io;
    t_io.start();
    Mesh m1, m2;
    std::ifstream in1(fileA);
    std::ifstream in2(fileB);
    if (!(in1 >> m1) || !(in2 >> m2)) {
        std::cerr << "Error: Could not read input OFF files." << std::endl;
        return 1;
    }
    t_io.stop();

    // [Step 0] Print Mesh Sizes
    // Added right after loading so you can correlate size with time
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "Mesh A: " << m1.number_of_faces() << " faces | " << m1.number_of_vertices() << " vertices" << std::endl;
    std::cout << "Mesh B: " << m2.number_of_faces() << " faces | " << m2.number_of_vertices() << " vertices" << std::endl;
    std::cout << "-----------------------------------" << std::endl;

    // [Step 3&4] BVH Build & Intersection Query
    CGAL::Real_timer t_query;
    t_query.start();
    std::vector<std::pair<face_descriptor, face_descriptor>> intersected_tris;
    
    CGAL::Polygon_mesh_processing::internal::compute_face_face_intersection(
        m1, 
        m2, 
        std::back_inserter(intersected_tris),
        CGAL::parameters::all_default(), 
        CGAL::parameters::all_default()  
    );
    t_query.stop();

    // [Step 8] Write Results to Disk
    CGAL::Real_timer t_write;
    t_write.start();
    std::ofstream out(outFile);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open output file " << outFile << std::endl;
        return 1;
    }
    out << "FaceIndexA FaceIndexB\n";
    for (const auto& pair : intersected_tris) {
        out << (size_t)pair.first << " " << (size_t)pair.second << "\n";
    }
    out.close();
    t_write.stop();

    total_timer.stop();

    // --- Output Stats ---
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "[Step 1&2] Load & Mesh Build:   " << t_io.time() * 1000.0 << " ms" << std::endl;
    std::cout << "[Step 3&4] Query (Math):        " << t_query.time() * 1000.0 << " ms" << std::endl;
    std::cout << "[Step 8]   Write CSV to Disk:   " << t_write.time() * 1000.0 << " ms" << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "Intersections found: " << intersected_tris.size() << std::endl;
    std::cout << "Total runtime:       " << total_timer.time() << " s" << std::endl;
    std::cout << "-----------------------------------" << std::endl;

    return 0;
}
