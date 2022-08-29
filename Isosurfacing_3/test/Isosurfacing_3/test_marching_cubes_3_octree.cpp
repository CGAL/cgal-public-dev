// Copyright (c) 2022  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ágoston Sipos

#include <CGAL/Octree.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Aff_transformation_3.h>

#include <CGAL/Marching_cubes_octree.h>
#include <CGAL/Minimal_area_triangulation.h>

#include <functional>
#include <algorithm>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef CGAL::Point_3<Kernel> Point;

typedef Octree_wrapper<Kernel> Octree;

typedef std::function<FT(Point)> ImplicitFunction;

// writes out polygons as returned by the algorithm to an OBJ file
void write_mesh(const std::vector<Point>& vertices, const std::vector<std::vector<size_t>>& faces, const std::string& filename) {
    std::ofstream mesh_out(filename.c_str());
    for(auto p : vertices) {
        mesh_out << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }
    for (auto f : faces) {
        mesh_out << "f ";
        for(auto p : f)
            mesh_out << (p+1) << " ";
        mesh_out << std::endl;
    }
    std::cout << "Result exported to " << filename << "\n";
}

// writes out edges of and octree to an OBJ file
void write_octree(const CGAL::Octree_domain<Kernel>& domain, const std::string& filename) {
    std::ofstream octree_out("octree.obj");
    int v = 1;
    auto ff = [&v,domain, &octree_out](const typename CGAL::Octree_domain<Kernel>::Edge_handle& eh){ auto ps = domain.edge_vertices(eh); auto p1 = domain.position(ps[0]), p2 = domain.position(ps[1]); octree_out << "v " << p1.x() << " " << p1.y() << " " << p1.z() << std::endl << "v " << p2.x() << " " << p2.y() << " " << p2.z() << std::endl << "l " << v++ << " " << v++ << std::endl; };
    domain.iterate_edges(ff);
}

// Custom refinement predicate, splits cell if mid point of it is "close" to isosurface
struct Split_by_closeness {
    ImplicitFunction func;
    Octree oct;
    Split_by_closeness(ImplicitFunction func, Octree oct) : func(func), oct(oct) {}

    bool operator()(const Octree::Node &n) const {
        CGAL::Bbox_3 b = oct.bbox(oct.root()); // maybe easier to store only this instead of octree
                                               // note: oct.bbox(n) does not give correct result here
        int d = n.depth();

        // calculating midpoint of box of node
        FT x = computeMiddle(b.xmin(), b.xmax(), n.global_coordinates()[0], n.depth());
        FT y = computeMiddle(b.ymin(), b.ymax(), n.global_coordinates()[1], n.depth());
        FT z = computeMiddle(b.zmin(), b.zmax(), n.global_coordinates()[2], n.depth());
        Point mid{x,y,z}; // note: oct.barycenter(n) does not give correct result here

        return (n.depth() <= 2 || func(mid) < 0.05) && n.depth() <= 4; // custom predicate, can be different
    }

private:
    FT computeMiddle(FT min, FT max, int c, int d) const {
        return min + (max - min) * (c + 0.5) / (1 << d);
    }
};

// creates an octree by splitting the bounding box and splitting one of the cells again
struct split_only_twice {
    bool operator()(const Octree::Node &n) const {
        return n.depth() <= 1 && n.global_coordinates() == std::array<uint32_t, 3>{0,0,0};
    }
};

void test1() {

    std::cout << "Test #1: one cell, two opposite corners are positive...\n";

    ImplicitFunction corners = [](Point p) {
        return std::abs(p.x()+p.y()+p.z()) - 2.5;
    };

    Octree octree(CGAL::Bbox_3(-1.2,-1.2,-1.2,1.2,1.2,1.2));
    octree.refine([](const Octree::Node &n){return false;}, corners);
    CGAL::Octree_domain domain(octree);

    std::vector<Point> vertices;
    std::vector<std::vector<size_t>> faces;

    make_polygon_mesh_using_marching_cubes_on_octree(domain, 0.0, vertices, faces);

    bool b = faces.size() == 2 && faces[0].size() == 3 && faces[1].size() == 3;
    std::cout << "\tExpected result: two triangles [" << (b ? "OK" : "Error!") << "]\n\n";
}

void test2() {

    std::cout << "Test #2: a cube's distance function on the border of bigger and smaller cells...\n";

    ImplicitFunction small_cube = [](Point p) {
        return std::max({std::abs(p.x() + 0.5), std::abs(p.y()), std::abs(p.z())}) - 0.4;
    };

    Octree octree(CGAL::Bbox_3(-1.2,-1.2,-1.2,1.2,1.2,1.2));
    octree.refine(split_only_twice(), small_cube);
    CGAL::Octree_domain domain(octree);

    std::vector<Point> vertices;
    std::vector<std::vector<size_t>> faces;

    make_polygon_mesh_using_marching_cubes_on_octree(domain, 0.0, vertices, faces);

    bool b = faces.size() == 4;
    for(auto it : faces) b = b && (it.size() == 3);
    std::cout << "\tExpected result: a tetrahedron [" << (b ? "OK" : "Error!") << "]\n\n";

}

int main(int argc, char** argv) {

    test1();

    test2();

    ImplicitFunction sphere = [](Point p) { return sqrt((p.x()-0.1)*(p.x()-0.1) + (p.y()-0.2)*(p.y()-0.2) + p.z()*p.z()) - 0.5; };
    ImplicitFunction cube = [](Point p) {
        return std::max({std::abs(p.x()), std::abs(p.y()), std::abs(p.z())}) - 0.5;
    };
    ImplicitFunction rotcube = [](Point p) {
        double cosa = cos(1), sina = sin(1), cosb = cos(0.5), sinb = sin(0.5);
        Kernel::Aff_transformation_3 rot1(
            1.0, 0.0, 0.0,
            0.0, cosa, -sina,
            0.0, sina, cosa);
        Kernel::Aff_transformation_3 rot2(
            cosb, -sinb, 0.0,
            sinb, cosb, 0.0,
            0.0, 0.0, 1.0);
        Point q = rot2(rot1(p));
        return std::max({std::abs(q.x()), std::abs(q.y()), std::abs(q.z())}) - 0.5;
    };
    ImplicitFunction torus = [](Point p) {
        return pow(sqrt(p.x()*p.x() + p.y()*p.y()) - 0.5, 2) + p.z()*p.z() - 0.2*0.2;
    };
    ImplicitFunction tanglecube = [](Point p) {
        double x = 2*p.x(), y = 2*p.y(), z = 2*p.z();
        double x2=x*x, y2=y*y, z2=z*z;
        double x4=x2*x2, y4=y2*y2, z4=z2*z2;
        return x4 - 5*x2 + y4 - 5*y2 + z4 - 5*z2 + 11.8;
    };

    ImplicitFunction func;
    if(argc >= 2 && std::string(argv[1]) == "--function") {
        if(std::string(argv[2]) == "sphere")
            func = sphere;
        else if(std::string(argv[2]) == "cube")
            func = cube;
        else if(std::string(argv[2]) == "rotcube")
            func = rotcube;
        else if(std::string(argv[2]) == "torus")
            func = torus;
        else if(std::string(argv[2]) == "tanglecube")
            func = tanglecube;
        else exit(1);
    }
    else
        func = sphere;

    std::string name = argc >= 2 ? argv[2] : "sphere";

    std::cout << "Using function " << name << "...\n";

    Octree octree(CGAL::Bbox_3(-1.2,-1.2,-1.2,1.2,1.2,1.2));
    octree.refine(Split_by_closeness(func, octree), func);

    CGAL::Octree_domain domain(octree);

    std::vector<Point> vertices;
    std::vector<std::vector<size_t>> faces;

    make_polygon_mesh_using_marching_cubes_on_octree(domain, 0.0, vertices, faces);

    std::cout << "# of polygons: " << faces.size() << std::endl;

    write_mesh(vertices, faces, name + ".obj");

    std::vector<std::vector<size_t>> triangles = minimal_area_triangulate(vertices, faces);

    std::cout << "# of triangles: " << triangles.size() << std::endl;

    write_mesh(vertices, triangles, name + "3.obj");

    return 0;
}