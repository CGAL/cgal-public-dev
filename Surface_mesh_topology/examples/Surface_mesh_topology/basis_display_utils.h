#ifndef CGAL_SURFACE_MESH_TOPOLOGY_EXAMPLES_BASIS_DISPLAY_UTILS_H
#define CGAL_SURFACE_MESH_TOPOLOGY_EXAMPLES_BASIS_DISPLAY_UTILS_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/squared_distance_3.h>

#include <iostream>
#include <fstream>

// Shared by minimal_length_homotopy_basis.cpp and
// minimal_length_homology_basis.cpp: both load a Surface_mesh, compute a
// basis of loops, then display/export it the same way -- only the basis
// computation itself differs between the two examples.
using Mesh            = CGAL::Surface_mesh<CGAL::Simple_cartesian<double>::Point_3>;
using Path_on_surface = CGAL::Surface_mesh_topology::Path_on_surface<Mesh>;

inline double cycle_length(const Mesh& mesh, const Path_on_surface& cycle)
{ // Compute the length of the given cycle.
  double res=0;
  for (std::size_t i=0; i<cycle.length(); ++i)
  { res+=std::sqrt
        (CGAL::squared_distance(mesh.point(mesh.vertex(mesh.edge(cycle[i]), 0)),
                                mesh.point(mesh.vertex(mesh.edge(cycle[i]), 1)))); }
  return res;
}

inline void display_cycle_info(const Mesh& mesh, const Path_on_surface& cycle, std::size_t i)
{ // Display information about the given cycle.
  if (cycle.is_empty())
  { std::cout<<"Loop "<<i<<": empty."<<std::endl; }
  else
  {
    std::cout<<"Loop "<<i<<": Root: "<<mesh.point(mesh.vertex(mesh.edge(cycle[0]), 0))<<"; "
             <<"Number of edges: "<<cycle.length()<<"; "
             <<"Length: "<<cycle_length(mesh, cycle)<<std::endl;
  }
}

inline void write_polylines_file(const Mesh& mesh, const std::vector<Path_on_surface>& basis,
                                  const std::string& filename)
{ // Export every loop of the basis as a polyline, in the CGAL .polylines.txt
  // format (one line per loop: "n x1 y1 z1 x2 y2 z2 ... xn yn zn"), so the
  // basis can be opened directly in CGALlab alongside the mesh.
  // Uses get_ith_real_dart(i), not cycle[i], since the latter ignores the
  // path's internal flip bit and can give points in the wrong order.
  std::ofstream out(filename);
  if (!out)
  {
    std::cerr<<"Cannot open '"<<filename<<"' for writing. Skipping polylines export."<<std::endl;
  }
  else
  {
    for (const Path_on_surface& cycle : basis)
    {
      std::size_t n=cycle.length();
      out<<(n+1); // +1: repeat the first point to close the loop visually
      for (std::size_t i=0; i<n; ++i)
      {
        auto p=mesh.point(CGAL::source(cycle.get_ith_real_dart(i), mesh));
        out<<" "<<p.x()<<" "<<p.y()<<" "<<p.z();
      }
      auto p0=mesh.point(CGAL::source(cycle.get_ith_real_dart(0), mesh));
      out<<" "<<p0.x()<<" "<<p0.y()<<" "<<p0.z()<<"\n";
    }
    std::cout<<"Wrote "<<basis.size()<<" polylines to "<<filename<<std::endl;
  }
}

#endif // CGAL_SURFACE_MESH_TOPOLOGY_EXAMPLES_BASIS_DISPLAY_UTILS_H
