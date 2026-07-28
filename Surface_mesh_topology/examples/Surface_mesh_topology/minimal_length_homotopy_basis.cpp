#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/draw_face_graph_with_paths.h>
#include <CGAL/IO/polygon_mesh_io.h>

#include <iostream>
#include <fstream>
#include <cstdlib>

using Mesh            = CGAL::Surface_mesh<CGAL::Simple_cartesian<double>::Point_3>;
using Path_on_surface = CGAL::Surface_mesh_topology::Path_on_surface<Mesh>;

double cycle_length(const Mesh& mesh, const Path_on_surface& cycle)
{ // Compute the length of the given cycle.
  double res=0;
  for (std::size_t i=0; i<cycle.length(); ++i)
  { res+=std::sqrt
        (CGAL::squared_distance(mesh.point(mesh.vertex(mesh.edge(cycle[i]), 0)),
                                mesh.point(mesh.vertex(mesh.edge(cycle[i]), 1)))); }
  return res;
}

void display_cycle_info(const Mesh& mesh, const Path_on_surface& cycle, std::size_t i)
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

void write_polylines_file(const Mesh& mesh, const std::vector<Path_on_surface>& basis,
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

int main(int argc, char* argv[])
{
  std::string filename=CGAL::data_file_path("meshes/3torus.off");
  bool draw=false;
  std::string dump_filename;
  for (int i=1; i<argc; ++i)
  {
    if (std::string(argv[i])=="-draw")
    { draw=true; }
    else if (std::string(argv[i])=="-dump" && i+1<argc)
    { dump_filename=argv[++i]; }
    else
    { filename=argv[i]; }
  }

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cout<<"Cannot read file '"<<filename<<"'. Exiting program"<<std::endl;
    return EXIT_FAILURE;
  }
  std::cout<<"File '"<<filename<<"' loaded. Computing minimal length homotopy basis..."<<std::endl;

  CGAL::Surface_mesh_topology::Curves_on_surface_topology<Mesh> cst(sm);
  auto root=*(sm.halfedges().begin()); // One dart of the mesh

  CGAL::Surface_mesh_topology::Euclidean_length_weight_functor<Mesh> wf(sm);
  auto basis=cst.compute_minimal_homotopy_basis_with_base_point(root, wf);

  std::cout<<"Basis has "<<basis.size()<<" loop(s)."<<std::endl;
  for (std::size_t i=0; i<basis.size(); ++i)
  { display_cycle_info(sm, basis[i], i); }

  if (!dump_filename.empty())
  { write_polylines_file(sm, basis, dump_filename); }

  if (draw)
  { CGAL::draw(sm, basis); }

  return EXIT_SUCCESS;
}
