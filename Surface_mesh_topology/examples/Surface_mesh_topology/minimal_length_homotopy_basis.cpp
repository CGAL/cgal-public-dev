#include "basis_display_utils.h"

#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/draw_face_graph_with_paths.h>
#include <CGAL/IO/polygon_mesh_io.h>

#include <cstdlib>

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
