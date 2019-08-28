#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;
using namespace CGAL::Surface_mesh_topology;

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  std::string file;
  for (unsigned i=1; i<argc; ++i)
  {
    LCC_3_cmap map;
    file = argv[i];
    std::ifstream in(file);
    if (!in.is_open())
    {
      std::cout<<"ERROR reading file " + file<<std::endl;
      exit(EXIT_FAILURE);
    }
    CGAL::load_off(map, file.c_str());
    std::vector<Path_on_surface<LCC_3_cmap> > paths;
    Curves_on_surface_topology<LCC_3_cmap> cst(map);
    std::cout << "handling " << file << "..." << std::endl;
    cst.compute_minimal_quadrangulation(true);
  
    std::cout<<"original_map characteristics : ";
    map.display_characteristics(std::cout) << std::endl;
    std::cout<<"quadrangulation characteristics : ";
    cst.get_minimal_quadrangulation().display_characteristics(std::cout) << std::endl << std::endl;
  }
  return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
