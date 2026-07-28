#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <iostream>
#include <fstream>
#include <cstdlib>

using Kernel = CGAL::Simple_cartesian<double>;
using Mesh   = CGAL::Surface_mesh<Kernel::Point_3>;

bool get_data(Mesh& mesh, const std::string& filename)
{
  std::ifstream in(filename);
  if (in.fail()) { return false; }
  in>>mesh;
  return true;
}

bool test_minimal_homotopy_basis(const std::string& filename, const std::string& label)
{
  Mesh mesh;
  if (!get_data(mesh, filename))
  {
    std::cerr<<"Fail "<<label<<": Cannot locate file "<<filename<<"\n";
    return false;
  }

  int V=mesh.number_of_vertices();
  int E=mesh.number_of_edges();
  int F=mesh.number_of_faces();
  int expected_nb_loops=2-(V-E+F);

  CGAL::Surface_mesh_topology::Curves_on_surface_topology<Mesh> cst(mesh);
  auto root=*(mesh.halfedges().begin());
  auto basis=cst.compute_minimal_homotopy_basis_with_base_point(root);

  if (static_cast<int>(basis.size())!=expected_nb_loops)
  {
    std::cerr<<"Fail "<<label<<": basis has "<<basis.size()<<" loop(s), expected "
             <<expected_nb_loops<<" (2-(V-E+F))\n";
    return false;
  }

  for (std::size_t i=0; i<basis.size(); ++i)
  {
    if (!basis[i].is_closed())
    {
      std::cerr<<"Fail "<<label<<": loop "<<i<<" is not closed\n";
      return false;
    }
    if (cst.is_contractible(basis[i]))
    {
      std::cerr<<"Fail "<<label<<": loop "<<i<<" is contractible (should not be)\n";
      return false;
    }
  }

  for (std::size_t i=0; i<basis.size(); ++i)
  {
    for (std::size_t j=i+1; j<basis.size(); ++j)
    {
      if (cst.are_freely_homotopic(basis[i], basis[j]))
      {
        std::cerr<<"Fail "<<label<<": loops "<<i<<" and "<<j<<" are freely homotopic "
                 <<"(should be independent)\n";
        return false;
      }
    }
  }

  return true;
}

int main()
{
  std::ifstream list_file("data/minimal_length_homotopy_basis_meshes.txt");
  if (!list_file)
  {
    std::cerr<<"Cannot open data/minimal_length_homotopy_basis_meshes.txt\n";
    return EXIT_FAILURE;
  }

  bool res=true;
  std::string mesh_name;
  while (std::getline(list_file, mesh_name))
  {
    if (!mesh_name.empty())
    { res=test_minimal_homotopy_basis(CGAL::data_file_path("meshes/"+mesh_name), mesh_name) && res; }
  }

  if (res)
  {
    std::cout<<"All tests passed\n";
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
