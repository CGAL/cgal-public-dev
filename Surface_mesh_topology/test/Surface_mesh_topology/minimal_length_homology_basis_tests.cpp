#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polygonal_schema.h>
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

bool test_minimal_homology_basis(const std::string& filename, const std::string& label)
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
  auto basis=cst.compute_minimal_homology_basis();

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

// n x m quad grid torus: vertex (i,j) at index i*m+j, at coordinates
// (i,j,0) -- indices wrap modulo n/m so the grid closes up into a torus
// (checked separately: is_closed() holds, Euler characteristic 0).
Mesh make_torus_grid(std::size_t n, std::size_t m)
{
  Mesh mesh;
  std::vector<Mesh::Vertex_index> v(n*m);
  for (std::size_t i=0; i<n; ++i)
  { for (std::size_t j=0; j<m; ++j)
    { v[i*m+j]=mesh.add_vertex(Kernel::Point_3(double(i), double(j), 0)); } }

  for (std::size_t i=0; i<n; ++i)
  {
    for (std::size_t j=0; j<m; ++j)
    {
      std::size_t i1=(i+1)%n, j1=(j+1)%m;
      mesh.add_face(v[i*m+j], v[i1*m+j], v[i1*m+j1], v[i*m+j1]);
    }
  }
  return mesh;
}

// Column i=0 (the "vertical" ring of edges connecting (0,j) to (0,j+1))
// costs 100 per edge; every other edge costs its Euclidean length -- so a
// generator wrapping in the j-direction should route through a
// *different* i column if the algorithm is genuinely weight-minimal, not
// just topologically correct.
struct Weighted_column_functor
{
  using Weight_t = double;
  Weighted_column_functor(const Mesh& mesh) : m_mesh(mesh) {}
  double operator()(Mesh::Halfedge_index he) const
  {
    auto p0=m_mesh.point(m_mesh.vertex(m_mesh.edge(he), 0));
    auto p1=m_mesh.point(m_mesh.vertex(m_mesh.edge(he), 1));
    if (p0.x()==0 && p1.x()==0) { return 100.0; }
    return std::sqrt(CGAL::to_double(CGAL::squared_distance(p0, p1)));
  }
private:
  const Mesh& m_mesh;
};

bool test_weighted_column_grid()
{
  const std::string label="weighted 3x3 torus grid";
  Mesh grid=make_torus_grid(3, 3);
  Weighted_column_functor wf(grid);

  CGAL::Surface_mesh_topology::Curves_on_surface_topology<Mesh> cst(grid);
  auto basis=cst.compute_minimal_homology_basis(wf);

  if (basis.size()!=2)
  {
    std::cerr<<"Fail "<<label<<": basis has "<<basis.size()<<" loop(s), expected 2\n";
    return false;
  }

  for (std::size_t i=0; i<basis.size(); ++i)
  {
    if (!basis[i].is_closed())
    {
      std::cerr<<"Fail "<<label<<": loop "<<i<<" is not closed\n";
      return false;
    }
    for (std::size_t k=0; k<basis[i].length(); ++k)
    {
      auto p0=grid.point(CGAL::source(basis[i].get_ith_real_dart(k), grid));
      auto p1=grid.point(CGAL::target(basis[i].get_ith_real_dart(k), grid));
      if (p0.x()==0 && p1.x()==0)
      {
        std::cerr<<"Fail "<<label<<": loop "<<i<<" uses the weight-100 column "
                 <<"(should route around it)\n";
        return false;
      }
    }
  }

  return true;
}

// Double torus (genus 2, 2g=4) built from an explicit octagon
// identification word, on a Mesh_ type other than Surface_mesh -- this is
// the exact fixture that caught the BGL-compatibility bug
// (num_vertices/num_edges/halfedges don't compile for
// Polygonal_schema_with_combinatorial_map, unlike Surface_mesh); kept as
// an official test so that specific regression can't come back silently.
bool test_bitorus_polygonal_schema()
{
  const std::string label="bitorus (Polygonal_schema)";
  using PS=CGAL::Surface_mesh_topology::Polygonal_schema_with_combinatorial_map<>;
  PS ps;
  ps.add_facet("a b -a -b c d -c -d");

  CGAL::Surface_mesh_topology::Curves_on_surface_topology<PS> cst(ps);
  auto basis=cst.compute_minimal_homology_basis();

  if (basis.size()!=4)
  {
    std::cerr<<"Fail "<<label<<": basis has "<<basis.size()<<" loop(s), expected 4\n";
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
  bool res=test_weighted_column_grid();
  res=test_bitorus_polygonal_schema() && res;

  std::ifstream list_file("data/minimal_length_homology_basis_meshes.txt");
  if (!list_file)
  {
    std::cerr<<"Cannot open data/minimal_length_homology_basis_meshes.txt\n";
    return EXIT_FAILURE;
  }

  std::string mesh_name;
  while (std::getline(list_file, mesh_name))
  {
    if (!mesh_name.empty())
    { res=test_minimal_homology_basis(CGAL::data_file_path("meshes/"+mesh_name), mesh_name) && res; }
  }

  if (res)
  {
    std::cout<<"All tests passed\n";
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
