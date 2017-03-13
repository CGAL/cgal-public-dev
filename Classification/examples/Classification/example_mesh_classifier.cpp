#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classifier.h>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Mesh_neighborhood.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>

#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;

typedef CGAL::Surface_mesh<Point> Mesh;
typedef typename Mesh::Face_range Face_range;
typedef typename boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef CGAL::Identity_property_map<face_descriptor> Face_map;

typedef CGAL::Classifier<Face_range, Face_map> Classifier;

typedef CGAL::Classification::Label_handle                                            Label_handle;
typedef CGAL::Classification::Feature_handle                                          Feature_handle;

typedef CGAL::Classification::Mesh_neighborhood<Mesh> Neighborhood;

///////////////////////////////////////////////////////////////////
//! [Analysis]

int main (int argc, char** argv)
{
  std::string filename (argc > 1 ? argv[1] : "data/example.off");
  std::ifstream in (filename.c_str());
  Mesh mesh;

  std::cerr << "Reading input" << std::endl;
  in >> mesh;

  Classifier classifier (mesh.faces(), Face_map());
  
  Neighborhood neighborhood (mesh);

  classifier.run_with_local_smoothing (neighborhood.one_ring_neighbor_query());

  
  return EXIT_SUCCESS;
}
