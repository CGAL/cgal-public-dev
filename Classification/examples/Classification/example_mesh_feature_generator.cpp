#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>

#include <CGAL/Surface_mesh/Surface_mesh.h>

#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;

typedef CGAL::Surface_mesh<Point> Mesh;
typedef typename Mesh::Face_range Face_range;
typedef typename boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef CGAL::Identity_property_map<face_descriptor> Face_map;

namespace Classif = CGAL::Classification;

typedef Classif::Sum_of_weighted_features_predicate Classification_predicate;

typedef Classif::Label_handle                                            Label_handle;
typedef Classif::Feature_handle                                          Feature_handle;
typedef Classif::Label_set                                               Label_set;
typedef Classif::Feature_set                                             Feature_set;

typedef Classif::Mesh_neighborhood<Mesh> Neighborhood;

template <typename FaceGraph, typename Point>
struct Face_graph_face_to_center_property_map
{
  typedef typename boost::graph_traits<Mesh>::face_descriptor key_type;
  typedef Point value_type;
  typedef Point reference;
  typedef boost::readable_property_map_tag category;

  typedef typename boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  
  const FaceGraph* mesh;

  Face_graph_face_to_center_property_map (const FaceGraph* mesh) : mesh (mesh) { }

  friend reference get (const Face_graph_face_to_center_property_map& map, key_type f)
  {
    std::vector<Point> points;
    BOOST_FOREACH(vertex_descriptor v, vertices_around_face(halfedge(f, *(map.mesh)), *(map.mesh)))
    {
      points.push_back (map.mesh->point(v));
    }
    return CGAL::centroid (points.begin(), points.end());
  }
};

typedef Face_graph_face_to_center_property_map<Mesh, Point> Face_center_map;

typedef Classif::Mesh_feature_generator<Kernel, Mesh, Face_center_map>             Generator;


int main (int argc, char** argv)
{
  std::string filename (argc > 1 ? argv[1] : "data/example.off");
  std::ifstream in (filename.c_str());
  Mesh mesh;

  std::cerr << "Reading input" << std::endl;
  in >> mesh;
  std::cerr << " * " << mesh.number_of_vertices() << " vertices" << std::endl;
  std::cerr << " * " << mesh.number_of_faces() << " faces" << std::endl;
  std::cerr << " * " << mesh.faces().size() << " faces" << std::endl;
  std::cerr << "Computing useful structures" << std::endl;

  Face_center_map fc_map (&mesh);
  Face_range faces = mesh.faces();

  
  std::cerr << "Computing features" << std::endl;
  Feature_set features;
  Generator generator (features, mesh, fc_map, 1);

  std::cerr << "Setting up labels" << std::endl;
  Label_set labels;
  Label_handle ground = labels.add ("ground");
  Label_handle vege = labels.add ("vegetation");
  Label_handle roof = labels.add ("roof");

  // Run classification
  std::cerr << "Classifying" << std::endl;
  
  std::vector<std::size_t> label_indices;

  Classification_predicate predicate (labels, features);
  
  CGAL::Real_timer t;
  t.start();
  Classif::classify_with_local_smoothing<CGAL::Parallel_tag>
    (faces, Face_map(), labels, predicate,
     generator.neighborhood().one_ring_neighbor_query(),
     label_indices);
  t.stop();
  std::cerr << "Classification with graphcut performed in " << t.time() << " second(s)" << std::endl;
  
  std::ofstream fout("out.off");
  fout << "COFF" << std::endl
       << mesh.number_of_vertices() << " " << mesh.number_of_faces() << " 0" << std::endl;
  BOOST_FOREACH (vertex_descriptor vd, vertices(mesh))
  {
    fout << mesh.point(vd) << std::endl;
  }
  BOOST_FOREACH (face_descriptor fd, mesh.faces())
  {
    fout << "3";
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(halfedge(fd, mesh), mesh))
    {
      fout << " " << std::size_t(vd);
    }

    Label_handle label = labels[label_indices[fd]];
    if (label == ground)
      fout << " 245 180 0" << std::endl;
    else if (label == vege)
      fout << " 0 255 27" << std::endl;
    else if (label == roof)
      fout << " 255 0 170" << std::endl;
    else
    {
      fout << " 0 0 0" << std::endl;
      std::cerr << "Error: unknown classification label" << std::endl;
    }
  }

  return EXIT_SUCCESS;
}
