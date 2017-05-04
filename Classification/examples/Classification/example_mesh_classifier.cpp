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

namespace Classification = CGAL::Classification;

typedef Classification::Sum_of_weighted_features_classifier Classifier;

typedef Classification::Label_handle                                            Label_handle;
typedef Classification::Feature_handle                                          Feature_handle;
typedef Classification::Label_set                                               Label_set;
typedef Classification::Feature_set                                             Feature_set;

typedef Classification::Mesh_neighborhood<Mesh> Neighborhood;

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

typedef Classification::Planimetric_grid<Kernel, Face_range, Face_center_map>             Planimetric_grid;
typedef Classification::Local_eigen_analysis                                              Local_eigen_analysis;

typedef Classification::Feature::Distance_to_plane<Face_range, Face_center_map>           Distance_to_plane;
typedef Classification::Feature::Linearity                                                Linearity;
typedef Classification::Feature::Omnivariance                                             Omnivariance;
typedef Classification::Feature::Planarity                                                Planarity;
typedef Classification::Feature::Surface_variation                                        Surface_variation;
typedef Classification::Feature::Elevation<Kernel, Face_range, Face_center_map>           Elevation;
typedef Classification::Feature::Vertical_dispersion<Kernel, Face_range, Face_center_map> Dispersion;


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
  
  Iso_cuboid_3 bbox = CGAL::bounding_box (mesh.points().begin(), mesh.points().end());
  
  double grid_resolution = 0.34;
  double radius_neighbors = 1.7;
  double radius_dtm = 15.0;

  Planimetric_grid grid (faces, fc_map, bbox, grid_resolution);
  Neighborhood neighborhood (mesh);

  Local_eigen_analysis eigen (mesh, neighborhood.one_ring_neighbor_query());
  
  std::cerr << "Computing features" << std::endl;
  Feature_set features;
  Feature_handle d2p = features.add<Distance_to_plane> (faces, fc_map, eigen);
  Feature_handle lin = features.add<Linearity> (faces, eigen);
  Feature_handle omni = features.add<Omnivariance> (faces, eigen);
  Feature_handle plan = features.add<Planarity> (faces, eigen);
  Feature_handle surf = features.add<Surface_variation> (faces, eigen);
  Feature_handle disp = features.add<Dispersion> (faces, fc_map, grid,
                                                  radius_neighbors);
  Feature_handle elev = features.add<Elevation> (faces, fc_map, grid,
                                                 radius_dtm);

  std::cerr << "Setting up labels" << std::endl;
  Label_set labels;
  Label_handle ground = labels.add ("ground");
  Label_handle vege = labels.add ("vegetation");
  Label_handle roof = labels.add ("roof");

  std::cerr << "Setting weights" << std::endl;
  Classifier classifier (labels, features);
  classifier.set_weight (d2p, 6.75e-2);
  classifier.set_weight (lin, 1.19);
  classifier.set_weight (omni, 1.34e-1);
  classifier.set_weight (plan, 7.32e-1);
  classifier.set_weight (surf, 1.36e-1);
  classifier.set_weight (disp, 5.45e-1);
  classifier.set_weight (elev, 1.47e1);
  
  std::cerr << "Setting effects" << std::endl;
  classifier.set_effect (ground, d2p, Classifier::NEUTRAL);
  classifier.set_effect (ground, lin,  Classifier::PENALIZING);
  classifier.set_effect (ground, omni, Classifier::NEUTRAL);
  classifier.set_effect (ground, plan, Classifier::FAVORING);
  classifier.set_effect (ground, surf, Classifier::PENALIZING);
  classifier.set_effect (ground, disp, Classifier::NEUTRAL);
  classifier.set_effect (ground, elev, Classifier::PENALIZING);
  
  classifier.set_effect (vege, d2p,  Classifier::FAVORING);
  classifier.set_effect (vege, lin,  Classifier::NEUTRAL);
  classifier.set_effect (vege, omni, Classifier::FAVORING);
  classifier.set_effect (vege, plan, Classifier::NEUTRAL);
  classifier.set_effect (vege, surf, Classifier::NEUTRAL);
  classifier.set_effect (vege, disp, Classifier::FAVORING);
  classifier.set_effect (vege, elev, Classifier::NEUTRAL);

  classifier.set_effect (roof, d2p,  Classifier::NEUTRAL);
  classifier.set_effect (roof, lin,  Classifier::PENALIZING);
  classifier.set_effect (roof, omni, Classifier::FAVORING);
  classifier.set_effect (roof, plan, Classifier::FAVORING);
  classifier.set_effect (roof, surf, Classifier::PENALIZING);
  classifier.set_effect (roof, disp, Classifier::NEUTRAL);
  classifier.set_effect (roof, elev, Classifier::FAVORING);

  // Run classification
  std::cerr << "Classifying" << std::endl;
  
  std::vector<std::size_t> label_indices;

  CGAL::Real_timer t;
  t.start();
  Classification::classify_with_local_smoothing<CGAL::Parallel_tag>
    (faces, Face_map(), labels, classifier,
     neighborhood.one_ring_neighbor_query(),
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
