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

typedef Classif::Planimetric_grid<Kernel, Face_range, Face_center_map>             Planimetric_grid;
typedef Classif::Local_eigen_analysis                                              Local_eigen_analysis;

typedef Classif::Feature::Distance_to_plane<Face_range, Face_center_map>           Distance_to_plane;
typedef Classif::Feature::Linearity                                                Linearity;
typedef Classif::Feature::Omnivariance                                             Omnivariance;
typedef Classif::Feature::Planarity                                                Planarity;
typedef Classif::Feature::Surface_variation                                        Surface_variation;
typedef Classif::Feature::Elevation<Kernel, Face_range, Face_center_map>           Elevation;
typedef Classif::Feature::Vertical_dispersion<Kernel, Face_range, Face_center_map> Dispersion;


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
                                                  grid_resolution,
                                                  radius_neighbors);
  Feature_handle elev = features.add<Elevation> (faces, fc_map, grid,
                                                 grid_resolution,
                                                 radius_dtm);

  std::cerr << "Setting up labels" << std::endl;
  Label_set labels;
  Label_handle ground = labels.add ("ground");
  Label_handle vege = labels.add ("vegetation");
  Label_handle roof = labels.add ("roof");

  std::cerr << "Setting weights" << std::endl;
  Classification_predicate predicate (labels, features);
  predicate.set_weight (d2p, 6.75e-2);
  predicate.set_weight (lin, 1.19);
  predicate.set_weight (omni, 1.34e-1);
  predicate.set_weight (plan, 7.32e-1);
  predicate.set_weight (surf, 1.36e-1);
  predicate.set_weight (disp, 5.45e-1);
  predicate.set_weight (elev, 1.47e1);
  
  std::cerr << "Setting effects" << std::endl;
  predicate.set_effect (ground, d2p, Classification_predicate::NEUTRAL);
  predicate.set_effect (ground, lin,  Classification_predicate::PENALIZING);
  predicate.set_effect (ground, omni, Classification_predicate::NEUTRAL);
  predicate.set_effect (ground, plan, Classification_predicate::FAVORING);
  predicate.set_effect (ground, surf, Classification_predicate::PENALIZING);
  predicate.set_effect (ground, disp, Classification_predicate::NEUTRAL);
  predicate.set_effect (ground, elev, Classification_predicate::PENALIZING);
  
  predicate.set_effect (vege, d2p,  Classification_predicate::FAVORING);
  predicate.set_effect (vege, lin,  Classification_predicate::NEUTRAL);
  predicate.set_effect (vege, omni, Classification_predicate::FAVORING);
  predicate.set_effect (vege, plan, Classification_predicate::NEUTRAL);
  predicate.set_effect (vege, surf, Classification_predicate::NEUTRAL);
  predicate.set_effect (vege, disp, Classification_predicate::FAVORING);
  predicate.set_effect (vege, elev, Classification_predicate::NEUTRAL);

  predicate.set_effect (roof, d2p,  Classification_predicate::NEUTRAL);
  predicate.set_effect (roof, lin,  Classification_predicate::PENALIZING);
  predicate.set_effect (roof, omni, Classification_predicate::FAVORING);
  predicate.set_effect (roof, plan, Classification_predicate::FAVORING);
  predicate.set_effect (roof, surf, Classification_predicate::PENALIZING);
  predicate.set_effect (roof, disp, Classification_predicate::NEUTRAL);
  predicate.set_effect (roof, elev, Classification_predicate::FAVORING);

  // Run classification
  std::cerr << "Classifying" << std::endl;
  
  std::vector<std::size_t> label_indices;

  CGAL::Real_timer t;
  t.start();
  Classif::classify_with_local_smoothing<CGAL::Parallel_tag>
    (faces, Face_map(), labels, predicate,
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
