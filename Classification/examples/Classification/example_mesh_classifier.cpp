#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classifier.h>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Mesh_neighborhood.h>
#include <CGAL/Classification/Planimetric_grid.h>
#include <CGAL/Classification/Feature/Eigen.h>
#include <CGAL/Classification/Feature/Distance_to_plane.h>
#include <CGAL/Classification/Feature/Vertical_dispersion.h>
#include <CGAL/Classification/Feature/Elevation.h>

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

typedef CGAL::Classifier<Face_range, Face_map> Classifier;

typedef CGAL::Classification::Label_handle                                            Label_handle;
typedef CGAL::Classification::Feature_handle                                          Feature_handle;

typedef CGAL::Classification::Mesh_neighborhood<Mesh> Neighborhood;

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

typedef CGAL::Classification::Planimetric_grid<Kernel, Face_range, Face_center_map>             Planimetric_grid;
typedef CGAL::Classification::Local_eigen_analysis<Kernel, Face_range, Face_center_map>         Local_eigen_analysis;

typedef CGAL::Classification::Feature::Distance_to_plane<Kernel, Face_range, Face_center_map>   Distance_to_plane;
typedef CGAL::Classification::Feature::Linearity<Kernel, Face_range, Face_center_map>           Linearity;
typedef CGAL::Classification::Feature::Omnivariance<Kernel, Face_range, Face_center_map>        Omnivariance;
typedef CGAL::Classification::Feature::Planarity<Kernel, Face_range, Face_center_map>           Planarity;
typedef CGAL::Classification::Feature::Surface_variation<Kernel, Face_range, Face_center_map>   Surface_variation;
typedef CGAL::Classification::Feature::Elevation<Kernel, Face_range, Face_center_map>           Elevation;
typedef CGAL::Classification::Feature::Vertical_dispersion<Kernel, Face_range, Face_center_map> Dispersion;


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

  Local_eigen_analysis eigen (faces, fc_map, Face_map(), neighborhood.one_ring_neighbor_query());


  Classifier classifier (faces, Face_map());
  
  std::cerr << "Computing features" << std::endl;
  Feature_handle d2p = classifier.add_feature<Distance_to_plane> (fc_map, eigen);
  Feature_handle lin = classifier.add_feature<Linearity> (eigen);
  Feature_handle omni = classifier.add_feature<Omnivariance> (eigen);
  Feature_handle plan = classifier.add_feature<Planarity> (eigen);
  Feature_handle surf = classifier.add_feature<Surface_variation> (eigen);
  Feature_handle disp = classifier.add_feature<Dispersion> (fc_map, grid,
                                                            grid_resolution,
                                                            radius_neighbors);
  Feature_handle elev = classifier.add_feature<Elevation> (fc_map, grid,
                                                           grid_resolution,
                                                           radius_dtm);

  std::cerr << "Setting weights" << std::endl;
  d2p->set_weight(6.75e-2);
  lin->set_weight(1.19);
  omni->set_weight(1.34e-1);
  plan->set_weight(7.32e-1);
  surf->set_weight(1.36e-1);
  disp->set_weight(5.45e-1);
  elev->set_weight(1.47e1);

  
  std::cerr << "Setting up labels" << std::endl;

  std::cerr << "Faces = " << faces.size() << std::endl;

  // Create label and define how features affect them
  Label_handle ground = classifier.add_label ("ground");
  ground->set_feature_effect (d2p,  CGAL::Classification::Feature::NEUTRAL);
  ground->set_feature_effect (lin,  CGAL::Classification::Feature::PENALIZING);
  ground->set_feature_effect (omni, CGAL::Classification::Feature::NEUTRAL);
  ground->set_feature_effect (plan, CGAL::Classification::Feature::FAVORING);
  ground->set_feature_effect (surf, CGAL::Classification::Feature::PENALIZING);
  ground->set_feature_effect (disp, CGAL::Classification::Feature::NEUTRAL);
  ground->set_feature_effect (elev, CGAL::Classification::Feature::PENALIZING);

  Label_handle vege = classifier.add_label ("vegetation");
  vege->set_feature_effect (d2p,  CGAL::Classification::Feature::FAVORING);
  vege->set_feature_effect (lin,  CGAL::Classification::Feature::NEUTRAL);
  vege->set_feature_effect (omni, CGAL::Classification::Feature::FAVORING);
  vege->set_feature_effect (plan, CGAL::Classification::Feature::NEUTRAL);
  vege->set_feature_effect (surf, CGAL::Classification::Feature::NEUTRAL);
  vege->set_feature_effect (disp, CGAL::Classification::Feature::FAVORING);
  vege->set_feature_effect (elev, CGAL::Classification::Feature::NEUTRAL);
  
  Label_handle roof = classifier.add_label ("roof");
  roof->set_feature_effect (d2p,  CGAL::Classification::Feature::NEUTRAL);
  roof->set_feature_effect (lin,  CGAL::Classification::Feature::PENALIZING);
  roof->set_feature_effect (omni, CGAL::Classification::Feature::FAVORING);
  roof->set_feature_effect (plan, CGAL::Classification::Feature::FAVORING);
  roof->set_feature_effect (surf, CGAL::Classification::Feature::PENALIZING);
  roof->set_feature_effect (disp, CGAL::Classification::Feature::NEUTRAL);
  roof->set_feature_effect (elev, CGAL::Classification::Feature::FAVORING);
  std::cerr << "Faces = " << faces.size() << std::endl;

  
  std::cerr << "Running with local smoothing" << std::endl;
  //  classifier.run_with_local_smoothing (neighborhood.one_ring_neighbor_query());
  classifier.run_with_one_graphcut (neighborhood.one_ring_neighbor_query(), 0.5);
  //  classifier.run();

  
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

      Label_handle label = classifier.label_of (fd);
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
