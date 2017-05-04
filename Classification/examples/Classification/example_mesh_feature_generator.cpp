#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#define CGAL_CLASSIFICATION_VERBOSE

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;

typedef CGAL::Surface_mesh<Point> Surface_mesh;

typedef CGAL::Polyhedron_3<Kernel,
                           CGAL::Polyhedron_items_with_id_3,
                           CGAL::HalfedgeDS_vector> Polyhedron;

template <typename FaceGraph, typename Point>
struct Face_graph_face_to_center_property_map
{
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor key_type;
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
      points.push_back (get(get(CGAL::vertex_point, *(map.mesh)), v));
    }
    return CGAL::centroid (points.begin(), points.end());
  }
};

template <typename Mesh>
void run (Mesh& mesh)
{
  typedef typename boost::graph_traits<Mesh>::face_iterator face_iterator;
  typedef typename boost::graph_traits<Mesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typedef CGAL::Iterator_range<face_iterator> Face_range;

  typedef CGAL::Identity_property_map<face_descriptor> Face_map;

  namespace Classif = CGAL::Classification;

  typedef Classif::Sum_of_weighted_features_predicate Classification_predicate;

  typedef Classif::Label_handle                                            Label_handle;
  typedef Classif::Label_set                                               Label_set;
  typedef Classif::Feature_set                                             Feature_set;

  typedef Face_graph_face_to_center_property_map<Mesh, Point> Face_center_map;

  typedef Classif::Mesh_feature_generator<Kernel, Mesh, Face_center_map>             Generator;

  CGAL::Real_timer ttot;
  ttot.start();
  
  Face_range range (faces(mesh));
  
  std::cerr << " * " << range.size() << " faces" << std::endl;
  std::cerr << "Computing useful structures" << std::endl;

  Face_center_map fc_map (&mesh);
  
  std::cerr << "Computing features" << std::endl;
  Feature_set features;
  Generator generator (features, mesh, fc_map, 5);

  std::cerr << "Setting up labels" << std::endl;
  Label_set labels;
  Label_handle ground = labels.add ("ground");
  Label_handle vege = labels.add ("vegetation");
  Label_handle roof = labels.add ("roof");

  ttot.stop();
  std::cerr << "Features computed in " << ttot.time() << "s" << std::endl;

  // Run classification
  std::cerr << "Classifying" << std::endl;
  
  std::vector<std::size_t> label_indices;

  Classification_predicate predicate (labels, features);
  
  CGAL::Real_timer t;
  t.start();
  Classif::classify_with_local_smoothing<CGAL::Parallel_tag>
    (range, Face_map(), labels, predicate,
     generator.neighborhood().one_ring_neighbor_query(),
     label_indices);
  t.stop();
  std::cerr << "Classification with graphcut performed in " << t.time() << " second(s)" << std::endl;
  
  std::ofstream fout("out.off");
  fout << "COFF" << std::endl
       << num_vertices(mesh) << " " << num_faces(mesh) << " 0" << std::endl;
  BOOST_FOREACH (vertex_descriptor vd, vertices(mesh))
  {
    fout << get(get(CGAL::vertex_point, mesh), vd) << std::endl;
  }
  BOOST_FOREACH (face_descriptor fd, range)
  {
    fout << "3";
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(halfedge(fd, mesh), mesh))
    {
      fout << " " << get(get(boost::vertex_index, mesh), vd);
    }

    Label_handle label = labels[label_indices[get(get(CGAL::face_index, mesh), fd)]];
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
}

int main (int argc, char** argv)
{
  std::string filename (argc > 1 ? argv[1] : "data/example.off");
  bool use_polyhedron = (argc > 2 && std::string(argv[2]) == std::string("-p"));
  
  if (use_polyhedron)
  {
    std::ifstream in (filename.c_str());
    Polyhedron mesh;
  
    std::cerr << "Reading polyhedron" << std::endl;
    in >> mesh;
    
    std::size_t id = 0;
    BOOST_FOREACH (typename boost::graph_traits<Polyhedron>::face_descriptor fd, faces(mesh))
      put (get(CGAL::face_index, mesh), fd, id ++);
    id = 0;
    BOOST_FOREACH (typename boost::graph_traits<Polyhedron>::vertex_descriptor vd, vertices(mesh))
      put (get(boost::vertex_index, mesh), vd, id ++);

    run (mesh);
  }
  else
  {
    std::ifstream in (filename.c_str());
    Surface_mesh mesh;
  
    std::cerr << "Reading surface mesh" << std::endl;
    in >> mesh;
    
    run (mesh);
  }

  return EXIT_SUCCESS;
}
