#include "Surface_mesh_item_classification.h"
#include "Color_ramp.h"

#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <CGAL/Three/Viewer_interface.h>

#include <set>
#include <stack>
#include <algorithm>
#include <boost/array.hpp>

Surface_mesh_item_classification::Surface_mesh_item_classification(Scene_surface_mesh_item* mesh)
  : m_mesh (mesh),
    m_selection (NULL),
    m_generator (NULL)
{
  m_index_color = 1;
  m_nb_scales = 5;
  m_nb_trials = 300;
  m_smoothing = 0.5;
  m_subdivisions = 16;

  backup_existing_colors_and_add_new();
  m_training = m_mesh->polyhedron()->add_property_map<face_descriptor, std::size_t>("f:training", std::size_t(-1)).first;
  m_classif = m_mesh->polyhedron()->add_property_map<face_descriptor, std::size_t>("f:label", std::size_t(-1)).first;

  m_labels.add("ground");
  m_labels.add("vegetation");
  m_labels.add("roof");
  m_labels.add("facade");
  
  m_label_colors.push_back (QColor(245, 180, 0));
  m_label_colors.push_back (QColor(0, 255, 27));
  m_label_colors.push_back (QColor(255, 0, 170));
  m_label_colors.push_back (QColor(100, 0, 255));
  
  m_sowf = new Sum_of_weighted_features (m_labels, m_features);
  m_random_forest = new Random_forest (m_labels, m_features);
}


Surface_mesh_item_classification::~Surface_mesh_item_classification()
{
  if (m_sowf != NULL)
    delete m_sowf;
  if (m_random_forest != NULL)
    delete m_random_forest;
  if (m_generator != NULL)
    delete m_generator;
}

void Surface_mesh_item_classification::backup_existing_colors_and_add_new()
{
  bool has_colors = false;
  boost::tie (m_color, has_colors) = m_mesh->polyhedron()->property_map<face_descriptor, CGAL::Color>("f:color");
  if (has_colors)
  {
    m_real_color
      = m_mesh->polyhedron()->add_property_map<face_descriptor, CGAL::Color>("f:real_color").first;
    BOOST_FOREACH(face_descriptor fd, faces(*(m_mesh->polyhedron())))
    {
      m_real_color[fd] = m_color[fd];
      m_color[fd] = CGAL::Color(128, 128, 128);
    }
  }
  else
    m_color =
      m_mesh->polyhedron()->add_property_map<face_descriptor, CGAL::Color>("f:color", CGAL::Color(128,128,128)).first;
}

bool Surface_mesh_item_classification::write_output(std::ostream& stream)
{
  // TODO
  return true;
}


void Surface_mesh_item_classification::change_color (int index)
{
  m_index_color = index;
  int index_color = index;
  if (index == 0 && m_real_color == Mesh::Property_map<face_descriptor, CGAL::Color>())
    index_color = -1;

  static Color_ramp ramp;
  ramp.build_red();

  if (index_color == -1) // item color
  {
    BOOST_FOREACH(face_descriptor fd, faces(*(m_mesh->polyhedron())))
      m_color[fd] = CGAL::Color(128,128,128);
  }
  else if (index_color == 0) // real colors
  {
    BOOST_FOREACH(face_descriptor fd, faces(*(m_mesh->polyhedron())))
      m_color[fd] = m_real_color[fd];
  }
  else if (index_color == 1) // classif
  {
    BOOST_FOREACH(face_descriptor fd, faces(*(m_mesh->polyhedron())))
    {
      QColor color (128, 128, 128);
      std::size_t c = m_classif[fd];
      
      if (c != std::size_t(-1))
        color = m_label_colors[c];

      m_color[fd] = CGAL::Color(color.red(), color.green(), color.blue());
    }
  }
  else if (index_color == 2) // training
  {
    BOOST_FOREACH(face_descriptor fd, faces(*(m_mesh->polyhedron())))
    {
      QColor color (128, 128, 128);
      std::size_t c = m_training[fd];
      std::size_t c2 = m_classif[fd];
          
      if (c != std::size_t(-1))
        color = m_label_colors[c];
      
      float div = 1;
      if (c != c2)
        div = 2;
      m_color[fd] = CGAL::Color(color.red() / div,
                                color.green() / div,
                                color.blue() / div);
    }
  }
  else
  {
    Feature_handle feature = m_features[index_color - 3];

    float max = 0.;
    BOOST_FOREACH(face_descriptor fd, faces(*(m_mesh->polyhedron())))
    {
      if (feature->value(fd) > max)
        max = feature->value(fd);
    }

    BOOST_FOREACH(face_descriptor fd, faces(*(m_mesh->polyhedron())))
    {
      float v = std::max (0.f, feature->value(fd) / max);
      m_color[fd] = CGAL::Color((unsigned char)(ramp.r(v) * 255),
                                (unsigned char)(ramp.g(v) * 255),
                                (unsigned char)(ramp.b(v) * 255));
    }
  }
}

void Surface_mesh_item_classification::compute_features ()
{
  std::cerr << "Computing features with " << m_nb_scales << " scale(s)" << std::endl;
  m_features.clear();

  if (m_generator != NULL)
    delete m_generator;

  Face_center_map fc_map (m_mesh->polyhedron());
  Face_range faces = m_mesh->polyhedron()->faces();
  
  m_generator = new Generator (m_features, *(m_mesh->polyhedron()), fc_map, m_nb_scales);
  
  delete m_sowf;
  m_sowf = new Sum_of_weighted_features (m_labels, m_features);
  delete m_random_forest;
  m_random_forest = new Random_forest (m_labels, m_features);
  std::cerr << "Features = " << m_features.size() << std::endl;
}

void Surface_mesh_item_classification::train(int predicate)
{
  if (m_features.size() == 0)
  {
    std::cerr << "Error: features not computed" << std::endl;
    return;
  }

  std::vector<std::size_t> indices (num_faces(*(m_mesh->polyhedron())), std::size_t(-1));
  BOOST_FOREACH(face_descriptor fd, faces(*(m_mesh->polyhedron())))
    indices[fd] = m_training[fd];
  
  if (predicate == 0)
  {
    m_sowf->train<Concurrency_tag>(indices, m_nb_trials);
    CGAL::Classification::classify<Concurrency_tag> (m_mesh->polyhedron()->faces(),
                                                     m_labels, *m_sowf,
                                                     indices);
  }
  else
  {
    m_random_forest->train (indices);
    CGAL::Classification::classify<Concurrency_tag> (m_mesh->polyhedron()->faces(),
                                                     m_labels, *m_random_forest,
                                                     indices);
  }

  BOOST_FOREACH(face_descriptor fd, faces(*(m_mesh->polyhedron())))
    m_classif[fd] = indices[fd];

  if (m_index_color == 1 || m_index_color == 2)
     change_color (m_index_color);
}

bool Surface_mesh_item_classification::run (int method, int predicate)
{
  if (m_features.size() == 0)
  {
    std::cerr << "Error: features not computed" << std::endl;
    return false;
  }

  if (predicate == 0)
    run (method, *m_sowf);
  else
    run (method, *m_random_forest);
  
  return true;
}
