#ifndef SURFACE_MESH_ITEM_CLASSIFICATION_H
#define SURFACE_MESH_ITEM_CLASSIFICATION_H

//#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#define CGAL_CLASSIFICATION_VERBOSE

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Classifier.h>
#include <CGAL/Classification/Trainer.h>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Feature/Vertical_dispersion.h>
#include <CGAL/Classification/Feature/Elevation.h>
#include <CGAL/Classification/Feature/Verticality.h>
#include <CGAL/Classification/Feature/Distance_to_plane.h>
#include <CGAL/Classification/Feature/Hsv.h>
#include <CGAL/Classification/Feature/Echo_scatter.h>
#include <CGAL/Classification/Feature/Eigen.h>

#include "Scene_surface_mesh_item.h"
#include "Item_classification_base.h"
#include "Kernel_type.h"

#include <iostream>

class Surface_mesh_item_classification : public Item_classification_base
{
public:

  typedef Scene_surface_mesh_item::SMesh Mesh;
  typedef Scene_surface_mesh_item::Point Point;
  typedef typename Mesh::Face_range Face_range;
  typedef typename boost::graph_traits<Mesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
  typedef CGAL::Identity_property_map<face_descriptor> Face_map;

  typedef CGAL::Classifier<Face_range, Face_map> Classifier;
  typedef CGAL::Classification::Mesh_neighborhood<Mesh> Neighborhood;
  typedef CGAL::Classification::Label_handle   Label_handle;
  typedef CGAL::Classification::Feature_handle Feature_handle;
  
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

public:
  
  Surface_mesh_item_classification(Scene_surface_mesh_item* mesh);
  ~Surface_mesh_item_classification();

  CGAL::Three::Scene_item* item() { return m_mesh; }
  void erase_item() { m_mesh = NULL; }

  void compute_features ();
  bool features_computed() const { return (m_classif->number_of_features() != 0); }
  std::size_t number_of_features() const { return m_classif->number_of_features(); }
  Feature_handle feature(std::size_t i) { return m_classif->feature(i); }

  void add_new_label (const char* name, const QColor& color)
  {
    m_labels.push_back (std::make_pair (m_classif->add_label(name),
                                       color));
  }
  void remove_label (const char* name)
  {
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_labels[i].first->name() == name)
        {
          m_classif->remove_label (m_labels[i].first);
          m_labels.erase (m_labels.begin() + i);
          break;
        }
  }
  
  void add_selection_to_training_set (const char* name)
  {
    // TODO
  }
  void reset_training_sets()
  {
    m_trainer->reset_inlier_sets();
  }
  void validate_selection ()
  {
    // TODO
  }
  void train();
  bool run (int method);
  
  void change_color (int index) = 0;
  void fill_display_combo_box (QComboBox* cb, QComboBox* cb1) const
  {
    // TODO
  }
  void generate_one_item_per_label(std::vector<CGAL::Three::Scene_item*>& items,
                                           const char* name) const
  {
    // TODO
  }

  bool write_output(std::ostream& out) = 0;
  void save_config(const char* filename)
  {
    // TODO
  }
  void load_config(const char* filename)
  {
    // TODO
  }


protected:

  Scene_surface_mesh_item* m_mesh;

  Classifier* m_classif;
  Trainer* m_trainer;
  
};




#endif // SURFACE_MESH_ITEM_CLASSIFICATION_H
