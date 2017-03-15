#ifndef POINT_SET_ITEM_CLASSIFICATION_H
#define POINT_SET_ITEM_CLASSIFICATION_H

//#define CGAL_DO_NOT_USE_BOYKOV_KOLMOGOROV_MAXFLOW_SOFTWARE
#define CGAL_CLASSIFICATION_VERBOSE

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Point_set_classifier.h>
#include <CGAL/Classification/Trainer.h>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Feature/Vertical_dispersion.h>
#include <CGAL/Classification/Feature/Elevation.h>
#include <CGAL/Classification/Feature/Verticality.h>
#include <CGAL/Classification/Feature/Distance_to_plane.h>
#include <CGAL/Classification/Feature/Hsv.h>
#include <CGAL/Classification/Feature/Echo_scatter.h>
#include <CGAL/Classification/Feature/Eigen.h>

#include "Scene_points_with_normal_item.h"
#include "Item_classification_base.h"
#include "Polyhedron_type_fwd.h"
#include "Kernel_type.h"
#include "Point_set_3.h"

#include <iostream>


// This class represents a point set in the OpenGL scene
class Point_set_item_classification : public Item_classification_base
{
 public:
  typedef Kernel::Point_3 Point_3;
  typedef Kernel::Vector_3 Vector_3;
  typedef CGAL::Classification::RGB_Color Color;
  
  typedef Point_set::Point_map Point_map;
  typedef Point_set::Vector_map Vector_map;

  typedef CGAL::Point_set_classifier<Kernel, Point_set, Point_map>                     PSC;
  typedef CGAL::Classification::Trainer<Point_set, Point_map>        Trainer;
  typedef CGAL::Classification::Label_handle                                               Label_handle;
  typedef CGAL::Classification::Feature_handle                                             Feature_handle;
  typedef CGAL::Classification::Feature::Vertical_dispersion<Kernel, Point_set, Point_map> Dispersion;
  typedef CGAL::Classification::Feature::Elevation<Kernel, Point_set, Point_map>           Elevation;
  typedef CGAL::Classification::Feature::Verticality<Kernel, Point_set, Point_map>         Verticality;
  typedef CGAL::Classification::Feature::Distance_to_plane<Kernel, Point_set, Point_map>   Distance_to_plane;


 public:
  
  Point_set_item_classification(Scene_points_with_normal_item* points);
  ~Point_set_item_classification();

  CGAL::Three::Scene_item* item() { return m_points; }
  void erase_item() { m_points = NULL; }

  void compute_features ();
  bool features_computed() const { return (m_psc->number_of_features() != 0); }
  std::size_t number_of_features() const { return m_psc->number_of_features(); }  
  Feature_handle feature(std::size_t i) { return m_psc->feature(i); }
  
  void add_new_label (const char* name, const QColor& color)
  {
    m_labels.push_back (std::make_pair (m_psc->add_label(name),
                                       color));
  }
  void remove_label (const char* name)
  {
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_labels[i].first->name() == name)
        {
          m_psc->remove_label (m_labels[i].first);
          m_labels.erase (m_labels.begin() + i);
          break;
        }
  }

  void add_selection_to_training_set (const char* name)
  {
    Label_handle label = get_label (name);

    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
      m_psc->set_label_of(*it, label);

    m_trainer->set_inliers(label,
                           boost::make_iterator_range(m_points->point_set()->first_selected(),
                                                      m_points->point_set()->end()));

    m_points->resetSelection();
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void reset_training_sets()
  {
    m_trainer->reset_inlier_sets();
  }
  void validate_selection ()
  {
    for (Point_set::const_iterator it = m_points->point_set()->first_selected();
         it != m_points->point_set()->end(); ++ it)
      {
        Label_handle t = m_psc->label_of(*it);
        m_trainer->set_inlier (t, *it);
      }

    m_points->resetSelection();
    if (m_index_color == 1 || m_index_color == 2)
      change_color (m_index_color);
  }
  void train();
  bool run (int method);

  void change_color (int index);
  void change_label_color (const char* name, const QColor& color)
  {
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_labels[i].first->name() == name)
        {
          m_labels[i].second = color;
          break;
        }
  }
  void fill_display_combo_box (QComboBox* cb, QComboBox* cb1) const
  {
    for (std::size_t i = 0; i < m_psc->number_of_features(); ++ i)
      {
        std::size_t scale = m_psc->scale_of_feature(m_psc->feature(i));
        std::ostringstream oss;
        oss << "Feature " << m_psc->feature(i)->name() << "_" << scale;
        cb->addItem (oss.str().c_str());
        cb1->addItem (oss.str().c_str());
      }
  }
  void generate_one_item_per_label(std::vector<CGAL::Three::Scene_item*>& items,
                                   const char* name) const
  {
    std::map<Label_handle, Scene_points_with_normal_item*> map_labels;
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      {
        Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item;
        new_item->setName (QString("%1 (%2)").arg(name).arg(m_labels[i].first->name().c_str()));
        new_item->setColor (m_labels[i].second);
        map_labels[m_labels[i].first] = new_item;
        items.push_back (new_item);
      }

    for (Point_set::const_iterator it = m_points->point_set()->begin();
         it != m_points->point_set()->end(); ++ it)
      {
        Label_handle c = m_psc->label_of(*it);
        if (c != Label_handle())
          map_labels[c]->point_set()->insert (m_points->point_set()->point(*it));
      }
  }
  
  bool write_output(std::ostream& out);
  void save_config(const char* filename)
  {
    if (m_psc->number_of_features() == 0)
      {
        std::cerr << "Error: features not computed" << std::endl;
        return;
      }

    std::ofstream f (filename);
    m_psc->save_configuration (f);
  }
  void load_config(const char* filename)
  {
    if (m_psc->number_of_features() != 0)
      m_psc->clear();

    reset_indices();
    
    bool normals = m_points->point_set()->has_normal_map();
    bool colors = (m_color != Point_set::Property_map<Color>());
    Point_set::Property_map<boost::uint8_t> echo_map;
    bool echo;
    boost::tie (echo_map, echo) = m_points->point_set()->template property_map<boost::uint8_t>("echo");

    std::ifstream f (filename);
    if (!normals && !colors && !echo)
      m_psc->load_configuration (f);
    else if (!normals && !colors && echo)
      m_psc->load_configuration (f, CGAL::Default(), CGAL::Default(), echo_map);
    else if (!normals && colors && !echo)
      m_psc->load_configuration (f, CGAL::Default(), m_color);
    else if (!normals && colors && echo)
      m_psc->load_configuration (f, CGAL::Default(), m_color, echo_map);
    else if (normals && !colors && !echo)
      m_psc->load_configuration (f, m_points->point_set()->normal_map());
    else if (normals && !colors && echo)
      m_psc->load_configuration (f, m_points->point_set()->normal_map(), CGAL::Default(), echo_map);
    else if (normals && colors && !echo)
      m_psc->load_configuration (f, m_points->point_set()->normal_map(), m_color);
    else
      m_psc->load_configuration (f, m_points->point_set()->normal_map(), m_color, echo_map);
    
    std::vector<std::pair<Label_handle, QColor> > new_labels;
    for (std::size_t i = 0; i < m_psc->number_of_labels(); ++ i)
      {
        Label_handle t = m_psc->label(i);
        QColor color (192 + rand() % 60,
                      192 + rand() % 60,
                      192 + rand() % 60);

        for (std::size_t j = 0; j < m_labels.size(); ++ j)
          if (t->name() == m_labels[j].first->name())
            {
              color = m_labels[j].second;
              break;
            }

        new_labels.push_back (std::make_pair (t, color));
      }
    m_labels.swap (new_labels);
  }


  int real_index_color() const;
  void reset_indices();
  void backup_existing_colors_and_add_new();
  void reset_colors();
  Label_handle get_label (const char* name)
  {
    for (std::size_t i = 0; i < m_labels.size(); ++ i)
      if (m_labels[i].first->name() == name)
        return m_labels[i].first;
    return Label_handle();
  }

 private:
  
  Scene_points_with_normal_item* m_points;

  Point_set::Property_map<unsigned char> m_red;
  Point_set::Property_map<unsigned char> m_green;
  Point_set::Property_map<unsigned char> m_blue;
  Point_set::Property_map<Color> m_color;
  
  PSC* m_psc;
  Trainer* m_trainer;
  
  int m_index_color;
  
}; // end class Point_set_item_classification




#endif // POINT_SET_ITEM_CLASSIFICATION_H
