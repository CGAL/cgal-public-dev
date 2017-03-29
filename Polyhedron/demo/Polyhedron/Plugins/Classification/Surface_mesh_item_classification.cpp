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
    m_classif (NULL),
    m_trainer (NULL)
{
  m_nb_scales = 5;
  m_nb_trials = 300;
  m_smoothing = 0.5;
  m_subdivisions = 16;

  m_classif = new PSC(m_mesh->polyhedron()->faces(), Face_map());
  m_trainer = new Trainer(*m_classif);

  Label_handle ground = m_classif->add_label("ground");
  Label_handle vegetation = m_classif->add_label("vegetation");
  Label_handle roof = m_classif->add_label("roof");
  Label_handle facade = m_classif->add_label("facade");
  m_labels.push_back (std::make_pair(ground, QColor(245, 180, 0)));
  m_labels.push_back (std::make_pair(vegetation, QColor(0, 255, 27)));
  m_labels.push_back (std::make_pair(roof, QColor(255, 0, 170)));
  m_labels.push_back (std::make_pair(facade, QColor(100, 0, 255)));
}


Surface_mesh_item_classification::~Surface_mesh_item_classification()
{
  if (m_classif != NULL)
    delete m_classif;
  if (m_trainer != NULL)
    delete m_trainer;
  if (m_mesh != NULL)
    {

    }
}


bool Surface_mesh_item_classification::write_output(std::ostream& stream)
{
  // TODO
  return true;
}


void Surface_mesh_item_classification::change_color (int index)
{
  // TODO
}

void Surface_mesh_item_classification::compute_features ()
{
  CGAL_assertion (m_classif != NULL);
  m_classif->clear_features();
  
  std::cerr << "Computing features with " << m_nb_scales << " scale(s)" << std::endl;
  if (m_classif->number_of_features() != 0)
    m_classif->clear();

  // TODO
}

void Surface_mesh_item_classification::train()
{
  if (m_classif->number_of_features() == 0)
    {
      std::cerr << "Error: features not computed" << std::endl;
      return;
    }
  reset_indices();
  
  m_trainer->train(m_nb_trials);
  m_classif->run();
  m_classif->info();
}

bool Surface_mesh_item_classification::run (int method)
{
  if (m_classif->number_of_features() == 0)
    {
      std::cerr << "Error: features not computed" << std::endl;
      return false;
    }

  if (method == 0)
    {
      m_classif->run();
      m_classif->info();
    }
  // TODO: other methods  
  
  return true;
}
