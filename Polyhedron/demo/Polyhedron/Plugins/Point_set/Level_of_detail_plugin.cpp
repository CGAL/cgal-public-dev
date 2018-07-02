#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polyhedron_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/Real_timer.h>
#include <CGAL/Memory_sizer.h>

#include <CGAL/Level_of_detail_include.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include "ui_Level_of_detail_plugin.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>
#include <QMessageBox>

using namespace CGAL::Three;

typedef CGAL::Level_of_detail::Level_of_detail<Kernel, Point_set, Point_set::Point_map> LOD;


class Polyhedron_demo_level_of_detail_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

private:
  QAction* actionLOD;
  
public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    scene = scene_interface;
    mw = mainWindow;
    actionLOD = new QAction(tr("Level Of Detail"), mainWindow);
    actionLOD->setObjectName("actionLOD");
    autoConnectActions();
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionLOD;
  }
  
  //! Applicable if the currently selected item is a
  //! points_with_normal_item.
  bool applicable(QAction*) const {
    Scene_points_with_normal_item* item =
      qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if (!item)
      return false;
    if (!item->point_set()->property_map<int>("label").second)
      return false;
    return true;
  }

public Q_SLOTS:
  void on_actionLOD_triggered();
}; // end Polyhedron_demo_level_of_detail_plugin

class Polyhedron_demo_lod_dialog : public QDialog, private Ui::LevelOfDetailDialog
{
  Q_OBJECT
public:
  
  Polyhedron_demo_lod_dialog(QWidget * /*parent*/ = 0)
  {
    setupUi(this);
  }

  double scale() const { return scaleDoubleSpinBox->value(); }
  double tolerance() const { return toleranceDoubleSpinBox->value(); }

  QString ground_indices() const { return groundLineEdit->text(); }
  QString boundary_indices() const { return buildingBoundariesLineEdit->text(); }
  QString interior_indices() const { return buildingInteriorLineEdit->text(); }
  
  QTextEdit* comment_section() { return textEdit; }
};

struct Semantic_map_from_labels
{
  typedef Point_set::Index key_type;
  typedef CGAL::Level_of_detail::Semantic_label value_type;
  typedef CGAL::Level_of_detail::Semantic_label reference;
  typedef boost::readable_property_map_tag category;

  Point_set* points;
  Point_set::Property_map<int> label_map;
  boost::shared_ptr<std::map<int, CGAL::Level_of_detail::Semantic_label> > map_l2sl;

  Semantic_map_from_labels (Point_set* points) : points (points)
                                               , map_l2sl (new std::map<int, CGAL::Level_of_detail::Semantic_label>())
  {
    label_map = points->property_map<int>("label").first;
  }

  friend value_type get (const Semantic_map_from_labels& map, const key_type& key)
  {
    int l = map.label_map[key];

    typename std::map<int, CGAL::Level_of_detail::Semantic_label>::const_iterator
      found = map.map_l2sl->find(l);
    if (found == map.map_l2sl->end())
      return CGAL::Level_of_detail::Semantic_label::UNASSIGNED;

    return found->second;
  }
};

void Polyhedron_demo_level_of_detail_plugin::on_actionLOD_triggered()
{
  
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    // Gets point set
    Point_set* points = item->point_set();
    if(points == NULL)
        return;
    Polyhedron_demo_lod_dialog dialog;
    dialog.comment_section()->setText(item->comments().c_str());
    if(!dialog.exec())
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);
    QApplication::processEvents();
    
    CGAL::Real_timer task_timer; task_timer.start();

    LOD::Parameters parameters;
    parameters.verbose() = true;
    parameters.scale() = dialog.scale();
    parameters.epsilon() = dialog.tolerance();
    parameters.update_dependent();
    
    LOD lod (*points, points->point_map(), parameters);
    Semantic_map_from_labels semantic_map (points);

    std::istringstream gi (dialog.ground_indices().toStdString());
    std::istringstream bi (dialog.boundary_indices().toStdString());
    std::istringstream ii (dialog.interior_indices().toStdString());
    int idx;
    while (gi >> idx)
    {
      std::cerr << idx << " is ground" << std::endl;
      semantic_map.map_l2sl->insert
        (std::make_pair (idx, CGAL::Level_of_detail::Semantic_label::GROUND));
    }
    while (bi >> idx)
    {
      std::cerr << idx << " is building boundary" << std::endl;
      semantic_map.map_l2sl->insert
        (std::make_pair (idx, CGAL::Level_of_detail::Semantic_label::BUILDING_BOUNDARY));
    }
    while (ii >> idx)
    {
      std::cerr << idx << " is building interior" << std::endl;
      semantic_map.map_l2sl->insert
        (std::make_pair (idx, CGAL::Level_of_detail::Semantic_label::BUILDING_INTERIOR));
    }
    
    typedef CGAL::Level_of_detail::Partition_element<Kernel, CGAL::Polygon_2<Kernel> > LOD_partition_face_2;
    CGAL::Level_of_detail::Visibility_from_classification_property_map_2<LOD_partition_face_2, Kernel,
                                                                         Point_set, Point_set::Point_map,
                                                                         Semantic_map_from_labels>
      visibility_map (*points, points->point_map(), semantic_map);
    
    lod.build (semantic_map, visibility_map);

    LOD::Lod_0 lod_0;
    lod.get_lod(lod_0);
    Scene_polyhedron_item* lod0_item = new Scene_polyhedron_item (Polyhedron());
    Polyhedron& lod0_poly = * const_cast<Polyhedron*>(lod0_item->polyhedron());
    CGAL::copy_face_graph(lod_0.mesh(), lod0_poly);
    lod0_item->setName(tr("%1 (LOD0)").arg(item->name()));
    scene->addItem(lod0_item);
    
    LOD::Lod_1 lod_1;
    lod.get_lod(lod_1);
    Scene_polyhedron_item* lod1_item = new Scene_polyhedron_item (Polyhedron());
    Polyhedron& lod1_poly = * const_cast<Polyhedron*>(lod1_item->polyhedron());
    CGAL::copy_face_graph(lod_1.mesh(), lod1_poly);
    lod1_item->setName(tr("%1 (LOD1)").arg(item->name()));
    scene->addItem(lod1_item);

    task_timer.stop();

    std::cerr << "Done in " << task_timer.time() << " second(s)" << std::endl;
    
    QApplication::restoreOverrideCursor();

  }
}


#include "Level_of_detail_plugin.moc"
