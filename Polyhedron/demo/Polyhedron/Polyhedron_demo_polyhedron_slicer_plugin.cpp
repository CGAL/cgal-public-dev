#include <QtCore/qglobal.h>
#include <CGAL/AABB_intersections.h>

#include "Messages_interface.h"
#include "Scene_item_with_display_list.h"
#include "Scene_plane_item.h"
#include "Scene_polyhedron_item.h"
#include "Scene_polylines_item.h"

#include "Polyhedron_demo_plugin_interface.h"
#include "ui_Polyhedron_slicer_widget.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/bounding_box.h> 
#include <CGAL/Polyhedron_slicer_3.h>

#include "Polyhedron_type.h"

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>

#include <vector>
#include <algorithm>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;

class Polyhedron_demo_polyhedron_slicer_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  bool applicable(QAction*) const { return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())); }
  void print_message(QString message) { messages->information(message);}

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m);
  QList<QAction*> actions() const;
  Scene_polyhedron_item* get_selected_item() {
    int item_id = scene->mainSelectionIndex();
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(item_id));
    if(!poly_item) {
      int counter = 0;
      for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end && counter < 2; ++i) {
        if(Scene_polyhedron_item* tmp = qobject_cast<Scene_polyhedron_item*>(scene->item(i))) { 
          poly_item = tmp;
          counter++; 
        }
      }
      if(counter != 1) { return NULL; }
    } 
    return poly_item;
  }
  bool get_base_1_2(double bases[6]) {
    bool oks[6];
    bases[0] = ui_widget->Base_1_x->text().toDouble(&oks[0]);
    bases[1] = ui_widget->Base_1_y->text().toDouble(&oks[1]);
    bases[2] = ui_widget->Base_1_z->text().toDouble(&oks[2]);

    bases[3] = ui_widget->Base_2_x->text().toDouble(&oks[3]);
    bases[4] = ui_widget->Base_2_y->text().toDouble(&oks[4]);
    bases[5] = ui_widget->Base_2_z->text().toDouble(&oks[5]);

    bool total_ok = true;
    for(int i = 0; i < 6; ++i && total_ok) { total_ok &= oks[i];}
    return total_ok;
  }
public Q_SLOTS:
  void slicer_widget_action();
  void on_Generate_button_clicked();
  bool on_Update_plane_button_clicked();
  void plane_manipulated_frame_modified();
  void plane_destroyed();

private:
  Scene_interface* scene;
  Messages_interface* messages;
  Scene_plane_item* plane_item;
  QAction* actionSlicerWidget;

  QDockWidget* dock_widget;
  Ui::Polyhedron_slicer_widget* ui_widget;

  void intersection_of_plane_Polyhedra_3_using_AABB_wrapper(Polyhedron& mesh, 
    const std::vector<Epic_kernel::Plane_3>& planes,
    const std::vector<qglviewer::Vec>& plane_positions,
    std::list<std::vector<Epic_kernel::Point_3> >& polylines);

}; // end Polyhedron_demo_polyhedron_slicer_plugin

void Polyhedron_demo_polyhedron_slicer_plugin::init(QMainWindow* mw,
                                      Scene_interface* scene_interface,
                                      Messages_interface* m)
{
  scene = scene_interface;
  messages = m;
  actionSlicerWidget = new QAction(tr("Polyhedron slicer"), mw);
  connect(actionSlicerWidget, SIGNAL(triggered()),
          this, SLOT(slicer_widget_action()));

  dock_widget = new QDockWidget("Polyhedron slicer parameters", mw);
  dock_widget->setVisible(false); // do not show at the beginning
  dock_widget->setObjectName("PolyhedronSlicerParametersDialog");
  ui_widget = new Ui::Polyhedron_slicer_widget();

  QWidget* qw =new QWidget();
  ui_widget->setupUi(qw);
  dock_widget->setWidget(qw);
  mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

  connect(ui_widget->Generate_button,  SIGNAL(clicked()), this, SLOT(on_Generate_button_clicked()));   
  connect(ui_widget->Update_plane_button,  SIGNAL(clicked()), this, SLOT(on_Update_plane_button_clicked())); 
}

QList<QAction*> Polyhedron_demo_polyhedron_slicer_plugin::actions() const {
  return QList<QAction*>() << actionSlicerWidget;
}

void Polyhedron_demo_polyhedron_slicer_plugin::slicer_widget_action(){
  if(dock_widget != NULL && !dock_widget->isVisible()) { 
    dock_widget->show(); 

    ///// from cut plugin /////
    plane_item = new Scene_plane_item(scene);
    const Scene_interface::Bbox& bbox = scene->bbox();
    plane_item->setPosition((bbox.xmin + bbox.xmax)/2.f,
      (bbox.ymin+bbox.ymax)/2.f,
      (bbox.zmin+bbox.zmax)/2.f);
    plane_item->setNormal(0., 0., 1.);
    plane_item->setManipulatable(true);
    plane_item->setClonable(false);
    plane_item->setColor(Qt::green);
    plane_item->setName(tr("Cutting plane"));
    connect(plane_item->manipulatedFrame(), SIGNAL(modified()),
      this, SLOT(plane_manipulated_frame_modified()));
    connect(plane_item, SIGNAL(destroyed()),
      this, SLOT(plane_destroyed()));
    scene->addItem(plane_item);

    // set distance_with_planes = bbox_diagona / 30
    double diagonal = std::sqrt(
      CGAL::squared_distanceC3( bbox.xmin, bbox.ymin, bbox.zmin, bbox.xmax, bbox.ymax, bbox.zmax) );
    ui_widget->Distance_with_planes->setText(QString::number(diagonal / 30.0));

    plane_manipulated_frame_modified(); // update text boxes
  }
}

// when manipulated frame of plane is modified, update line-edits
void Polyhedron_demo_polyhedron_slicer_plugin::plane_manipulated_frame_modified() {
  qglviewer::ManipulatedFrame* mf = plane_item->manipulatedFrame();
  const qglviewer::Vec& pos = mf->position();
  ui_widget->Center_x->setText(QString::number(pos.x));
  ui_widget->Center_y->setText(QString::number(pos.y));
  ui_widget->Center_z->setText(QString::number(pos.z));

  const qglviewer::Vec& base_1 = mf->inverseTransformOf(qglviewer::Vec(1., 0., 0.));
  const qglviewer::Vec& base_2 = mf->inverseTransformOf(qglviewer::Vec(0., 1., 0.));

  ui_widget->Base_1_x->setText(QString::number(base_1.x));
  ui_widget->Base_1_y->setText(QString::number(base_1.y));
  ui_widget->Base_1_z->setText(QString::number(base_1.z));

  ui_widget->Base_2_x->setText(QString::number(base_2.x));
  ui_widget->Base_2_y->setText(QString::number(base_2.y));
  ui_widget->Base_2_z->setText(QString::number(base_2.z));
}

// when Update Plane button is clicked, update manipulated frame of plane with line-edits
bool Polyhedron_demo_polyhedron_slicer_plugin::on_Update_plane_button_clicked() {
  qglviewer::ManipulatedFrame* mf = plane_item->manipulatedFrame();
  // get center
  bool ok_1 = true, ok_2 = true, ok_3 = true;
  double center_x = ui_widget->Center_x->text().toDouble(&ok_1);
  double center_y = ui_widget->Center_y->text().toDouble(&ok_2);
  double center_z = ui_widget->Center_z->text().toDouble(&ok_3);
  if(!ok_1 || !ok_2 || !ok_3) 
  { print_message("Error: center coordinates not convertible to double."); return false; }

  // set center
  bool oldState = mf->blockSignals(true); // dont let it signal, it will invoke plane_manipulated_frame_modified otherwise
  mf->setPosition(center_x, center_y, center_z);
  mf->blockSignals(oldState);

  // get base 1 and base 2
  double bases[6];
  if(!get_base_1_2(bases)) 
  { print_message("Error: Base-1, Base-2 coordinates not convertible to double."); return false; }

  // compute other axis
  qglviewer::Vec base_1(bases[0], bases[1], bases[2]);
  qglviewer::Vec base_2(bases[3], bases[4], bases[5]);
  qglviewer::Vec other = cross(base_1, base_2);
  if(other.norm() == 0.0) { print_message("Error: collinear base vectors are not accepted!"); return false; }
  
  // set orientation
  qglviewer::Quaternion orientation_from_bases;
  orientation_from_bases.setFromRotatedBasis(base_1, base_2, other);

  oldState = mf->blockSignals(true); // dont let it signal, it will invoke plane_manipulated_frame_modified otherwise
  mf->setOrientation(orientation_from_bases);
  mf->blockSignals(oldState);

  scene->itemChanged(plane_item); // redraw
  return true;
}

// generate multiple cuts, until any cut does not intersect with bbox
void Polyhedron_demo_polyhedron_slicer_plugin::on_Generate_button_clicked()
{
  Scene_polyhedron_item* item = get_selected_item();
  if(!item) { 
    print_message("Error: There is no selected Scene_polyhedron_item!");
    return; 
  }

  if(!on_Update_plane_button_clicked()) { return; }

  // get plane position and normal
  qglviewer::ManipulatedFrame* mf = plane_item->manipulatedFrame();
  const qglviewer::Vec& pos = mf->position();
  // WARNING: due to fp arithmetic (setting quaternion based orientation from base vectors then getting plane normal back from this orientation)
  // for base vectors like: 1,0,0 - 0,1,0 we might not have exact corresponding normal vector.
  // So not using below normal but construct plane directly from bases from text boxes
  const qglviewer::Vec& n = mf->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));

  // get bases
  double bases[6];
  get_base_1_2(bases); // no need to check since we call on_Update_plane_button_clicked
  Epic_kernel::Vector_3 base_1(bases[0], bases[1], bases[2]);
  Epic_kernel::Vector_3 base_2(bases[3], bases[4], bases[5]);

  // get distance between planes
  bool to_double_ok = true;
  double distance_with_planes = ui_widget->Distance_with_planes->text().toDouble(&to_double_ok);
  if(!to_double_ok) { 
    print_message("Error: Set Distance_with_planes text box!");
    return; 
  }

  // construct a bbox for selected polyhedron
  const Scene_interface::Bbox& bbox = item->bbox();
  CGAL::Bbox_3 cgal_bbox(bbox.xmin, bbox.ymin, bbox.zmin,
    bbox.xmax, bbox.ymax, bbox.zmax);
  Polyhedron* poly = item->polyhedron();

  // continue generating planes while inside bbox
  std::vector<Epic_kernel::Plane_3> planes;
  std::vector<qglviewer::Vec> plane_positions;

  for(int dir = 1, step = 0; /* */ ; ++step) 
  {
    double distance_norm = (dir * step) * distance_with_planes;
    qglviewer::Vec new_pos = pos + (n*distance_norm);

    //Epic_kernel::Plane_3 plane(n[0], n[1],  n[2], - n * new_pos);
    Epic_kernel::Point_3 new_pos_cgal(new_pos[0], new_pos[1], new_pos[2]);
    Epic_kernel::Plane_3 plane(new_pos_cgal, new_pos_cgal + base_1, new_pos_cgal + base_2);

    if(!CGAL::do_intersect(cgal_bbox, plane)) { 
      if(dir == -1) { break; }
      std::reverse(planes.begin(), planes.end());
      std::reverse(plane_positions.begin(), plane_positions.end());
      dir = -1; // reverse direction
      step = 0; // we should skip the plane itself, and we will when continue cause ++step
      continue;
    }
    planes.push_back(plane);
    plane_positions.push_back(new_pos); 
  }
  print_message(QString("Created %1 cuts inside bbox...").arg(planes.size()));

  bool new_polyline_item_for_polylines = ui_widget->newPolylineItemCheckBox->checkState() == Qt::Checked;
  if(!new_polyline_item_for_polylines) 
  {
    Scene_polylines_item* new_polylines_item = new Scene_polylines_item();
    QTime time; time.start();
    // call algorithm and fill polylines in polylines_item
    intersection_of_plane_Polyhedra_3_using_AABB_wrapper(*poly, planes, plane_positions, new_polylines_item->polylines);
    // set names etc and print timing
    print_message( QString("Done: processed %1 cuts - generated %2 polylines in %3 ms!").
      arg(planes.size()).arg(new_polylines_item->polylines.size()).arg(time.elapsed()) );

    new_polylines_item->setName(QString("%1 with %2 cuts").
      arg(item->name()).arg(planes.size()) );
    new_polylines_item->setColor(Qt::green);
    new_polylines_item->setRenderingMode(Wireframe);
    scene->addItem(new_polylines_item);
  }
  else {
    QTime time; time.start();
    std::list<std::vector<Epic_kernel::Point_3> > polylines;
    // call algorithm and fill polylines in polylines_item
    intersection_of_plane_Polyhedra_3_using_AABB_wrapper(*poly, planes, plane_positions, polylines);
    // set names etc and print timing
    print_message( QString("Done: processed %1 cuts - generated %2 polylines in %3 ms!").
      arg(planes.size()).arg(polylines.size()).arg(time.elapsed()) );

    int counter = 0;
    for(std::list<std::vector<Epic_kernel::Point_3> >::iterator it = polylines.begin(); it != polylines.end(); ++it, ++counter) {
      Scene_polylines_item* new_polylines_item = new Scene_polylines_item();
      new_polylines_item->polylines.push_back(*it);
      new_polylines_item->setName(QString("%1 with %2 cuts %3").
        arg(item->name()).arg(planes.size()).arg(counter) );
      new_polylines_item->setColor(Qt::green);
      new_polylines_item->setRenderingMode(Wireframe);
      scene->addItem(new_polylines_item);
    }
  }
}

void Polyhedron_demo_polyhedron_slicer_plugin::plane_destroyed() {
  dock_widget->hide(); 
}

// this function assumes 'planes' are parallel
void Polyhedron_demo_polyhedron_slicer_plugin::intersection_of_plane_Polyhedra_3_using_AABB_wrapper(
  Polyhedron& poly, 
  const std::vector<Epic_kernel::Plane_3>& planes,
  const std::vector<qglviewer::Vec>& plane_positions,
  std::list<std::vector<Epic_kernel::Point_3> >& polylines) 
{
  typedef std::list<std::vector<Epic_kernel::Point_3> >::iterator Polyline_it;

  std::size_t nb_projection = 0;
  CGAL::Polyhedron_slicer_3<Polyhedron, Epic_kernel> slicer(poly);
  std::vector<qglviewer::Vec>::const_iterator plane_position_it = plane_positions.begin();
  for(std::vector<Epic_kernel::Plane_3>::const_iterator plane_it = planes.begin(); plane_it != planes.end(); ++plane_it, ++plane_position_it) 
  {
    Polyline_it last_processed_polyline = polylines.begin();
    slicer(*plane_it, std::front_inserter(polylines));
    
    double a = planes.front().a(); double b = planes.front().b(); double c = planes.front().c();
    // std::cout << "plane a, b, c: " << a << " " << b << " " << c << std::endl;
    int on_axis = (a == 0.0 && b == 0.0) ? 2 :
      (a == 0.0 && c == 0.0) ? 1 :
      (b == 0.0 && c == 0.0) ? 0 : -1;
    
    // continue if planes are not axis oriented
    if( on_axis == -1) { continue; }
    ++nb_projection;

    for(Polyline_it polyline_it = polylines.begin(); polyline_it != last_processed_polyline; ++polyline_it) 
    {
      for(std::vector<Epic_kernel::Point_3>::iterator point_it = polyline_it->begin();
        point_it != polyline_it->end(); ++ point_it) 
      {
        *point_it = Epic_kernel::Point_3(
          on_axis != 0 ? point_it->x() : plane_position_it->x,
          on_axis != 1 ? point_it->y() : plane_position_it->y,
          on_axis != 2 ? point_it->z() : plane_position_it->z);
      }
    }
  }
  print_message(QString("%1 axis aligned planes are found, and points are projected...").arg(nb_projection));
}
Q_EXPORT_PLUGIN2(Polyhedron_demo_polyhedron_slicer_plugin, Polyhedron_demo_polyhedron_slicer_plugin)

#include "Polyhedron_demo_polyhedron_slicer_plugin.moc"
