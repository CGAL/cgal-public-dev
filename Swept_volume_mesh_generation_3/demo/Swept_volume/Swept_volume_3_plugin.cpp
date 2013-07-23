#define CGAL_MESH_3_VERBOSE 1

#include "config.h"
#include "Scene_swept_volume_3_item_config.h"

#include <CGAL_demo/Plugin_helper.h>
#include <CGAL_demo/Plugin_interface.h>
#include <CGAL_demo/Messages_interface.h>
#include "ui_Meshing_dialog.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QMessageBox>
#include <QInputDialog>
#include <QFileDialog>
#include <QReal_timer>

#include <Scene_swept_volume_3_item.h>
#include <Scene_c3t3_item.h>
#include <Meshing_thread.h>
#include <Swept_volume_3_mesh_function.h>
#include <SV/Mesh_domain_3.h>

#include <iostream>
#include <fstream>
#include <math.h>

// Constants
const QColor default_mesh_color(45,169,70);


inline 
Meshing_thread* cgal_code_mesh_3(Swept_volume_3* p_swept_volume_3,
                                 const double facet_angle,
                                 const double facet_sizing,
                                 const double facet_approx,
                                 const double tets_sizing,
                                 const double tets_shape){

  std::cerr << "BEGIN  cgal_code_mesh_3 Swept_volume_3"  << std::endl;
  typedef SV::Mesh_domain_3<Swept_volume_3> Mesh_domain; 
  typedef Swept_volume_3_mesh_function<Mesh_domain> Swept_volume_3_mesh_function;
  
  if( NULL == p_swept_volume_3 ) { return NULL; }
  

  Mesh_domain*     p_domain   = new Mesh_domain(*p_swept_volume_3);
  Scene_c3t3_item* p_new_item = new Scene_c3t3_item();
  
  Mesh_parameters param;
  param.facet_angle = facet_angle;
  param.facet_sizing = facet_sizing;
  param.facet_approx = facet_approx;
  param.tet_sizing = tets_sizing;
  param.tet_shape = tets_shape;
  
  Swept_volume_3_mesh_function* p_mesh_function = new Swept_volume_3_mesh_function(p_new_item->c3t3(), p_domain, param);

  Meshing_thread* thread = new Meshing_thread(p_mesh_function, p_new_item);
  return thread;
}

inline double get_approximate(double d, int precision, int& decimals){
  if ( d<0 ) { return 0; }
  
  double i = std::pow(10.,precision-1);
  
  decimals = 0;
  while ( d > i*10 ) { d = d/10.; ++decimals; }
  while ( d < i ) { d = d*10.; --decimals; }
  
  return std::floor(d)*std::pow(10.,decimals);
}



class Swept_volume_3_plugin : 
  public QObject,
  protected Plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Plugin_interface);
public:
  Swept_volume_3_plugin();
  
  virtual void init(QMainWindow* mainWindow, 
                    Scene_interface* scene_interface,
                    Messages_interface* msg_interface)
  {
    this->scene = scene_interface;
    
    this->mw = mainWindow;
    //actionSwept_volume_3 = this->getActionFromMainWindow(mw, "Mesh Swept_volume_3");
    actionSwept_volume_3 = new QAction(tr("Mesh Swept_volume_3"), mw);
    if(actionSwept_volume_3)
    {
      connect(actionSwept_volume_3, SIGNAL(triggered()), this, SLOT(mesh_3()));
    }
    this->msg = msg_interface;
  }

  virtual QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionSwept_volume_3;
  }
  
public slots:
  void mesh_3();
  void meshing_done(Meshing_thread* t);
  void status_report(QString str);
  
private:
  void launch_thread(Meshing_thread* mesh_thread);
  void treat_result(Scene_item& source_item, Scene_c3t3_item& result_item) const;

private:
  QAction* actionSwept_volume_3;
  Messages_interface* msg;
  QMessageBox* message_box_;
  Scene_item* source_item_;
  
}; // end class Swept_volume_3_plugin

Swept_volume_3_plugin::
Swept_volume_3_plugin()
  : actionSwept_volume_3(NULL)
  , msg(NULL)
  , message_box_(NULL)
  , source_item_(NULL)
{
}


void Swept_volume_3_plugin::mesh_3()
{
  std::cerr << "BEGIN Swept_volume_3_plugin::mesh_3 " << std::endl;  
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  // -----------------------------------
  // Check if selected item is meshable
  // -----------------------------------
  Scene_swept_volume_3_item* sv3_item = 
    qobject_cast<Scene_swept_volume_3_item*>(scene->item(index));

  // Get item
  Scene_item* item = NULL;
  if( NULL != sv3_item ) { item = sv3_item; }

  if ( NULL == item )
  {
    QMessageBox::warning(mw,tr(""),
                          tr("Selected object can't be meshed")); 
    return;
  }
  
  // -----------------------------------
  // Create Mesh dialog
  // -----------------------------------
  QDialog dialog(mw);
  Ui::Meshing_dialog ui;
  ui.setupUi(&dialog);
  connect(ui.buttonBox, SIGNAL(accepted()),
          &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()),
          &dialog, SLOT(reject()));

  // Connect checkboxes to spinboxes
  connect(ui.noApprox, SIGNAL(toggled(bool)),
          ui.approx,   SLOT(setEnabled(bool)));

  connect(ui.noFacetSizing, SIGNAL(toggled(bool)),
          ui.facetSizing,   SLOT(setEnabled(bool)));

  connect(ui.noAngle,    SIGNAL(toggled(bool)),
          ui.facetAngle, SLOT(setEnabled(bool)));

  connect(ui.noTetSizing, SIGNAL(toggled(bool)),
          ui.tetSizing,   SLOT(setEnabled(bool)));

  connect(ui.noTetShape, SIGNAL(toggled(bool)),
          ui.tetShape,   SLOT(setEnabled(bool)));

  // Set default parameters
  Scene_interface::Bbox bbox = item->bbox();
  ui.objectName->setText(item->name());
  ui.objectNameSize->setText(tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
                             .arg(bbox.width(),0,'g',3)
                             .arg(bbox.height(),0,'g',3)
                             .arg(bbox.depth(),0,'g',3) );
  
  double diag = bbox.diagonal_length();
  int decimals = 0;
  double sizing_default = get_approximate(diag * 0.05, 2, decimals);
  ui.facetSizing->setDecimals(-decimals+2);
  ui.facetSizing->setSingleStep(std::pow(10.,decimals));
  ui.facetSizing->setRange(diag * 10e-6, // min
                           diag); // max
  ui.facetSizing->setValue(sizing_default); // default value
  
  ui.tetSizing->setDecimals(-decimals+2);
  ui.tetSizing->setSingleStep(std::pow(10.,decimals));
  ui.tetSizing->setRange(diag * 10e-6, // min
                         diag); // max
  ui.tetSizing->setValue(sizing_default); // default value
  
  double approx_default = get_approximate(diag * 0.005, 2, decimals);
  ui.approx->setDecimals(-decimals+2);
  ui.approx->setSingleStep(std::pow(10.,decimals));
  ui.approx->setRange(diag * 10e-7, // min
                      diag); // max
  ui.approx->setValue(approx_default);

  // -----------------------------------
  // Get values
  // -----------------------------------
  int i = dialog.exec();
  if( i == QDialog::Rejected ) { return; }
  
  // 0 means parameter is not considered
  const double angle = !ui.noAngle->isChecked() ? 0 : ui.facetAngle->value();
  const double approx = !ui.noApprox->isChecked() ? 0 : ui.approx->value();
  const double facet_sizing = !ui.noFacetSizing->isChecked() ? 0 : ui.facetSizing->value();
  const double radius_edge = !ui.noTetShape->isChecked() ? 0 : ui.tetShape->value();
  const double tet_sizing = !ui.noTetSizing->isChecked() ? 0  : ui.tetSizing->value();
    

  // -----------------------------------
  // Dispatch mesh process
  // -----------------------------------
  QApplication::setOverrideCursor(Qt::WaitCursor);
  
  Meshing_thread* thread = NULL;
  
  // Swept_volume_3
  if ( NULL != sv3_item )
  {
    Swept_volume_3* p_swept_volume_3 = sv3_item->swept_volume_3();
    if( NULL == p_swept_volume_3 )
      {
        QMessageBox::critical(mw,tr(""),tr("ERROR: no data in selected item")); 
        return;
      }
    
    thread = cgal_code_mesh_3(p_swept_volume_3,angle, facet_sizing, approx, tet_sizing, radius_edge);
  }

  if ( NULL == thread )
  {
    QMessageBox::critical(mw,tr(""),tr("ERROR: no thread created")); 
    return;
  }
  
  std::cerr << "GOING TO LAUNCH THREAD NOW" << std::endl;
  // Launch thread
  source_item_ = item;
  launch_thread(thread);
  QApplication::restoreOverrideCursor();
}


void
Swept_volume_3_plugin::
launch_thread(Meshing_thread* mesh_thread)
{
  // -----------------------------------
  // Create message box with stop button
  // -----------------------------------
  message_box_ = new QMessageBox(QMessageBox::NoIcon,
                                 "Meshing",
                                 "Mesh generation in progress...",
                                 QMessageBox::Cancel,
                                 mw);
  
  message_box_->setDefaultButton(QMessageBox::Cancel);
  QAbstractButton* cancelButton = message_box_->button(QMessageBox::Cancel);
  cancelButton->setText(tr("Stop"));
  
  QObject::connect(cancelButton, SIGNAL(clicked()),
                   mesh_thread,  SLOT(stop()));
  
  message_box_->show();
  
  // -----------------------------------
  // Connect main thread to meshing thread
  // -----------------------------------
  QObject::connect(mesh_thread, SIGNAL(done(Meshing_thread*)),
                   this,        SLOT(meshing_done(Meshing_thread*)));
  
  QObject::connect(mesh_thread, SIGNAL(status_report(QString)),
                   this,        SLOT(status_report(QString)));
  
  // -----------------------------------
  // Launch mesher
  // -----------------------------------
  std::cerr << "mesh_thread->start(); "  << std::endl;
  mesh_thread->start();
}


void
Swept_volume_3_plugin::
status_report(QString str)
{
  if ( NULL == message_box_ ) { return; }
  
  message_box_->setInformativeText(str);
}


void
Swept_volume_3_plugin::
meshing_done(Meshing_thread* thread)
{
  // Print message in console
  QString str = QString("Meshing of \"%1\" done in %2s<br>")
    .arg(source_item_->name())
    .arg(thread->time());
  
  Q_FOREACH( QString param, thread->parameters_log() )
  {
    str.append(QString("( %1 )<br>").arg(param));
  }
  
  msg->information(qPrintable(str));
  
  // Treat new c3t3 item
  Scene_c3t3_item* result_item = thread->item();
  treat_result(*source_item_, *result_item);
  
  // close message box
  message_box_->close();
  message_box_ = NULL;
  
  // free memory
  // TODO: maybe there is another way to do that
  delete thread;
}


void
Swept_volume_3_plugin::
treat_result(Scene_item& source_item,
             Scene_c3t3_item& result_item) const
{
  result_item.setName(tr("%1 [3D Mesh]").arg(source_item.name()));

  result_item.c3t3_changed();
  
  const Scene_item::Bbox& bbox = result_item.bbox();
  result_item.setPosition((bbox.xmin + bbox.xmax)/2.f,
                          (bbox.ymin + bbox.ymax)/2.f,
                          (bbox.zmin + bbox.zmax)/2.f);
    
  result_item.setColor(default_mesh_color);
  result_item.setRenderingMode(source_item.renderingMode());
  result_item.set_data_item(&source_item);
  
  source_item.setVisible(false);
    
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  scene->itemChanged(index);
    
  Scene_interface::Item_id new_item_id = scene->addItem(&result_item);
  scene->setSelectedItem(new_item_id);
}


Q_EXPORT_PLUGIN2(Swept_volume_3_plugin, Swept_volume_3_plugin);

#include "Swept_volume_3_plugin.moc"
