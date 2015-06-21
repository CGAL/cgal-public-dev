#include <QTime>
#include <QApplication>
#include <QAction>
#include <QList>
#include <QMainWindow>
#include <QObject>

#include "Scene_polygon_soup_item.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type.h"

#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include "ui_Polyhedron_demo_scale_space_reconstruction_plugin.h"


class Polyhedron_demo_scale_space_reconstruction_plugin_dialog : public QDialog, private Ui::ScaleSpaceOptionsDialog
{
  Q_OBJECT
  public:
    Polyhedron_demo_scale_space_reconstruction_plugin_dialog(QWidget* /*parent*/ = 0)
    {
      setupUi(this);
    }

    double neighbors() const { return m_neighbors->value(); }
    double iterations() const { return m_iterations->value(); }
    double samples() const { return m_samples->value(); }
    bool generate_smoothed() const { return m_genSmooth->isChecked(); }
};

#include <CGAL/Scale_space_surface_reconstruction_3.h>

class Polyhedron_demo_scale_space_reconstruction_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  QAction* actionScaleSpaceReconstruction;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {

    actionScaleSpaceReconstruction = new QAction(tr("Scale-space surface reconstruction"), mainWindow);
    actionScaleSpaceReconstruction->setObjectName("actionScaleSpaceReconstruction");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  //! Applicate for Point_sets with normals.
  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionScaleSpaceReconstruction;
  }

public Q_SLOTS:
  void on_actionScaleSpaceReconstruction_triggered();
}; // end class Polyhedron_scale_space_reconstruction_plugin


void Polyhedron_demo_scale_space_reconstruction_plugin::on_actionScaleSpaceReconstruction_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* pts_item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(pts_item)
  {
    //generate the dialog box to set the options
    Polyhedron_demo_scale_space_reconstruction_plugin_dialog dialog;
    if(!dialog.exec())
      return;

    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);

    QTime time;
    time.start();
    std::cout << "Scale scape surface reconstruction...";

    typedef CGAL::Scale_space_surface_reconstruction_3<Kernel> Recontructor;
    Recontructor reconstruct( dialog.neighbors(), dialog.samples() );
    reconstruct.reconstruct_surface(
      pts_item->point_set()->begin(),
      pts_item->point_set()->end(),
      dialog.iterations()
    );
    std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

    for( unsigned int sh = 0; sh < reconstruct.number_of_shells(); ++sh ) {
        // collect the number of triples.
        std::ptrdiff_t num = std::distance( reconstruct.shell_begin( sh ),
                                            reconstruct.shell_end( sh ) );

        //create item for the reconstruction output with input point set
        Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item();
        new_item->init_polygon_soup(pts_item->point_set()->size(),
                                    num );

        typedef Point_set::iterator Point_iterator;

        for(Point_iterator it = pts_item->point_set()->begin(),
                           end = pts_item->point_set()->end(); it!=end; ++it)
        {
          new_item->new_vertex(it->x(), it->y(), it->z());
        }

        for (Recontructor::Triple_iterator it=reconstruct.shell_begin( sh ),
                                           end=reconstruct.shell_end( sh );it!=end;++it)
        {
          new_item->new_triangle( (*it)[0], (*it)[1], (*it)[2] );
        }

        new_item->finalize_polygon_soup();

        new_item->setName(tr("%1-shell %2 (ss reconstruction)").arg(scene->item(index)->name()).arg(sh+1));
        new_item->setColor(Qt::magenta);
        new_item->setRenderingMode(FlatPlusEdges);
        scene->addItem(new_item);

        if ( dialog.generate_smoothed() ){
          //create item for the reconstruction output with input point set smoothed
          Scene_polygon_soup_item *new_item_smoothed = new Scene_polygon_soup_item();

          new_item_smoothed->init_polygon_soup(pts_item->point_set()->size(),
                                               num );

          typedef Recontructor::Point_iterator SS_point_iterator;
          for(SS_point_iterator it = reconstruct.points_begin(),
                                end = reconstruct.points_end(); it!=end; ++it)
          {
            new_item_smoothed->new_vertex(it->x(), it->y(), it->z());
          }

          for (Recontructor::Triple_iterator it=reconstruct.shell_begin( sh ),
                                             end=reconstruct.shell_end( sh );it!=end;++it)
          {
            new_item_smoothed->new_triangle( (*it)[0], (*it)[1], (*it)[2] );
          }

          new_item_smoothed->finalize_polygon_soup();

          new_item_smoothed->setName(tr("%1-shell %2 (ss smoothed reconstruction)").arg(scene->item(index)->name()).arg(sh+1));
          new_item_smoothed->setColor(Qt::magenta);
          new_item_smoothed->setRenderingMode(FlatPlusEdges);
          scene->addItem(new_item_smoothed);
        }
    }

    // default cursor
    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_scale_space_reconstruction_plugin, Polyhedron_demo_scale_space_reconstruction_plugin)

#include "Polyhedron_demo_scale_space_reconstruction_plugin.moc"
