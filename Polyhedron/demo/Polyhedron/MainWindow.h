#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include "config.h"

#include <QtOpenGL/qgl.h>
#include <CGAL/Qt/DemosMainWindow.h>
#ifdef QT_SCRIPT_LIB
#  include  <QScriptEngine>
#endif

#include <QVector>
#include <QList>
#include <QFileInfo>
#include <QStringList>
#include <QSet>

class Scene;
class Viewer;
class QTreeView;
class QMenu;
class Polyhedron_demo_io_plugin_interface;
class Polyhedron_demo_plugin_interface;
class Scene_item;
class QSortFilterProxyModel;

namespace Ui {
  class MainWindow;
}

#include "Polyhedron_type_fwd.h"

#include "Messages_interface.h"

class MainWindow : 
  public CGAL::Qt::DemosMainWindow,
  public Messages_interface
{
  Q_OBJECT
  Q_INTERFACES(Messages_interface)
public:
  MainWindow(QWidget* parent = 0);
  ~MainWindow();

  /// Find an IO plugin.
  /// @throws `std::invalid_argument` if no loader with that argument can be found
  /// @returns the IO plugin associated with `loader_name`
  Polyhedron_demo_io_plugin_interface* find_loader(const QString& loader_name) const;
  
  /// Load an item with a given loader.
  ///
  /// @throws `std::logic_error` if loading does not succeed or
  /// `std::invalid_argument` if `fileinfo` specifies an invalid file
  Scene_item* load_item(QFileInfo fileinfo, Polyhedron_demo_io_plugin_interface*);

public Q_SLOTS:
  void updateViewerBBox();
  void open(QString);

  /// given a file extension file, returns true if `filename` matches the filter
  bool file_matches_filter(const QString& filters, const QString& filename);

  /// Open a file with a given loader, and return true iff it was successful.
  ///
  /// This slot is for use by scripts.
  bool open(QString filename, QString loader_name);

  /// Reloads an item. Expects to be called by a QAction with the
  /// index of the item to be reloaded as data attached to the action.
  /// The index must identify a valid `Scene_item`.
  void reload_item();
  
  bool load_script(QString filename);
  bool load_script(QFileInfo);

  void setFocusToQuickSearch();

  void selectSceneItem(int i);
  void showSelectedPoint(double, double, double);
  void unSelectSceneItem(int i);
  void selectAll();
  void addSceneItemInSelection(int i);
  void removeSceneItemFromSelection(int i); // same as unSelectSceneItem

  void setAddKeyFrameKeyboardModifiers(Qt::KeyboardModifiers);

  void clearMenu(QMenu*);
  void addAction(QAction*);
  void addAction(QString actionName,
                 QString actionText,
                 QString menuName);
  void viewerShow(float, float, float);
  void viewerShow(float, float, float, float, float, float);
  void viewerShowObject();

  void information(QString);
  void warning(QString);
  void error(QString);
  void message(QString, QString, QString = QString("normal"));

  bool hasPlugin(const QString&) const;
  void enableScriptDebugger(bool = true);

protected Q_SLOTS:
  void selectionChanged();

  void contextMenuRequested(const QPoint& global_pos);
  void showSceneContextMenu(int selectedItemIndex,
                            const QPoint& global_pos);
  void showSceneContextMenu(const QPoint& local_pos_of_treeview);

  void updateInfo();
  void updateDisplayInfo();
  void removeManipulatedFrame(Scene_item*);

  // settings
  void quit();
  void readSettings();
  void writeSettings();

  // load, erase, duplicate
  void on_actionEraseAll_triggered();
  void on_actionLoad_triggered();
  bool on_actionErase_triggered();
  void on_actionDuplicate_triggered();
  void on_actionLoad_Script_triggered();

  // Show/Hide
  void on_actionShowHide_triggered();

  // Select A/B
  void on_actionSetPolyhedronA_triggered();
  void on_actionSetPolyhedronB_triggered();

  //Preferences edition
  void on_actionPreferences_triggered();
  // save as...
  void on_actionSaveAs_triggered(); 
  void save(QString filename, Scene_item* item);

  void on_actionSetBackgroundColor_triggered();

  void on_action_Look_at_triggered();

  QString camera_string() const;
  void on_actionDumpCamera_triggered();
  void on_action_Copy_camera_triggered();
  void on_action_Paste_camera_triggered();

  void filterOperations();

  void on_actionRecenterScene_triggered();
protected:
  void loadPlugins();
  bool initPlugin(QObject*);
  bool initIOPlugin(QObject*);

  void closeEvent(QCloseEvent *event);

  bool onePolygonIsSelected() const;
  int getSelectedSceneItemIndex() const;
  QList<int> getSelectedSceneItemIndices() const;

private:
  QString strippedName(const QString &fullFileName);

  /// plugin black-list
  QSet<QString> plugin_blacklist;

  Scene* scene;
  Viewer* viewer;
  QSortFilterProxyModel* proxyModel;
  QTreeView* sceneView;
  Ui::MainWindow* ui;
  QVector<Polyhedron_demo_io_plugin_interface*> io_plugins;
  QMap<QString,QString> default_plugin_selection;
  // typedef to make Q_FOREACH work
  typedef QPair<Polyhedron_demo_plugin_interface*, QString> PluginNamePair;
  QVector<PluginNamePair > plugins;
#ifdef QT_SCRIPT_LIB
  QScriptEngine* script_engine;
public:
  void evaluate_script(QString script, 
                       const QString & fileName = QString(),
                       const bool quiet = false);
  void evaluate_script_quiet(QString script, 
                             const QString & fileName = QString());
#endif
};

#endif // ifndef MAINWINDOW_H
