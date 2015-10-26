#ifndef WINDOW_
#define WINDOW_

// SLT
#include <list>

// Qt
#include <QWidget>
#include <QString>

// local
#include "scene.h"
#include "ui_pwsrec.h"

class MainWindow : public QMainWindow, public Ui_MainWindow
{
  Q_OBJECT

private:
  Scene* m_scene;

  // main options
  int  m_verbose;
  int  m_mchoice;
  bool m_use_flip;
  double m_alpha;
  double m_ghost;
  double m_percent;
  double m_norm_tol;
  double m_tang_tol;
  double m_relocation;

  unsigned int maxNumRecentFiles;
  QAction* recentFilesSeparator;
  QVector<QAction*> recentFileActs;

public:
  MainWindow();
  ~MainWindow();

  // Parameters
  void set_scene_options();

  double percentage() const { return m_percent / 100.0; }

  double norm_tol() const { return m_norm_tol / 100.0; }

  double tang_tol() const { return m_tang_tol / 100.0; }

  double alpha() const { return m_alpha / 100.0; }

  double min_mass() const
  {
    double value = double(min_mass_slider->value());
    double max_value = double(min_mass_slider->maximum());
    return value / max_value;
  }

  int ignore_edges() const
  {
    return discard_spinbox->value();
  }

  protected Q_SLOTS:
  // drag and drop
  void dropEvent(QDropEvent *event);
  void closeEvent(QCloseEvent *event);
  void dragEnterEvent(QDragEnterEvent *event);

  // recent files
  void openRecentFile_aux();
  void updateRecentFileActions();
  void addToRecentFiles(QString fileName);
  void addRecentFiles(QMenu* menu, QAction* insertBefore = 0);
  unsigned int maxNumberOfRecentFiles() const {return maxNumRecentFiles;}

  // io
  void open(const QString& file);
  void save(const QString& file);

  public Q_SLOTS:
  // render
  void update();
  void on_actionRecenter_triggered();

  // io
  void on_actionClear_triggered();
  void on_actionLoadPoints_triggered();
  void on_actionSave_triggered();
  void on_actionInsertPoint_toggled();
  void on_actionSnapshot_triggered();
  void on_actionInvert_mass_triggered();
  void on_actionClamp_mass_triggered();
  void on_actionSubdivide_triggered();
  void on_actionDecimate_triggered();
  void on_actionKeep_one_point_out_of_n_triggered();

  // data
  void on_actionStar_triggered();
  void on_actionBox_triggered();
  void on_actionLine_triggered();
  void on_actionStair_triggered();
  void on_actionBoxes_triggered();
  void on_actionNoise_triggered();
  void on_actionSpiral_triggered();
  void on_actionCircle_triggered();
  void on_actionSkyline_triggered();
  void on_actionHalf_circle_triggered();
  void on_actionAdd_outliers_triggered();
  void on_actionParallel_lines_triggered();
  void on_actionBox_with_boundaries_triggered();
  void on_actionBox_with_missing_corners_triggered();
  void on_actionIncreasingly_sharp_angles_triggered();
  void on_actionWidely_variable_sampling_triggered();

  // reconstruction
  void on_actionSet_MChoice_toggled();
  void on_actionSet_options_triggered();
  void on_actionReconstruction_one_step_triggered();
  void on_actionReconstruction_10_steps_triggered();
  void on_actionReconstruction_100_steps_triggered();
  void on_actionReconstruction_1000_steps_triggered();
  void on_actionReconstruction_until_triggered();
  void on_actionRelocate_vertices_triggered();

  // view
  void on_actionView_points_toggled();
  void on_actionView_vertices_toggled();
  void on_actionView_edges_toggled();
  void on_actionView_ghost_toggled();
  void on_actionView_edge_cost_toggled();
  void on_actionView_edge_priority_toggled();
  void on_actionView_relevance_toggled();

  void on_actionView_bins_toggled();
  void on_actionView_foot_points_toggled();
  void on_actionView_relocation_toggled();
  void on_actionView_tolerance_toggled();

  Q_SIGNALS:
  void openRecentFile(QString filename);
};

#endif // WINDOW_
