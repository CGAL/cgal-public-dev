#ifndef POLYGON_OFFSET_2_MAIN_WINDOW_H
#define POLYGON_OFFSET_2_MAIN_WINDOW_H

#include "EditParametersDialog.h"
#include "typedefs.h"
#include "colors.h"

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>

// for viewportsBbox
#include <CGAL/Qt/utility.h>
// the two base classes
#include "ui_Polygon_offset_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Polygon_offset_2
{

Q_OBJECT

public:
  MainWindow();

  typedef EditParametersDialog::Parameters Parameters;
  typedef Minsum_with_circle_2 MSWC;
  typedef MSWC::Data_type Data_type;
  typedef MSWC::Input_rational Input_rational;

public slots:

  void processInput(CGAL::Object o);

  // File menu
  void on_actionLoadPolygon_triggered();
  void on_actionSavePolygon_triggered();
  void on_actionClear_triggered();

  // Edit menu
  void on_actionEditInputPolygon_toggled(bool checked);
  void on_actionEditParameters_triggered();
  void changedParameters(const Parameters& parameters);
  void on_actionConstructionType_toggled(bool checked);

  // View menu
  void on_actionViewPolygon_toggled(bool checked);

  void on_actionViewKgon_toggled(bool checked);
  void on_actionConstructKgonInducedCircles_toggled(bool checked);
  void on_actionShowKgonSkippedCircles_toggled(bool checked);

  void on_actionRecenter_triggered();

  // Constructions menu
  void update(const Parameters& parameters);

  // forward engineering
  void on_actionConstructKgonSum_toggled(bool checked);
  void on_actionConstructKgonOffset_toggled(bool checked);
  void on_actionConstructExactOffset_toggled(bool checked);
  void on_actionConstructAppOffset_toggled(bool checked);

  // straight skeleton
  void on_actionConstructInnerSkeleton_toggled(bool checked);
  void on_actionConstructOuterSkeleton_toggled(bool checked);

  // reverse engineering
  void on_actionConstructOuterCore_toggled(bool checked);
  void on_actionConstructInnerCore_toggled(bool checked);

  void on_actionCheckApproximability_toggled(bool checked);

  void on_actionConstructInnerTolerance_toggled(bool checked);
  void on_actionConstructOuterTolerance_toggled(bool checked);
  void on_actionConstructOuterCoreApp_toggled(bool checked);
  void on_actionConstructInnerCoreApp_toggled(bool checked);
  void on_actionConstructOuterCoreSum_toggled(bool checked);
  void on_actionConstructInnerCoreSum_toggled(bool checked);
  void on_actionConstructInnerOffsetApp_toggled(bool checked);
  void on_actionConstructOuterOffsetApp_toggled(bool checked);  
  void on_actionConstructInnerOffset_toggled(bool checked);
  void on_actionConstructOuterOffset_toggled(bool checked);
  void on_actionConstructOuterDiffApp_toggled(bool checked);
  void on_actionConstructInnerDiffApp_toggled(bool checked);
  void on_actionYesApproximability_toggled(bool checked);
  void on_actionNoApproximability_toggled(bool checked);
  void on_actionCheckApproximabilityApp_toggled(bool checked);

  // search
  void on_actionSearchMinEpsilon_triggered();  
 
  // Reverse Engineering menu
  void compute_lines();
  void compute_circles();

  void on_actionShow_all_triggered();
  void on_actionHide_all_triggered();
  
  void on_actionLines_show_all_triggered();
  void on_actionLines_hide_all_triggered();

  void on_actionLines_eps_toggled(bool checked);
  void on_actionLines_r_minus_eps_toggled(bool checked);
  void on_actionLines_r_toggled(bool checked);
  void on_actionLines_r_plus_eps_toggled(bool checked);

  void on_actionCircles_show_all_triggered();
  void on_actionCircles_hide_all_triggered();

  void on_actionCircles_eps_toggled(bool checked);
  void on_actionCircles_r_minus_eps_toggled(bool checked);
  void on_actionCircles_r_toggled(bool checked);
  void on_actionCircles_r_plus_eps_toggled(bool checked);

/*
    QAction *actionInsertPoint;
    QAction *actionKgon;
*/

signals:
  void changed();
  
  void radiusVal(QString);
  void epsHatVal(QString);
  void pgonSize(QString);
  void kgonSize(QString);

private:  

  struct LineStyle
  {
    const QColor& color;
    Qt::PenStyle style;
    bool is_bold;
    double width;
    static double pen_width;

    LineStyle(const QColor& l_color, Qt::PenStyle l_style = Qt::SolidLine, bool l_bold = false):
      color(l_color), style(l_style), is_bold(l_bold), width(is_bold? 2*pen_width: pen_width) {}

    // Qt::SolidLine, Qt::DashLine, Qt::DotLine
    
    QPen getPen() const
    {
      return QPen(color, width, style, Qt::RoundCap, Qt::RoundJoin);
    }
  };

  std::string rat_to_str(const Input_rational& i_rat);
  
  template <typename GI>
  void setupGI(const Data_type& i_data, GI* i_gi, QAction* i_action, const LineStyle& line_style);

  template <typename GI>
  void setupGI(const Data_type& i_data, GI* i_fwd_gi, GI* i_bwd_gi, 
      QAction* i_action, const LineStyle& line_style);

  template <typename GI>
  void setupGI(GI* i_gi, QAction* i_action, const LineStyle& line_style);

  template <typename GI>
  void set_color(const QColor& color, GI* i_gi);

  template <typename GI>
  void set_width(const double& width, GI* i_gi);

  template <typename GI>
  void set_style(const Qt::PenStyle& style, GI* i_gi);

  // checked: true/activate, false/deactivate
  void toggle_action_constructions();
  void toggle_forward_constructions(bool checked);
  void toggle_reverse_constructions(bool checked);

  void enable(bool status, const Data_type& i_data);
  void check(bool status, const Data_type& i_data);
  void enable_and_hide(const Data_type& i_data);
  void disable_and_hide(const Data_type& i_data);
  void enable_and_show(const Data_type& i_data);

  void update_if_needed();
  void update_if_needed(const Data_type& i_data);
  void set_visible(bool status, const Data_type& i_data);
  void update_approximate_approximability_color();
  void update_parameter_indicators();

  QGraphicsScene scene;  
  MSWC mswc;
  GI* gi[MSWC::Data_type_size][2]; // graphic item - can have fwd and bwd actions
  QAction* action[MSWC::Data_type_size];

  // the gui input
  GVPolylineInput * gv_polyline_input;

  // the polygon
  PolygonGI * polygon_gi;

  // the kgon
  PolygonGI * kgon_gi;
  
  // the kgon sum polygon
  SlopedPolygonGI * sloped_kgon_sum_gi;
  HoledPolygonGI * holed_kgon_sum_gi;

  CirclesGI * kgon_circles_gi;
  CirclesGI * skipped_kgon_circles_gi;
  CirclesGI * intersecting_kgon_circles_gi;

  ApproximateHoledPolygonGI * kgon_offset_gi;

  // the exact offset
  ExactHoledPolygonGI * exact_offset_gi;
//  ExactPolygonGI * exact_offset_contour_gi;

  // the approximate offset
  ApproximateHoledPolygonGI * app_offset_gi;
//  ApproximatePolygonGI * app_offset_contour_gi;

  // the skeletons
  StraightSkeletonGI* inner_straight_skeleton_fwd_gi;
  StraightSkeletonGI* inner_straight_skeleton_bwd_gi;
  StraightSkeletonGI* outer_straight_skeleton_fwd_gi;
  StraightSkeletonGI* outer_straight_skeleton_bwd_gi;


  // the inner and outer core (first contour) inset r +/- eps
  ExactPolygonGI * inner_core_contour_fwd_gi;  
  ExactPolygonGI * inner_core_contour_bwd_gi;  
  ExactPolygonGI * outer_core_contour_fwd_gi;  
  ExactPolygonGI * outer_core_contour_bwd_gi;  
  
  // the approximability
  SegmentsGI * approximability_fwd_pos_gi;
  SegmentsGI * approximability_fwd_neg_gi;
  SegmentsGI * approximability_bwd_pos_gi;
  SegmentsGI * approximability_bwd_neg_gi;

  // the approximate inner and outer inflation
  // inner/outer eps-offset
  HoledPolygonGI * inner_eps_app_fwd_gi;
  HoledPolygonGI * inner_eps_app_bwd_gi;
  HoledPolygonGI * outer_eps_app_fwd_gi;
  HoledPolygonGI * outer_eps_app_bwd_gi;
  
  // the approximate inner and outer core (first contour)
  // inner/outer eps-offset with outer/inner r-inset
  PolygonGI * inner_core_app_fwd_gi;
  PolygonGI * inner_core_app_bwd_gi;
  PolygonGI * outer_core_app_fwd_gi;
  PolygonGI * outer_core_app_bwd_gi;

  // the approximate offset of inner and outer core (first contour)
  // inner/outer r+eps-offset
  HoledPolygonGI * inner_core_ms_fwd_gi;
  HoledPolygonGI * inner_core_ms_bwd_gi;
  HoledPolygonGI * outer_core_ms_fwd_gi;
  HoledPolygonGI * outer_core_ms_bwd_gi;

  // the approximate offset of inner and outer core (first contour)
  // inner/outer r-offset
  HoledPolygonGI * inner_core_msr_fwd_gi;
  HoledPolygonGI * inner_core_msr_bwd_gi;
  HoledPolygonGI * outer_core_msr_fwd_gi;
  HoledPolygonGI * outer_core_msr_bwd_gi;

  // the exact offset of inner and outer core (first contour)
  // inner/outer r-offset
  ExactPolygonGI * inner_core_offr_fwd_gi;
  ExactPolygonGI * inner_core_offr_bwd_gi;
  ExactPolygonGI * outer_core_offr_fwd_gi;
  ExactPolygonGI * outer_core_offr_bwd_gi;
  
  // the approximability
  HoledPolygonGI * inner_approximability_fwd_diff_gi;
  HoledPolygonGI * inner_approximability_bwd_diff_gi;
  HoledPolygonGI * outer_approximability_fwd_diff_gi;
  HoledPolygonGI * outer_approximability_bwd_diff_gi;

  PolygonGI * polygon_bwd_gi;
  PolygonGI * polygon_fwd_gi;

  int m_mode; // 0 - forward, 1 - backward
  // CircularArcGraphicsItem

  std::list< Line_2 > _m_eps_lines[2];
  std::list< Line_2 > _m_rme_lines[2];
  std::list< Line_2 > _m_rad_lines[2];
  std::list< Line_2 > _m_rpe_lines[2];

  std::list< Circle_2 > _m_eps_circles;
  std::list< Circle_2 > _m_rme_circles;
  std::list< Circle_2 > _m_rad_circles;
  std::list< Circle_2 > _m_rpe_circles;

  LinesGI * epslgi[2];
  LinesGI * rmelgi[2];
  LinesGI * radlgi[2];
  LinesGI * rpelgi[2];

  CirclesGI * epscgi;
  CirclesGI * rmecgi;
  CirclesGI * radcgi;
  CirclesGI * rpecgi;

};

#endif //POLYGON_OFFSET_2_MAIN_WINDOW_H
