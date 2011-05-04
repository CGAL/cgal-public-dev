#include "include/MainWindow.h"

//#define SET_WIDTH // uncomment to modify the width of drawn lines (for taking pictures)
double MainWindow::LineStyle::pen_width;

template <typename GI>
void 
MainWindow::setupGI(const Data_type& i_data, GI* i_gi, QAction* i_action, const LineStyle& line_style)
//                    const QColor& color, bool is_solid, bool is_bold)
{
  setupGI(i_data, i_gi, i_gi, i_action, line_style);
}

template <typename GI>
void 
MainWindow::setupGI(const Data_type& i_data, GI* i_fwd_gi, GI* i_bwd_gi,  
                    QAction* i_action, const LineStyle& line_style)
{
  static int FWD = 0;
  static int BWD = 1;
  gi[i_data][FWD] = i_fwd_gi;
  gi[i_data][BWD] = i_bwd_gi;
  action[i_data] = i_action;
//  LOG_DEBUG << "setup GI[" << i_data << "] with " << i_fwd_gi
//    << (i_fwd_gi == i_bwd_gi? "": " and ") << (i_fwd_gi == i_bwd_gi? NULL: i_bwd_gi) << std::endl;
//  LOG_DEBUG << "setup action[" << i_data << "] with " << i_action << std::endl;

  QObject::connect(this, SIGNAL(changed()), i_fwd_gi, SLOT(modelChanged()));
  if(i_fwd_gi != i_bwd_gi)
    QObject::connect(this, SIGNAL(changed()), i_bwd_gi, SLOT(modelChanged()));

  QPen epen = line_style.getPen();

  i_fwd_gi->setEdgesPen(epen);
  i_fwd_gi->setVerticesPen(epen);
  scene.addItem(i_fwd_gi);
  i_fwd_gi->hide();

  if(i_fwd_gi != i_bwd_gi)
  {
    i_bwd_gi->setEdgesPen(epen);
    i_bwd_gi->setVerticesPen(epen);
    scene.addItem(i_bwd_gi);
    i_bwd_gi->hide();
  }
}

template <typename GI>
void 
MainWindow::setupGI(GI* i_gi, QAction* i_action, const LineStyle& line_style)
{
  QObject::connect(this, SIGNAL(changed()), i_gi, SLOT(modelChanged()));

  QPen epen = line_style.getPen();

  i_gi->setEdgesPen(epen);
  i_gi->setVerticesPen(epen);
  scene.addItem(i_gi);
  i_gi->hide();
}

template <>
void 
MainWindow::setupGI<CirclesGI>(const Data_type& i_data, CirclesGI* i_gi,
                               QAction* i_action, const LineStyle& line_style)
{
  static int FWD = 0;
  static int BWD = 1;
  gi[i_data][FWD] = i_gi;
  gi[i_data][BWD] = i_gi;
  action[i_data] = i_action;
//  LOG_DEBUG << "setup GI[" << i_data << "] with " << i_gi << std::endl;
//  LOG_DEBUG << "setup action[" << i_data << "] with " << i_action << std::endl;

  QObject::connect(this, SIGNAL(changed()), i_gi, SLOT(modelChanged()));

  QPen epen = line_style.getPen();
  i_gi->setPen(epen);

  scene.addItem(i_gi);
  i_gi->hide();
}

template <>
void 
MainWindow::setupGI<CirclesGI>(CirclesGI* i_gi, QAction* i_action, const LineStyle& line_style)
{
  QObject::connect(this, SIGNAL(changed()), i_gi, SLOT(modelChanged()));

  QPen epen = line_style.getPen();
  i_gi->setPen(epen);

  scene.addItem(i_gi);
  i_gi->hide();
}

template <typename GI>
void
MainWindow::set_color(const QColor& color, GI* i_gi)
{
  QPen epen = i_gi->edgesPen();
  epen.setColor(color);
  i_gi->setEdgesPen(epen);
}

template <typename GI>
void
MainWindow::set_width(const double& width, GI* i_gi)
{
  QPen epen = i_gi->edgesPen();
  epen.setWidthF(width);
  i_gi->setEdgesPen(epen);
}

template <typename GI>
void
MainWindow::set_style(const Qt::PenStyle& style, GI* i_gi)
{
  QPen epen = i_gi->edgesPen();
  epen.setStyle(style);
  i_gi->setEdgesPen(epen);
}

MainWindow::MainWindow(): DemosMainWindow()
{
  LOG_DEBUG << "Construct MainWindow " << std::endl;

#ifdef SET_WIDTH
//  LineStyle::pen_width = 0.035;
//  LineStyle::pen_width = 0.002; // kaz
  LineStyle::pen_width = 0.05; // stars_part
#else
  LineStyle::pen_width = 0;
#endif // SET_WIDTH

  m_mode = static_cast<int>(!mswc.forward_construction());


/*
 *  Qt Automatic Connections
 *  http://doc.trolltech.com/4.4/designer-using-a-component.html#automatic-connections
 *
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */

  setupUi(this);

  // Add a GraphicItem for the Polygon
  polygon_gi = new PolygonGI(&mswc.polygon());
  setupGI(MSWC::Data_polygon, polygon_gi, actionViewPolygon, LineStyle(polygon_color.polygon, Qt::SolidLine, true));

  // Add a GraphicItem for the Kgon 
  kgon_gi = new PolygonGI(&mswc.kgon());
  setupGI(MSWC::Data_kgon, kgon_gi, actionViewKgon, LineStyle(polygon_color.kgon));

  // Add a GraphicItem for the Kgon sum with Polygon
  // TODO replace outer boundary, by all boundaries
  sloped_kgon_sum_gi = new SlopedPolygonGI(
    &mswc.kgon_sum_outer_boundary(), &mswc.circle_app(), polygon_color.kgon_circles);
  setupGI(sloped_kgon_sum_gi, actionConstructKgonSum, LineStyle(polygon_color.kgon_sum));

  // Add a GraphicItem for the Kgon sum with Polygon 
  // with all boundary (including holes)
  //Circle_app_col_2 circle_app_col(mswc.circle_app(), polygon_color.kgon_circles, polygon_color.kgon_sum);
  holed_kgon_sum_gi = new HoledPolygonGI(&mswc.kgon_sum());
  setupGI(MSWC::Data_kgon_sum, holed_kgon_sum_gi, actionConstructKgonSum, LineStyle(polygon_color.kgon_sum));

  kgon_circles_gi = new CirclesGI(&mswc.kgon_circles());
  setupGI(MSWC::Data_kgon_induced_circles, kgon_circles_gi, actionConstructKgonInducedCircles, 
    LineStyle(polygon_color.kgon_circles, Qt::DashLine));

  skipped_kgon_circles_gi = new CirclesGI(&mswc.kgon_skipped_circles());
  setupGI(skipped_kgon_circles_gi, actionShowKgonSkippedCircles, 
    LineStyle(polygon_color.kgon_skipped_circles, Qt::DotLine));

  intersecting_kgon_circles_gi = new CirclesGI(&mswc.kgon_intersecting_circles());
  setupGI(intersecting_kgon_circles_gi, actionShowKgonSkippedCircles,
    LineStyle(polygon_color.core, Qt::DotLine));

  kgon_offset_gi = new ApproximateHoledPolygonGI(&mswc.kgon_offset_polygon());
  setupGI(MSWC::Data_kgon_offset_polygon, kgon_offset_gi, actionConstructKgonOffset, 
    LineStyle(polygon_color.kgon_offset, Qt::SolidLine));

  exact_offset_gi = new ExactHoledPolygonGI(&mswc.exact_offset_polygon());
  setupGI(MSWC::Data_exact_offset_polygon, exact_offset_gi, actionConstructExactOffset, 
    LineStyle(polygon_color.exact_offset));
  
  app_offset_gi = new ApproximateHoledPolygonGI(&mswc.approximate_offset_polygon());
  setupGI(MSWC::Data_approximate_offset_polygon, app_offset_gi, actionConstructAppOffset, 
    LineStyle(polygon_color.approximate_offset, Qt::DashLine));
  
  inner_straight_skeleton_fwd_gi = new StraightSkeletonGI(mswc.inner_skeleton_ptr(true));
  inner_straight_skeleton_bwd_gi = new StraightSkeletonGI(mswc.inner_skeleton_ptr(false));
  setupGI(MSWC::Data_inner_skeleton, inner_straight_skeleton_fwd_gi, inner_straight_skeleton_bwd_gi,
    actionConstructInnerSkeleton, LineStyle(polygon_color.skeleton));

  outer_straight_skeleton_fwd_gi = new StraightSkeletonGI(mswc.outer_skeleton_ptr(true));
  outer_straight_skeleton_bwd_gi = new StraightSkeletonGI(mswc.outer_skeleton_ptr(false));
  setupGI(MSWC::Data_outer_skeleton, outer_straight_skeleton_fwd_gi, outer_straight_skeleton_bwd_gi,
    actionConstructOuterSkeleton, LineStyle(polygon_color.skeleton));

  inner_core_contour_fwd_gi = new ExactPolygonGI(&mswc.inner_core(true));
  inner_core_contour_bwd_gi = new ExactPolygonGI(&mswc.inner_core(false));
  setupGI(MSWC::Data_inner_core, inner_core_contour_fwd_gi, inner_core_contour_bwd_gi,
    actionConstructInnerCore, LineStyle(polygon_color.core));

  outer_core_contour_fwd_gi = new ExactPolygonGI(&mswc.outer_core(true));
  outer_core_contour_bwd_gi = new ExactPolygonGI(&mswc.outer_core(false));
  setupGI(MSWC::Data_outer_core, outer_core_contour_fwd_gi, outer_core_contour_bwd_gi,
    actionConstructOuterCore, LineStyle(polygon_color.core));

  approximability_fwd_pos_gi = new SegmentsGI(&mswc.approximable_data(true));
  approximability_fwd_neg_gi = new SegmentsGI(&mswc.non_approximable_data(true));
  approximability_bwd_pos_gi = new SegmentsGI(&mswc.approximable_data(false));
  approximability_bwd_neg_gi = new SegmentsGI(&mswc.non_approximable_data(false));

  setupGI(MSWC::Data_approximability, approximability_fwd_pos_gi, approximability_bwd_pos_gi,
    actionCheckApproximability, LineStyle(polygon_color.approximability_positive));
  setupGI(approximability_fwd_neg_gi, actionCheckApproximability,
          LineStyle(polygon_color.approximability_negative));
  setupGI(approximability_bwd_neg_gi, actionCheckApproximability,
          LineStyle(polygon_color.approximability_negative));
    
  inner_core_app_fwd_gi = new PolygonGI(&mswc.inner_core_app(true));
  inner_core_app_bwd_gi = new PolygonGI(&mswc.inner_core_app(false));
  setupGI(MSWC::Data_approximate_inner_core, inner_core_app_fwd_gi, inner_core_app_bwd_gi,
    actionConstructInnerCoreApp, LineStyle(polygon_color.green, Qt::SolidLine));

  outer_core_app_fwd_gi = new PolygonGI(&mswc.outer_core_app(true));
  outer_core_app_bwd_gi = new PolygonGI(&mswc.outer_core_app(false));
  setupGI(MSWC::Data_approximate_outer_core, outer_core_app_fwd_gi, outer_core_app_bwd_gi,
    actionConstructOuterCoreApp, LineStyle(polygon_color.red, Qt::SolidLine)); // Qt::DashLine));

  inner_core_ms_fwd_gi = new HoledPolygonGI(&mswc.inner_core_ms(true));
  inner_core_ms_bwd_gi = new HoledPolygonGI(&mswc.inner_core_ms(false));
  setupGI(MSWC::Data_approximate_inner_kgon_sum, inner_core_ms_fwd_gi, inner_core_ms_bwd_gi,
    actionConstructInnerCoreSum, LineStyle(polygon_color.cyan, Qt::SolidLine));

  outer_core_ms_fwd_gi = new HoledPolygonGI(&mswc.outer_core_ms(true));
  outer_core_ms_bwd_gi = new HoledPolygonGI(&mswc.outer_core_ms(false));
  setupGI(MSWC::Data_approximate_outer_kgon_sum, outer_core_ms_fwd_gi, outer_core_ms_bwd_gi,
    actionConstructOuterCoreSum, LineStyle(polygon_color.magenta, Qt::SolidLine)); // Qt::DotLine)); 
  
  inner_approximability_fwd_diff_gi = new HoledPolygonGI(&mswc.inner_core_diff(true));
  inner_approximability_bwd_diff_gi = new HoledPolygonGI(&mswc.inner_core_diff(false));
  setupGI(MSWC::Data_approximate_inner_decision, inner_approximability_fwd_diff_gi, inner_approximability_bwd_diff_gi,
    actionConstructInnerDiffApp, LineStyle(polygon_color.orange));
  
  outer_approximability_fwd_diff_gi = new HoledPolygonGI(&mswc.outer_core_diff(true));
  outer_approximability_bwd_diff_gi = new HoledPolygonGI(&mswc.outer_core_diff(false));
  setupGI(MSWC::Data_approximate_outer_decision, outer_approximability_fwd_diff_gi, outer_approximability_bwd_diff_gi,
    actionConstructOuterDiffApp, LineStyle(polygon_color.red));

  polygon_fwd_gi = new PolygonGI(&mswc.kgon_sum_outer_boundary());
  polygon_bwd_gi = new PolygonGI(&mswc.polygon());

  setupGI(MSWC::Data_approximate_approximability, polygon_fwd_gi, polygon_bwd_gi,
    actionCheckApproximabilityApp,  LineStyle(polygon_color.orange, Qt::SolidLine, true)); // Qt::DashLine, true));

//  core_offset_contour_fwd_gi = new ExactPolygonGI(&mswc.core_exact_outer_boundary(true));
//  core_offset_contour_bwd_gi = new ExactPolygonGI(&mswc.core_exact_outer_boundary(false));
//  setupGI(MSWC::Data_core_offset, core_offset_contour_fwd_gi, core_offset_contour_bwd_gi,
//    actionConstructCoreOffset, LineStyle(polygon_color.exact_offset));
    
  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the polygon offset with 
  // the signal/slot mechanism    
  gv_polyline_input = new GVPolylineInput(this, &scene, 0, true); 
// emits points
  
  actionEditInputPolygon->setChecked(false);
  QObject::connect(gv_polyline_input, SIGNAL(generate(CGAL::Object)),
       this, SLOT(processInput(CGAL::Object)));

  // 
  // Manual handling of actions
  //
  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
       this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
//  QActionGroup* ag = new QActionGroup(this);
//  ag->addAction(this->actionInsertPoint);

  // Check two actions 
//  this->actionInsertPoint->setChecked(true);
//  this->actionShowRegular->setChecked(true);


  // reverse engineering lines
  
  for (int i = 0; i < 2; i++) {
    // Add a GraphicItem for the eps lines
    
    epslgi[i] = new LinesGI(&(_m_eps_lines[i]));
    
    QObject::connect(this, SIGNAL(changed()),
                     epslgi[i], SLOT(modelChanged()));
    
    (epslgi[i])->setPen(QPen(line_color.eps, 0, Qt::SolidLine, 
                             Qt::RoundCap, Qt::RoundJoin));
    
    scene.addItem(epslgi[i]);
    

    rmelgi[i] = new LinesGI(&(_m_rme_lines[i]));
    
    QObject::connect(this, SIGNAL(changed()),
                     rmelgi[i], SLOT(modelChanged()));
    
    (rmelgi[i])->setPen(QPen(line_color.rme, 0, Qt::SolidLine, 
                             Qt::RoundCap, Qt::RoundJoin));
    
    scene.addItem(rmelgi[i]);
    
    
    radlgi[i] = new LinesGI(&(_m_rad_lines[i]));
    
    QObject::connect(this, SIGNAL(changed()),
                     radlgi[i], SLOT(modelChanged()));
    
    (radlgi[i])->setPen(QPen(line_color.rad, 0, Qt::SolidLine, 
                             Qt::RoundCap, Qt::RoundJoin));
    
    scene.addItem(radlgi[i]);
  
    
    rpelgi[i] = new LinesGI(&(_m_rpe_lines[i]));
    
    QObject::connect(this, SIGNAL(changed()),
                     rpelgi[i], SLOT(modelChanged()));
    
    (rpelgi[i])->setPen(QPen(line_color.rpe, 0, Qt::SolidLine, 
                             Qt::RoundCap, Qt::RoundJoin));
    
    scene.addItem(rpelgi[i]);
       
  }
  
  
  // circles

  epscgi = new CirclesGI(&_m_eps_circles);
  
  QObject::connect(this, SIGNAL(changed()), 
                   epscgi, SLOT(modelChanged()));
  
  epscgi->setPen(QPen(circle_color.eps, 0, Qt::DashLine, 
                      Qt::RoundCap, Qt::RoundJoin));
  
  scene.addItem(epscgi);
  

  rmecgi = new CirclesGI(&_m_rme_circles);
  
  QObject::connect(this, SIGNAL(changed()), 
                   rmecgi, SLOT(modelChanged()));
    
  rmecgi->setPen(QPen(circle_color.rme, 0, Qt::DashLine, 
                      Qt::RoundCap, Qt::RoundJoin));
  
  scene.addItem(rmecgi);
    
  
  radcgi = new CirclesGI(&_m_rad_circles);
  
  QObject::connect(this, SIGNAL(changed()), 
                   radcgi, SLOT(modelChanged()));
    
  radcgi->setPen(QPen(circle_color.rad, 0, Qt::DashLine, 
                      Qt::RoundCap, Qt::RoundJoin));
  
  scene.addItem(radcgi);


  rpecgi = new CirclesGI(&_m_rpe_circles);
  
  QObject::connect(this, SIGNAL(changed()), 
                   rpecgi, SLOT(modelChanged()));
    
  rpecgi->setPen(QPen(circle_color.rpe, 0, Qt::DashLine, 
                      Qt::RoundCap, Qt::RoundJoin));
  
  scene.addItem(rpecgi);
    
  //
  // Setup the scene and the view
  //
  
  // setup parameter indicators
  QObject::connect(this, SIGNAL(radiusVal(QString)),
                   radius_lab, SLOT(setText(QString)));
  QObject::connect(this, SIGNAL(epsHatVal(QString)),
                   eps_hat_lab, SLOT(setText(QString)));
  QObject::connect(this, SIGNAL(pgonSize(QString)),
                   pgon_size_lab, SLOT(setText(QString)));
  QObject::connect(this, SIGNAL(kgonSize(QString)),
                   kgon_size_lab, SLOT(setText(QString)));
  
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(-100, -100, 100, 100);
  
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMouseTracking(true);

  // Turn the vertical axis upside down
  //this->graphicsView->matrix().scale(1, -1);

  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  //this->addAboutDemo(":/cgal/help/about_Polygon_offset_2.html");
  //this->addAboutCGAL();

  actionClear->trigger();

}


void MainWindow::enable(bool status, const Data_type& i_data)
{
    action[i_data]->setEnabled(status);
}

void MainWindow::check(bool status, const Data_type& i_data)
{
    action[i_data]->setChecked(status);

    // update data before visualization if needed
    set_visible(status, i_data);
}

void MainWindow::enable_and_hide(const Data_type& i_data)
{
    enable(true, i_data);
    check(false, i_data);
}

void MainWindow::enable_and_show(const Data_type& i_data)
{
    enable(true, i_data);
    check(true, i_data);
}

void MainWindow::disable_and_hide(const Data_type& i_data)
{
    enable(false, i_data);
    check(false, i_data);
}

void MainWindow::set_visible(bool status, const Data_type& i_data)
{
  if(status)
  {
    update_if_needed(i_data);
    gi[i_data][(m_mode + 1) % 2]->setVisible(false); // always hide the gi of another mode
    gi[i_data][m_mode]->modelChanged();
  }
  gi[i_data][m_mode]->setVisible(status);
}

// Edit Polygon actions

void 
MainWindow::on_actionLoadPolygon_triggered()
{
  LOG_DEBUG << "LoadPolygon " << std::endl;
  
  QString fileName = 
    QFileDialog::getOpenFileName(this, tr("Open Polygon file"), ".");

  if (!fileName.isEmpty()) {
    
    actionClear->trigger();
    actionClear->setEnabled(true);

    Polygon_2 r_polygon;
    std::ifstream ifs(qPrintable(fileName));
    ifs >> r_polygon;
    LOG_DEBUG << "polygon of size " << r_polygon.size() << std::endl;

    std::list<Point_2> upoints;
    unique_copy(r_polygon.vertices_begin(), r_polygon.vertices_end(), std::back_inserter(upoints));

    // check that last point is not the first point
    if(upoints.size() >= 2)
    {
      if(*(upoints.begin()) == *(upoints.rbegin()))
        upoints.pop_back();
    }
    
    MSWC::Polygon_2 polygon;
    // update polygon
    polygon.insert(polygon.vertices_begin(), upoints.begin(), upoints.end());
    LOG_DEBUG << "polygon of size " << polygon.size() << std::endl;

    mswc.polygon(polygon);

    LOG_DEBUG << "polygon is simple " << polygon.is_simple() << std::endl;
   
    actionEditInputPolygon->setChecked(false);
    actionSavePolygon->setEnabled(true);
    enable_and_show(MSWC::Data_polygon);
    actionRecenter->trigger();

    toggle_forward_constructions(actionConstructionType->isChecked());
    toggle_reverse_constructions(!actionConstructionType->isChecked());

    if (!mswc.polygon().is_empty()) {

      // reverse engineering
      compute_lines();
      compute_circles();
    }

    // update parameter indicators
    update_parameter_indicators();
    
    emit(changed());
    toggle_action_constructions();
  }
}

void
MainWindow::on_actionEditInputPolygon_toggled(bool checked)
{
  LOG_DEBUG << "actionEditInputPolygon " << checked << std::endl;

  actionEditInputPolygon->setChecked(checked);

  if(checked)
  {
//    disable_and_hide(MSWC::Data_polygon);
    scene.installEventFilter(gv_polyline_input);
  }
  else
  {
    scene.removeEventFilter(gv_polyline_input);
    //actionSavePolygon->setEnabled(true);
 //   enable_and_show(MSWC::Data_polygon);
 //   enable(false, MSWC::Data_kgon_induced_circles);
 //   actionShowKgonSkippedCircles->setEnabled(false);
 //   actionRecenter->trigger();
  }

  actionSavePolygon->setEnabled(!checked);

 // toggle_forward_constructions(actionConstructionType->isChecked());
 // toggle_reverse_constructions(!actionConstructionType->isChecked());

 // emit(changed());

}


void
MainWindow::processInput(CGAL::Object o)
{
  LOG_DEBUG << "processInput " << std::endl;

  std::list<Point_2> points;

  if(CGAL::assign(points, o)){
    actionEditInputPolygon->setChecked(false);
    actionClear->trigger();

    LOG_DEBUG << "processInput of size " << points.size() << std::endl;

    // accidentally pressing 2 times at the same point
    // one can create polygon with 0 size edges that cause crashes
    // remove such accidental duplications

    // points without duplications
    std::list<Point_2> upoints;
    unique_copy(points.begin(), points.end(), std::back_inserter(upoints));

    // check that last point is not the first point
    if(upoints.size() >= 2)
    {
      if(*(upoints.begin()) == *(upoints.rbegin()))
        upoints.pop_back();
    }

    MSWC::Polygon_2 polygon;
    // update polygon
    polygon.insert(polygon.vertices_begin(), upoints.begin(), upoints.end());
    LOG_DEBUG << "polygon of size " << polygon.size() << std::endl;

    mswc.polygon(polygon);

    LOG_DEBUG << "polygon is simple " << polygon.is_simple() << std::endl;

    // update parameter indicators
    update_parameter_indicators();

    //actionRecenter->trigger();
    enable_and_show(MSWC::Data_polygon);
    actionSavePolygon->setEnabled(true);
    actionClear->setEnabled(true);

    emit(changed());
    toggle_action_constructions();
  }

}

/*
  void
  MainWindow::on_actionInsertPoint_toggled(bool checked)
  {
  if(checked){
  scene.installEventFilter(pi);
    scene.installEventFilter(trv);
    } else {
    scene.removeEventFilter(pi);
    scene.removeEventFilter(trv);
    }
    }
*/

/*
void
MainWindow::on_actionInsertRandomPoints_triggered()
{
  QRectF rect = CGAL::Qt::viewportsBbox(&scene);
  CGAL::Qt::Converter<K> convert;
  Iso_rectangle_2 isor = convert(rect);
  CGAL::Random_points_in_iso_rectangle_2<Point_2> pg(isor.min(), isor.max());
  const int number_of_points =
    QInputDialog::getInteger(this,
                             tr("Number of random points"),
                             tr("Enter number of random points"), 100, 0);

  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);
  std::vector<Point_2> points;
  points.reserve(number_of_points);
  for(int i = 0; i < number_of_points; ++i){
    points.push_back(*pg++);
  }
  polygon.insert(points.begin(), points.end());
  // default cursor
  QApplication::setOverrideCursor(Qt::ArrowCursor);
  emit(changed());
}
*/

void
MainWindow::on_actionClear_triggered()
{
  LOG_DEBUG << "Clear " << std::endl;

  epscgi->hide();
  rmecgi->hide();
  radcgi->hide();
  rpecgi->hide();

  on_actionHide_all_triggered();

  (epslgi[0])->hide();
  (epslgi[1])->hide();
  (rmelgi[0])->hide();
  (rmelgi[1])->hide();
  (radlgi[0])->hide();
  (radlgi[1])->hide();
  (rpelgi[0])->hide();
  (rpelgi[1])->hide();

  _m_eps_lines[0].clear();
  _m_eps_lines[1].clear();
  _m_rme_lines[0].clear();
  _m_rme_lines[1].clear();
  _m_rad_lines[0].clear();
  _m_rad_lines[1].clear();
  _m_rpe_lines[0].clear();
  _m_rpe_lines[1].clear();

  _m_eps_circles.clear();
  _m_rme_circles.clear();
  _m_rad_circles.clear();
  _m_rpe_circles.clear();

  mswc.clear();

  actionClear->setEnabled(false);

  actionEditInputPolygon->setChecked(false);
  actionSavePolygon->setEnabled(false);
  disable_and_hide(MSWC::Data_polygon);

  toggle_forward_constructions(actionConstructionType->isChecked());
  toggle_reverse_constructions(!actionConstructionType->isChecked());

  update_parameter_indicators();

  emit(changed());
  toggle_action_constructions();
}


// Edit Parameters actions

void
MainWindow::on_actionEditParameters_triggered()
{
  LOG_DEBUG << "actionEditParameters " << std::endl;

  EditParametersDialog dialog(this);
  dialog.exec();
}

void MainWindow::changedParameters(const Parameters& parameters)
{
  LOG_DEBUG << "changedParameters: type(" << parameters.type
            << ") size(" << parameters.size 
            << ") radius(" << parameters.radius << ")"
            << " epsilon(" << parameters.epsilon << ") "
            << " delta(" << parameters.delta << ") " << std::endl;
     
  // update data from parameters
  // and graphic items from data when needed
  update(parameters);

  if (!mswc.polygon().is_empty()) {

    // reverse engineering
    compute_lines();
    compute_circles();
  }
}


void MainWindow::update(const Parameters& parameters)
{
    // copy parameters values from UI to the application
    
    Minsum_with_circle_2::Input_rational offset_radius(parameters.radius);
    mswc.offset(offset_radius);

    mswc.kgon_type(Minsum_with_circle_2::Kgon_type_2(parameters.type));

//    // changed size can override _unchanged_ epsilon
//    // changed epsilon can override size (epsilon "beats" size)
    Minsum_with_circle_2::Input_rational epsilon(parameters.epsilon);
    mswc.epsilon(epsilon, false);
    /*
    bool changed_size = (parameters.size != mswc.kgon_size());
    bool changed_epsilon = (epsilon != mswc.epsilon());

    if (changed_size)
    {
      // update eps from size only if it was not changed
      mswc.kgon_size(parameters.size, !changed_epsilon);
      if(!changed_epsilon)
      {
        // update UI epsilon from size via mswc
        std::ostringstream epsilon_stream;
        epsilon_stream << mswc.epsilon();
        std::string epsilon_string = epsilon_stream.str();
        EditParametersDialog::saved_epsilon(epsilon_string);
      }
    }

    if(changed_epsilon)
    {
      // always update size from changed epsilon
      mswc.epsilon(epsilon, true);
      EditParametersDialog::saved_size(mswc.kgon_size());
    }
    */

    Minsum_with_circle_2::Input_rational delta(parameters.delta);
    mswc.delta(delta);

    // update all graphic items that could have been affected
    // by recent change in data

    update_if_needed();
}

void MainWindow::update_if_needed()
{
  for(int data_type = MSWC::Data_polygon; data_type != MSWC::Data_type_size; ++data_type)
  {
    // update data for activated graphics items only
    if(action[data_type]->isChecked())
    {
      bool has_changed = mswc.update(static_cast<Data_type>(data_type));
      if(has_changed && data_type == MSWC::Data_approximate_approximability)
      {
          update_approximate_approximability_color();
      }
    }
  }

  // update parameter indicators
  update_parameter_indicators();

  emit(changed());
  toggle_action_constructions();  
}

void MainWindow::update_approximate_approximability_color()
{
  set_color(mswc.approximability_indicator(),
    static_cast<PolygonGI*>(gi[MSWC::Data_approximate_approximability][m_mode]));
}

void MainWindow::update_if_needed(const Data_type& i_data)
{
    bool has_changed = mswc.update(i_data);
    
    if(has_changed)
    {
      if(i_data == MSWC::Data_approximate_approximability)
      {
         update_approximate_approximability_color();
      }

      update_parameter_indicators();
       
      emit(changed());
      toggle_action_constructions();
      /*
      gi[i_data][m_mode]->modelChanged();
      // just in case induced update got missed by UI
      switch(i_data)
      {
      case MSWC::Data_kgon_induced_circles:
        if(skipped_kgon_circles_gi->isVisible())
          skipped_kgon_circles_gi->modelChanged();
      case MSWC::Data_kgon_offset_polygon:
      case MSWC::Data_outer_core:
      case MSWC::Data_inner_core:
      case MSWC::Data_approximability:
        if(gi[MSWC::Data_kgon_sum][m_mode]->isVisible())
          gi[MSWC::Data_kgon_sum][m_mode]->modelChanged();
        // fall thru
      case MSWC::Data_kgon_sum:
        if(sloped_kgon_sum_gi->isVisible())
          sloped_kgon_sum_gi->modelChanged();
        if(gi[MSWC::Data_kgon][m_mode]->isVisible())
          gi[MSWC::Data_kgon][m_mode]->modelChanged();
        // fall thru
      case MSWC::Data_inner_skeleton:
      case MSWC::Data_outer_skeleton:
      case MSWC::Data_kgon:
      case MSWC::Data_exact_offset_polygon:
      case MSWC::Data_approximate_offset_polygon: 
        if(gi[MSWC::Data_polygon][m_mode]->isVisible())
          gi[MSWC::Data_polygon][m_mode]->modelChanged();
        break;
        
      };
      */
    }
}

void MainWindow::toggle_action_constructions()
{
  // ensure that changed() signal does not make hidden things visible
  for(int data_type = MSWC::Data_kgon; data_type != MSWC::Data_type_size; ++data_type)
  {
    gi[data_type][(m_mode + 1) % 2]->setVisible(false); // always hide the gi of another mode
    bool checked = action[data_type]->isChecked();
    gi[data_type][m_mode]->setVisible(checked); // always show construction of checked action
  }
  sloped_kgon_sum_gi->setVisible(action[MSWC::Data_kgon_sum]->isEnabled() && action[MSWC::Data_kgon_sum]->isChecked());
}


void
MainWindow::update_parameter_indicators()
{
  EditParametersDialog::saved_polygon_size(mswc.polygon().size());

  double bboxSize = 1;
  if(mswc.polygon().size() != 0)
  {
    CGAL::Bbox_2 bbox = mswc.polygon().bbox();
    double bboxXDim = (bbox.xmax() - bbox.xmin());
    double bboxYDim = (bbox.ymax() - bbox.ymin());
    bboxSize = bboxXDim > bboxYDim? bboxXDim: bboxYDim;
    LOG_DEBUG << "bbox of size " << bboxSize << std::endl;
  }
  Minsum_with_circle_2::Input_rational bboxRatSize(ceil(bboxSize));
  std::ostringstream bbox_stream;
  bbox_stream << bboxRatSize;
  std::string bbox_string = bbox_stream.str();
  EditParametersDialog::saved_polygon_bbox(bbox_string);

  QString pgon_size_str = QString(" polygon size = ") + QString::number(uint(mswc.polygon().size())) + " " +
    QString(" bbox = ") + QString::number(uint(ceil(bboxSize))) + " ";
  emit pgonSize(pgon_size_str);

  QString radius_str = QString(" radius = ") +
    QString::number(to_double(mswc.offset()),'g', 3) + " ";
  emit radiusVal(radius_str);

  QString eps_hat_str = QString(" epsilon/radius = ") +
    QString::number(to_double(mswc.epsilon()/mswc.offset()),'g', 6) + " " +
    QString(" delta/epsilon = ") +
    QString::number(to_double(mswc.delta()/mswc.epsilon()),'g', 6) + " ";
  emit epsHatVal(eps_hat_str);

  EditParametersDialog::saved_size(mswc.kgon().size());
  QString kgon_size_str = QString(" kgon size = ") + QString::number(uint(mswc.kgon().size())) + " ";
  emit kgonSize(kgon_size_str);
}




#include "MainWindow.moc"
