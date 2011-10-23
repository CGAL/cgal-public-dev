#include "include/MainWindow.h"


void MainWindow::toggle_reverse_constructions(bool checked)
{
  for(int data_type = MSWC::Data_inner_skeleton; data_type != MSWC::Data_type_size; ++data_type)
  {
    check(false, static_cast<Data_type>(data_type));

    mswc.needs_update(static_cast<Data_type>(data_type));
  }
}

// skeleton constructions
void
MainWindow::on_actionConstructInnerSkeleton_toggled(bool checked)
{
  std::clog << "actionConstructInnerSkeleton " << checked << std::endl;

  set_visible(checked, MSWC::Data_inner_skeleton);
}

void
MainWindow::on_actionConstructOuterSkeleton_toggled(bool checked)
{
  std::clog << "actionConstructOuterSkeleton " << checked << std::endl;

  set_visible(checked, MSWC::Data_outer_skeleton);
}

// exact approximability constructions
void
MainWindow::on_actionConstructOuterCore_toggled(bool checked)
{
  std::clog << "actionConstructOuterCore " << checked << std::endl;

  set_visible(checked, MSWC::Data_outer_core);
}

void
MainWindow::on_actionConstructInnerCore_toggled(bool checked)
{
  std::clog << "actionConstructInnerCore " << checked << std::endl;

  set_visible(checked, MSWC::Data_inner_core);
}


void
MainWindow::on_actionCheckApproximability_toggled(bool checked)
{
  std::clog << "actionCheckApproximability " << checked << std::endl;

  set_visible(checked, MSWC::Data_approximability);

  //if(!mswc.approximability())
  // show non_approximable data
  {
    GI* app_neg_gi = mswc.forward_construction()? approximability_fwd_neg_gi: approximability_bwd_neg_gi;
    GI* other_neg_gi = mswc.forward_construction()? approximability_bwd_neg_gi: approximability_fwd_neg_gi;
    app_neg_gi->setVisible(checked);
    if(checked)
    {
      app_neg_gi->modelChanged();
      other_neg_gi->setVisible(!checked);
    }
  }
}


// rational approximability constructions

void
MainWindow::on_actionConstructInnerTolerance_toggled(bool checked)
{
  std::clog << "actionConstructInnerTolerance " << checked << std::endl;

  set_visible(checked, MSWC::Data_approximate_inner_eps);
}

void
MainWindow::on_actionConstructOuterTolerance_toggled(bool checked)
{
  std::clog << "actionConstructOuterTolerance " << checked << std::endl;

  set_visible(checked, MSWC::Data_approximate_outer_eps);
}

void
MainWindow::on_actionConstructOuterCoreApp_toggled(bool checked)
{
  std::clog << "actionConstructOuterCoreApp " << checked << std::endl;

  set_visible(checked, MSWC::Data_approximate_outer_core);
}

void
MainWindow::on_actionConstructInnerCoreApp_toggled(bool checked)
{
  std::clog << "actionConstructInnerCoreApp " << checked << std::endl;

  set_visible(checked, MSWC::Data_approximate_inner_core);
}

void
MainWindow::on_actionConstructOuterCoreSum_toggled(bool checked)
{
  std::clog << "actionConstructOuterCoreSum " << checked << std::endl;

  set_visible(checked, MSWC::Data_approximate_outer_kgon_sum);
}

void
MainWindow::on_actionConstructInnerCoreSum_toggled(bool checked)
{
  std::clog << "actionConstructInnerCoreAppSum " << checked << std::endl;

  set_visible(checked, MSWC::Data_approximate_inner_kgon_sum);
}

void
MainWindow::on_actionConstructInnerOffsetApp_toggled(bool checked)
{
  std::clog << "actionConstructInnerOffsetApp " << checked << std::endl;

  set_visible(checked, MSWC::Data_approximate_inner_offset);
}

void
MainWindow::on_actionConstructOuterOffsetApp_toggled(bool checked)
{
  std::clog << "actionConstructOuterOffsetApp " << checked << std::endl;

  set_visible(checked, MSWC::Data_approximate_outer_offset);
}

void
MainWindow::on_actionConstructInnerOffset_toggled(bool checked)
{
  std::clog << "actionConstructInnerOffset " << checked << std::endl;

//  set_visible(checked, MSWC::Data_approximate_inner_offset_exact);
}

void
MainWindow::on_actionConstructOuterOffset_toggled(bool checked)
{
  std::clog << "actionConstructOuterOffset " << checked << std::endl;

//  set_visible(checked, MSWC::Data_approximate_outer_offset_exact);
}

void
MainWindow::on_actionConstructInnerDiffApp_toggled(bool checked)
{
  std::clog << "actionConstructInnerDiff " << checked << std::endl;

  set_visible(checked, MSWC::Data_approximate_inner_decision);
}

void
MainWindow::on_actionConstructOuterDiffApp_toggled(bool checked)
{
  std::clog << "actionConstructOuterDiff " << checked << std::endl;

  set_visible(checked, MSWC::Data_approximate_outer_decision);
}

void
MainWindow::on_actionYesApproximability_toggled(bool checked)
{
  std::clog << "actionYesApproximability " << checked << std::endl;

//  set_visible(checked, MSWC::Data_approximate_inner_decision);
}

void
MainWindow::on_actionNoApproximability_toggled(bool checked)
{
  std::clog << "actionNoApproximability " << checked << std::endl;

//  set_visible(checked, MSWC::Data_approximate_outer_decision);
}

void
MainWindow::on_actionCheckApproximabilityApp_toggled(bool checked)
{
  std::clog << "actionCheckApproximabilityApp " << checked << std::endl;

  set_visible(checked, MSWC::Data_approximate_approximability);
}


// search activation
void
MainWindow::on_actionSearchMinEpsilon_triggered()
{
  std::clog << "actionSearchMinEpsilon " << std::endl;

  bool update_params = true;
  mswc.compute_search_epsilon(update_params);
  //mswc.update(MSWC::Data_search_epsilon);

  if(update_params)
    update_if_needed();
}


// auxilary constructions

void
MainWindow::compute_lines()
{
  std::clog << "RE: Compute lines " << std::endl;

  _m_eps_lines[0].clear();
  _m_eps_lines[1].clear();
  _m_rme_lines[0].clear();
  _m_rme_lines[1].clear();
  _m_rad_lines[0].clear();
  _m_rad_lines[1].clear();
  _m_rpe_lines[0].clear();
  _m_rpe_lines[1].clear();

  const Polygon_2& polygon = mswc.polygon();

  for (Polygon_2::Edge_const_iterator eit = polygon.edges_begin();
       eit != polygon.edges_end();
       ++eit) {

    Line_2 line = eit->supporting_line();
    NT a = line.a();
    NT b = line.b();
    NT c = line.c();
    
    NT rme = mswc.offset() - mswc.epsilon();
    NT rpe = mswc.offset() + mswc.epsilon();

    Kernel::Vector_2 vec = line.to_vector();
    
    // TODO remove double
    double angle = atan2(CGAL::to_double(vec.y()), CGAL::to_double(vec.x()));

    NT cosalpha(std::cos(angle));

    // move eps away!
    _m_eps_lines[0].push_back(Line_2(a,b,c+b*(mswc.epsilon()/cosalpha)));
    _m_eps_lines[1].push_back(Line_2(a,b,c-b*(mswc.epsilon()/cosalpha)));

    // move rme away!
    _m_rme_lines[0].push_back(Line_2(a,b,c+b*(rme/cosalpha)));
    _m_rme_lines[1].push_back(Line_2(a,b,c-b*(rme/cosalpha)));

    // move rad away!
    _m_rad_lines[0].push_back(Line_2(a,b,c+b*(mswc.offset()/cosalpha)));
    _m_rad_lines[1].push_back(Line_2(a,b,c-b*(mswc.offset()/cosalpha)));
    
    // move rpe away!
    _m_rpe_lines[0].push_back(Line_2(a,b,c+b*(rpe/cosalpha)));
    _m_rpe_lines[1].push_back(Line_2(a,b,c-b*(rpe/cosalpha)));

  }
}

void
MainWindow::compute_circles()
{
  _m_eps_circles.clear();
  _m_rme_circles.clear();
  _m_rad_circles.clear();
  _m_rpe_circles.clear();

  std::clog << "RE: Compute circles " << std::endl;

  NT eps2 = mswc.epsilon() * mswc.epsilon();
  NT rme2 = (mswc.offset() - mswc.epsilon());
  rme2 = rme2 * rme2;
  NT rad2 = mswc.offset() *mswc.offset();
  NT rpe2 = (mswc.offset() + mswc.epsilon());
  rpe2 = rpe2 * rpe2;

  const Polygon_2& polygon = mswc.polygon();
  for (Polygon_2::Vertex_const_iterator vit = polygon.vertices_begin();
      vit != polygon.vertices_end();
      ++vit) {
    _m_eps_circles.push_back(Circle_2(*vit, eps2));
    _m_rme_circles.push_back(Circle_2(*vit, rme2));
    _m_rad_circles.push_back(Circle_2(*vit, rad2));
    _m_rpe_circles.push_back(Circle_2(*vit, rpe2));
  }
}

void
MainWindow::on_actionShow_all_triggered()
{
  std::clog << "RE: Show all " << std::endl;

  actionLines_show_all->trigger();
  actionCircles_show_all->trigger();

}

void
MainWindow::on_actionHide_all_triggered()
{
  std::clog << "RE: Hide all " << std::endl;

  actionLines_hide_all->trigger();
  actionCircles_hide_all->trigger();
  
}

void
MainWindow::on_actionLines_show_all_triggered()
{
  std::clog << "RE: Lines show all " << std::endl;

  actionLines_eps->setChecked(true);
  actionLines_r_minus_eps->setChecked(true);
  actionLines_r->setChecked(true);
  actionLines_r_plus_eps->setChecked(true);
}

void
MainWindow::on_actionLines_hide_all_triggered()
{
  std::clog << "RE: Lines hide all " << std::endl;

  actionLines_eps->setChecked(false);
  actionLines_r_minus_eps->setChecked(false);
  actionLines_r->setChecked(false);
  actionLines_r_plus_eps->setChecked(false);
}

void
MainWindow::on_actionLines_eps_toggled(bool checked)
{
  std::clog << "RE: Lines_eps " << checked << std::endl;

  (epslgi[0])->setVisible(checked);
  (epslgi[1])->setVisible(checked);
}

void
MainWindow::on_actionLines_r_minus_eps_toggled(bool checked)
{
  std::clog << "RE: Lines_r_minus_eps " << checked << std::endl;

  (rmelgi[0])->setVisible(checked);
  (rmelgi[1])->setVisible(checked);
}

void
MainWindow::on_actionLines_r_toggled(bool checked)
{
  std::clog << "RE: Lines_r " << checked << std::endl;

  (radlgi[0])->setVisible(checked);
  (radlgi[1])->setVisible(checked);
}

void
MainWindow::on_actionLines_r_plus_eps_toggled(bool checked)
{
  std::clog << "RE: Lines_r_plus_eps " << checked << std::endl;

  (rpelgi[0])->setVisible(checked);
  (rpelgi[1])->setVisible(checked);
}

void
MainWindow::on_actionCircles_show_all_triggered()
{
  std::clog << "RE: Circles show all " << std::endl;

  actionCircles_eps->setChecked(true);
  actionCircles_r_minus_eps->setChecked(true);
  actionCircles_r->setChecked(true);
  actionCircles_r_plus_eps->setChecked(true);

}

void
MainWindow::on_actionCircles_hide_all_triggered()
{
  std::clog << "RE: Circles hide all " << std::endl;

  actionCircles_eps->setChecked(false);
  actionCircles_r_minus_eps->setChecked(false);
  actionCircles_r->setChecked(false);
  actionCircles_r_plus_eps->setChecked(false);

}

void
MainWindow::on_actionCircles_eps_toggled(bool checked)
{
  std::clog << "RE: Circles_eps " << checked << std::endl;

  epscgi->setVisible(checked);
}

void
MainWindow::on_actionCircles_r_minus_eps_toggled(bool checked)
{
  std::clog << "RE: Circles_r_minus_eps " << checked << std::endl;

  rmecgi->setVisible(checked);
}

void
MainWindow::on_actionCircles_r_toggled(bool checked)
{
  std::clog << "RE: Circles_r " << checked << std::endl;

  radcgi->setVisible(checked);
}

void
MainWindow::on_actionCircles_r_plus_eps_toggled(bool checked)
{
  std::clog << "RE: Circles_r_plus_eps " << checked << std::endl;

  rpecgi->setVisible(checked);
}

