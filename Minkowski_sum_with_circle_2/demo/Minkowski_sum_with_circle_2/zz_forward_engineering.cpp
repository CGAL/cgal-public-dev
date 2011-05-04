#include "include/MainWindow.h"
#include "include/EditParametersDialog.h"


void MainWindow::toggle_forward_constructions(bool checked)
{
  for(int data_type = MSWC::Data_kgon; data_type != MSWC::Data_inner_skeleton; ++data_type)
  {
    if(checked && data_type != MSWC::Data_kgon_induced_circles)
      enable_and_hide(static_cast<Data_type>(data_type));
    else
      disable_and_hide(static_cast<Data_type>(data_type));
  }

  sloped_kgon_sum_gi->setVisible(false);

  actionShowKgonSkippedCircles->setEnabled(false);
  actionShowKgonSkippedCircles->setChecked(false);
  skipped_kgon_circles_gi->setVisible(false);

}
  
void
MainWindow::on_actionViewKgon_toggled(bool checked)
{
  std::clog << "actionViewKgon " << checked << std::endl;

  set_visible(checked, MSWC::Data_kgon);
}

void
MainWindow::on_actionConstructKgonInducedCircles_toggled(bool checked)
{
  // ConstructKgonInducedCircles is checkable (and checked) only if minkowski sum was computed,
  // that is if actionConstructKgonSum or actionConstructKgonOffset are active
  std::clog << "actionConstructKgonInducedCircles " << checked << std::endl;

  set_visible(checked, MSWC::Data_kgon_induced_circles);
}

void
MainWindow::on_actionShowKgonSkippedCircles_toggled(bool checked)
{
  // ConstructKgonInducedCircles is checkable only if ConstructKgonInducedCircles is
  std::clog << "actionShowKgonSkippedCircles " << checked << std::endl;

  if(checked)
  {
    update_if_needed(MSWC::Data_kgon_induced_circles);
    update_if_needed(MSWC::Data_kgon_offset_polygon);
  }

  skipped_kgon_circles_gi->setVisible(checked);
  if(checked) skipped_kgon_circles_gi->modelChanged();

  intersecting_kgon_circles_gi->setVisible(checked);
  if(checked) intersecting_kgon_circles_gi->modelChanged();
}


void
MainWindow::on_actionConstructKgonSum_toggled(bool checked)
{
  std::clog << "actionConstructKgonSum " << checked << std::endl;

  set_visible(checked, MSWC::Data_kgon_sum);

  sloped_kgon_sum_gi->setVisible(checked);
  if(checked) sloped_kgon_sum_gi->modelChanged();

  // ConstructKgonInducedCircles is checkable (and checked) only if minkowski sum was computed,
  // that is if actionConstructKgonSum or actionConstructKgonOffset are active
  bool enableInducedCircles = checked || actionConstructKgonOffset->isChecked();
  if (!enableInducedCircles) {
    actionConstructKgonInducedCircles->setChecked(false);
    actionShowKgonSkippedCircles->setChecked(false);
  }
  actionConstructKgonInducedCircles->setEnabled(enableInducedCircles);
  actionShowKgonSkippedCircles->setEnabled(enableInducedCircles);
}

void
MainWindow::on_actionConstructKgonOffset_toggled(bool checked)
{
  std::clog << "actionConstructKgonOffset " << checked << std::endl;

  set_visible(checked, MSWC::Data_kgon_offset_polygon);

  // ConstructKgonInducedCircles is checkable (and checked) only if minkowski sum was computed,
  // that is if actionConstructKgonSum or actionConstructKgonOffset are active
  bool enableInducedCircles = checked || actionConstructKgonSum->isChecked();
  if (!enableInducedCircles) {
    actionConstructKgonInducedCircles->setChecked(false);
    actionShowKgonSkippedCircles->setChecked(false);
  }
  actionConstructKgonInducedCircles->setEnabled(enableInducedCircles);
  actionShowKgonSkippedCircles->setEnabled(enableInducedCircles);
}


void
MainWindow::on_actionConstructExactOffset_toggled(bool checked)
{
  std::clog << "actionConstructExactOffset " << checked << std::endl;

  set_visible(checked, MSWC::Data_exact_offset_polygon);
}

void
MainWindow::on_actionConstructAppOffset_toggled(bool checked)
{
  std::clog << "actionConstructAppOffset " << checked << std::endl;

  set_visible(checked, MSWC::Data_approximate_offset_polygon);
}

