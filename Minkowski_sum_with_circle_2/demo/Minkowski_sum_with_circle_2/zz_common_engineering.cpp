#include "include/MainWindow.h"


void
MainWindow::on_actionViewPolygon_toggled(bool checked)
{
  LOG_DEBUG << "actionViewPolygon " << checked << std::endl;

  set_visible(checked, MSWC::Data_polygon);
}


void
MainWindow::on_actionSavePolygon_triggered()
{
  LOG_DEBUG << "SavePolygon " << std::endl;

  QString fileName = QFileDialog::getSaveFileName(this,
              tr("Save polygon file"),
              ".");
  if(! fileName.isEmpty()){
    std::ofstream ofs(qPrintable(fileName));
    ofs << mswc.polygon();
  }
}

void
MainWindow::on_actionRecenter_triggered()
{
  LOG_DEBUG << "Recenter " << std::endl;
  QRectF bbox = polygon_gi->boundingRect();

  //bbox = QRectF(-4.4,-4.4,13.8,13.8);

  double margin = to_double(mswc.offset() + mswc.epsilon());
  margin *= 4;
  bbox.adjust(-margin, -margin, margin, margin);

  this->graphicsView->setSceneRect(bbox);
  this->graphicsView->fitInView(bbox, Qt::KeepAspectRatio);
}



void MainWindow::on_actionConstructionType_toggled(bool checked)
{
  LOG_DEBUG << "actionConstructionType " << checked << std::endl;

  // reset actions and constructions of the other mode
  m_mode = static_cast<int>(checked);
  toggle_forward_constructions(checked);
  toggle_reverse_constructions(!checked);

  // set current mode
  m_mode = static_cast<int>(!checked);
  mswc.forward_construction(checked);

  // remove or add back "forward-only" mode toolbars
//   if(checked)
//   {
//     this->addToolBar(Qt::TopToolBarArea, compareToolBar);
//     this->addToolBar(Qt::TopToolBarArea, constructionToolBar);
//     retranslateUi(this);
//   }
//   else
//   {
//     this->removeToolBar(compareToolBar);
//     this->removeToolBar(constructionToolBar);
//   }

  // TODO: when on - can change kgon type and size
  // TODO: when off - can't change kgon type and size
}
