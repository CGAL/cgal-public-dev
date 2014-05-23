#include "BezierExampleWindow.h"

BezierExampleWindow::BezierExampleWindow(QWidget* parent):
  CGAL::Qt::DemosMainWindow( parent ),
  ui( new Ui::BezierExampleWindow )
{
  setupUi( );
}

void BezierExampleWindow::setupUi( )
{
  ui->setupUi( this );
}

void BezierExampleWindow::on_actionQuit_triggered( )
{
  qApp->exit( );
}
