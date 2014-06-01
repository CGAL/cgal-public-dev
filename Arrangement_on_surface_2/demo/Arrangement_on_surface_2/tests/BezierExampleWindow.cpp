#include "BezierExampleWindow.h"

#include <QGraphicsScene>

BezierExampleWindow::BezierExampleWindow(QWidget* parent):
  CGAL::Qt::DemosMainWindow( parent ),
  m_scene( new QGraphicsScene ),
  m_arr( NULL ),
  ui( new Ui::BezierExampleWindow )
{
  setupUi( );
}

BezierExampleWindow::~BezierExampleWindow( )
{
    if ( m_arr )
        delete m_arr;
}

void BezierExampleWindow::setArrangement( Arrangement_2* arr )
{
    m_arr = arr;
}

Arrangement_2* BezierExampleWindow::arr( )
{
    return m_arr;
}

void BezierExampleWindow::clear( )
{
    m_arr = NULL;
}

void BezierExampleWindow::setupUi( )
{
  ui->setupUi( this );
  ui->graphicsView->setScene( m_scene );
  m_scene->addEllipse( 0, 0, 100, 100 );
}

void BezierExampleWindow::on_actionQuit_triggered( )
{
  qApp->exit( );
}
