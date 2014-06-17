#include "BezierExampleWindow.h"

#include <QGraphicsScene>

BezierExampleWindow::BezierExampleWindow(QWidget* parent):
  CGAL::Qt::DemosMainWindow( parent ),
  m_scene( new QGraphicsScene ),
  m_arrangement( NULL ),
  m_arrangementGraphicsItem( NULL ),
  ui( new Ui::BezierExampleWindow )
{
  setupUi( );
}

BezierExampleWindow::~BezierExampleWindow( )
{
  if ( m_arrangement )
    delete m_arrangement;
}

void BezierExampleWindow::setArrangement( Arrangement_2* arr )
{
  m_arrangement = arr;
  if ( m_arrangementGraphicsItem )
  {
    // TODO: Remove it from the scene
    delete m_arrangementGraphicsItem;
  }
  m_arrangementGraphicsItem =
    new BezierArrangementGraphicsItem( m_arrangement );

  // TODO: Add it to the scene
  m_scene->addItem( m_arrangementGraphicsItem );
}

Arrangement_2* BezierExampleWindow::arrangement( )
{
  return m_arrangement;
}

void BezierExampleWindow::clear( )
{
  // TODO: Remove it from the scene
}

void BezierExampleWindow::setupUi( )
{
  ui->setupUi( this );
  this->addNavigation( ui->graphicsView );
  ui->graphicsView->setScene( m_scene );
  m_scene->addEllipse( 0, 0, 100, 100 );
}

void BezierExampleWindow::on_actionQuit_triggered( )
{
  qApp->exit( );
}
