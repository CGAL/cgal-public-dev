#include "BezierTabWindow.h"
#include "Utils.h"
#include "ArrangementDemoTab.h"

BezierTabWindow::BezierTabWindow(QWidget* parent):
  CGAL::Qt::DemosMainWindow( parent ),
  m_arr( new ArrangementType ),
  m_tab( new TabType( m_arr ) ),
  ui( new Ui::BezierTabWindow )
{
  setupUi( );
}

BezierTabWindow::~BezierTabWindow( )
{

}

void BezierTabWindow::load( const std::string& filename )
{
  LoadArrFromFile< DemoTraitsType > load_arr_from_file;
  load_arr_from_file( filename, m_arr );
}

void BezierTabWindow::setupUi( )
{
  ui->setupUi( this );
  ui->tabWidget->addTab( m_tab, "Bezier" );

  this->addNavigation( m_tab->getView( ) );
}
