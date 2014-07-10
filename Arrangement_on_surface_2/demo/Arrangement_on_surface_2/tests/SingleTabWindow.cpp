#include <QObject>
#include "SingleTabWindow.h"
#include "Utils.h"
#include "ArrangementDemoTab.h"

SingleTabWindow::SingleTabWindow(QWidget* parent):
  CGAL::Qt::DemosMainWindow( parent ),
  m_arr( new ArrangementType ),
  m_tab( new TabType( m_arr ) ),
  ui( new Ui::SingleTabWindow )
{
  setupUi( );
}

SingleTabWindow::~SingleTabWindow( )
{

}

void SingleTabWindow::load( const std::string& filename )
{
  LoadArrFromFile< DemoTraitsType > load_arr_from_file;
  load_arr_from_file( filename, m_arr );
}

void SingleTabWindow::fitInView( const QRectF& rect )
{
  m_tab->getView( )->fitInView( rect, ::Qt::KeepAspectRatio );
}

void SingleTabWindow::setupUi( )
{
  ui->setupUi( this );
  ui->tabWidget->addTab( m_tab, DemoTraitsType::Name.c_str() );

  this->addNavigation( m_tab->getView( ) );
}
