#include "BezierTabWindow.h"

BezierTabWindow::BezierTabWindow(QWidget* parent):
  CGAL::Qt::DemosMainWindow( parent ),
  m_arr( new Bezier_arrangement_2 ),
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
  std::ifstream inputFile( filename.c_str( ) );
  // Read the curves from the input file.
  unsigned int               n_curves;
  std::list<Bezier_curve_2>  curves;
  Bezier_curve_2             B;
  unsigned int               k;

  inputFile >> n_curves;
  for (k = 0; k < n_curves; k++) {
    // Read the current curve (specified by its control points).
    inputFile >> B;
    curves.push_back (B);

    std::cout << "B = {" << B << "}" << std::endl;
  }

  // Construct the arrangement.
  insert (*m_arr, curves.begin(), curves.end());

  // Print the arrangement size.
  std::cout << "The arrangement size:" << std::endl
            << "   V = " << m_arr->number_of_vertices()
            << ",  E = " << m_arr->number_of_edges()
            << ",  F = " << m_arr->number_of_faces() << std::endl;
  //m_tab->getArrangementGraphicsItem( )->update( );
  m_tab->setArrangement( m_arr );
}

void BezierTabWindow::setupUi( )
{
  ui->setupUi( this );
  ui->tabWidget->addTab( m_tab, "Bezier" );

  this->addNavigation( m_tab->getView( ) );
}
