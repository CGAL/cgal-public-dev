#include <QApplication>
#include "BezierTabWindow.h"

int main( int argc, char *argv[] )
{
  QApplication app( argc, argv );

  BezierTabWindow window;
  window.load( "Bezier.dat" );
  window.show( );

  return app.exec( );
}
