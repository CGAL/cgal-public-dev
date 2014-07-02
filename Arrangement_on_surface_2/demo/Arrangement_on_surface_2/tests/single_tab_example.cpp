#include <QApplication>
#include "SingleTabWindow.h"

int main( int argc, char *argv[] )
{
  QApplication app( argc, argv );

  SingleTabWindow window;
  //window.load( "Bezier.dat" );
  window.show( );

  return app.exec( );
}
