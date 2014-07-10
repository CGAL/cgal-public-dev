#include <QApplication>
#include "SingleTabWindow.h"

int main( int argc, char *argv[] )
{
  QApplication app( argc, argv );
  if ( argc < 2 )
  {
    std::cout << "give me data\n";
    return 0;
  }

  SingleTabWindow window;
  window.fitInView( QRectF( -20, -20, 50, 50 ) );
  window.load( argv[1] );
  window.show( );

  return app.exec( );
}
