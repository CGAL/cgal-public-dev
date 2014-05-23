#include "BezierExampleWindow.h"
#include <QApplication>

int main( int argc, char* argv[] )
{
  QApplication app( argc, argv );

  BezierExampleWindow window;
  window.show( );

  return app.exec( );
}
