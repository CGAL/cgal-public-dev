#include "window.h"
#include <QApplication>

#include <CGAL/Qt/resources.h>
#include <QMimeData>

int main(int argc, char **argv)
{
	srand(0);
  QApplication app(argc, argv);
  app.setApplicationName("Smooth Tests");

#if (QT_VERSION >= QT_VERSION_CHECK(5, 3, 0))
  app.setAttribute(Qt::AA_UseDesktopOpenGL);
#endif

  // Import resources from libCGALQt (Qt5).
  CGAL_QT_INIT_RESOURCES;

  MainWindow mainWindow;
  mainWindow.show();

  return app.exec();
}
