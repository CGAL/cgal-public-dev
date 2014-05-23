#ifndef BEZIER_EXAMPLE_WINDOW_H
#define BEZIER_EXAMPLE_WINDOW_H
#include <CGAL/Qt/DemosMainWindow.h>
#include "ui_BezierExampleWindow.h"

class BezierExampleWindow : public CGAL::Qt::DemosMainWindow
{
  Q_OBJECT
public:
  BezierExampleWindow(QWidget* parent = 0);

public slots:
  void on_actionQuit_triggered( );

protected:
  void setupUi( );

protected:
  Ui::BezierExampleWindow* ui;

};
#endif // BEZIER_EXAMPLE_WINDOW_H
