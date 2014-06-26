#ifndef BEZIER_TAB_WINDOW_H
#define BEZIER_TAB_WINDOW_H
#include <CGAL/Qt/DemosMainWindow.h>
#include "BezierTabTypes.h"
#include "ui_BezierTabWindow.h"
#include "ArrangementDemoTab.h"

class BezierTabWindow : public CGAL::Qt::DemosMainWindow
{
  Q_OBJECT

  typedef ArrangementDemoTab< Bezier_arrangement_2 > TabType;

public:
  BezierTabWindow(QWidget* parent = 0);
  virtual ~BezierTabWindow( );

  void load( const std::string& filename );

protected:
  void setupUi( );

  Bezier_arrangement_2* m_arr;
  TabType* m_tab;

  Ui::BezierTabWindow* ui;
};
#endif // BEZIER_TAB_WINDOW_H
