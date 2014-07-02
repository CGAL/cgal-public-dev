#ifndef BEZIER_TAB_WINDOW_H
#define BEZIER_TAB_WINDOW_H
#include <CGAL/Qt/DemosMainWindow.h>
#include "BezierDemoTraits.h"
#include "ui_BezierTabWindow.h"

template < class Arr_ >
class ArrangementDemoTab;

class BezierTabWindow : public CGAL::Qt::DemosMainWindow
{
  Q_OBJECT

public:
  typedef BezierDemoTraits DemoTraitsType;
  typedef DemoTraitsType::ArrTraitsType ArrTraitsType;
  typedef DemoTraitsType::ArrangementType ArrangementType;
  typedef ArrTraitsType::Curve_2 Curve_2;
  typedef ArrangementDemoTab< ArrangementType > TabType;

protected:
  ArrangementType* m_arr;
  TabType* m_tab;
  Ui::BezierTabWindow* ui;

public:
  BezierTabWindow(QWidget* parent = 0);
  virtual ~BezierTabWindow( );

  void load( const std::string& filename );

protected:
  void setupUi( );
};
#endif // BEZIER_TAB_WINDOW_H
