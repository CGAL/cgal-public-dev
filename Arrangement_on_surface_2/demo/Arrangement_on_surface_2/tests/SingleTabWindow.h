#ifndef SINGLE_TAB_WINDOW_H
#define SINGLE_TAB_WINDOW_H
#include <QObject>
#include "PolynomialParser.h"
#include <CGAL/Qt/DemosMainWindow.h>
#include "BezierDemoTraits.h"
#include "AlgebraicDemoTraits.h"
#include "ui_SingleTabWindow.h"

template < class Arr_ >
class ArrangementDemoTab;

class SingleTabWindow : public CGAL::Qt::DemosMainWindow
{
  Q_OBJECT

public:
  //typedef BezierDemoTraits DemoTraitsType;
  typedef AlgebraicDemoTraits DemoTraitsType;
  typedef DemoTraitsType::ArrTraitsType ArrTraitsType;
  typedef DemoTraitsType::ArrangementType ArrangementType;
  typedef ArrTraitsType::Curve_2 Curve_2;
  typedef ArrangementDemoTab< ArrangementType > TabType;

protected:
  ArrangementType* m_arr;
  TabType* m_tab;
  Ui::SingleTabWindow* ui;

public:
  SingleTabWindow(QWidget* parent = 0);
  virtual ~SingleTabWindow( );

  void load( const std::string& filename );
  void fitInView( const QRectF& rect );

protected:
  void setupUi( );
};
#endif // SINGLE_TAB_WINDOW_H
