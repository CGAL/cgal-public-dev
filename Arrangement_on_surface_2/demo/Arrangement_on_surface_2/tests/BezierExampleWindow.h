#ifndef BEZIER_EXAMPLE_WINDOW_H
#define BEZIER_EXAMPLE_WINDOW_H
#include <CGAL/Qt/DemosMainWindow.h>
#include "ui_BezierExampleWindow.h"
#include "BezierExampleTypes.h"
#include "PointSetGraphicsItem.h"

class QGraphicsScene;

class BezierExampleWindow : public CGAL::Qt::DemosMainWindow
{
  Q_OBJECT
public:
  BezierExampleWindow(QWidget* parent = 0);
  virtual ~BezierExampleWindow( );

  /**
  Display the arrangement. Passes ownership of the item to the window.
  */
  void setArrangement(
    Arrangement_2* item );

  /**
  Get a pointer to the arrangement. Don't delete it -- it belongs to us.
  */
  Arrangement_2* arrangement( );

  /**
  Clear the window. Takes ownership of any attached arrangement.
  */
  void clear( );

  /**
  Iterators have value-type of std::pair< double, double >
  */
  template < class InputIterator >
  void drawXYPairs( InputIterator begin, InputIterator end )
  {
    std::cout << "TODO: draw points\n";
    std::vector< QPointF > points;
    for (InputIterator it = begin; it != end; ++it )
    {
      points.push_back( QPointF( it->first, it->second ) );
      ++it;
    }
    PointSetGraphicsItem *point_set = new PointSetGraphicsItem( points.begin(), points.end() );
    m_scene->addItem( point_set );
  }

public slots:
  void on_actionQuit_triggered( );

protected:
  void setupUi( );

protected:
  QGraphicsScene* m_scene;
  Arrangement_2* m_arrangement;

  // internal
  BezierArrangementGraphicsItem* m_arrangementGraphicsItem;

  Ui::BezierExampleWindow* ui;
};
#endif // BEZIER_EXAMPLE_WINDOW_H
