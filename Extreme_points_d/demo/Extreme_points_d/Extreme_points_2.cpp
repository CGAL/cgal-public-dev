#include <fstream>
#include <vector>
#include <iterator>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Extreme_points_d.h>
#include <CGAL/Extreme_points_traits_d.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QFileDialog>
#include <QGraphicsRectItem>
#include <QGraphicsLineItem>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/PointsGraphicsItem.h>
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/Qt/utility.h>
  
// the two base classes
#include "ui_Extreme_points_2.h"
#include <CGAL/Qt/DemosMainWindow.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point_2;
typedef K::Vector_2                                         Vector_2;
typedef K::Iso_rectangle_2                                  Iso_rectangle_2;
typedef CGAL::Extreme_points_traits_d<Point_2>              EP_Traits_2;

class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Extreme_points_2
{
  Q_OBJECT
  
private:

  Iso_rectangle_2 square;
  CGAL::Extreme_points_d<EP_Traits_2> ep;
  CGAL::Qt::Converter<K> convert;
  std::vector<Point_2> points;
  std::vector<Point_2> extreme_points;
  QGraphicsScene scene;

  CGAL::Qt::PointsGraphicsItem<std::vector<Point_2> > * pgi;
  CGAL::Qt::PointsGraphicsItem<std::vector<Point_2> > * epgi;
  CGAL::Qt::GraphicsViewPolylineInput<K> * pi;

  template <typename G>
  void
  on_actionGenerate_triggered()
  {
    Iso_rectangle_2 isor = square;
    Point_2 center = CGAL::midpoint(isor[0], isor[2]);
    Vector_2 offset = center - CGAL::ORIGIN;
    double w = isor.xmax() - isor.xmin();
    double h = isor.ymax() - isor.ymin();
    double radius = (w<h) ? w/2 : h/2;

    G pg(radius);
    bool ok = false;
    const int number_of_points = 
      QInputDialog::getInteger(this, 
                               tr("Number of random points"),
                               tr("Enter number of random points"),
                               100,
                               0,
                               std::numeric_limits<int>::max(),
                               1,
                               &ok);

    if(!ok) {
      return;
    }

    // wait cursor
    QApplication::setOverrideCursor(Qt::WaitCursor);

    points.reserve(points.size() + number_of_points);
    for(int i = 0; i < number_of_points; ++i) {
      Point_2 p = *pg + offset;
      points.push_back(p);
      ep.insert(p);
      ++pg;
    }
    
    // default cursor
    QApplication::restoreOverrideCursor();
    emit(changed());
  }

public:
  MainWindow();

public slots:

  void on_actionClear_triggered();

  void processInput(CGAL::Object);

  void on_actionRecenter_triggered();
  void on_actionGeneratePointsOnCircle_triggered();
  void on_actionGeneratePointsInSquare_triggered();
  void on_actionGeneratePointsInDisc_triggered();
  void on_actionGeneratePointsOnSquare_triggered();
  void clear();

  void update_extreme_points();

signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow(), square(Point_2(-1, -1), Point_2(1,1)), ep(2)
{
  setupUi(this);

  // Add a GraphicItem for the point set
  pgi = new CGAL::Qt::PointsGraphicsItem<std::vector<Point_2> >(&points);
  epgi = new CGAL::Qt::PointsGraphicsItem<std::vector<Point_2> >(&extreme_points);

 
  QObject::connect(this, SIGNAL(changed()),
		   this, SLOT(update_extreme_points()));
  
  QObject::connect(this, SIGNAL(changed()),
                   pgi, SLOT(modelChanged()));
  
  QObject::connect(this, SIGNAL(changed()),
                   epgi, SLOT(modelChanged()));
  
  pgi->setVerticesPen(QPen(Qt::red, 2, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  epgi->setVerticesPen(QPen(Qt::black, 5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(pgi);
  scene.addItem(epgi);

  // 
  // Manual handling of actions
  //
  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   this, SLOT(close()));

 
  pi = new CGAL::Qt::GraphicsViewPolylineInput<K>(this, &scene, 1, false); // inputs a list with one point
  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));
   
  scene.installEventFilter(pi);
  //
  // Setup the scene and the view
  //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(-1, -1, 1, 1);
  this->graphicsView->setScene(&scene);

  // Uncomment the following line to get antialiasing by default.
//   actionUse_Antialiasing->setChecked(true);

  // Turn the vertical axis upside down
  this->graphicsView->scale(1, -1);
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Extreme_points_2.html");
  this->addAboutCGAL();


}


/* 
 *  Qt Automatic Connections
 *  http://doc.trolltech.com/4.4/designer-using-a-component.html#automatic-connections
 * 
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */



void
MainWindow::processInput(CGAL::Object o)
{
  std::list<Point_2> input;
  if(CGAL::assign(input, o)){
    if(input.size() == 1) {
      Point_2 p = input.front();
      points.push_back(p);
      ep.insert(p);
    }
    emit(changed());
  }

}

void
MainWindow::update_extreme_points()
{
    extreme_points.clear();
    ep.extreme_points(std::back_inserter(extreme_points));
}


void
MainWindow::on_actionClear_triggered()
{
  clear();
  emit(changed());
}

void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(convert(square));
  this->graphicsView->fitInView(convert(square), Qt::KeepAspectRatio);  
}

void
MainWindow::on_actionGeneratePointsOnCircle_triggered()
{
  typedef CGAL::Random_points_on_circle_2<Point_2> Generator;
  on_actionGenerate_triggered<Generator>();
}


void
MainWindow::on_actionGeneratePointsInSquare_triggered()
{
  typedef CGAL::Random_points_in_square_2<Point_2> Generator;
  on_actionGenerate_triggered<Generator>();
}


void
MainWindow::on_actionGeneratePointsInDisc_triggered()
{
  typedef CGAL::Random_points_in_disc_2<Point_2> Generator;
  on_actionGenerate_triggered<Generator>();
}

void
MainWindow::on_actionGeneratePointsOnSquare_triggered()
{
    typedef CGAL::Random_points_on_square_2<Point_2> Generator;
    on_actionGenerate_triggered<Generator>();
}


void
MainWindow::clear()
{
  points.clear();
  extreme_points.clear();
  ep.clear();
}


#include "Extreme_points_2.moc"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

  app.setOrganizationDomain("ethz.ch");
  app.setOrganizationName("ETH Zürich");
  app.setApplicationName("Extreme_points_2 demo");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Extreme_points_2);
  Q_INIT_RESOURCE(Input);
  Q_INIT_RESOURCE(CGAL);

  MainWindow mainWindow;
  mainWindow.show();
  mainWindow.on_actionRecenter_triggered();
  return app.exec();
}
