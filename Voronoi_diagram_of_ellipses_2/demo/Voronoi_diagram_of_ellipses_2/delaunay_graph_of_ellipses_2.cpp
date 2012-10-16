//    (c) 2007-2009 National and Kapodistrian University of Athens
//    (c) 2009-2011 INRIA Nancy
//    (c) 2011-2012 National and Kapodistrian University of Athens
//
//    Author:
//        George M. Tzoumas <geortz@gmail.com>
//
//  THIS SOURCE CODE IS CONSIDERED EXPERIMENTAL AND PROVIDED WITHOUT ANY WARRANTY
//  YOU ARE FREE TO USE THIS CODE FOR NON-COMMERCIAL PURPOSES
//  BUT PROPER CREDITS MUST BE GIVEN
//  MORE DETAILS ON THE LICENCE WHEN OFFICIALLY RELEASED

//#define CGAL_DISABLE_ROUNDING_MATH_CHECK
//#define CGAL_NO_ASSERTIONS
//#define NO_VALIDITY_CHECK

//#define CGAL_OBJECT_CACHE_DISABLE

//#define CGAL_GMPFR_NO_REFCOUNT
//#define WITH_BOOST_INTERVAL

//#define DGET_VERBOSE
//#define VERBOSE 3
//#define VERBOSE 2
//#define COMPARE

//#define WITH_PSEUDOCIRCLES_PATCH

#ifdef COMPARE
#define CGAL_VORELL_CIRCLIFY
#endif


#include <fstream>
#include<boost/shared_ptr.hpp>
// CGAL headers

#include <CGAL/basic.h>
#include <CGAL/Gmpq.h>
#ifdef WITH_BOOST_INTERVAL
#include <CGAL/Boost_interval_Gmpfr.h>
#else
#include <CGAL/Gmpfi.h>
#endif

//#if defined (__APPLE__)
//#define __USE_ISOC99
//#endif

#include <CGAL/Ellipse_traits.h>
#include <CGAL/Ellipse_2.h>
#include <CGAL/Bitangent.h>

#include <CGAL/Delaunay_graph_of_ellipses_2.h>

#ifdef COMPARE

#include <CGAL/Apollonius_vs_ellipses_traits_2.h>

#else

#include <CGAL/Delaunay_graph_of_ellipses_traits_2.h>
#include <CGAL/Delaunay_graph_of_ellipses_profiled_traits_2.h>

#endif

#include <CGAL/Real_timer.h>
#include <CGAL/Random.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QFileDialog>
#include <QDir>
#include <QGraphicsLineItem>
#include <QPrinter>
#include <QTimer>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/EllipseGraphicsItem.h>
#include <CGAL/Qt/EllipseBisectorGraphicsItem.h>
#include <CGAL/Qt/VoronoiEllipsesGraphicsItem.h>
#include <CGAL/Qt/DelaunayEllipsesGraphicsItem.h>
#include <CGAL/Qt/GraphicsViewEllipseInput.h>
  
// the two base classes
#include "ui_Delaunay_graph_of_ellipses_2.h"
#include <CGAL/Qt/DemosMainWindow.h>


typedef CGAL::Gmpq QQ;

#ifdef WITH_BOOST_INTERVAL
typedef CGAL::Boost_interval_Gmpfr INT;
#else
typedef CGAL::Gmpfi INT;
#endif

typedef CGAL::Ellipse_traits<QQ, INT> ET;
typedef CGAL::Ellipse_2<ET> Ellipse_2;
typedef CGAL::VORELL::Bitangent<ET> Bitangent;
typedef CGAL::Ellipse_bisector_2<ET> Ellipse_bisector_2;

#ifdef COMPARE
typedef CGAL::Apollonius_vs_ellipses_traits_2<ET> GT;
#else
typedef CGAL::Delaunay_graph_of_ellipses_traits_2<ET> DGT;
typedef CGAL::Delaunay_graph_of_ellipses_profiled_traits_2<DGT> GT;
#endif
typedef CGAL::Delaunay_graph_of_ellipses_2< GT > DG;

typedef CGAL::Qt::EllipseGraphicsItem<ET> EGI;


class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Delaunay_graph_of_ellipses_2
{
  Q_OBJECT
  
private:

  QGraphicsScene *scene;  
  CGAL::Qt::GraphicsViewEllipseInput<ET> *pi;
  std::vector<Ellipse_2> ells;
  QGraphicsItemGroup *ellg;

  CGAL::Qt::VoronoiEllipsesGraphicsItem<DG> *vor;
  CGAL::Qt::DelaunayEllipsesGraphicsItem<DG> *del;
  QGraphicsPixmapItem *bgpix;

  DG graph;

  int argc;
  char **argv;

  QString fname;
  bool validity_check;
  bool no_intersecting;
  double ell_stroke, bis_stroke;
  bool grayscale;
  int ell_size_factor;
  int skip, limit;
  CGAL::Real_timer vtimer;
  CGAL::Real_timer gtimer;
  CGAL::Random rng;
  
public:          
  bool benchmark;
  
private:

  void parse_command_line() {
      int i;

      validity_check = true;
      ell_stroke = bis_stroke = 0.0;
      ell_size_factor = 5;
      no_intersecting = false;
      benchmark = false;
      grayscale = false;
      skip = 0; limit = -1;

#define ARGIS(x) strcmp(argv[i], x) == 0

      for (i = 1; i < argc; i++) {
          if (ARGIS("-i")) {
              if (i+1 < argc) fname = QString(argv[++i]);
              else std::cerr << "-i requires an argument\n";
          } else if (ARGIS("--no-validity-check")) validity_check = false;
          else if (ARGIS("--no-graph")) actionCompute_graph->setChecked(false);
          else if (ARGIS("--ellipse-stroke")) {
              if (i+1 < argc) ell_stroke = atof(argv[++i]);
              else std::cerr << "--ellipse-stroke requires an argument\n";
          } else if (ARGIS("--bisector-stroke")) {
              if (i+1 < argc) bis_stroke = atof(argv[++i]);
              else std::cerr << "--bisector-stroke requires an argument\n";
          } else if (ARGIS("--ell-size-factor")) {
              if (i+1 < argc) bis_stroke = atoi(argv[++i]);
              else std::cerr << "--ell-size-factor requires an argument\n";
          } else if (ARGIS("--no-intersecting")) {
              no_intersecting = true;
          } else if (ARGIS("--grayscale")) {
              grayscale = true;
          } else if (ARGIS("--skip")) {
              if (i+1 < argc) skip = atoi(argv[++i]);
              else std::cerr << "--skip requires an argument\n";
          } else if (ARGIS("--limit")) {
              if (i+1 < argc) limit = atoi(argv[++i]);
              else std::cerr << "--limit requires an argument\n";
          } else if (ARGIS("--benchmark")) {
              benchmark = true;
          } else if (ARGIS("--precision")) {
              if (i+1 < argc) CGAL::Gmpfr::set_default_precision(atoi(argv[++i]));
              else std::cerr << "--skip requires an argument\n";
          } else if (ARGIS("--help")) {
              std::cout << "options: \n"
                           "  -i <inputfile>          read input data\n"
                           "  --no-validity-check     do not check graph validity after insertion\n"
                           "  --no-graph              do not compute graph by default (insert only ellipses)\n"
                           "  --ellipse-stroke <f>    ellipse stroke size (pen)\n"
                           "  --bisector-stroke <f>   bisector stroke size\n"
                           "  --ell-size-factor <n>   factor for random ellipse size\n"
                           "  --no-intersecting       do not allow pseudo-circles\n"
                           "  --grayscale             grayscale colors\n"
                           "  --skip <n>              skip first <n> ellipses from input file\n"
                           "  --limit <n>             read only <n> ellipses from input file\n"
                           " --benchmark             create only the graph and exit\n";

          } else {
              std::cerr << "unknown option " << argv[i] << std::endl;
          }
      }
//      std::cerr << "valid = " << validity_check << " test = " << test_only << " estroke = "
//              << ell_stroke << " bstroke = " << bis_stroke << std::endl;
#ifdef CGAL_HAS_THREADS
      std::cerr << "CGAL is MT\n";
#endif
#ifdef _OPENMP
      std::cerr << "OpenMP enabled\n";
#endif
      CGAL::set_pretty_mode (std::cerr);
      CGAL::set_pretty_mode (std::cout);
  }

  EGI *make_graphics(const Ellipse_2& e);

  void newScene() {
      if (!scene) {
          scene = new QGraphicsScene();
          scene->setItemIndexMethod(QGraphicsScene::NoIndex);
          graphicsView->setScene(scene);
          if (pi) delete pi;
          pi = new CGAL::Qt::GraphicsViewEllipseInput<ET>(this, scene);
          QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
                       this, SLOT(processInput(CGAL::Object)));
          vor = 0; del = 0; bgpix = 0;
          ellg = new QGraphicsItemGroup();
          scene->addItem(ellg);
      } else {
          remove_vor_from_scene();
          remove_del_from_scene();
          remove_ells_from_scene();
          remove_bg_from_scene();
      }
      ells.clear();
      graph.clear();
      scene->setSceneRect(-100, -100, 100, 100);
  }

//  void destroyScene(QGraphicsScene *scene) {
//      QList<QGraphicsItem *>::const_iterator it;
//      for (it = scene->items().constBegin(); it != scene->items().constEnd();
//           it++) {
//          scene->remove(*it);
//          delete(*it);
//      }
//  }

  void insertNewEllipse(const Ellipse_2 &e);
  void recompute_graph();
  void update_graph();
  void openFile(QString fname);
  void saveFile(QString fname);
  void openBackground(QString fname);

  void remove_vor_from_scene();
  void remove_del_from_scene();
  void remove_ells_from_scene();
  void remove_bg_from_scene();

public:
  MainWindow(int argc_, char **argv_);

public slots:

   bool processInput(CGAL::Object o);

   void on_actionRecenter_triggered();
   void on_actionI_points_triggered();
   void on_actionIds_triggered();
   void on_actionBisectors_triggered();
   void on_actionDelaunay_edges_triggered();
   void on_actionCompute_graph_triggered();
   void on_actionExport_triggered();
   void on_actionRandom_ellipse_triggered();

   void on_actionNew_triggered();
   void on_action_Open_triggered();
   void on_action_Save_triggered();
   void on_actionBackground_triggered();

   void on_actionInsert_ellipse_toggled(bool checked);

//signals:
//  void changed();
};


EGI *MainWindow::make_graphics(const Ellipse_2& e)
{
    EGI *gi;
    gi = new EGI();
    if (grayscale)
        gi->setPen(QPen(Qt::red, ell_stroke, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    else
        gi->setPen(QPen(Qt::black, ell_stroke, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    gi->setEllipse(e);
    gi->setIdVisible(actionIds->isChecked());
    gi->setIptVisible(actionI_points->isChecked());
    gi->show();
    return gi;
}

MainWindow::MainWindow(int argc_, char **argv_)
  : DemosMainWindow()
{
  setupUi(this);

//  QObject::connect(this, SIGNAL(changed()),
//		   gi, SLOT(modelChanged()));
  
  // Setup input handlers. They get events before the scene gets them

  // 
  // Manual handling of actions
  //
  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   this, SLOT(close()));

 
  argc = argc_; argv = argv_;
  parse_command_line();
  //
  // Setup the scene and the view
  //

  scene = 0; pi = 0;
  newScene();

  // Uncomment the following line to get antialiasing by default.
//   actionUse_Antialiasing->setChecked(true);

  // Turn the vertical axis upside down
  this->graphicsView->scale(1, -1);
  this->graphicsView->setMouseTracking(true);
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo("about_delaunay_graph_of_ellipses_2.html");
  this->addAboutCGAL();
//  actionRecenter->trigger();
  actionInsert_ellipse->toggle();
  if (!fname.isEmpty()) openFile(fname);

//  connect(this->actionBackground, SIGNAL(triggered()), this, SLOT(on_actionBackground_triggered()));
}

bool
MainWindow::processInput(CGAL::Object o)
{
  std::pair<ET::Kernel_scr::Point_2, ET::Kernel_scr::FT > center_sqr;
  Ellipse_2 e;
  if(CGAL::assign(center_sqr, o)){
    std::cerr << center_sqr.first << ',' << center_sqr.second << std::endl;
  } else if (CGAL::assign(e, o)) {
        int i;
        for (i = 0; i < ells.size(); i++) {
            if (e == ells[i]) {
                std::cerr << "ERROR: new ellipse is identical to " << ells[i].get_id() << std::endl;
                return false;
            }
            Bitangent bt = Bitangent(e, ells[i]);
//          std::cerr << "checking relpos with " << ells[i].get_id() << ": " << bt.relative_position() << std::endl;
            if (no_intersecting) {
                int p = bt.relative_position();
                if (p != CGAL::VORELL::SEPARATED) {
                  std::cerr << "ERROR: new ellipse is intersecting with ellipse "
                            << ells[i].get_id() << std::endl;
                  return false;
                }
            } else {
                if (bt.relative_position() == CGAL::VORELL::NOT_PC) {
                  std::cerr << "ERROR: new ellipse is not forming pseudo-circles with ellipse "
                            << ells[i].get_id() << std::endl;
                  return false;
                }
            }
        }
        if (i == ells.size()) {
          insertNewEllipse(e);
          update_graph();
        }
    }
  return true;
  //emit(changed());
}


// *  Qt Automatic Connections
// *  http://doc.trolltech.com/4.4/designer-using-a-component.html#automatic-connections
// *
// *  setupUi(this) generates connections to the slots named
// *  "on_<action_name>_<signal_name>"




void
MainWindow::on_actionRecenter_triggered()
{
  QRectF b;
  QList<QGraphicsItem *>::const_iterator it;
  for (it = scene->items().constBegin(); it != scene->items().constEnd(); it++) {
//      std::cerr << "object type = " << (*it)->type() << std::endl;
      if ((*it)->type() == EGI::Type) {
          EGI *e = static_cast<EGI *>(*it);
          b = b | e->boundingRect();
      }
  }
  if (bgpix && bgpix->isVisible()) b = b | bgpix->boundingRect();
  if (b.isValid()) {
      this->graphicsView->setSceneRect(b);
      this->graphicsView->fitInView(b, Qt::KeepAspectRatio);
  }
}

void MainWindow::on_actionI_points_triggered()
{
    QRectF b;
    QList<QGraphicsItem *>::const_iterator it;
    for (it = scene->items().constBegin(); it != scene->items().constEnd(); it++) {
        if ((*it)->type() == EGI::Type) {
            EGI *e = static_cast<EGI *>(*it);
            e->setIptVisible(actionI_points->isChecked());
            b = b | e->boundingRect();
        }
    }
    if (b.isValid()) {
        scene->update(b);
    }
}

void MainWindow::on_actionIds_triggered()
{
    QRectF b;
    QList<QGraphicsItem *>::const_iterator it;
    for (it = scene->items().constBegin(); it != scene->items().constEnd(); it++) {
        if ((*it)->type() == EGI::Type) {
            EGI *e = static_cast<EGI *>(*it);
            e->setIdVisible(actionIds->isChecked());
            b = b | e->boundingRect();
        }
    }
    if (b.isValid()) {
        scene->update(b);
    }
}

void MainWindow::on_actionBisectors_triggered()
{
    if (vor) vor->setVisible(actionBisectors->isChecked());
}

void MainWindow::on_actionDelaunay_edges_triggered()
{
    if (del) del->setVisible(actionDelaunay_edges->isChecked());
}

void MainWindow::on_actionCompute_graph_triggered()
{
    recompute_graph();
    update_graph();
}

void MainWindow::on_actionExport_triggered()
{
    QPrinter printer;
    QString fname = QFileDialog::getSaveFileName(this, "Export current view");
    if (fname.isEmpty()) return;
    printer.setOutputFileName(fname);
    printer.setFullPage(true);
    QRectF r = graphicsView->sceneRect();
    //std::cerr << "rw = " << r.width() << " rh = " << r.height() << std::endl;
    printer.setPaperSize(QSizeF(r.width(), r.height()),QPrinter::Point);
    printer.setPageMargins(0, 0, 0, 0, QPrinter::Point);
    QPainter painter;
    if (painter.begin(&printer)) {
        //painter.setTransform(graphicsView->viewportTransform());
        graphicsView->render(&painter);
        //scene->render(&painter, printer.pageRect(QPrinter::Point), r);
    } else {
        QMessageBox::critical(this, "Export failure", "Could not export file");
    }
}

void MainWindow::on_actionRandom_ellipse_triggered()
{
    Ellipse_2 e;
    do {
        QQ a,b,w,xc,yc;
        double sw = graphicsView->sceneRect().width();
        double sh = graphicsView->sceneRect().height();
        double left = graphicsView->sceneRect().left();
        double top = graphicsView->sceneRect().top();
        int gres = 30;
        double xs = sw/gres, ys = sh/gres;
        a = CGAL::Qt::rationalize<QQ>(rng.get_double()*xs*ell_size_factor+0.1,0.01);
        b = CGAL::Qt::rationalize<QQ>(rng.get_double()*ys*0.75*ell_size_factor+0.1,0.01);
        if (a < b) { w = a; a = b; b = w; }
        w = CGAL::Qt::rationalize<QQ>(rng.get_double()*2-1,0.01);
        xc = CGAL::Qt::rationalize<QQ>(double(rng.get_int(1,gres))*xs+left,0.01);
        yc = CGAL::Qt::rationalize<QQ>(double(rng.get_int(1,gres))*ys+top,0.01);
        e = Ellipse_2(a, b, w, xc, yc);
        //std::cerr << e << std::endl;
    //    std::cerr << graphicsView->sceneRect().left()  << ',' << graphicsView->sceneRect().top()
    //            << ',' << graphicsView->sceneRect().right() << ',' << graphicsView->sceneRect().bottom() << std::endl;
    } while (!processInput(CGAL::make_object(e)));
    std::cerr << "total ellipses = " << ells.size() << std::endl;
}

void MainWindow::recompute_graph()
{
    if (!actionCompute_graph->isChecked()) return;
    graph.clear();
    for (int i = 0; i < ells.size(); i++) {
        gtimer.start();
        graph.insert(ells[i]);
	gtimer.stop();
        if (validity_check) {
            vtimer.start();
            std::cerr << "Checking validity: [ ";
            graph.is_valid(true, 1);
            std::cerr << " ]\n";
            vtimer.stop();
        }
    }
}

void MainWindow::insertNewEllipse(const Ellipse_2 &e)
{
    std::cout << e << std::endl;
    ells.push_back(e);
    ellg->addToGroup(make_graphics(e));

    if (actionCompute_graph->isChecked()) {
        gtimer.start();
        graph.insert(e);
	gtimer.stop();
	if (benchmark) {
            std::cerr << "ellipse " << e.get_id() << " @time " << gtimer.time() << std::endl;
	}
        if (validity_check) {
            vtimer.start();
            std::cerr << "Checking validity: [ ";
            graph.is_valid(true, 1);
            std::cerr << " ]\n";
            vtimer.stop();
        }
    }
}

void MainWindow::remove_vor_from_scene()
{
    QList<QGraphicsItem *>::const_iterator it;
    if (vor) {
        for (it = vor->childItems().constBegin(); it != vor->childItems().constEnd(); it++) {
           scene->removeItem(*it);
           delete *it;
        }
        scene->destroyItemGroup(vor);
        vor = 0;
    }
}

void MainWindow::remove_del_from_scene()
{
    QList<QGraphicsItem *>::const_iterator it;
    if (del) {
        for (it = del->childItems().constBegin(); it != del->childItems().constEnd(); it++) {
            scene->removeItem(*it);
            delete *it;
        }
        scene->destroyItemGroup(del);
        del = 0;
    }
}

void MainWindow::remove_ells_from_scene()
{
    QList<QGraphicsItem *>::const_iterator it;
    for (it = ellg->childItems().constBegin(); it != ellg->childItems().constEnd(); it++) {
        scene->removeItem(*it);
        delete *it;
    }
}

void MainWindow::remove_bg_from_scene()
{
    if (bgpix) {
        scene->removeItem(bgpix);
        delete bgpix;
        bgpix = 0;
    }
}

void MainWindow::update_graph()
{
    remove_vor_from_scene();
    remove_del_from_scene();

    if (!actionCompute_graph->isChecked()) return;

    CGAL::Real_timer tm;
    tm.start();
    if (grayscale)
        vor = new CGAL::Qt::VoronoiEllipsesGraphicsItem<DG>(graph,
                QPen(Qt::gray, bis_stroke, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    else
        vor = new CGAL::Qt::VoronoiEllipsesGraphicsItem<DG>(graph,
                QPen(Qt::blue, bis_stroke, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    tm.stop();
    std::cerr << "Bisectors drawn in " << tm.time() << " seconds" << std::endl;
    on_actionBisectors_triggered();
    scene->addItem(vor);

    del = new CGAL::Qt::DelaunayEllipsesGraphicsItem<DG>(graph,
            QPen(Qt::green, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    on_actionDelaunay_edges_triggered();
    scene->addItem(del);

}

void MainWindow::openFile(QString fname)
{
    if (fname.isEmpty()) return;
    newScene();
    std::ifstream fin(fname.toUtf8().data());

    while (fin) {
        Ellipse_2 e;
        fin >> e;
        if (!fin) break;
        if (skip-- > 0) continue;
        if (limit-- == 0) break;
        std::cerr << "Inserting ellipse: " << e.get_id() << std::endl;
        insertNewEllipse(e);
    }
    fin.close();
    std::cerr << "Inserted " << gtimer.intervals() << " ellipses in " << gtimer.time() << " seconds, " <<
            gtimer.time()/gtimer.intervals() << " s/e" << std::endl;
    double tot = gtimer.time() + vtimer.time();
    std::cerr << "Validity check took " << vtimer.time() << " seconds, " << vtimer.time()*100.0/tot << "% of time\n";
    std::cerr << "Total time: " << gtimer.intervals() << " ellipses in " << tot << " seconds, " <<
            tot/gtimer.intervals() << " s/e" << std::endl;

#ifdef DELAUNAY_GRAPH_OF_ELLIPSES_PROFILED_TRAITS_2_H
    GT::profile_report();
#endif

    if (!benchmark) {
        update_graph();
//      TODO: unknown BUG here: segfaults/doesn't run at all when called instantly!
//        actionRecenter->trigger();
//      nasty workaround!
        QTimer::singleShot(50, this, SLOT(on_actionRecenter_triggered()));
    } else {
        QTimer::singleShot(50, this, SLOT(close()));
    }
}

void MainWindow::saveFile(QString fname)
{
    if (fname.isEmpty()) return;
    std::ofstream fout(fname.toUtf8().data());
    //int ce = 0;
    for (int i = 0; i < ells.size(); i++) {
        fout << ells[i] << std::endl;
    }
    fout.close();
    std::cerr << "Saved " << ells.size() << " ellipses\n";
}

void MainWindow::openBackground(QString fname)
{
    if (fname.isEmpty()) return;
    remove_bg_from_scene();
    bgpix = new QGraphicsPixmapItem(QPixmap(fname));
    bgpix->setZValue(-100);
    bgpix->scale(1,-1);
    scene->addItem(bgpix);
    //actionRecenter->trigger();
}

void MainWindow::on_actionNew_triggered()
{
    newScene();
}

void MainWindow::on_action_Open_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this, "Open File");
  openFile(fileName);
}

void MainWindow::on_action_Save_triggered()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Save File");
  saveFile(fileName);
}

void MainWindow::on_actionBackground_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this, "Background Image File");
  openBackground(fileName);
}

void
MainWindow::on_actionInsert_ellipse_toggled(bool checked)
{
  if(checked){
    scene->installEventFilter(pi);
  } else {
    scene->removeEventFilter(pi);
  }
}

#include "delaunay_graph_of_ellipses_2.moc"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);

//  mpfr_set_default_prec(128);

  app.setOrganizationDomain("di.uoa.gr");
  app.setOrganizationName("NKUA");
  app.setApplicationName("2D Delaunay graph of ellipses demo");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Delaunay_graph_of_ellipses_2);
  Q_INIT_RESOURCE(Input);
  Q_INIT_RESOURCE(CGAL);

  MainWindow *mainWindow = new MainWindow(argc, argv);
  if (!mainWindow->benchmark) mainWindow->show();
  // possible BUG: doesn't draw properly, some geometry change maybe
  // mainWindow->actionRecenter->trigger();
  return app.exec();
}
