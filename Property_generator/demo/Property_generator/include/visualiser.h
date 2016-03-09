/******************************************************************************/
#ifndef VISUALISER_H
#define VISUALISER_H
/******************************************************************************/

// Qt Components.
#include <QMainWindow>
#include <QtGui>
#include <QAbstractTableModel>

// CGAL Componenets.
#include <CGAL/Qt/Converter.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <CGAL/point_generators_2.h>

// Local Includes.
#include "statistics_plotter.h"
#include "walk.h"
#include "point_generator.h"
#include "triangulation_viewer.h"

// Property generating classes.
#include <CGAL/Property_generator.h>
#include <CGAL/property_generator_traits/triangulation_2_property_traits.h>
#include "walk_traits.h"

/******************************************************************************/
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
/******************************************************************************/

typedef Kernel::Point_2                              Point_2;
typedef Kernel::Vector_2		                     Vector_2;
typedef Kernel::Point_3     	                     Point_3;
typedef Kernel::Plane_3      					     Plane;

typedef CGAL::Triangulation_2<Kernel>                Triangulation_2;
typedef CGAL::Delaunay_triangulation_2<Kernel>       Delaunay;
typedef CGAL::Creator_uniform_2<double,Point_2>      Creator;

typedef Delaunay::Face_handle                        Face_handle;
typedef Delaunay::Vertex_handle                      Vertex_handle;
typedef Delaunay::Face                               Face;
typedef Delaunay::Edge                               Edge;
typedef Delaunay::Segment                            Segment;
typedef Delaunay::Face_circulator                    Face_circulator;
typedef Delaunay::Vertex_circulator                  Vertex_circulator;
typedef Delaunay::Line_face_circulator               Line_face_circulator;
typedef Delaunay::Finite_faces_iterator              Finite_faces_iterator;
typedef Delaunay::Finite_edges_iterator              Finite_edges_iterator;
typedef Delaunay::Finite_vertices_iterator           Finite_vertices_iterator;

// The property generator traits classes.
typedef CGAL::Triangulation_2_edge_property_traits<Delaunay>
                                                     Dt_edge_traits;
typedef CGAL::Triangulation_2_vertex_property_traits<Delaunay>
                                                     Dt_vertex_traits;
typedef CGAL::Triangulation_2_face_property_traits<Delaunay>
                                                     Dt_face_traits;
typedef CGAL::Walk_2_face_property_traits<Delaunay>  Walk_traits;

// The property generators.
typedef CGAL::Property_generator<Dt_face_traits>     Face_stats;
typedef CGAL::Property_generator<Dt_vertex_traits>   Vertex_stats;
typedef CGAL::Property_generator<Dt_edge_traits>     Edge_stats;
typedef CGAL::Property_generator<Walk_traits>        Walk_stats;

/*******************************************************************************
* Worker thread for doing long running tasks.
*******************************************************************************/

class Worker : public QObject {
    Q_OBJECT

public slots:
    void build_triangulation         ( std::list<Point_2>* points,
                                       Delaunay*    dt,
                                       int          size        );

    void build_triangulation_fractal ( std::list<Point_2>* points,
                                       Delaunay*    dt,
                                       int          start_size,
                                       int          iterations  );

    void build_triangulation         ( std::list<Point_2>* points,
                                       Delaunay*    dt,
                                       QString      filename    );

    void build_triangulation_poisson ( std::list<Point_2>* points,
                                       Delaunay*    dt,
                                       double       rate,
                                       double       side        );

private:
    void load          ( const QString& filename, std::list<Point_2>* points );
    void load_xy_file  ( const QString& fileName, std::list<Point_2>* points );
    void load_xyz_file ( const QString& fileName, std::list<Point_2>* points );

signals:
     void finished();
     void error(QString err);

};

/******************************************************************************
* Visualiser class
******************************************************************************/

namespace Ui {
class Visualiser;
}

class Visualiser : public QMainWindow
{
    Q_OBJECT

public:
    explicit Visualiser(QWidget *parent = 0);
    ~Visualiser();

protected:
    bool eventFilter(QObject *obj, QEvent *event);
    void resizeEvent (QResizeEvent * event);

private slots:

    // Functions written for Qt to auto connect.
    void on_chk_enable_sampling_stateChanged( int state );
    void on_actionLoad_Pointset_triggered();

    void spin_density_change();
    void update_display_details();

    // Slots relating to walk.
    // TODO: This in a less stupid way.
    void new_walk();
    void set_straight_walk_state(int b);
    void set_visibility_walk_state(int b);
    void set_pivot_walk_state(int b);
    void reseed_generator();
    void clear_walks();

    // Slots relating to scene.
    void updateScene();
    void resizeScene();

    // Stastics plotting stuff.
    void select_plot_type_walk();
    void select_stats_object(QString type);
    void draw_stats();
    void draw_walk_scatter_plot();
    void draw_walk_histogram();

    // Slots for involved in building new triangulations.
    void show_new_point_dialog();
    void build_triangulation(int size);
    void build_triangulation(QString filename);
    void build_triangulation_poisson(double rate, double side);
    void build_triangulation_fractal(int start_size, int iterations);
    void triangulation_finished();

    // Deal with selections.
    void selection_update_object_type();
    void selection_select();
    void selection_select_previous();
    void selection_select_next();
    void selection_select_first();
    void selection_select_last();
    void selection_clear();

signals:

    // Signalling mechanism to communicate with worker thread.
    // We use these to pass the object references to the worker thread
    // from the main thread.
    void do_load_triangulation          ( std::list<Point_2>* points,
                                          Delaunay*      dt,
                                          QString        filename     );

    void do_build_triangulation         ( std::list<Point_2>* points,
                                          Delaunay*      dt,
                                          int            size         );

    void do_build_triangulation_fractal ( std::list<Point_2>* points,
                                          Delaunay*      dt,
                                          int            start_size,
                                          int            iterations   );

    void do_build_triangulation_poisson ( std::list<Point_2>* points,
                                          Delaunay*      dt,
                                          double         rate,
                                          double         side         );

private:
    // Qt generated UI.
    Ui::Visualiser*                 ui;

    // Thread for dealing with triangulation constrution.
    QThread*                        thread;
    Worker*                         worker;

    // Dialog for creating new pointsets.
    Point_Generator*           		point_generator;

    // Current random seed.
    long                            rseed;

    // Disable interaction with the interface at certain points.
    bool 							DISABLE_INTERFACE;

    // Which walks to draw.
    bool                            drawPivotWalk;
    bool                            drawStraightWalk;
    bool                            drawVisibilityWalk;

    std::string                     selected_statistic_text;
    QString							selected_statistic_sorting_text;

    // When doing a selection, we set this as the first point.
    bool                            doing_selection;
    QPointF                         mouse_selection_point;
    QGraphicsRectItem*              rectItem;

    // If an item is selected, one of these will be
    // non negative, and relate to the vectors below.
    int                             selected_face_index;
    int                             selected_edge_index;
    int                             selected_vertex_index;

    std::vector<Face_handle>*       selected_faces;
    std::vector<Vertex_handle>*     selected_vertices;
    std::vector<Edge>* 		        selected_edges;


    // Point location status indicator and values.
    int                             point_location_status;
    QPointF                         point_location_points[2];

    // Label for current details.
    // Sometimes removed from UI so kept here.
    QLabel*                         label_display_details;

    // Progress bar in status bar.
    QProgressBar *progress_bar;

    // Objects.
    CGAL::Qt::Converter<Kernel>     c;
    Delaunay*                       dt;
    std::list<Point_2>*  			points;
    QGraphicsScene*                 scene;
    Triangulation_Viewer<Delaunay>* triangulation_viewer;
    QList<QGraphicsItem*>           walk_items;
    VisibilityWalk<Delaunay>*       visibility_walk;
    PivotWalk<Delaunay>*            pivot_walk;
    StraightWalk<Delaunay>*         straight_walk;

    // Functions //

    Vertex_handle closest_vertex(Delaunay *dt, Vertex_handle v);
    void init_worker_thread();
    void triangulation_starting();
    void view_selection();
    void selection_update_details_text();
    void initialise_walk_tools();
    void initialise_selection_tools();
    void initialise_statistics_tools();
    void connect_signals_and_slots();
    void draw_scatter_plot       ( std::vector<double>* v1,
                                   std::vector<double>* v2,
                                   std::string          s1,
                                   std::string          s2  );

    void draw_scatter_plot_color ( std::vector<double>* v1,
                                   std::vector<double>* v2,
                                   std::vector<double>* v3,
                                   std::string          s1,
                                   std::string          s2,
                                   std::string          s3  );

    void draw_histogram(           std::vector<double>* v1,
                                   std::string          s1  );


};

/*****************************************************************************/
#endif // VISUALISER_H
/*****************************************************************************/
