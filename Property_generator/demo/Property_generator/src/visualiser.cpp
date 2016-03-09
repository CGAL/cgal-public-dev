/******************************************************************************/

#include <QGLWidget>
#include <CGAL/Random.h>
#include <boost/random/poisson_distribution.hpp>

// Local.
#include "visualiser.h"
#include "ui_visualiser.h"
#include "point_generator.h"
#include "walk.h"
#include "triangulation_viewer.h"
#include "vertex_walk.h"


/*******************************************************************************
* Worker Members
*******************************************************************************/

// The fractal version of the build triangulation routine.
void Worker::build_triangulation_fractal( std::list<Point_2>* points,
                                          Delaunay*           dt,
                                          int                 start_size,
                                          int                 iterations  )
{
    int initial_size=start_size;

    std::cout << "Creating new fractal triangulation with "
              << iterations
               << " points."
              << endl;

    CGAL::Random_points_in_square_2<Point_2,Creator>   g(0.5);

    // The original size of the triangulation (consider what happens
    // when we set this to one?).
    CGAL::copy_n( g, initial_size, std::back_inserter(*points) );
    dt->insert(points->begin(), points->end());

    // Count the number of times the process fails.
    int fails=0;

    CGAL::Random random;

    // Now loop inserting random points int uniformly chosen triangles
    // in the triangulatino tht we have up until now:
    for (int i=0; i<iterations; i++)
    {
        // Choose a random triangle.
        int random_index = random.get_int (0, dt->number_of_faces());

        // TODO: do this more efficiently.
        Face_handle f;

        Finite_faces_iterator it;
        for (it=dt->finite_faces_begin(); it != dt->finite_faces_end(); ++it)
        {
            if (random_index == 0)
            {
                f = it;
                break;
            }
            random_index--;
        }

        // insert a point uniformly into that triangle.
        Point_2 cc = dt->circumcenter(f);

        // Choose any vertex of this face.
        const Point_2 & v = f->vertex(0)->point();

        // The circum radius is the distance from the circum
        // center to any vertex.
        double r2 = ( cc - v).squared_length();

        // Get circumradius
        double r = sqrt(r2);
        CGAL::Random_points_in_disc_2<Point_2, Creator> g2(r);

        // Give up after this number of attempts.
        int max_tries = 10e3;
        while (true)
        {
            // give up.
            if (max_tries == 0)
            {
                fails++;
                break;
            }
            max_tries--;

            Point_2 p = *(++g2);

            // Move point to this triangle.
            p = p + Vector_2(cc.x(), cc.y());

            // Is point inside triangle?
            bool inside = true;
            for (int i=0; i<3; i++)
            {
                const Point_2 & p0 = f->vertex(i)->point();
                const Point_2 & p1 = f->vertex(f->cw(i))->point();

                // If we have found a face that can see the point.
                if ( orientation(p0,p1,p) == CGAL::POSITIVE )
                {
                    inside = false;
                    break;
                }
            }

            if (inside)
            {
                // Insert into triangulation
                dt->insert(p);

                // Add to point list
                points->push_back(p);

                break;
            }
        }
    }

    std::cout << "Failed: " << fails << " times" << endl;
    emit finished();
}

/******************************************************************************/

void Worker::build_triangulation(std::list<Point_2>* points, Delaunay* dt, int size)
{
    std::cout << "Creating new triangulation with "<< size << " points." <<endl;

    CGAL::Random_points_in_square_2<Point_2,Creator> g(0.5);
    CGAL::copy_n( g, size, std::back_inserter(*points) );
    dt->insert(points->begin(), points->end());

    std::cout << "Done." << endl;
    emit finished();
}

/******************************************************************************/

void Worker::build_triangulation_poisson( std::list<Point_2>* points,
                                          Delaunay*           dt,
                                          double              rate,
                                          double              side    )
{
    std::cout << "Creating new poisson triangulation of rate "
         << rate << endl;

    // Generate a Poisson r.v using boost.    
    boost::mt19937 rng;    
    rng.seed(static_cast<unsigned int>(std::time(0)));
    boost::poisson_distribution<> gen(side*side*rate);

    int lambda = gen(rng);

    CGAL::Random_points_in_square_2<Point_2,Creator> g(side/2);
    CGAL::copy_n( g, lambda, std::back_inserter(*points) );
    dt->insert(points->begin(), points->end());

    std::cout << "Done." << endl;

    emit finished();
}

/******************************************************************************/

void Worker::build_triangulation( std::list<Point_2>* points,
                                  Delaunay*           dt,
                                  QString             filename )
{

    std::cout << "Creating new triangulation from file." << endl;

    load(filename, points);
    dt->insert(points->begin(), points->end());
    std::cout << "Done." << endl;

    emit finished();
}

/******************************************************************************/
// The following functions borrowed from Pierre Alliez.

void Worker::load(const QString& filename, std::list<Point_2> *points)
{
    if (filename.contains (".xyz", Qt::CaseInsensitive))
    {
        load_xyz_file (filename, points);
    }
        else if (filename.contains(".xy", Qt::CaseInsensitive))
    {
        load_xy_file(filename, points);
        return;
    }
}

/******************************************************************************/

void Worker::load_xy_file(const QString& fileName,std::list<Point_2> *points)
{
    std::ifstream ifs(qPrintable(fileName));
    std::cerr << "read xy...";
    Point_2 point;
    while (ifs >> point)
        points->push_back(point);
    std::cerr << "done (" << points->size() << " samples)" << std::endl;
    ifs.close();
}

/******************************************************************************/

void Worker::load_xyz_file(const QString& fileName,std::list<Point_2> *points)
{
    std::ifstream ifs(qPrintable(fileName));
    std::cerr << "read xyz...";

    std::vector<Point_3> temp;
    Point_3 point;
    while (ifs >> point)
        temp.push_back (point);

    Plane p;
    linear_least_squares_fitting_3
    (temp.begin (), temp.end (), p, CGAL::Dimension_tag<0>());

    for (uint i = 0; i < temp.size (); i ++)
        points->push_back(p.to_2d (temp[i]));

    ifs.close();
}

/******************************************************************************
* Visualiser Members
******************************************************************************/

void Visualiser::connect_signals_and_slots()
{
    // Connect up the UI using signals and slots.
    connect( ui->btn_walk_histogram,  SIGNAL( released()                     ),
             this,                    SLOT  ( draw_stats()                   ));
    connect( ui->combo_object_type,   SIGNAL( currentIndexChanged(QString)   ),
             this,                    SLOT  ( select_stats_object(QString)   ));
    connect( ui->radio_draw_hist,     SIGNAL( toggled(bool)                  ),
             this,                    SLOT  ( select_plot_type_walk()        ));
    connect( ui->actionNew_Pointset,  SIGNAL( triggered()                    ),
             this,                    SLOT  ( show_new_point_dialog()        ));
    connect( ui->actionNew_Walk,      SIGNAL( triggered()                    ),
             this,                    SLOT  ( new_walk()                     ));
    connect( ui->btn_choose_points,   SIGNAL( released()                     ),
             this,                    SLOT  ( new_walk()                     ));
    connect( ui->btn_new_seed,        SIGNAL( released()                     ),
             this,                    SLOT  ( reseed_generator()             ));
    connect( ui->btn_clear,           SIGNAL( released()                     ),
             this,                    SLOT  ( clear_walks()                  ));
    connect( ui->cb_pivot,            SIGNAL( stateChanged(int)              ),
             this,                    SLOT  ( set_pivot_walk_state(int)      ));
    connect( ui->cb_visibility,       SIGNAL( stateChanged(int)              ),
             this,                    SLOT  ( set_visibility_walk_state(int) ));
    connect( ui->cb_straight,         SIGNAL( stateChanged(int)              ),
             this,                    SLOT  ( set_straight_walk_state(int)   ));
    connect( ui->actionNew_Pointset,  SIGNAL( triggered()                    ),
             this,                    SLOT  ( show_new_point_dialog()        ));
    connect( ui->combo_selection_object_type,
                                      SIGNAL( currentIndexChanged (int)      ),
             this,                    SLOT  ( selection_update_object_type() ));
    connect( ui->btn_select,          SIGNAL( released()                     ),
             this,                    SLOT  ( selection_select()             ));
    connect( ui->btn_select_previous, SIGNAL( released()                    ),
             this,                    SLOT  ( selection_select_previous()    ));
    connect( ui->btn_select_next,     SIGNAL( released()                     ),
             this,                    SLOT  ( selection_select_next()        ));
    connect( ui->btn_select_last,     SIGNAL( released()                     ),
             this,                    SLOT  ( selection_select_last()        ));
    connect( ui->btn_select_first,    SIGNAL( released()                     ),
             this,                    SLOT  ( selection_select_first()       ));
    connect( ui->btn_clear_selection, SIGNAL( released()                     ),
             this,                    SLOT  ( selection_clear()              ));
    connect( ui->spin_density,        SIGNAL( valueChanged(int)              ),
             this,                    SLOT  ( spin_density_change()          ));
    connect( scene,                   SIGNAL( changed(QList<QRectF>)         ),
             this,                    SLOT  ( update_display_details()       ));
}

/******************************************************************************/
// Find closest vertex to a vertex in a triangulation.

Vertex_handle Visualiser::closest_vertex(Delaunay *dt, Vertex_handle v)
{
    double dist = std::numeric_limits<double>::infinity();
    double l;

    Vertex_handle nn;

    Vertex_circulator i = dt->incident_vertices(v), done(i);
    if (i!=0)
    {
        do
        {
            if (!dt->is_infinite(i))
            {
                l = (v->point() - i->point()).squared_length();
                if ( l < dist )
                {
                    dist = l;
                    nn   = i;
                }
            }
        } while ( ++i != done);
    }
    return nn;
}


/******************************************************************************/
// Create a dialog for loading pointsets.

void Visualiser::on_actionLoad_Pointset_triggered()
{
    QString filename = QFileDialog::getOpenFileName(this,
          tr("Load Pointset"), "",
          tr("2D Pointsets (*.xy);;3D Pointsets (*.xyz);;All Files (*)"));

    // Tell the worker thread what we want using the signals and slots
    // mechanism.
    build_triangulation(filename);

}

/******************************************************************************/
// View currently selected item.

void Visualiser::view_selection()
{
    // Clear previous selection.
    triangulation_viewer->clear_highlighted_faces();
    triangulation_viewer->clear_highlighted_edges();
    triangulation_viewer->clear_highlighted_vertices();

    // Center and feature size of the selection.
    Point_2 center;
    double d;

    if (selected_face_index != -1)
    {
        Face_handle f = selected_faces->at(selected_face_index);

        // Find shortest edge;
        Point_2  p1 = f->vertex(0)->point();
        Point_2  p2 = f->vertex(1)->point();
        Point_2  p3 = f->vertex(2)->point();

        double d1 = (p1-p2).squared_length();
        double d2 = (p2-p3).squared_length();
        double d3 = (p3-p1).squared_length();

        d = (d1 < d2) ? d1 : d2;
        d = (d3 < d ) ? d3 : d;
        d = sqrt(d);

        // Barycenter of triangle.
        center = CGAL::ORIGIN + (((p1-CGAL::ORIGIN) +
                                  (p2-CGAL::ORIGIN) +
                                  (p3-CGAL::ORIGIN)   ) /3);

        triangulation_viewer->add_highlighted_face(f);
    }

    if (selected_edge_index != -1)
    {
        Edge e = selected_edges->at(selected_edge_index);

        Segment s  = dt->segment(e);
        Point_2   p1 = s.point(0);
        Point_2   p2 = s.point(1);

        d      = sqrt( (p1-p2).squared_length() );
        center = CGAL::ORIGIN + (( (p1-CGAL::ORIGIN) +
                                   (p2-CGAL::ORIGIN)   ) /2);

        triangulation_viewer->add_highlighted_edge(s);
    }

    if (selected_vertex_index != -1)
    {
        Vertex_handle v, v1;

        v  = selected_vertices->at(selected_vertex_index);
        v1 = closest_vertex(dt,v);

        // Local feature size.
        d      = sqrt((v->point()-v1->point()).squared_length());
        center = v->point();

        Face_circulator fit = dt->incident_faces(v), done(fit);
        if (fit!=0)
        {
            do
            {
                if (!dt->is_infinite(fit))
                    triangulation_viewer->add_highlighted_face(fit);
            } while ( ++fit != done);
        }
    }

    double x = center.x() - d;
    double y = center.y() - d;
    double w = 2*d;
    double h = 2*d;

    if (ui->chk_focus_on_current->isChecked())
    {
        // Set the visible region.
        ui->graphicsView->fitInView(QRectF(x,y,w,h) , Qt::KeepAspectRatio );
        triangulation_viewer->update_view_rectangle(ui->graphicsView, true);
    }

    triangulation_viewer->update();
}

/******************************************************************************/

void Visualiser::initialise_selection_tools()
{

    selected_face_index   = -1;
    selected_edge_index   = -1;
    selected_vertex_index = -1;

    selected_faces        = NULL;
    selected_vertices     = NULL;
    selected_edges        = NULL;

    // Hide the 'currently selected' view from the tab bar.
    ui->tab_bar->removeTab(4);
    selection_update_object_type();
}

/******************************************************************************/
// This event is triggered when the user selects some items using the
// select tab in the side pane.

void Visualiser::selection_select()
{
    if (DISABLE_INTERFACE)
        return;

    // Get information from the ui.
    std::string stat = ui->combo_select_property->currentText().toStdString();

    // Remember the current statistic string so we can display it.
    selected_statistic_text = stat;

    selected_statistic_sorting_text = ui->combo_select_ordering->currentText();

    // Should we use all points in the triangulation, or should we
    // limit ourselves to a region?
    bool use_filter = ui->radio_selection_selected_region->isChecked();

    // Now, create the selection list.
    bool sorting = TRUE;
    if (ui->combo_select_ordering->currentText().compare("None") == 0)
        sorting = FALSE;

    bool ascending = TRUE;
    if (ui->combo_select_ordering->currentText().compare("Ascending") == 0)
        ascending = FALSE;

    // Figure out what kind of object the user is selecting, and
    // decide what to do depending on the result
    if (ui->combo_selection_object_type->currentIndex() == 1)
    {
        if (use_filter)
            selected_vertices = triangulation_viewer->get_selected_vertices();
        else
        {
            selected_vertices = new std::vector<Vertex_handle>();
            Finite_vertices_iterator i;
            for (i=dt->finite_vertices_begin();i!=dt->finite_vertices_end();++i)
                selected_vertices->push_back(i);
        }
        selected_vertex_index = 0;

        if (sorting)
        {


            // sort depending on statistic.
            Vertex_stats::Property s;
            s = Vertex_stats::string_to_property(stat);
            sort(selected_vertices->begin(), selected_vertices->end(),
                  CGAL::Property_sorter<Vertex_stats>(s, *dt) );
        }
    }

    if (ui->combo_selection_object_type->currentIndex() == 2)
    {
        if (use_filter)
            selected_edges = triangulation_viewer->get_selected_edges();
        else
        {
            selected_edges = new std::vector<Edge>();
            Finite_edges_iterator i;
            for (i=dt->finite_edges_begin(); i!=dt->finite_edges_end(); ++i)
                selected_edges->push_back(*i);
        }
        selected_edge_index = 0;

        if (sorting)
        {
            // sort depending on statistic.
            Edge_stats::Property s;
            s = Edge_stats::string_to_property(stat);
            sort(selected_edges->begin(), selected_edges->end(),
                  CGAL::Property_sorter<Edge_stats>(s,*dt) );
        }
    }

    if (ui->combo_selection_object_type->currentIndex() == 0)    {
        if (use_filter)
            selected_faces = triangulation_viewer->get_selected_faces();
        else
        {
            selected_faces = new std::vector<Face_handle>();
            Finite_faces_iterator i;
            for (i=dt->finite_faces_begin(); i!=dt->finite_faces_end(); ++i)
                selected_faces->push_back(i);
        }

        selected_face_index = 0;

        if (sorting)
        {
            // sort depending on statistic.
            Face_stats::Property s;
            s = Face_stats::string_to_property(stat);
            sort(selected_faces->begin(), selected_faces->end(),
                  CGAL::Property_sorter<Face_stats>(s,*dt) );
        }
    }


    // Show the 'current selection' tab and populate it.
    ui->tab_bar->removeTab(3);
    ui->tab_bar->insertTab(3, ui->tab_selection, "Selection");
    ui->tab_bar->setCurrentIndex(3);

    selection_update_details_text();
    view_selection();
}

/******************************************************************************/

void Visualiser::selection_select_previous()
{
    if (selected_face_index != -1 && selected_face_index != 0)
        selected_face_index--;

    if (selected_vertex_index != -1 && selected_vertex_index != 0)
        selected_vertex_index--;

    if (selected_edge_index != -1 && selected_edge_index != 0)
        selected_edge_index--;

    selection_update_details_text();
    view_selection();
}

/******************************************************************************/

void Visualiser::selection_select_next()
{
    if (selected_face_index != -1
        && selected_face_index != selected_faces->size()-1)
        selected_face_index++;

    if (selected_edge_index != -1
        && selected_edge_index != selected_edges->size()-1)
        selected_edge_index++;

    if (selected_vertex_index != -1
        && selected_vertex_index != selected_vertices->size()-1)
        selected_vertex_index++;

    selection_update_details_text();
    view_selection();
}

/******************************************************************************/

void Visualiser::selection_select_first()
{
    if (selected_face_index != -1)
        selected_face_index = 0;

    if (selected_edge_index != -1)
        selected_edge_index = 0;

    if (selected_vertex_index != -1)
        selected_vertex_index = 0;

    selection_update_details_text();
    view_selection();
}

/******************************************************************************/

void Visualiser::selection_select_last()
{
    if (selected_face_index != -1)
        selected_face_index = selected_faces->size()-1;

    if (selected_edge_index != -1)
        selected_edge_index = selected_edges->size()-1;

    if (selected_vertex_index != -1)
        selected_vertex_index = selected_vertices->size()-1;

    selection_update_details_text();
    view_selection();
}

/******************************************************************************/

void Visualiser::selection_clear()
{
    int page = ui->tab_bar->currentIndex();

    triangulation_viewer->clear_highlighted_faces();
    triangulation_viewer->clear_highlighted_edges();
    triangulation_viewer->clear_highlighted_vertices();

    // CLEAR OLD DATA
    if (selected_edge_index!=-1)
        delete selected_edges;
    if (selected_face_index!=-1)
        delete selected_faces;
    if (selected_vertex_index!=-1)
        delete selected_vertices;

    selected_edge_index   = -1;
    selected_face_index   = -1;
    selected_vertex_index = -1;

    ui->tab_bar->removeTab(3);
    ui->tab_bar->insertTab(3, ui->tab_select, "Selection");
    ui->tab_bar->setCurrentIndex(page);

    triangulation_viewer->update();
}


/******************************************************************************/

void Visualiser::spin_density_change()
{
    triangulation_viewer->set_density(ui->spin_density->value());
}

/******************************************************************************/

void Visualiser::on_chk_enable_sampling_stateChanged( int state )
{
    bool enabled =  (state == Qt::Checked);

    ui->spin_density->setEnabled(enabled);

    if (!enabled)
        triangulation_viewer->set_density(0);
    else
        triangulation_viewer->set_density(ui->spin_density->value());

}

/******************************************************************************/
// Populate the display details label with information about the
// currently visualised triangulation.

void Visualiser::update_display_details()
{
    QString details;

    details += "<strong>Points loaded: </strong>";
    details +=QString::number(triangulation_viewer->number_of_points_in_view());
    details += ",  <strong>Vertices: </strong>";
    details += QString::number(dt->number_of_vertices ());
    details += ",  <strong>Faces: </strong>";
    details += QString::number(dt->number_of_faces ());

    label_display_details->setText(details);

    if (triangulation_viewer->sampling_enabled())
    {
        ui->chk_enable_sampling->setCheckState(Qt::Checked);
        ui->spin_density->setEnabled(true);
    }
    else
    {
        ui->chk_enable_sampling->setCheckState(Qt::Unchecked);
        ui->spin_density->setEnabled(false);
    }

    if (triangulation_viewer->sampling_enabled())
        ui->spin_density->setValue(triangulation_viewer->get_density());
}

/******************************************************************************/

void Visualiser::selection_update_details_text()
{
    int         current;
    int         max; 
    std::string object_type;

    if (selected_face_index != -1)
    {
        object_type = "Triangulation Faces";
        current = selected_face_index;
        max     = selected_faces->size();
    }

    if (selected_vertex_index != -1)
    {
        object_type = "Trianguation Vertices";
        current = selected_vertex_index;
        max     = selected_vertices->size();
    }

    if (selected_edge_index != -1)
    {
        object_type = "Triangulation Edges";
        current = selected_edge_index;
        max     = selected_edges->size();
    }

    // all the statistics for this object.
    std::vector< std::pair<std::string,double> > *all_stats;

    if (selected_face_index != -1)
    {
        Face_handle f = selected_faces->at(selected_face_index);
        all_stats = Face_stats::get_all_properties(*dt, f);
    }

    if (selected_vertex_index != -1)
    {
        Vertex_handle v=selected_vertices->at(selected_vertex_index);
        all_stats = Vertex_stats::get_all_properties(*dt, v);
    }

    if (selected_edge_index != -1)
    {
        Edge e = selected_edges->at(selected_edge_index);
        all_stats = Edge_stats::get_all_properties(*dt, e);
    }

    ui->table_selection_properties->setRowCount(4+all_stats->size());

    QString sort_text = "None";
    if (selected_statistic_sorting_text.compare("None") != 0)
    {
        sort_text = selected_statistic_sorting_text;
        sort_text += QString(" (");
        sort_text += QString::fromStdString(selected_statistic_text);
        sort_text += QString(")");
    }

    QTableWidgetItem* tw_object_name =
        new QTableWidgetItem(QString::fromStdString(object_type));
    QTableWidgetItem* tw_sorting     =
        new QTableWidgetItem(sort_text);
    QTableWidgetItem* tw_current_val =
        new QTableWidgetItem(QString::number(current+1));
    QTableWidgetItem* tw_max_val     =
        new QTableWidgetItem(QString::number(max));

    ui->table_selection_properties->setItem(0,0,
        new QTableWidgetItem("Object") );
    ui->table_selection_properties->setItem(1,0,
        new QTableWidgetItem("Sorting") );
    ui->table_selection_properties->setItem(2,0,
        new QTableWidgetItem("Selected Item") );
    ui->table_selection_properties->setItem(3,0,
        new QTableWidgetItem("Selection Size") );

    ui->table_selection_properties->setItem(0,1, tw_object_name);
    ui->table_selection_properties->setItem(1,1, tw_sorting);
    ui->table_selection_properties->setItem(2,1, tw_current_val);
    ui->table_selection_properties->setItem(3,1, tw_max_val);

    int j=4;
    // Display all the statistics for the selected item.
    std::vector< std::pair<std::string,double> >::iterator i;
    for (i=all_stats->begin(); i!=all_stats->end(); ++i, ++j)
    {
        QTableWidgetItem* property;
        QTableWidgetItem* value;

        property = new QTableWidgetItem(QString::fromStdString((*i).first));
        value    = new QTableWidgetItem(QString::number((*i).second));

        ui->table_selection_properties->setItem(j,0, property);
        ui->table_selection_properties->setItem(j,1, value);
    }
}

/******************************************************************************/

void Visualiser::selection_update_object_type()
{
    // Use the statistics class to give us a list of possible things
    // to select.

    std::vector<std::string>*    properties;

    if (ui->combo_selection_object_type->currentIndex() == 1)
    {
        // Add the properties that we can do to points.
        properties = Vertex_stats::get_capabilities();
    }

    if (ui->combo_selection_object_type->currentIndex() == 2)
    {
        // Add the properties we can find about edges.
        properties = Edge_stats::get_capabilities();
    }

    if (ui->combo_selection_object_type->currentIndex() == 0)
    {
        // Add the properties that we can find about triangles.
        properties = Face_stats::get_capabilities();
    }

    ui->combo_select_property->clear();
    std::vector<std::string>::iterator i;
    for (i=properties->begin(); i!=properties->end(); ++i)
    {
        ui->combo_select_property->addItem(QString::fromStdString(*i));
    }
}

/******************************************************************************/

void Visualiser::initialise_walk_tools()
{
    straight_walk      = new StraightWalk<Delaunay>(dt);
    visibility_walk    = new VisibilityWalk<Delaunay>(dt);
    pivot_walk         = new PivotWalk<Delaunay>(dt);

    // Init. bools to determine whether or not we draw given walks.
    drawStraightWalk   = FALSE;
    drawVisibilityWalk = FALSE;
    drawPivotWalk      = FALSE;

    // Point location status.The default state is to not take new input points.
    point_location_status = -1;
}

/******************************************************************************/

void Visualiser::initialise_statistics_tools()
{
    select_stats_object(ui->combo_object_type->currentText());
}

/******************************************************************************/

Visualiser::Visualiser(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::Visualiser)
{
    // Allow interaction with the interface.
    // This is to give thread safety when constructing triangulations.
    DISABLE_INTERFACE=false;

    // "Hacky" but effective way to generate a pseudorandom random seed
    // on user request.
    srand(time(NULL));

    QMessageBox msgBox;
    msgBox.setText("WARNING! This is a pre-pre-alpha build and not for     \
                    redistribution. Please also note the following:<br><br>\
                    - Vertex Walk is not yet finished at all,              \
                      and so is very buggy.                                \
                    - To select a region, hold shift whilst dragging a box \
                      around points you want to select (This allows you    \
                      to restrict statistics and 'selections' (better name \
                      pending) to a subset of the points).<br>             \
                    - Note that you can view very large triangulations     \
                      (tested up to 20,000,000). But the plotting          \
                      and selection tools are not yet ready to deal with   \
                      that many points. <br>                               \
                    - Straight walk is just the built in line-walk function\
                      from CGAL so has no statistics etc. It's just there  \
                      for comparison.<br>                                  \
                    - For updates, bug reports, and feature requests,      \
                      e-mail rlhemsley@gmail.com.<br><br>                  \
                    Thanks, Ross.                                          ");
    msgBox.exec();

    // Initialise user interface and connections.
    ui->setupUi(this);
    init_worker_thread();

    // Build dialog which we use to create new pointsets on start up.
    point_generator             = new Point_Generator(this);

    // Connect the point generator to the worker thread.
    connect( point_generator, SIGNAL( build_triangulation(int) ),
             this,            SLOT  ( build_triangulation(int) ));

    connect( point_generator, SIGNAL( build_triangulation_fractal(int, int) ),
             this,            SLOT  ( build_triangulation_fractal(int, int) ));

    connect( point_generator, SIGNAL( build_triangulation(QString) ),
             this,            SLOT  ( build_triangulation(QString) ));

    connect( point_generator, SIGNAL(
                                build_triangulation_poisson(double, double) ),
             this,            SLOT  (
                                build_triangulation_poisson(double, double)  ));

    // Initialise the main graphics items and data structures.
    dt                          = new Delaunay();
    points                      = new std::list<Point_2>();
    scene                       = new QGraphicsScene();

    label_display_details       = new QLabel();
    ui->statusbar->addWidget(label_display_details, 0.2);

    // Initialise the progress bar.
    progress_bar                = new QProgressBar();
    ui->statusbar->addWidget(progress_bar,1);
    ui->statusbar->setContentsMargins(10,0,0,0);
    progress_bar->hide();

    // Do other initialisations.
    connect_signals_and_slots();
    initialise_selection_tools();
    initialise_walk_tools();
    initialise_statistics_tools();

    // Connect up the layers of the UI.
    ui->graphicsView->setScene(scene);
    ui->graphicsView->setViewport(new QGLWidget(QGLFormat(QGL::SampleBuffers)));
    ui->graphicsView->setRenderHint(QPainter::Antialiasing);

    // Initliase the triangulation viewer.
    triangulation_viewer = new Triangulation_Viewer<Delaunay>(dt);
    ui->graphicsView->installEventFilter(triangulation_viewer);
    ui->graphicsView->viewport()->installEventFilter(triangulation_viewer);

    // Event filters for dealing with mouse input.
    // These are attached to both the GrahpicsView and GraphicsScene.
    ui->graphicsView->installEventFilter(this);
    ui->graphicsView->viewport()->installEventFilter(this);
    ui->graphicsView->viewport()->setMouseTracking(true);

    scene->addItem(triangulation_viewer);

    // Create a random triangulation on start up.
    build_triangulation(200);

    // TODO: disentangle selecting in the UI with selecting properties...
    // This is set when the user starts dragging a region.
    doing_selection=false;
    // This is the rectangle we use for selections.
    rectItem = new QGraphicsRectItem;
    QColor rect_color(250, 221, 0);
    rect_color.setAlpha(60);
    rectItem->setBrush(rect_color);
    rect_color.setAlpha(255);
    rectItem->setPen(Qt::NoPen);
    rectItem->hide();
    rectItem->setZValue(10000);
    scene->addItem(rectItem);
}

/******************************************************************************/
// This thread deals with triangulation generation.

void Visualiser::init_worker_thread()
{
    // Create a worker thread to build the triangulation.
    thread = new QThread;
    worker = new Worker();
    worker->moveToThread(thread);

    // connect(worker, SIGNAL(finished()), thread, SLOT(quit()));

    connect(worker, SIGNAL(finished()), this, SLOT(triangulation_finished() ));


    connect(this,   SIGNAL( do_load_triangulation(          std::list<Point_2>*,
                                                            Delaunay*,
                                                            QString         )),
            worker, SLOT  ( build_triangulation(            std::list<Point_2>*,
                                                            Delaunay*,
                                                            QString         )));


    connect(this,   SIGNAL( do_build_triangulation(         std::list<Point_2>*,
                                                            Delaunay*, int  )),
            worker, SLOT  ( build_triangulation(            std::list<Point_2>*,
                                                            Delaunay*, int  )));


    connect(this,   SIGNAL( do_build_triangulation_fractal( std::list<Point_2>*,
                                                            Delaunay*,
                                                            int,
                                                            int             )),
            worker, SLOT  ( build_triangulation_fractal(    std::list<Point_2>*,
                                                            Delaunay*,
                                                            int,
                                                            int             )));


    connect(this,   SIGNAL( do_build_triangulation_poisson( std::list<Point_2>*,
                                                            Delaunay*,
                                                            double,
                                                            double          )),
            worker, SLOT  ( build_triangulation_poisson(    std::list<Point_2>*,
                                                            Delaunay*,
                                                            double,
                                                            double          )));

    thread->start();
}


/******************************************************************************/

void Visualiser::select_plot_type_walk()
{
    ui->combo_walk_properties_2->setEnabled(!ui->radio_draw_hist->isChecked());
    ui->combo_walk_properties_3->setEnabled(!ui->radio_draw_hist->isChecked());
}


/******************************************************************************/
// Take logs of all the values in this vector.

void take_logs(std::vector<double>* v)
{
    std::vector<double>::iterator i;
    for (i=v->begin(); i!=v->end(); ++i)
        *i = std::log(*i);        
}

/******************************************************************************/
// TODO: Clean up this mess of a function.

void Visualiser::draw_stats()
{
    if (DISABLE_INTERFACE)
        return;

    // Should we limit the statistics to the subset defined by
    // the user's current selection?
    bool use_filter = ui->radio_plotting_selected_region->isChecked();

    if (use_filter)
        std::cout << "USING FILTER"<<endl;

    std::string stat_1 
        = ui->combo_walk_properties_1->currentText().toStdString();
    std::string stat_2 
        = ui->combo_walk_properties_2->currentText().toStdString();
    std::string stat_3 
        = ui->combo_walk_properties_3->currentText().toStdString();

    std::vector<double> *values_1;
    std::vector<double> *values_2;
    std::vector<double> *values_3;

    QString type = ui->combo_object_type->currentText();

    if(ui->radio_draw_hist->isChecked())
    {
        if (type.compare("Walk (faces)") == 0)
        {
            draw_walk_histogram();
            return;
        }
        else if (type.compare("Triangulation (faces)") == 0)
        {
            Face_stats::Property p = Face_stats::string_to_property(stat_1);
            if (use_filter)
            {
                std::vector<Face_handle>* faces;
                faces = triangulation_viewer->get_selected_faces();
                values_1 = Face_stats::compute_property(*dt, p, faces);
            } else {
                values_1 = Face_stats::compute_property(*dt, p);
            }
        }
        else if (type.compare("Triangulation (edges)") == 0)
        {
            Edge_stats::Property p = Edge_stats::string_to_property(stat_1);            
            if (use_filter)
            {
                std::vector<Edge>* edges;
                edges = triangulation_viewer->get_selected_edges();
                values_1 = Edge_stats::compute_property(*dt, p, edges);
            } else {
                values_1 = Edge_stats::compute_property(*dt, p);
            }
        }
        else if (type.compare("Triangulation (vertices)") == 0)
        {
            Vertex_stats::Property p = Vertex_stats::string_to_property(stat_1);            
            if (use_filter)
            {
                std::vector<Vertex_handle>* vertices;
                vertices = triangulation_viewer->get_selected_vertices();
                values_1 = Vertex_stats::compute_property_log(*dt, p, vertices);
            } else {
                values_1 = Vertex_stats::compute_property_log(*dt, p);
            }
        }
        else
        {
            std::cerr << "Unknown object type" << endl;
            exit(1);
        }
        draw_histogram(values_1, stat_1);
    }
    else if (ui->combo_walk_properties_3->currentText().compare("Fixed") == 0)
    {
        if (type.compare("Walk (faces)") == 0)
        {
            draw_walk_scatter_plot();
            return;
        }
        else if (type.compare("Triangulation (faces)") == 0)
        {
            Face_stats::Property p1 = Face_stats::string_to_property(stat_1);
            Face_stats::Property p2 = Face_stats::string_to_property(stat_2);            
            if (use_filter)
            {
                std::vector<Face_handle>* faces;
                faces = triangulation_viewer->get_selected_faces();
                values_1 = Face_stats::compute_property(*dt, p1, faces);
                values_2 = Face_stats::compute_property(*dt, p2, faces);
            } else {
                values_1 = Face_stats::compute_property(*dt, p1);
                values_2 = Face_stats::compute_property(*dt, p2);
            }
        }
        else if (type.compare("Triangulation (edges)") == 0)
        {
            Edge_stats::Property p1 = Edge_stats::string_to_property(stat_1);
            Edge_stats::Property p2 = Edge_stats::string_to_property(stat_2);                  
            if (use_filter)
            {
                std::vector<Edge>* edges;
                edges = triangulation_viewer->get_selected_edges();
                values_1 = Edge_stats::compute_property(*dt, p1, edges);
                values_2 = Edge_stats::compute_property(*dt, p2, edges);
            } else {
                values_1 = Edge_stats::compute_property(*dt, p1);
                values_2 = Edge_stats::compute_property(*dt, p2);
            }

        }
        else if (type.compare("Triangulation (vertices)") == 0)
        {
            Vertex_stats::Property p1 = Vertex_stats::string_to_property(stat_1);
            Vertex_stats::Property p2 = Vertex_stats::string_to_property(stat_2);                 
            if (use_filter)
            {
                std::vector<Vertex_handle>* vertices;
                vertices = triangulation_viewer->get_selected_vertices();
                values_1 = Vertex_stats::compute_property(*dt, p1, vertices);
                values_2 = Vertex_stats::compute_property(*dt, p2, vertices);
            } else {
                values_1 = Vertex_stats::compute_property(*dt, p1);
                values_2 = Vertex_stats::compute_property(*dt, p2);
            }
        }
        else
        {
            std::cerr << "Unknown object type" << endl;
            exit(1);
        }

        if (ui->chk_plot_log_x->checkState() == Qt::Checked)
        {
            stat_1 = stat_1 + " (log)";
            take_logs(values_1);            
        }
            
        if (ui->chk_plot_log_y->checkState() == Qt::Checked)
        {
            stat_2 = stat_2 + " (log)";
            take_logs(values_2);            
        }
            
        draw_scatter_plot(values_1, values_2, stat_1, stat_2);

    } else {
        // Draw a colored plot.
        if (type.compare("Triangulation (faces)") == 0)
        {
            Face_stats::Property p1 = Face_stats::string_to_property(stat_1);
            Face_stats::Property p2 = Face_stats::string_to_property(stat_2);                  
            Face_stats::Property p3 = Face_stats::string_to_property(stat_3);                  
            if (use_filter)
            {
                std::vector<Face_handle>* faces;
                faces = triangulation_viewer->get_selected_faces();
                values_1 = Face_stats::compute_property(*dt, p1, faces);
                values_2 = Face_stats::compute_property(*dt, p2, faces);
                values_3 = Face_stats::compute_property(*dt, p3, faces);
            } else {
                values_1 = Face_stats::compute_property_log(*dt, p1);
                values_2 = Face_stats::compute_property_log(*dt, p2);
                values_3 = Face_stats::compute_property_log(*dt, p3);
            }
        }
        else if (type.compare("Triangulation (edges)") == 0)
        {
            Edge_stats::Property p1 = Edge_stats::string_to_property(stat_1);
            Edge_stats::Property p2 = Edge_stats::string_to_property(stat_2);                  
            Edge_stats::Property p3 = Edge_stats::string_to_property(stat_3);      
            if (use_filter)
            {
                std::vector<Edge>* edges;
                edges = triangulation_viewer->get_selected_edges();
                values_1 = Edge_stats::compute_property_log(*dt, p1, edges);
                values_2 = Edge_stats::compute_property_log(*dt, p2, edges);
                values_3 = Edge_stats::compute_property_log(*dt, p3, edges);
            } else {
                values_1 = Edge_stats::compute_property_log(*dt, p1);
                values_2 = Edge_stats::compute_property_log(*dt, p2);
                values_3 = Edge_stats::compute_property_log(*dt, p3);
            }
        }
        else if (type.compare("Triangulation (vertices)") == 0)
        {
            Vertex_stats::Property p1 = Vertex_stats::string_to_property(stat_1);
            Vertex_stats::Property p2 = Vertex_stats::string_to_property(stat_2);                  
            Vertex_stats::Property p3 = Vertex_stats::string_to_property(stat_3);                          
            if (use_filter)
            {
                std::vector<Vertex_handle>* vertices;
                vertices = triangulation_viewer->get_selected_vertices();
                values_1 = Vertex_stats::compute_property_log(*dt, p1, vertices);
                values_2 = Vertex_stats::compute_property_log(*dt, p2, vertices);
                values_3 = Vertex_stats::compute_property_log(*dt, p3, vertices);
            } else {
                values_1 = Vertex_stats::compute_property_log(*dt, p1);
                values_2 = Vertex_stats::compute_property_log(*dt, p2);
                values_3 = Vertex_stats::compute_property_log(*dt, p3);
            }
        }
        else
        {
            // TODO: Add 3D walk plots.
            QMessageBox msgBox;
            msgBox.setText("Not yet added coloured walk plot \
                           (set 'Ticks' to 'Fixed' and try again).");
            msgBox.exec();
            return;
        }

        if (ui->chk_plot_log_x->checkState() == Qt::Checked)
        {
            stat_1 = stat_1 + " (log)";
            take_logs(values_1);
        }
        
        if (ui->chk_plot_log_y->checkState() == Qt::Checked)
        {
            stat_2 = stat_2 + " (log)";
            take_logs(values_2);            
        }

        if (ui->chk_plot_log_z->checkState() == Qt::Checked)
        {
            stat_3 = stat_3 + " (log)";
            take_logs(values_3);
        }

        draw_scatter_plot_color( values_1, values_2, values_3,
                                 stat_1,   stat_2,   stat_3    );
    }
}

/******************************************************************************/

void Visualiser::draw_scatter_plot_color( std::vector<double>* v1,
                                          std::vector<double>* v2,
                                          std::vector<double>* v3,
                                          std::string          s1,
                                          std::string          s2,
                                          std::string          s3  )
{
    Statistics_plotter* plot = new Statistics_plotter(this);

    if (v1->size()==0)
    {
        QMessageBox::information( this, "Delaunay Visualiser",
                                  "No data available to plot.\n" );
        return;
    }

    Scatter_Plot_Color* scatter_plot_color;

    scatter_plot_color = new Scatter_Plot_Color("Title");
    scatter_plot_color->setValues(v1, v2,v3);
    scatter_plot_color->attach(plot);


    plot->add_color_bar( QString::fromStdString(s3),
                         scatter_plot_color->getInterval() );

    plot->resize(600,400);
    plot->set_titles(QString::fromStdString(s1), QString::fromStdString(s2));
    plot->show();
}

/******************************************************************************/

void Visualiser::draw_scatter_plot( std::vector<double>* v1,
                                    std::vector<double>* v2,
                                    std::string          s1,
                                    std::string          s2  )
{
    Statistics_plotter* plot = new Statistics_plotter(this);

    if (v1->size()==0)
    {
        QMessageBox::information( this, "Delaunay Visualiser",
                                  "No data available to plot.\n" );
        return;
    }

    Scatter_Plot* scatter_plot;

    scatter_plot = new Scatter_Plot("Title", Qt::red);
    scatter_plot->setValues(v1, v2);
    scatter_plot->attach(plot);

    plot->resize(600,400);
    plot->set_titles(QString::fromStdString(s1), QString::fromStdString(s2));
    plot->show();

}

/******************************************************************************/

void Visualiser::draw_histogram(std::vector<double>* v1, std::string s1)
{
    Statistics_plotter* plot = new Statistics_plotter(this);

    if (v1->size()==0)
    {
        QMessageBox::information( this, "Delaunay Visualiser",
                                  "No data available to plot.\n" );
        return;
    }

    Histogram* hist;

    hist = new Histogram("Title", Qt::red);
    hist->setValues(v1);
    hist->attach(plot);

    plot->resize(600,400);
    plot->set_titles(QString::fromStdString(s1), QString("Frequency"));

    // plot->fix_legend();
    plot->show();
}

/******************************************************************************/
// Walks are a special case because we have multiple data sets
// for one object.

void Visualiser::draw_walk_scatter_plot()
{
    // Local variables defined here to save horizontal space.
    std::vector<double>*  visibility_v1;
    std::vector<double>*  visibility_v2;
    std::vector<double>*  pivot_v1;
    std::vector<double>*  pivot_v2;
    std::vector<double>*  straight_v1;
    std::vector<double>*  straight_v2;
    Scatter_Plot*         scatter_visibility;
    Scatter_Plot*         scatter_pivot;
    Scatter_Plot*         scatter_straight;

    // The plot window we will plot to.
    Statistics_plotter* plot = new Statistics_plotter(this);

    // Find the enum value of the statistics we are interested in.
    Walk_stats::Property s1 = Walk_stats::string_to_property(
        ui->combo_walk_properties_1->currentText().toStdString());
    Walk_stats::Property s2 = Walk_stats::string_to_property(
        ui->combo_walk_properties_2->currentText().toStdString());

    // Get value series for each plot.
    visibility_v1 = Walk_stats::compute_property(*visibility_walk, s1);
    visibility_v2 = Walk_stats::compute_property(*visibility_walk, s2);
    pivot_v1      = Walk_stats::compute_property(*pivot_walk,      s1);
    pivot_v2      = Walk_stats::compute_property(*pivot_walk,      s2);
    straight_v1   = Walk_stats::compute_property(*straight_walk,   s1);
    straight_v2   = Walk_stats::compute_property(*straight_walk,   s2);

    if ( visibility_v1 ->size() == 0 &&
         pivot_v1      ->size() == 0 &&
         straight_v1   ->size() == 0      )
    {
        QMessageBox::information( this, "Delaunay Visualiser",
                                  "No walks available to plot.\n" );

        delete plot;
        delete visibility_v1;
        delete visibility_v2;
        delete pivot_v1;
        delete pivot_v2;
        delete straight_v1;
        delete straight_v2;

        return;
    }

    if (visibility_v1->size() != 0)
    {
        scatter_visibility = new Scatter_Plot("Visibility Walk", Qt::blue);
        scatter_visibility->setValues(visibility_v1, visibility_v2);
        scatter_visibility->attach(plot);
    }

    if (pivot_v1->size() != 0)
    {
        scatter_pivot = new Scatter_Plot("Pivot Walk", Qt::red);
        scatter_pivot->setValues(pivot_v1, pivot_v2);
        scatter_pivot->attach(plot);
    }

    if (straight_v1->size() != 0)
    {
        scatter_straight = new Scatter_Plot("Straight Walk", Qt::green);
        scatter_straight->setValues(straight_v1, straight_v2);
        scatter_straight->attach(plot);
    }

    plot->resize(600,400);
    plot->set_titles(ui->combo_walk_properties_1->currentText(), 
                     ui->combo_walk_properties_2->currentText());
    plot->add_legend();
    plot->show();

}

/******************************************************************************/
// TODO: Ensure we free values allocated when histogram is destroyed.

void Visualiser::draw_walk_histogram()
{
    Statistics_plotter* plot = new Statistics_plotter(this);
    plot->resize(600,400);

    QString _s = ui->combo_walk_properties_1->currentText();
    Walk_stats::Property stat;
    stat = Walk_stats::string_to_property(_s.toStdString());

    std::vector<double>* visibility_values;
    std::vector<double>* straight_values;
    std::vector<double>* pivot_values;

    visibility_values = Walk_stats::compute_property(*visibility_walk, stat);
    pivot_values      = Walk_stats::compute_property(*pivot_walk, stat);
    straight_values   = Walk_stats::compute_property(*straight_walk, stat);

    if ( visibility_values->size() == 0 &&
         pivot_values     ->size() == 0 &&
         straight_values  ->size() == 0      )
    {
        QMessageBox::information( this, "Delaunay Visualiser",
                                  "No walks available to plot.\n" );

        delete plot;
        delete visibility_values;
        delete straight_values;
        delete pivot_values;

        return;
    }

    if (visibility_values->size() != 0)
    {
        Histogram* hist_visibility = new Histogram("Visibility Walk",Qt::blue);
        hist_visibility->setValues( visibility_values );
        hist_visibility->attach(plot);
    }

    if (pivot_values->size() != 0)
    {
        Histogram* hist_pivot = new Histogram("Pivot Walk", Qt::red);
        hist_pivot->setValues(pivot_values);
        hist_pivot->attach(plot);
    }

    if (straight_values->size() != 0)
    {
        Histogram* hist_straight = new Histogram("Straight Walk", Qt::green);
        hist_straight->setValues(straight_values);
        hist_straight->attach(plot);
    }

    plot->add_legend();
    plot->set_titles(_s, QString("Frequency"));
    plot->show();
}

/******************************************************************************/

void Visualiser::select_stats_object(QString type)
{
    std::vector<std::string>* stats=NULL;

    if (type.compare("Walk (faces)") == 0)
        stats = Walk_stats::get_capabilities();
    else if (type.compare("Triangulation (faces)") == 0)
        stats = Face_stats::get_capabilities();
    else if (type.compare("Triangulation (edges)") == 0)
        stats = Edge_stats::get_capabilities();
    else if (type.compare("Triangulation (vertices)") == 0)
        stats = Vertex_stats::get_capabilities();
    else
    {
        std::cerr << "Unknown object type" << endl;
        exit(1);
    }

    ui->combo_walk_properties_1->clear();
    ui->combo_walk_properties_2->clear();
    ui->combo_walk_properties_3->clear();
    ui->combo_walk_properties_3->addItem("Fixed");

    std::vector<std::string>::iterator i;
    for (i=stats->begin(); i!=stats->end(); ++i)
    {
        ui->combo_walk_properties_1->addItem(QString::fromStdString(*i));
        ui->combo_walk_properties_2->addItem(QString::fromStdString(*i));
        ui->combo_walk_properties_3->addItem(QString::fromStdString(*i));
    }
}

/******************************************************************************/

void Visualiser::reseed_generator()
{
    rseed = rand();
    updateScene();
}

/******************************************************************************/
// Identical functions to flip drawing state of the different walks.

void Visualiser::set_straight_walk_state(int b)
{
    drawStraightWalk = b;
    updateScene();
}

/******************************************************************************/
void Visualiser::set_visibility_walk_state(int b)
{
    drawVisibilityWalk = b;
    updateScene();
}

/******************************************************************************/

void Visualiser::set_pivot_walk_state(int b)
{
    drawPivotWalk = b;
    updateScene();
}

/******************************************************************************/
// When the user wants to create a new pointset this is activated.

void Visualiser::show_new_point_dialog()
{
    point_generator->show();
    updateScene();
}

/******************************************************************************/

void Visualiser::clear_walks()
{
    point_location_status = -1;
    updateScene();
}

/******************************************************************************/
// Main event handler, attached to scene and graphics view.

bool Visualiser::eventFilter(QObject *obj, QEvent *event)
{
    // In order to correctly resize the scene, we need to wait until
    // the graphicsview is fully constructed.
    if( obj == ui->graphicsView && event->type() == QEvent::Show )
        resizeScene();

    QMouseEvent* mouseEvent= static_cast<QMouseEvent*>(event);

    // Ignore non-mouse events.
    if (mouseEvent == NULL) return false;

    // Mouse position in scene coordiantes.
    QPointF pos = ui->graphicsView->mapToScene(mouseEvent->pos());

    switch( event->type() )
    {
        case QEvent::MouseButtonPress:
        {
            if (mouseEvent->modifiers() & Qt::ShiftModifier)
            {
                // Start doing a selection.
                if (!doing_selection)
                {
                    doing_selection = true;
                    mouse_selection_point =
                        ui->graphicsView->mapToScene(mouseEvent->pos());
                    rectItem->show();
                    rectItem->setRect(
                        QRectF(mouse_selection_point, QSize(0,0)));
                }
                return true;
            }
            return false;
        }

        // ** MOUSE BUTTON RELEASE ************************************* //
        case QEvent::MouseButtonRelease:
        {
            if (doing_selection)
            {
                doing_selection=false;
                triangulation_viewer->setSelection(rectItem->rect());
                rectItem->hide();
            }

            if ( mouseEvent->button() == Qt::LeftButton )
            {
                // If we are selecting the first or second point.
                if (    point_location_status == 0
                     || point_location_status == 1  )
                {
                    // Set this point, and increment status counter.
                    point_location_points[point_location_status++] = pos;
                }

                // If we just finished inserting the points.
                if ( point_location_status == 2 )
                {
                    ui->graphicsView->setCursor(Qt::ArrowCursor);
                    updateScene();
                    return false;
                }
            }

            // Allow other objects to receive mouse release events.
            return false;
        }

        // ** MOUSE MOVED ********************************************** //
        case QEvent::MouseMove:
        {

            // If we're currently inserting the first or second point.
            if ( point_location_status == 0 || point_location_status == 1 )
            {
                // Temporarily set the first point to be the current
                // position.
                point_location_points[ point_location_status ] = pos;
                updateScene();
                return true;
            }

            // Draw a rectangle in the scene to represent selection.
            if(rectItem->isVisible()) {
                QPointF size = ui->graphicsView->mapToScene(mouseEvent->pos());
                size = size - mouse_selection_point;

                QPointF TL = mouse_selection_point;

                if (size.x() < 0)
                    TL.setX( TL.x() + size.x() );
                if (size.y() < 0)
                    TL.setY( TL.y() + size.y() );

                rectItem->setRect(
                    TL.x(),
                    TL.y(),
                    abs(size.x()),
                    abs(size.y()) );
                return false;
            }
            return false;
        }

        // ** DEFAULT ************************************************** //
        default:
            return false;
    }

    // pass the event on to the parent class
    return false;
}

/******************************************************************************/

void Visualiser::new_walk()
{

    if (DISABLE_INTERFACE)
        return;

    // We are going to start taking point inputs.
    // This says that we are currently learning point 1.
    point_location_status = 0;

    // Set the mouse to a crosshair when moving over the graphics view.
    ui->graphicsView->setCursor(Qt::CrossCursor);

    // Small hack to stop any faces being displayed before any mouse over
    // events occur.
    point_location_points[0] = QPointF(-1000,-1000);

    // Clear any old walk graphics from the scene.
    updateScene();
}

/******************************************************************************/
// Automatically refresh sizing of triangulation when window is resized.

void Visualiser::resizeEvent (QResizeEvent * event)
{
    triangulation_viewer->update_view_rectangle(ui->graphicsView, true);
}

/******************************************************************************/
// Call this when we are about to create a new triangulation.

void Visualiser::triangulation_starting()
{
    //qDebug() << "Preparing to build..." <<endl;
    progress_bar->setMaximum ( 0 );
    progress_bar->setMinimum ( 0 );
    progress_bar->reset ();
    progress_bar->show();

    // Don't interact with the triangulation.
    DISABLE_INTERFACE=true;

    // Clear anything that may have been selected
    selection_clear();

    // Do not interact with the triangulation in this thread
    triangulation_viewer->setValid(false);

    // Claer old triangulation in this thread.
    dt->clear();
    points->clear();

    // Clear the walks.
    while (! walk_items.isEmpty() )
        scene->removeItem(walk_items.takeFirst());

    visibility_walk->clear();
    pivot_walk->clear();
    straight_walk->clear();

    label_display_details->setText("");
}

/******************************************************************************/

void Visualiser::triangulation_finished()
{
    triangulation_viewer->setValid(true);
    triangulation_viewer->modelChanged();

    progress_bar->hide();

    resizeScene();

    // Clear old walk.
    point_location_status = -1;
    updateScene();
    update_display_details();

    triangulation_viewer->update_view_rectangle(ui->graphicsView, true);

    DISABLE_INTERFACE=false;
}

/******************************************************************************/
// Ensure the triangulation fits neatly into the viewer.

void Visualiser::resizeScene()
{
    // Bounding rectangle for the triangulation.
    // TODO: Remove magic numbers.
    ui->graphicsView->setSceneRect( QRectF(-1000,-1000,2000,2000) );
    ui->graphicsView->fitInView( triangulation_viewer->boundingRect(),
                                                         Qt::KeepAspectRatio );
}

/******************************************************************************/
// Re-draw the walks and things which may change on mouse events.

void Visualiser::updateScene()
{
    if (DISABLE_INTERFACE)
        return;

    // Style for points.
    QPen   pen(Qt::black);
    QBrush brush(Qt::blue);

    // This is where we store details about each walk to be displayed.
    QString details;

    // Remove all the previous walkItems, that is, the previous triangles
    // drawn as part of the walk.
    while (! walk_items.isEmpty() )
        scene->removeItem(walk_items.takeFirst());

    visibility_walk->clear();
    pivot_walk->clear();
    straight_walk->clear();

    // If we are drawing walks.
    if ( point_location_status >= 0 )
    {
        // Find the face we are hovering over.
        Face_handle f = dt->locate(c(point_location_points[0]));

        // Check the face is finite, and then draw it.
        if (! dt->is_infinite(f) )
        {
            // Draw a shaded in face where the user's mouse is.
            QGraphicsItem *tr;
            tr = Walk<Delaunay>::drawTriangle(f, QPen(),
                                                 QColor("#D2D2D2"));
            scene->addItem(tr);
            walk_items.append(tr);
        }
    }

    // If we have enough data to draw a walk, then do so.
    if (point_location_status > 0)
    {
        // Start and end points.
        Point_2 p = c( point_location_points[0] );
        Point_2 q = c( point_location_points[1] );

        // Start and end faces.
        Face_handle f = dt->locate(p);
        Face_handle g = dt->locate(q);

        if ( !dt->is_infinite(f) && !dt->is_infinite(g) )
        {
            if (drawStraightWalk)
            {
                straight_walk->do_walk(q,f);

                QGraphicsItem* walkGraphics;
                walkGraphics = straight_walk->getGraphics(QPen(),
                                                          QColor("#EBEBD2"));

                walk_items.append(walkGraphics);
                scene->addItem(walkGraphics);

                int o_count = straight_walk->getNumOrientationsPerformed();
                int t_count = straight_walk->getNumTrianglesVisited();

                // Details about this walk.
                details += "<b>Straight Walk</b><br>";
                details += "Orientations: ";
                details += QString::number(o_count);
                details += "<br>Triangles Visited: ";
                details += QString::number(t_count);
                details += "<br><br>";
            }

            if (drawVisibilityWalk)
            {
                visibility_walk->do_walk(q, f, rseed);

                QGraphicsItem* walkGraphics;
                walkGraphics = visibility_walk->getGraphics(QPen(),
                                                            QColor("#D2D2EB"));

                walk_items.push_back(walkGraphics);
                scene->addItem(walkGraphics);

                int o_count = visibility_walk->getNumOrientationsPerformed();
                int t_count = visibility_walk->getNumTrianglesVisited();

                // Details about this walk.
                details += "<b>Visibility Walk</b><br>";
                details += "Orientations: ";
                details += QString::number(o_count);
                details += "<br>Triangles Visited: ";
                details += QString::number(t_count);
                details += "<br><br>";
            }

            if (drawPivotWalk)
            {
                pivot_walk->do_walk(q, f, rseed);

                QGraphicsItem* walkGraphics;
                walkGraphics = pivot_walk->getGraphics(QPen(),
                                                       QColor("#EBD2D2"));

                walk_items.append(walkGraphics);
                scene->addItem(walkGraphics);

                // Details about this walk.
                details += "<b>Pivot Walk</b><br>";
                details += "Orientations: ";
                details += QString::number(
                               pivot_walk->getNumOrientationsPerformed());
                details += "<br>Triangles Visited: ";
                details += QString::number(
                               pivot_walk->getNumTrianglesVisited());
                details += "<br><br>";
            }

            // TODO: Quite a lot of re-factoring here.
            // It seems that we should implement all walks
            // as single classes with options rather than
            // have so many different inhereted classes
            // for different strategies.

            if (ui->chk_vertex_walk->isChecked())
            {
                Vertex_Walk<Delaunay> vertex_walk(dt);
                vertex_walk.do_walk(q,f, Vertex_Walk<Delaunay>::STANDARD);
                QGraphicsItem* walkGraphics = vertex_walk.getGraphics();
                walk_items.append(walkGraphics);
                scene->addItem(walkGraphics);
            }
            
            if (ui->chk_monotone_vertex_walk->isChecked())
            {
                Vertex_Walk<Delaunay> vertex_walk(dt);
                vertex_walk.do_walk(q,f, Vertex_Walk<Delaunay>::MONOTONE);
                QGraphicsItem* walkGraphics = vertex_walk.getGraphics();
                walk_items.append(walkGraphics);
                scene->addItem(walkGraphics);
            }                    

            if (ui->chk_restricted_vertex_walk->isChecked())
            {
                Vertex_Walk<Delaunay> vertex_walk(dt);
                vertex_walk.do_walk(q,f, Vertex_Walk<Delaunay>::RESTRICTED);
                QGraphicsItem* walkGraphics = vertex_walk.getGraphics();
                walk_items.append(walkGraphics);
                scene->addItem(walkGraphics);
            }
        }

        // If we've finished choosing points for the walk, then draw a point
        // to mark where the user selected the end point to be.
        if (point_location_status == 2)
        {
            QGraphicsEllipseItem *e;
            QPointF p = point_location_points[1];

            e = scene->addEllipse( QRectF( QPointF(-3.5,-3.5),
                                           QSizeF(7,7)), pen, brush );

            // Don't resize this point, but do move it.
            e->setFlag(QGraphicsItem::ItemIgnoresTransformations);
            e->setPos(p);

            walk_items.append(e);
        }
    }

    // update_view_rectangle(v, false);
    // triangulation_viewer->update(  );

    // Add the walk details to the UI.
    ui->label_walk_status->setText(details);
}

/******************************************************************************/

void Visualiser::build_triangulation(int size)
{
    if (DISABLE_INTERFACE)
        return;

    triangulation_starting();
    emit do_build_triangulation( points, dt, size );
}

/******************************************************************************/

void Visualiser::build_triangulation_fractal(int start_size, int iterations)
{
    if (DISABLE_INTERFACE)
        return;

    triangulation_starting();
    emit do_build_triangulation_fractal( points, dt, start_size, iterations);
}

/******************************************************************************/

void Visualiser::build_triangulation_poisson(double rate, double side)
{
    if (DISABLE_INTERFACE)
        return;

    triangulation_starting();
    emit do_build_triangulation_poisson( points, dt, rate, side);
}

/******************************************************************************/

void Visualiser::build_triangulation(QString filename)
{
    if (DISABLE_INTERFACE)
        return;

    triangulation_starting();
    emit do_load_triangulation(points, dt, filename );
}

/******************************************************************************/

Visualiser::~Visualiser()
{
    delete ui;
}

/******************************************************************************/
