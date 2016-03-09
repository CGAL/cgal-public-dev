/******************************************************************************
* Written by Ross Hemsley for INRIA.fr.
* A simple application to simulate different walks on Delaunay Triangulations.
******************************************************************************/

#include "simulate.h"
#include "walk.h"
#include "walk_3d.h"
#include "vertex_walk.h"    
#include <fstream>
#include <math.h>
#include <CGAL/Timer.h>
#include <CGAL/Real_timer.h>

using namespace std;

/*****************************************************************************/
// Display the current progress using a loadbar.

static inline void loadBar(int x, int n, int r, int w)
{
    // Only update r times.
    if ( x % (n/r) != 0 ) return;

    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;

    // Show the percentage complete.§
    printf("%3d%% [", (int)(ratio*100) );

    // Show the load bar.
    for (int x=0; x<c; x++)
       printf("=");

    for (int x=c; x<w; x++)
       printf(" ");

    // ANSI Control codes to go back to the
    // previous line and clear it.
    printf("]\n\033[F\033[J");
}

/*****************************************************************************/

// Run some statistics gathering on the visibility walk.
void runStats_timing()
{
    double side = 2000;
    int tests   = 10e5;

    std::list<Point_2>   points;
    Delaunay_2           dt;

    boost::mt19937 rng;    
    rng.seed(static_cast<unsigned int>(std::time(0)));
    boost::poisson_distribution<> gen(side*side);

    int lambda = gen(rng);

    std::cout << "Lambda: " << lambda << endl;
    CGAL::Random_points_in_square_2<Point_2,Creator_2> g(side/2);

    // For generating points well within the border.
    CGAL::Random_points_in_square_2<Point_2,Creator_2> g2(200/2);

    CGAL::copy_n(g, lambda, std::back_inserter(points));
    dt.insert(points.begin(), points.end());

    // The start and end points we are going to use.
    std::vector<Face_handle_2> start_face;
    std::vector<Point_2>       end_point;

    for (int i=0; i<tests; i++)
    {
        start_face.push_back(dt.locate( *(++g2) ));
        end_point.push_back(*(++g));
    }

    // The walking strategies.
    VisibilityWalk<Delaunay_2> w_vis(&dt);
    PivotWalk<Delaunay_2>      w_piv(&dt);
    Vertex_Walk<Delaunay_2>    w_vertex(&dt);

    CGAL::Real_timer t;

    int piv_orientations = 0;
    int piv_index        = 0;
    int vis_index        = 0;

    cout << "Doing Vertex " <<  endl;
    
    double average_x_progress        = 0;
    double average_angle             = 0;
    double average_radius            = 0;
    double average_distance_per_step = 0;
    double average_substeps          = 0;
    double average_path_length       = 0;
    double average_path_steps        = 0;
    double average_step_length       = 0;
    double average_neighbour_count   = 0;
    double num_steps                 = 0;
    
    int count=0;
    for (int i=0; i<tests; i++)
    {
        loadBar(i, tests, 100, 50);

        Face_handle_2 f = start_face[i];
        Point_2       p = end_point[i];

        w_vertex.do_walk(p,f, Vertex_Walk<Delaunay_2>::RESTRICTED);

        // If NaN.
        if (w_vertex.progress != w_vertex.progress)
            continue;

        count++;
        num_steps                 += w_vertex.step_count;
        average_distance_per_step += w_vertex.distance_per_step;
        average_x_progress        += w_vertex.progress;
        average_angle             += w_vertex.angle;
        average_radius            += w_vertex.radius;
        average_substeps          += w_vertex.substeps;
        average_step_length       += w_vertex.step_length;
        average_neighbour_count   += w_vertex.neighbour_count;
        average_path_steps        += w_vertex.path_steps;
        average_path_length       += w_vertex.path_length;
    }
    
    count = num_steps;
    cout << "dist. per step:      " << average_distance_per_step/count <<endl;
    cout << "x progress:          " << average_x_progress/count        <<endl;
    cout << "angle:               " << average_angle/count             <<endl;
    cout << "radius:              " << average_radius/count            <<endl;
    cout << "average substeps:    " << average_substeps/count          <<endl;
    cout << "average neighbour:   " << average_neighbour_count/count   <<endl;
    cout << "average length:      " << average_step_length/count       <<endl;
    cout << "Average path steps:  " << average_path_steps/count        <<endl;
    cout << "Average path length: " << average_path_length/count       <<endl;
    
    exit(0);

    cout << "Doing pivot" << endl;

    t.start();
    for (int i=0; i<tests; i++)
    {
        Face_handle_2 f = start_face[i];
        Point_2       p = end_point[i];

        w_piv.do_walk(p,f);

        piv_index+=w_piv.index_count;
        piv_orientations+= w_piv.getNumOrientationsPerformed();

    }
    t.stop();

    cout << "Piv: " << t.time() << endl;

    int vis_orientations=0;

    cout << "Doing vis" << endl;
    t.reset();
    t.start();
    for (int i=0; i<tests; i++)
    {
        Face_handle_2 f = start_face[i];
        Point_2       p = end_point[i];
        w_vis.do_walk(p,f);

        vis_index        += w_vis.index_count;
        vis_orientations += w_vis.getNumOrientationsPerformed();
    }
    t.stop();

    cout << "Vis: " << t.time() << endl;


    cout << "Vis did " << vis_orientations << " orientations" << endl;
    cout << "Piv did " << piv_orientations << " orientations" << endl;


    cout << "Vis did " << vis_index << " index tests" << endl;
    cout << "Piv did " << piv_index << " index tests" << endl;

}



/*****************************************************************************/

// Run some statistics gathering on the visibility walk.
void runStats_comparisons()
{

    int n       = 10e3;
    int tests   = 10e4;
    int buckets = 30;

    std::list<Point_2>   points;


    // Use a seperate triangulation for this.
    Delaunay_2 *dt  = new Delaunay_2();

    // Create 1000 random points in a square of side sqrt(n).
    CGAL::Random_points_in_square_2<Point_2,Creator_2> g(sqrt(n)/2);
    CGAL::copy_n( g, n, std::back_inserter(points) );
    dt->insert(points.begin(), points.end());

    // Variables to store statistics.
    float orientations_per_triangle_vis = 0;
    float orientations_per_triangle_piv = 0;
    float triangles_per_pivot           = 0;
    float improvement_ratio             = 0;
    int   orientation_difference        = 0;

    // Maximum possible distance between two points.
    float max_dist = sqrt(2)*sqrt(n);

    // Bucket the number of orientations for each distance.
    int piv_buckets[buckets];
    int vis_buckets[buckets];
    int bucket_size[buckets];

    // Zero the bucket arrays.
    for (int i=0; i<buckets; ++i)
    {
        piv_buckets[i] = 0;
        vis_buckets[i] = 0;
        bucket_size[i] = 0;
    }

    VisibilityWalk<Delaunay_2> w_vis(dt);
    PivotWalk<Delaunay_2>      w_piv(dt);


    for (int i=0; i<tests; ++i)
    {
        // Show current progress.
        loadBar(i, tests, 100, 50);

        // A new random point, which will be the destination for this walk.
        Point_2 p = *(++g);

        // The start point of this walk.
        Point_2 q = *(++g);

        // Locate the face containing the start point using CGAL's built-in
        // walk.
        Face_handle_2 f = dt->locate(q);

        Face_handle_2 f_vis = w_vis.do_walk(p,f);
        Face_handle_2 f_piv = w_piv.do_walk(p,f);

        Face_handle_2 f_actual = dt->locate(p);

        bool fail=false;
        if (f_vis != f_actual)
        {

            if (! (dt->is_infinite(f_vis) && dt->is_infinite(f_actual)))
            {
                fail=true;
                cout << "Visibility walk failed" << endl;
            }
        }
        if (f_piv != f_actual)
        {
            if (! (dt->is_infinite(f_piv) && dt->is_infinite(f_actual)) )
            {
                fail=true;
                cout << "Pivot walk failed" << endl;
            }
        }

        if (fail)
            exit(1);

        // Number of orientations for each walk.
        int o_vis = w_vis.getNumOrientationsPerformed();
        int o_piv = w_piv.getNumOrientationsPerformed();

        // Number of triangles visited for each walk.
        int t_vis = w_vis.getNumTrianglesVisited();
        int t_piv = w_piv.getNumTrianglesVisited();

        // Ratio of orientations per triangle for each walk.
        orientations_per_triangle_vis += o_vis / (float)t_vis;
        orientations_per_triangle_piv += o_piv / (float)t_piv;

        // Ratio of triangles per pivot for each walk.
        if (w_piv.getNumPivots() != 0)
            triangles_per_pivot  += t_piv /((float)w_piv.getNumPivots());

        // difference between number of orientations.
        orientation_difference += o_vis - o_piv;

        // Relative improvement of pivot relative to visibility
        improvement_ratio      += (o_vis - o_piv)/(float)o_vis;

        // Inter-point distance.
        float dist_pq = std::sqrt( (p - q).squared_length() );

        // The appropriate bucket for this walk.
        int this_bucket = floor( dist_pq/max_dist * buckets );

        // We'll need to average these later, keep track of the size of
        // each bucket.
        piv_buckets[this_bucket] += w_piv.getNumOrientationsPerformed();
        vis_buckets[this_bucket] += w_vis.getNumOrientationsPerformed();
        bucket_size[this_bucket] ++;
    }


    float count = tests;

    cout <<                                             endl
         << "Average orientations saved:  "
         << orientation_difference/count                 << endl

         << "Average saving:              "
         << improvement_ratio/count*100<< "%"        << endl << endl

         << "Orientations/triangle (vis): "
         << orientations_per_triangle_vis/count      << endl
         << "Orientations/triangle (piv): "
         << orientations_per_triangle_piv/count      << endl << endl

         << "Triangles per pivot:         "
         << triangles_per_pivot/count                << endl;


    // Write the bucket information to a file.
    ofstream f;
    f.open ("stats.txt");

    for (int i=0; i<buckets; ++i)
    {
        f << i << "  "
          << (i)/(float)buckets * max_dist          << "  "
          << piv_buckets[i] / (float)bucket_size[i] << "  "
          << vis_buckets[i] / (float)bucket_size[i] << "  "
          << endl;
    }

    f.close();
}

/*****************************************************************************/

void runStats_3d()
{

    int n       = 10e4;
    int tests   = 10e3;
    int buckets = 30;

    std::list<Point_3>   points;

    Delaunay_3 *dt = new Delaunay_3();

    // Insert n random points into the 3D Delaunay tessellation.
    CGAL::Random_points_in_sphere_3<Point_3> g;
    CGAL::copy_n( g, n, std::back_inserter(points) );
    dt->insert(points.begin(), points.end());

    // Do a 3D visibility walk.
    PivotWalk_3d<Delaunay_3>      w_piv(dt);
    VisibilityWalk_3d<Delaunay_3> w_vis(dt);

    // Variables to store statistics.
    float orientations_per_triangle_vis = 0;
    float orientations_per_triangle_piv = 0;
    float improvement_ratio             = 0;
    int   orientation_difference        = 0;

    // Maximum possible distance between two points.
    float max_dist = sqrt(2)*sqrt(n);

    // Bucket the number of orientations for each distance.
    int piv_buckets[buckets];
    int vis_buckets[buckets];
    int bucket_size[buckets];

    // Zero the bucket arrays.
    for (int i=0; i<buckets; ++i)
    {
        piv_buckets[i] = 0;
        vis_buckets[i] = 0;
        bucket_size[i] = 0;
    }

    float count=0;

    for (int i=0; i<tests; ++i)
    {
        // Show current progress.
        loadBar(i, tests, 100, 50);

        // A new random point, which will be the destination for this walk.
        Point_3 p = *(++g);

        // The start point of this walk.
        Point_3 q = *(++g);

        // Locate the face containing the start point using CGAL's built-in
        // walk.
        Cell_handle_3 f = dt->locate(q);

       if (w_vis.do_walk(p, f) != 1)
       {
           cout << "Failed test (vis). " << endl;
           continue;
       }

        if (w_piv.do_walk(p, f) != 1)
        {
            cout << "Failed test (piv.) " << endl;
            continue;
        }


        // Number of orientations for each walk.
        int o_vis = w_vis.getNumOrientationsPerformed();
        int o_piv = w_piv.getNumOrientationsPerformed();

        // Number of triangles visited for each walk.
        int t_vis = w_vis.getNumTrianglesVisited();
        int t_piv = w_piv.getNumTrianglesVisited();

        // Ratio of orientations per triangle for each walk.
        orientations_per_triangle_vis += o_vis / (float)t_vis;
        orientations_per_triangle_piv += o_piv / (float)t_piv;

        // difference between number of orientations.
        orientation_difference += o_vis - o_piv;

        // Relative improvement of pivot relative to visibility
        improvement_ratio      += (o_vis - o_piv)/(float)o_vis;

        // Inter-point distance.
        float dist_pq = std::sqrt( (p - q).squared_length() );

        // The appropriate bucket for this walk.
        int this_bucket = floor( dist_pq/max_dist * buckets );

        // We'll need to average these later, keep track of the size of
        // each bucket.
        piv_buckets[this_bucket] += w_piv.getNumOrientationsPerformed();
        vis_buckets[this_bucket] += w_vis.getNumOrientationsPerformed();
        bucket_size[this_bucket] ++;

        count ++;
    }

    cout <<                                             endl
         << "Average orientations saved:  "
         << orientation_difference/count                 << endl

         << "Average saving:              "
         << improvement_ratio/count*100<< "%"          << endl << endl

         << "Orientations/triangle (vis): "
         << orientations_per_triangle_vis/count      << endl
         << "Orientations/triangle (piv): "
         << orientations_per_triangle_piv/count      << endl << endl;

}

/*****************************************************************************/

int main(int argc, char **argv)
{
    // Do this every time to make sure that the walks are working correctly.
    //runStats_comparisons();

    // Do some timing.
    runStats_timing();

    return 0;
}

/*****************************************************************************/
