// (c) 2007-2009 George Tzoumas <geortz@gmail.com>

//#define DGET_VERBOSE
//#define VERBOSE 1

#define CGAL_GMPFR_NO_REFCOUNT

#include <iostream>

#include <CGAL/basic.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfi.h>

#include <CGAL/Ellipse_traits.h>
//#include <CGAL/Apollonius_graph_2.h>
//#include <CGAL/Apollonius_ellipses_traits_2.h>
#include <CGAL/Delaunay_graph_of_ellipses_2.h>
#include <CGAL/Delaunay_graph_of_ellipses_traits_2.h>


typedef CGAL::Gmpq QQ;
typedef CGAL::Gmpfi INT;
typedef CGAL::Ellipse_traits<QQ, INT> ET;
typedef CGAL::Ellipse_2<ET> Ellipse_2;
typedef CGAL::Delaunay_graph_of_ellipses_traits_2<ET> GT;
//typedef CGAL::Apollonius_ellipses_traits_2<ET> GT;
typedef CGAL::Delaunay_graph_of_ellipses_2< GT > DG;

//using namespace std;

void compute_dg()
{
    DG graph;
    Ellipse_2 e;
    int ce = 0;
    
    while (std::cin) {
        std::cin >> e;
        if (!std::cin) break;
        ce++;
        std::cerr << "Inserting ellipse: " << ce << std::endl;
        graph.insert(e);
        std::cerr << "Checking validity: [ ";
        graph.is_valid(true, 1);
        std::cerr << " ]\n";
    }
}

int main(int argc, char **argv)
{
    CGAL::set_pretty_mode (std::cerr);
    compute_dg();

    return 0;
}

