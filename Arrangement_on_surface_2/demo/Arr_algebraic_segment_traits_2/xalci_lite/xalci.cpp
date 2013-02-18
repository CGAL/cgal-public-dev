// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : AlciX
// File          : demos/xalci/xalci.cpp
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.43 $
// Revision_date : $Date: 2009-06-30 13:14:58 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-inf.mpg.de>
//                 Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#define NDEBUG 1

//#define CGAL_NO_LEDA

#include <CGAL/basic.h>
#include <CGAL/Timer.h>

#include <set>
#include <ctime>
#include <fstream>

#define CGAL_POLYNOMIAL_USE_NTL_MUL

#include "include/misc.h"
#include "include/xalci.h"

#define XALCI_USE_FLAT_COLOR_SCHEME 0

#ifdef CGAL_ACK_BENCHMARK_RES
namespace CGAL {
CGAL::Timer res_tm;
}
#endif

static const Kernel_2& kernel_2 = CKvA_2::instance().kernel();

bool check_testsuite(int argc, char** argv) {
    if (argc > 1) {
        for (int a = 0; a < argc; a++) {
            if (std::strcmp( argv[a], "--test-suite") == 0) {
                std::cerr << "This interactive program terminates "
                          << "immediately for the test-suite." << std::endl;
                return true;
            }
        }
    }
    return false;
}

#if !CGAL_HAS_QT3

int main (int argc, char** argv) {

    if(check_testsuite(argc, argv)) 
        return 0;
    
    std::cerr << "This demo requires Qt!" << std::endl;
    return 0;
}
#else // !CGAL_USE_QT

extern QColor rasterize_colors[];
extern int n_rast_colors;

extern Object_vector oc_objects;
extern Object_vector* curr_objects;
extern bool subdiv_layer_changed;
 
extern Graphic_layer *subdiv_layer;
extern Layers oc_layers;

extern CGAL::Timer timer;

typedef std::map<std::size_t, int> Color_map;

// picks up a color index based on id parameter
static int pick_color(std::size_t id, Color_map& color_map) {

    int cindex = color_map.size();
    std::pair<Color_map::iterator, bool> ret =
        color_map.insert(std::make_pair(id, cindex));
    if(!ret.second) // already exists in the map
        cindex = ret.first->second;
    return cindex;
}

void xAlci_main_window::oc_analyse_click()
{
    CGAL::set_pretty_mode(std::cout);
    CGAL::set_pretty_mode(std::cerr);
    
    arr_curves.clear();
    if(!input_poly(arr_curves, oc_input->text().ascii()))
        return;

    timer.reset();
    timer.start();

    arr_compute_arrangement();
}

void xAlci_main_window::arr_compute_arrangement() {

    oc_seg_list->clear();
    oc_objects.clear();
    for(Layers::iterator it= oc_layers.begin(); it != oc_layers.end();
            it++) {
        widget->detach(*it);
        delete (*it);
    } 
    oc_layers.clear();
    
    timer.reset();
    timer.start();
    
    Kernel_2::Construct_curve_2 cc_2 = kernel_2.construct_curve_2_object();
    typedef std::vector<Curve_analysis_2> Curve_vector;

    Curve_vector curves;
    int n = static_cast<int>(arr_curves.size());
    bool one_curve = (n == 1);
    
    for(int j = 0; j < n; j++) {

        Curve_analysis_2 tmp = cc_2(make_square_free(arr_curves[j]));
        curves.push_back(tmp);
    }

        typedef CGAL::Arrangement_2<CKvA_2> Arrangement;
        Arrangement arr;

        CGAL::insert(arr, curves.begin(), curves.end(), boost::false_type());
                  
        timer.stop();
        std::cout << "\n\nAnalyse elapsed time: " << timer.time() << std::endl;

        Color_map color_map;
        
        int i = 0, cindex;
        for(Arrangement::Vertex_const_iterator vit = arr.vertices_begin();
                vit != arr.vertices_end(); vit++) {

            std::ostringstream os;
            const Point_2& pt = vit->point();
            print_point(pt, os);

            if(!vit->is_isolated()) {
                continue;
            }
            
            os << "; isolated";
            oc_seg_list->insertItem(os.str());

            if(one_curve)
                cindex = i;            
            else
                cindex = pick_color(pt.curve().id(), color_map);

            oc_layers.push_back(new Graphic_layer(widget, i, cindex));
            oc_layers[i]->deactivate();
            oc_objects.push_back(CGAL::make_object(pt));
            i++;
        }

        for(Arrangement::Edge_const_iterator eit = arr.edges_begin();
                eit != arr.edges_end(); eit++, i++) {

            std::ostringstream os;
            const Arc_2& arc = eit->curve();
            print_arc(arc, os);
            oc_seg_list->insertItem(os.str());

            if(one_curve)
                cindex = i;            
            else
                cindex = pick_color(arc.curve().id(), color_map);

            oc_layers.push_back(new Graphic_layer(widget, i, cindex));
            oc_layers[i]->deactivate();
            oc_objects.push_back(CGAL::make_object(arc));
       }

    oc_rasterize_btn->setEnabled(true);
    widget->redraw();
}

Poly_int2 xAlci_main_window::make_square_free(const Poly_int2& poly) {
#if !AcX_SQRT_EXTENSION

    return kernel_2.decompose_2_object()(poly);
#else
    return poly; // no gcds available for sqrt_extensions
#endif
}

void xAlci_main_window::visualize()
{
    show();
    widget->zoom(1);
}

void xAlci_main_window::tab_changed(QWidget*) {

    curr_objects = &oc_objects;
}

#include "xalci.moc"

int main(int argc, char **argv)
{
    if (check_testsuite(argc, argv)) {
        return (0);
    }
    QApplication app(argc, argv);
    xAlci_main_window widget(1024, 700); // physical window size
    app.setMainWidget(&widget);
    widget.setCaption("Curve rasterizer demo");
    widget.setMouseTracking(TRUE);
    widget.visualize();
    return app.exec();
}

#endif // CGAL_HAS_QT3
