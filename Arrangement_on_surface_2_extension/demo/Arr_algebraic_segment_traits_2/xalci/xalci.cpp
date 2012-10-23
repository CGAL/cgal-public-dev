// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
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

#include "include/misc.h"
#include "include/xalci.h"

#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>

#define XALCI_USE_FLAT_COLOR_SCHEME 0

// static const Kernel_2& kernel_2 = CKvA_2::instance().kernel();

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

extern Object_vector objects;
extern Object_vector* curr_objects;
extern bool subdiv_layer_changed;
 
extern Graphic_layer *subdiv_layer;
extern Layers cad_layers, layers, arr_layers;

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

void xAlci_main_window::switch_method(int index)
{
    if(index != cur_method) {
        cur_method = index;
        if(cur_method == 0) { // segment renderer
            subdiv_layer->deactivate();
            activate_layers();
            analyse_btn->setEnabled(true);
            seg_list->setEnabled(true);
            complete_check->setEnabled(true);
        } else {            // space subdivision
            deactivate_layers();
            subdiv_layer->activate();
            analyse_btn->setEnabled(false);
                        seg_list->setEnabled(false);
            complete_check->setEnabled(false);
            rasterize_btn->setEnabled(true);
        }
        widget->redraw();
    }
}

Lite_CA_2 lca;

void xAlci_main_window::analyse_click()
{
    int i;
    seg_list->clear();
    objects.clear();
    for(Layers::iterator it= layers.begin(); it != layers.end(); it++) {
        widget->detach(*it);
        delete (*it);
    }
    
    layers.clear();
    rasterize_btn->setEnabled(false);
    
    Poly_int2 f;
    if(!input_poly(f, input->text().ascii()))
        return;

    CGAL::Polynomial_traits_d< Poly_int2 >::Differentiate diff;

    Poly_int2 fx = diff(f, 0), fxx = diff(f, 0), fxy = diff(fx, 1),
        fy = diff(f, 1), fyy = diff(fy, 1);
//    Poly_int2 ress = fxx*fy*fy - ((fx*fy*fxy)*Poly_int1(2,0)) + fyy*fx*fx;
    Poly_int2 res1 = f*fx, res2 = f*fy;
  
    CGAL::set_pretty_mode(std::cout);
//     std::cout << "curv:\n " << ress << "\n\n";
//     std::cout << "fx:\n " << fx << "\n\n";
//     std::cout << "fy:\n " << fy << "\n\n";
    std::cout << "f*fx:\n " << res1 << "\n\n";
    std::cout << "f*fy:\n " << res2 << "\n\n";

    CGAL::set_pretty_mode(std::cout);
    std::cout << "f:\n " << f << "\n\n";
    
    CGAL::Timer tm;
    tm.start();

    objects.clear();
    lca = Lite_CA_2(f);

    std::vector< Arc_2 > arcs;
    lca.compute(arcs);
    tm.stop();

    std::cout << arcs.size() << " arcs found..\n";
    std::cout << "##### Analysis time elapsed: " << tm.time() << "\n";

    objects = arcs;

    {
#if 0
    Kernel_2::Construct_curve_2 cc_2 = kernel_2.construct_curve_2_object();
    
Lbegin:
    try {
        Curve_analysis_2 curve = cc_2(f);
        CKvA_2::instance().make_x_monotone_2_object()(curve,
            std::back_inserter(objects));

    } catch(...) {
        std::cout << "Seemingly non-square free.. restarting..\n";
        f = make_square_free(f);

        std::cout << "Square-free part: " << f << "\n";
        goto Lbegin;
    }

        tm.stop();
        std::cout << "\n\nAnalyse elapsed time: " << tm.time() << std::endl;
//         std::cout << "\nResultant time: " << CGAL::res_tm.time() << "\n";
        std::cout << objects.size() << " arcs found (incl isolated points)"
             << std::endl;
        
        Arc_2 arc;
        Point_2 pt;
#endif
        Object_vector::const_iterator oit;
        for(oit = objects.begin(), i = 0; oit != objects.end();
                oit++, i++) {

            std::ostringstream os;
//             if(CGAL::assign(arc, *oit))
                print_arc(*oit, os);
//             else if(CGAL::assign(pt, *oit)) {
//                 print_point(pt, os);
//                 os << "; isolated";
//             } else
//                 CGAL_error_msg("Malformed object found..\n");

            seg_list->insertItem(os.str());
            int cindex;
#if !XALCI_USE_FLAT_COLOR_SCHEME
            cindex = i;
#else
            cindex = 1;
#endif
            layers.push_back(new Graphic_layer(widget, i, cindex));
            layers[i]->deactivate();
        }
    } 
    subdiv_layer_changed = true;
    //layers.push_back(new Graphic_layer(-1,widget));
    rasterize_btn->setEnabled(complete_check->isChecked());
    widget->redraw();
}

Poly_int2 xAlci_main_window::make_square_free(const Poly_int2& poly) {
#if !AcX_SQRT_EXTENSION

    return poly;//kernel_2.decompose_2_object()(poly);
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
