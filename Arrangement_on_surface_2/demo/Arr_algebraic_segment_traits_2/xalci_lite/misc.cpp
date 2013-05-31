// ============================================================================
//
// Copyright (c) 2001-2010 Max-Planck-Institut Saarbruecken (Germany).
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
// File          : demos/xalci/misc.cpp
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.12 $
// Revision_date : $Date: 2009-06-30 13:14:58 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-inf.mpg.de>
//                 
// ============================================================================

/*!\brief
 * declares miscellaneous routines in a separate file (to speed-up
 * compilation)
 */
   
#include <CGAL/basic.h>
#include <CGAL/Timer.h>

#include <set>
#include <ctime>
#include <fstream>

#define CGAL_CKvA_LINE_THICKNESS 2
#define CGAL_CKvA_NO_AXES 1

#include "include/misc.h"
#include "include/xalci.h"

#include "axis.xpm"

#include <CGAL/Curved_kernel_via_analysis_2/Qt_widget_Curve_renderer_2.h>

#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>

#include <CGAL/Polynomial_type_generator.h>

CGAL::Timer refine_timer;

QColor rasterize_colors[] = {
    QColor(139, 0, 0),
    QColor(138, 43, 226),
    QColor(95, 158, 160),
    QColor(0, 0, 139),
    QColor(205, 149, 12),
    QColor(0, 100, 0),
    QColor(139, 0, 139),
    QColor(85, 107, 47),
    QColor(255, 127, 0),
    QColor(0, 206, 209),
    QColor(238, 18, 137),
    QColor(238, 99, 99),
    QColor(205, 55, 0)
};
int n_rast_colors = 13;

static CGAL::Bbox_2 bbox(0.0, 0.0, 0.0, 0.0);
//static SoX::Subdiv_renderer subdiv_renderer;

Object_vector oc_objects;
Object_vector* curr_objects;

Graphic_layer *subdiv_layer;
Layers oc_layers;

QPixmap *subdiv_plot, *arcs_plot;
bool subdiv_layer_changed = true;

CGAL::Timer timer;

void draw_topology_graph(CGAL::Qt_widget* widget);

void Graphic_layer::draw()
{
    QPainter *ppnt = &widget->get_painter();
    QPen old_pen = ppnt->pen();

    if(index == -1) { // this layer is dedicated to axis drawing

#if !CGAL_CKvA_NO_AXES
        RasterOp old_raster = widget->rasterOp();
        widget->setRasterOp(XorROP);
        ppnt->setPen(QPen(QColor(150,150,0),1,Qt::DashDotDotLine));
        ppnt->moveTo(0,widget->y_pixel(0));           
        ppnt->lineTo(widget->width(),widget->y_pixel(0));
        ppnt->moveTo(widget->x_pixel(0),0);           
        ppnt->lineTo(widget->x_pixel(0),widget->height());
        widget->setRasterOp(old_raster);
#endif

    } else {
        widget->setRasterOp(CopyROP);
        ppnt->setPen(QPen(rasterize_colors[color_index % n_rast_colors],
             CGAL_CKvA_LINE_THICKNESS,  Qt::SolidLine, Qt::RoundCap,
                 Qt::MiterJoin));

        const CGAL::Object obj = (*curr_objects)[index];
        Arc_2 arc;
        Point_2 pt;

        QBrush b1(Qt::black);
        ppnt->setBrush(b1);
        
        timer.start();
        if(CGAL::assign(arc, obj)) 
           *widget << arc;
        else if(CGAL::assign(pt, obj))
           *widget << pt;
        else
            CGAL_error_msg("Malformed object found..\n");
        timer.stop();
    }
    ppnt->setPen(old_pen);
}


void xAlci_main_window::oc_activate_layers()
{
    Layers::iterator it = oc_layers.begin();
    int i = 0;
    while(it != oc_layers.end()) {
        if(oc_seg_list->isSelected(i) || oc_complete_check->isChecked()) {
            (*it)->activate();
        } else if((*it)->is_active()) {
            (*it)->deactivate();
        }
        it++; i++;
    }
}

void xAlci_main_window::oc_rasterize_click()
{
    oc_activate_layers();
    
    refine_timer.reset();
    timer.reset();
   
    widget->redraw();
    refine_timer.stop();

    std::cout << "\n\nrefine end-points: " << refine_timer.time() <<
        std::endl;
    std::cout << "\n\nRasterize elapsed time: " << timer.time() << std::endl;
}
   
bool xAlci_main_window::input_poly(
    std::vector< Poly_int2 >& polys, const char *ascii) {

    if(ascii == NULL || ascii == "")
    {
        ascii = "y^7 + (-3)*y^6 + (2*x^2 + (-1)*x + 2)*y^5 + (x^3 + (-6)*x^2 + x + 2)*y^4 + (x^4 + (-2)*x^3 + 2*x^2 + x + (-3))*y^3 + (2*x^5 + (-3)*x^4 + x^3 + 10*x^2 + (-1)*x + 1)*y^2 + ((-1)*x^5 + 3*x^4 + 4*x^3 + (-12)*x^2)*y + (x^7 + (-3)*x^5 + (-1)*x^4 + (-4)*x^3 + 4*x^2)";
    }
        
    typedef CGAL::Polynomial_type_generator< Rational, 2 >::Type Poly_rat_2;

    CGAL::Polynomial_parser_d< Poly_rat_2,
        CGAL::Mixed_floating_point_parser_policy< Poly_rat_2 > > parser;
    std::string str(ascii);


    typedef std::set< Poly_rat_2 > Polyset;
    Polyset poly_set;
    // tokenize input delimited by commas:

    char *saveptr, *token;
    int n_polys = 0;
    
    token = strtok_r((char *)ascii, ",", &saveptr);
    try {

    Poly_rat_2 tmp_rat;
    while(token != NULL) {
        std::string str(token);
        if(!parser(str, tmp_rat)) {
            std::cerr << "Syntax error while parsing " <<
                    n_polys << " polynomial\n";
            return false;
        }
        if(tmp_rat.is_zero()) {
            std::cerr << n_polys << ": zero polynomial\n";
            return false;
        }

        if(!poly_set.insert(tmp_rat).second)
            std::cerr << n_polys << "duplicate polynomial\n";
        n_polys++;
        token = strtok_r(NULL, ",", &saveptr);
    }
    } catch(...) {
        std::cerr << "Invalid polynomial\n";
        return false;
    }

    if(poly_set.size() == 0) { // this could happen ?
        std::cerr << "no polynomials found\n";
        return false;
    }

    typedef CGAL::Fraction_traits< Poly_rat_2 > FTraits;
    FTraits::Denominator_type det(1);
    FTraits::Decompose decompose;

    Polyset::const_iterator pit;
    for(pit = poly_set.begin(); pit != poly_set.end(); pit++) {

        Poly_int2 tmp_int;
        decompose(*pit, tmp_int, det);
        polys.push_back(tmp_int);
        
        std::cout << tmp_int << std::endl;
    }
    return true;
}

void xAlci_main_window::print_arc(const Arc_2& arc, std::ostream& os) {

    print_endpoint(arc, CGAL::ARR_MIN_END, os);
    os << " - ";
    print_endpoint(arc, CGAL::ARR_MAX_END, os);

    if(arc.is_vertical())
        os << "; vertical";
    else
        os << "; arcno: " << arc.arcno();
}

void xAlci_main_window::print_endpoint(const Arc_2& arc,
    CGAL::Arr_curve_end end, std::ostream& os) {

    CGAL::Arr_parameter_space loc = arc.location(end);

    os << "(x: ";
    if(loc == CGAL::ARR_LEFT_BOUNDARY)
        os << "-oo";
    else if(loc == CGAL::ARR_RIGHT_BOUNDARY)
        os << "+oo";
    else 
        os << CGAL::to_double(arc.curve_end_x(end));

    switch(loc) {
    case CGAL::ARR_BOTTOM_BOUNDARY:
        os << "; y: -oo)";
        break;

    case CGAL::ARR_TOP_BOUNDARY:
        os << "; y: +oo)";
        break;
        
    case CGAL::ARR_LEFT_BOUNDARY:
    case CGAL::ARR_RIGHT_BOUNDARY: {

        CGAL::Object obj =
            arc.curve().asymptotic_value_of_arc(loc, arc.arcno());
            
        CGAL::Arr_parameter_space loc2;
        if(CGAL::assign(loc2, obj)) 
            os << (loc2 == CGAL::ARR_BOTTOM_BOUNDARY ? "; y: -oo)" :
                "; y: +oo)");
        else {
            Kernel_2::Coordinate_1 y;
            if(CGAL::assign(y, obj)) 
                os << "; y: " << CGAL::to_double(y) << " (asym))";
            else
                CGAL_error_msg("Ill-typed object returned..\n");
        }
        break;
    }    
    default:
        os << "; arcno: " << arc.curve_end(end).arcno() << ')';
    }
}

void xAlci_main_window::print_point(const Point_2& pt,
    std::ostream& os) {

    os << "(x: " << CGAL::to_double(pt.x()) <<
        "; arcno: " << pt.arcno() << ")";
}

void xAlci_main_window::oc_deactivate_layers()
{
    Layers::iterator it = oc_layers.begin();
    while(it != oc_layers.end())
        (*it++)->deactivate();
}

void xAlci_main_window::oc_complete_toggle(bool on) {
    if(on && !oc_objects.empty())
        oc_rasterize_btn->setEnabled(true);
}

void xAlci_main_window::oc_seg_list_click()
{
    oc_rasterize_btn->setEnabled(true);
}

void xAlci_main_window::axis_toggle()
{
    widget->lock();
    axis->draw();
    widget->unlock();
    widget->repaint(false);
}

void xAlci_main_window::setup(int w, int h)
{
    err_msg_dialog = new QErrorMessage(this);

    central_widget = new QWidget(this);
    setCentralWidget(central_widget);

    subdiv_plot = new QPixmap(w, h);
    arcs_plot = new QPixmap(w, h);
    subdiv_plot->fill();
    arcs_plot->fill();
    
    subdiv_layer = NULL;
    QBoxLayout *hbox = new QHBoxLayout(central_widget, 10, 10);
    widget = new CGAL::Qt_widget(central_widget);
    hbox->addWidget(widget,8);
       
    tab_widget = new QTabWidget(central_widget);
    hbox->addWidget(tab_widget,4);
    one_curve_tab = new QFrame(tab_widget,"curve arrangement");

    curr_objects=&oc_objects;
    tab_widget->addTab(one_curve_tab,"curve arrangement");

    *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor(CGAL::WHITE);
    resize(w,h);
    double ratio = 1.0;//(double)h/w;
    widget->set_window(-1, 1, -ratio, ratio, true);
        widget->setMouseTracking(TRUE);
    subdiv_layer = new Graphic_layer(widget, -2, 0);
    subdiv_layer->deactivate();
    axis = new Graphic_layer(widget, -1, 0);

    // ONE CURVE TAB
    QBoxLayout* oc_vbox = new QVBoxLayout(one_curve_tab,10,10);
    oc_vbox->addWidget(new QLabel("<b>Input polynomial:</b>",one_curve_tab));
    
    oc_input = new QTextEdit("", QString::null,one_curve_tab);
    oc_vbox->addWidget(oc_input,6);
    oc_analyse_btn = new QPushButton("Analyse",one_curve_tab);
    oc_vbox->addWidget(oc_analyse_btn);
    oc_vbox->addWidget(new QLabel("<b>Curve arcs:</b>",one_curve_tab));
               
    oc_seg_list = new QListBox(one_curve_tab);
    oc_seg_list->setSelectionMode(QListBox::Multi);
    oc_vbox->addWidget(oc_seg_list,6);

    oc_complete_check =
        new QCheckBox("rasterize complete curve",one_curve_tab);
    oc_complete_check->setChecked(true);
    oc_vbox->addWidget(oc_complete_check);
    
    oc_method_box = new QComboBox("Rasterization method", one_curve_tab);
    oc_method_box->insertItem("Segment Renderer");
    oc_method_box->setEditable(false);
        
    oc_vbox->addWidget(oc_method_box);
    oc_rasterize_btn = new QPushButton("Rasterize",one_curve_tab);
    oc_vbox->addWidget(oc_rasterize_btn);

            
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), 
        CTRL+Key_Q );   
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), 0);
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );
 
    QToolBar *layers_toolbar = 
    new QToolBar("Tools", this, QMainWindow::DockTop, TRUE, "Tools");
    QToolButton *axis_button = new QToolButton(QPixmap(axis_xpm),
               "Show axis", 0, this, SLOT(axis_toggle()), layers_toolbar);
    axis_button->setToggleButton(true);
    axis_button->toggle();

    connect(tab_widget,SIGNAL(currentChanged(QWidget*)),
        SLOT(tab_changed(QWidget*)));

    //connect(widget, SIGNAL(rangesChanged()), SLOT(rasterize()));

    connect(oc_analyse_btn, SIGNAL(clicked()), SLOT(oc_analyse_click()));
    connect(oc_rasterize_btn, SIGNAL(clicked()), SLOT(oc_rasterize_click()));
    connect(oc_complete_check, SIGNAL(toggled(bool)), 
        SLOT(oc_complete_toggle(bool)));
    connect(oc_seg_list, SIGNAL(selectionChanged()),
        SLOT(oc_seg_list_click()));

    stoolbar = new CGAL::Qt_widget_standard_toolbar(widget, this, "ST");
    oc_rasterize_btn->setEnabled(false);

    tab_widget->showPage(one_curve_tab);

    QRect frect = frameGeometry();
    frect.moveCenter(QDesktopWidget().availableGeometry().center());
    move(frect.topLeft());

}


#include "misc.moc"
