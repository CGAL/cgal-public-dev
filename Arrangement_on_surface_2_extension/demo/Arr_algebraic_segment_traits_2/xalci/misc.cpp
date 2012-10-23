// ============================================================================
//
// Copyright (c) 2001-2010 Max-Planck-Institut Saarbruecken (Germany).
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
   
#define NDEBUG 1

#include <CGAL/basic.h>
#include <CGAL/Timer.h>

#include <set>
#include <ctime>
#include <fstream>

#define CGAL_CKVA_CR_TIMING
#define CGAL_CKvA_LINE_THICKNESS 4
#define CGAL_CKvA_NO_AXES 1

#include "include/misc.h"
#include "include/xalci.h"

#include "axis.xpm"

#include <CGAL/Curved_kernel_via_analysis_2/Qt_widget_Curve_renderer_2.h>

#include <CGAL/Algebraic_kernel_d/Polynomial_parser_d.h>

#include <CGAL/Polynomial_type_generator.h>

CGAL::Timer tm_refine, tm_complete, tm_clip_pts;

/*QColor rasterize_colors[] = {
    QColor(0, 206, 209),
    QColor(50, 230, 0),
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
};*/
QColor rasterize_colors[] = {
    QColor(QRgb(0xcc0000)),
    QColor(QRgb(0xff3300)),
    QColor(QRgb(0xb3b300)),
    QColor(QRgb(0x00cc00)),
    QColor(QRgb(0x009999)),
    QColor(QRgb(0x0099ff)),
    QColor(QRgb(0x0000ff)),
    QColor(QRgb(0x9900cc)),
    QColor(QRgb(0xff0099))
};

int n_rast_colors = sizeof(rasterize_colors) / sizeof(QColor);

static CGAL::Bbox_2 bbox(0.0, 0.0, 0.0, 0.0);

Object_vector objects;
Object_vector* curr_objects;

Graphic_layer *subdiv_layer;
Layers layers;

QPixmap *subdiv_plot, *arcs_plot;
bool subdiv_layer_changed = true;

void Graphic_layer::draw()
{
#if 1
    QPainter *ppnt = &widget->get_painter();
    QPen old_pen = ppnt->pen();
    
    if(index == -2) {
        CGAL::Bbox_2 new_box(widget->x_min(), widget->y_min(),
                             widget->x_max(), widget->y_max());
        if(bbox != new_box) {
            bbox = new_box;
    
            //subdiv_renderer.setup(bbox,
                // widget->width(), widget->height());
            subdiv_layer_changed = true;
        }
        
        if(subdiv_layer_changed) {
            QPainter offscreen(subdiv_plot);
            //offscreen.setPen(QPen(Qt::black, 2, Qt::SolidLine, Qt::SquareCap,
              //  Qt::MiterJoin));
            subdiv_plot->fill();

            subdiv_layer_changed = false;
            typedef std::pair<int, int> Int_pair;
            std::list<Int_pair> points;
           
            //subdiv_renderer.draw(std::back_inserter(points));
            if(!points.empty()) {
                
                offscreen.setPen(QPen(Qt::black ,1));
                std::list<Int_pair>::iterator it = points.begin();
                while(it != points.end()) {
                    offscreen.drawPoint(it->first, it->second);
                    it++;
                }
            }
            subdiv_layer_changed = false;
        }
        ppnt->drawPixmap(0,0,*subdiv_plot);
        
    } else if(index == -1) { // this layer is dedicated to axis drawing

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

        const Arc_2& obj = (*curr_objects)[index];
//         Point_2 pt;

//#if !XALCI_USE_FLAT_COLOR_SCHEME
         QBrush b1(Qt::black);
/*#else
        QBrush b1(Qt::NoBrush);
#endif*/
        ppnt->setBrush(b1);
        
        tm_complete.start();
//         if(CGAL::assign(arc, obj)) 
        *widget << obj;
/*        else if(CGAL::assign(pt, obj))
           *widget << pt;
        else
            CGAL_error_msg("Malformed object found..\n");*/
        tm_complete.stop();
    }
    ppnt->setPen(old_pen);
#endif
}

void xAlci_main_window::activate_layers()
{
    Layers::iterator it = layers.begin();
    int i = 0;
    while(it != layers.end()) {
        if(seg_list->isSelected(i) || complete_check->isChecked()) {
            (*it)->activate();
        } else if((*it)->is_active()) {
            (*it)->deactivate();
        }
        it++; i++;
    }
}

void xAlci_main_window::rasterize_click()
{
    if(cur_method == 0) {
        activate_layers();
    } else {
//         CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
        ::CGAL::set_pretty_mode(std::cout);
        Poly_int2 f;

        if(!input_poly(f, input->text().ascii()))
            return;
        
        //subdiv_renderer.set_polynomial(f);
        subdiv_layer_changed = true;
    }

    tm_refine.reset();
    tm_clip_pts.reset();
    tm_complete.reset();
   
    widget->redraw();

    std::cout << "\nprocess end-points total: " << tm_refine.time() <<
        std::endl;
    std::cout << "\ncompute clip-points: " << tm_clip_pts.time() <<
        std::endl;
    std::cout << "\nrasterize complete: " << tm_complete.time() << std::endl;
    tm_refine.reset();
    tm_clip_pts.reset();
    tm_complete.reset();
}
   
bool xAlci_main_window::input_poly(Poly_int2& p, const char *ascii) {

    if(ascii == NULL)
        return false;

    typedef CGAL::Polynomial_type_generator< Rational, 2 >::Type Poly_rat_2;

    CGAL::Polynomial_parser_d< Poly_rat_2, CGAL::Mixed_floating_point_parser_policy< Poly_rat_2 > > parser;
    std::string str(ascii);

    if(str.length() == 0) {
        str = "y^7 + (-3)*y^6 + (2*x^2 + (-1)*x + 2)*y^5 + (x^3 + (-6)*x^2 + x + 2)*y^4 + (x^4 + (-2)*x^3 + 2*x^2 + x + (-3))*y^3 + (2*x^5 + (-3)*x^4 + x^3 + 10*x^2 + (-1)*x + 1)*y^2 + ((-1)*x^5 + 3*x^4 + 4*x^3 + (-12)*x^2)*y + (x^7 + (-3)*x^5 + (-1)*x^4 + (-4)*x^3 + 4*x^2)";
    }

    Poly_rat_2 prat;
    if(parser(str, prat)) {
        
        typedef CGAL::Fraction_traits< Poly_rat_2 > FTraits;
        FTraits::Denominator_type det(1);
        FTraits::Decompose decompose;
        decompose(prat, p, det);

    } else {
        std::cerr << "Parser error, trying another format..\n";
        try {
            std::stringstream ss(str);
            ss >> p;
        }
        catch(...) {
            std::cerr << "Invalid format of polynomial..\n";
            return false;
        }
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
    os << ")";
}

void xAlci_main_window::deactivate_layers()
{
    Layers::iterator it = layers.begin();
    while(it != layers.end())
        (*it++)->deactivate();
}

void xAlci_main_window::complete_toggle(bool on) {
    if(on && !objects.empty())
        rasterize_btn->setEnabled(true);
}


void xAlci_main_window::seg_list_click()
{
    rasterize_btn->setEnabled(true);
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
        //QBoxLayout *vbox = new QVBoxLayout(0,0,5);    
        //hbox->addLayout(vbox);
       
    tab_widget = new QTabWidget(central_widget);
    hbox->addWidget(tab_widget,4);
    one_curve_tab = new QFrame(tab_widget,"one_curve");

    curr_objects=&objects;
    tab_widget->addTab(one_curve_tab,"one_curve");

    *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor(CGAL::WHITE);
    resize(w,h);
    double ratio = 1.0;//(double)h/w;
    widget->set_window(-1, 1, -ratio, ratio, true);
        widget->setMouseTracking(TRUE);
    subdiv_layer = new Graphic_layer(widget, -2, 0);
    subdiv_layer->deactivate();
    axis = new Graphic_layer(widget, -1, 0);

        // ONE CURVE TAB
    QBoxLayout* vbox = new QVBoxLayout(one_curve_tab,10,10);
    vbox->addWidget(new QLabel("<b>Input polynomial:</b>",one_curve_tab));
    
    input = new QTextEdit("", QString::null,one_curve_tab);
    vbox->addWidget(input,6);
    analyse_btn = new QPushButton("Analyse",one_curve_tab);
    vbox->addWidget(analyse_btn);
    vbox->addWidget(new QLabel("<b>Curve segments:</b>",one_curve_tab));
               
    seg_list = new QListBox(one_curve_tab);
    seg_list->setSelectionMode(QListBox::Multi);
    vbox->addWidget(seg_list,6);

    complete_check =
        new QCheckBox("rasterize complete curve",one_curve_tab);
    vbox->addWidget(complete_check);
    
    method_box = new QComboBox("Rasterization method", one_curve_tab);
    method_box->insertItem("Segment Renderer");
    method_box->insertItem("Space Subdivision");
    method_box->setEditable(false);
        
    vbox->addWidget(method_box);
    rasterize_btn = new QPushButton("Rasterize",one_curve_tab);
    vbox->addWidget(rasterize_btn);

            
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

    connect(analyse_btn, SIGNAL(clicked()), SLOT(analyse_click()));
    connect(rasterize_btn, SIGNAL(clicked()), SLOT(rasterize_click()));
    connect(complete_check, SIGNAL(toggled(bool)),
        SLOT(complete_toggle(bool)));
    connect(seg_list, SIGNAL(selectionChanged()),
        SLOT(seg_list_click()));
    connect(method_box, SIGNAL(activated(int)), this,
        SLOT(switch_method(int)));

    stoolbar = new CGAL::Qt_widget_standard_toolbar(widget, this, "ST");
    rasterize_btn->setEnabled(false);

    QRect frect = frameGeometry();
    frect.moveCenter(QDesktopWidget().availableGeometry().center());
    move(frect.topLeft());

//     CORE::CORE_init(2);
//     CORE::setDefaultPrecision(150, CORE::extLong::getPosInfty());
}


#include "misc.moc"
