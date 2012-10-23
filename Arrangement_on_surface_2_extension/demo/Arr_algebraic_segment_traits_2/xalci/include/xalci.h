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
// File          : demos/xalci/include/xalci.h
// SoX_release   : $Name:  $
// Revision      : $Revision: 1.37 $
// Revision_date : $Date: 2009-06-30 13:14:59 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef XALCI_H
#define XALCI_H

#include <iostream>
#include <vector>

#include <qplatinumstyle.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qerrormessage.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qtimer.h>
#include <qlistbox.h>
#include <qlayout.h>
#include <qpushbutton.h>
#include <qlabel.h>
#include <qtextedit.h>
#include <qcheckbox.h>
#include <qcombobox.h> 
#include <qpixmap.h>
#include <qpainter.h>
#include <qtabwidget.h>
#include <qhbuttongroup.h>
#include <qradiobutton.h>
#include <qvbox.h>
#include <qfiledialog.h>
#include <qprogressdialog.h>

// #define CGAL_RS_OLD_INCLUDES

// #if CGAL_BISOLVE_USE_CGAL_SYMBOLIC
#include <CGAL/symbolic_exports.h>
// #endif

#define CGAL_MODULAR_FILTER_OFF     // do not use modular filter

//! HACK HACK HACK: \c CGAL_HAS_THREADS is undefined
//! in Installation/include/CGAL/config.h !!!
#include "misc.h"

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/flags.h>
#include <CGAL/Algebraic_kernel_d_1_generator.h>

// #include <CGAL/Arrangement_2.h>
// #include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h>

#define CGAL_LITE_CA_VERBOSE 1
#include <CGAL/Curved_kernel_via_analysis_2/Lite_curve_analysis_2.h>

#include <CGAL/Bbox_2.h>

#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>

// typedef CGAL::CORE_arithmetic_kernel AK;
typedef CGAL::GMP_arithmetic_kernel AK;
typedef AK::Rational Rational;
typedef AK::Integer Integer;

typedef Integer Coefficient;


// #if !CGAL_BISOLVE_USE_RS_ISOLATOR

// typedef CGAL::Algebraic_kernel_d_1_generator<Coefficient, Rational>
//     ::Algebraic_kernel_with_qir_and_bitstream_1 Algebraic_kernel_1;

// typedef CGAL::Algebraic_kernel_d_1_generator<Coefficient, Rational>
//     ::Algebraic_kernel_with_qir_and_descartes_1 Algebraic_kernel_1;

// typedef CGAL::Algebraic_kernel_d_1_generator<Coefficient, Rational>
//     ::Algebraic_kernel_with_bisection_and_descartes_1 Algebraic_kernel_1;

typedef CGAL::Algebraic_kernel_d_1_generator<Coefficient, Rational>
    ::Algebraic_kernel_with_qir_and_rs_1 Algebraic_kernel_1;

typedef CGAL::Lite_curve_analysis_2< Algebraic_kernel_1 > Lite_CA_2;

// types of supporting polynomials
typedef Lite_CA_2::Polynomial_2 Poly_int2;
typedef Lite_CA_2::Polynomial_1 Poly_int1;

// typedef CGAL::Curved_kernel_via_analysis_2<Kernel_2> CKvA_2;
// typedef CKvA_2::X_monotone_curve_2 Arc_2;
// typedef CKvA_2::Point_2 Point_2;

typedef Lite_CA_2::Arc_2 Arc_2;

typedef std::vector< Arc_2 > Object_vector;

class xAlci_main_window : public QMainWindow
{
    Q_OBJECT

public:

    xAlci_main_window(int w, int h) : bbox(0.0, 0.0, 0.0, 0.0), cur_method(0) {
        setup(w, h);
    }

    void visualize();

public slots:

    void analyse_click();
    void rasterize_click();
    void seg_list_click();
    void complete_toggle(bool);
    void switch_method(int);

    void axis_toggle();

    void new_instance()
    {
        widget->lock();
        widget->clear_history();
        widget->set_window(-1.1, 1.1, -1.1, 1.1);
        // set the Visible Area to the Interval
        widget->unlock();
    }

protected slots:

    void about()
    {
        QMessageBox::about(this, "About",
            "This program demonstrates the use of Segment renderern"
            "and space subdivision to visualize algebraic curvesn");
    }
    void aboutQt()
    {
        QMessageBox::aboutQt(this, "About Qt");
    }
    void howto()
    {
        CGAL::Qt_help_window *help =
          new CGAL::Qt_help_window("help/index.html", ".", 0, "help viewer");
        help->resize(400, 400);
        help->setCaption("Demo HowTo");
        help->show();
    }
    void new_window()
    {
        xAlci_main_window *ed = new xAlci_main_window(500, 500);
        ed->setCaption("Layer");
        ed->widget->clear_history();
        ed->widget->set_window(-1.1, 1.1, -1.1, 1.1);
        ed->show();
    }

    void tab_changed(QWidget*);
    void setup(int, int);

    void print_arc(const Arc_2& arc, std::ostream& os);
    void print_endpoint(const Arc_2& arc,
        CGAL::Arr_curve_end end, std::ostream& os);
//     void print_point(const Point_2& pt, std::ostream& os);
        
    void activate_layers();
    void deactivate_layers();
protected:

    Poly_int2 make_square_free(const Poly_int2& poly);

    bool input_poly(Poly_int2& p, const char *ascii);

    CGAL::Bbox_2 bbox;
    QErrorMessage* err_msg_dialog;

    QPushButton *analyse_btn, *rasterize_btn;
          
    QCheckBox *cad_complete_check,*complete_check, *arr_complete_check;
    QListBox *seg_list;
          
    QTabWidget *tab_widget;
    QFrame* one_curve_tab;
    
    QTextEdit *input;
    QComboBox *method_box;
    Graphic_layer *axis;
        
    QWidget *central_widget;
    CGAL::Qt_widget *widget;
    CGAL::Qt_widget_standard_toolbar *stoolbar;

    int cur_method;
};

#endif


