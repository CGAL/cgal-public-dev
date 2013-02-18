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

#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qerrormessage.h>
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

#define CGAL_BISOLVE_ENABLE_ARCAVOID 1 // default TODO?
#define CGAL_ACK_BITSTREAM_USES_E08_TREE 1 // do not change
#define CGAL_BISOLVE_USE_ADJUSTABLE_PRECISION 1 // 1 is default

#define CGAL_ACK_CURVE_ANALYSES_BISOLVE_USE_TEISSIER 1
#define CGAL_ACK_CURVE_ANALYSES_USE_BISOLVE 1

#define CGAL_BISOLVE_DEBUG 0
#define CGAL_BISOLVE_VERBOSE 0

#define CGAL_BISOLVE_USE_RS_AK 0
#define CGAL_BISOLVE_USE_RS_ISOLATOR 1 // 1 is default

#define CGAL_ACK_DEBUG_FLAG 0
#define CGAL_ACK_DEBUG_PRINT std::cout

#define CGAL_BISOLVE_USE_RESULTANT_COFACTORS 1
#define CGAL_BISOLVE_ARRANGEMENTS 1 // 0 is default

#define CGAL_BISOLVE_USE_GMP  1
#define CGAL_BISOLVE_USE_CORE 0

#define CGAL_AK_USE_OLD_BITSTREAM_DESCARTES 0
#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 1

#define CGAL_MODULAR_FILTER_OFF

#include <CGAL/Arithmetic_kernel.h>

#if CGAL_BISOLVE_USE_GMP
#include <CGAL/GMP_arithmetic_kernel.h>
typedef CGAL::GMP_arithmetic_kernel AK;
#endif

#if CGAL_BISOLVE_USE_CORE
#include <CGAL/CORE_arithmetic_kernel.h>
typedef CGAL::CORE_arithmetic_kernel AK;
#endif

/** **************************************************************************/

#define CGAL_BISOLVE_USE_BIGCD 1
#define CGAL_BIGCD_USE_SHIFT 0
#define CGAL_BIGCD_CHECK_SANITY 1

#define CGAL_BISOLVE_USE_GPU_RESULTANTS 1 // default?
#define CGAL_BISOLVE_CHECK_GPU_RESULTANTS_SANITY 1 // default 0

#define CGAL_BISOLVE_USE_GPU_GCDS 1  // default?
#define CGAL_BISOLVE_CHECK_GPU_GCDS_SANITY 1 // default 1

#define CGAL_BISOLVE_USE_NTL  1 // default 1 ??

#include <CGAL/symbolic_standalone.h>

/** **************************************************************************/


#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

#include <CGAL/Algebraic_kernel_2/Rounding_ak_d_1.h>

#include <CGAL/Algebraic_kernel_d/Generic_isolator.h>

#if CGAL_BISOLVE_USE_RS_AK
#include <CGAL/Algebraic_kernel_d/Float_traits.h>
#include <CGAL/Algebraic_kernel_rs_gmpz_d_1.h>
#else
#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d_1_generator.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_quadratic_refinement_rep_bfi.h>
#endif

#include <CGAL/Algebraic_kernel_d_2.h>
// #include <CGAL/Arcavoid_root_isolator.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h>

#include <CGAL/Bbox_2.h>

#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>

#include "misc.h"

typedef AK::Rational Rational;
typedef AK::Integer Integer;

typedef Integer Coefficient;

#if CGAL_BISOLVE_USE_RS_AK
  typedef Algebraic_kernel_rs_gmpz_d_1 Algebraic_kernel_d_1;
#else

#if CGAL_BISOLVE_USE_RS_ISOLATOR
#warning using RS root isolator !
   typedef CGAL::Algebraic_kernel_d_1_generator< Integer, Rational >
    ::Algebraic_kernel_with_qir_and_rs_1 Internal_ak_1;

    typedef CGAL::internal::Rounding_ak_d_1< Internal_ak_1 > AK_1;

#else
#warning using bitstream root isolator !
    typedef CGAL::Algebraic_kernel_d_1
    < Integer, Rational,
      CGAL::internal::Algebraic_real_quadratic_refinement_rep_bfi
           < Integer, Rational >,
      CGAL::internal::Bitstream_descartes
        < CGAL::internal::Bitstream_descartes_rndl_tree_traits
            < CGAL::internal::Bitstream_coefficient_kernel< Integer > >
        >
    > AK_1;
#endif
#endif // !CGAL_BISOLVE_USE_RS_AK

typedef CGAL::Algebraic_curve_kernel_2< AK_1 > Kernel_2;
typedef Kernel_2::Curve_analysis_2 Curve_analysis_2;

typedef Curve_analysis_2::Polynomial_2 Poly_int2;
typedef Poly_int2::NT Poly_int1;

typedef CGAL::Curved_kernel_via_analysis_2<Kernel_2> CKvA_2;
typedef CKvA_2::X_monotone_curve_2 Arc_2;
typedef CKvA_2::Point_2 Point_2;

// typedef Curve_analysis_2::Status_line_1 Status_line_1;
// typedef Kernel_2::Algebraic_real_1 X_coordinate_1;
// typedef Kernel_2::Algebraic_real_2 Xy_coordinate_2;

typedef std::vector<CGAL::Object> Object_vector;

class xAlci_main_window : public QMainWindow
{
    Q_OBJECT

public:

    xAlci_main_window(int w, int h) : bbox(0.0, 0.0, 0.0, 0.0), cur_method(0) {
        setup(w, h);
    }

    void visualize();

public slots:

    void oc_analyse_click();
    void oc_rasterize_click();
    void oc_seg_list_click();
    void oc_complete_toggle(bool);

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
    void print_point(const Point_2& pt, std::ostream& os);
        
    void oc_activate_layers();

    void oc_deactivate_layers();

    
protected:

    Poly_int2 make_square_free(const Poly_int2& poly);

    bool input_poly(std::vector< Poly_int2 >& p, const char *ascii);

    void arr_compute_arrangement();

    CGAL::Bbox_2 bbox;
  
    std::vector< Poly_int2 > arr_curves;
    Curve_selection_dialog* curve_selection_dialog;

    QErrorMessage* err_msg_dialog;

    QPushButton *cad_analyse_btn, *cad_rasterize_btn, *cad_file_search,
          *cad_partial_selection,
          *arr_analyse_btn, *arr_rasterize_btn, *arr_file_search, 
          *arr_partial_selection,
          *oc_analyse_btn, *oc_rasterize_btn;
          
//     QHButtonGroup *arr_method;
//     QRadioButton *arr_cgal, *arr_leda;
    QCheckBox *oc_complete_check;
    QListBox *oc_seg_list;
    //*arr_edge_list, *arr_node_list;
          
    QTabWidget *tab_widget;
    QFrame* one_curve_tab, *cad_tab, *arr_tab;
//     QLabel *arr_node_label,*arr_edge_label;
    
    QTextEdit *oc_input;
    QComboBox *oc_method_box;
    Graphic_layer *axis;
        
    QWidget *central_widget;
    CGAL::Qt_widget *widget;
    CGAL::Qt_widget_standard_toolbar *stoolbar;

    int cur_method;
};

#endif


