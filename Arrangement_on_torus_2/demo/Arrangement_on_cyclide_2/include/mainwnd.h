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
// Library       : QdX
// File          : demos/xsurface/include/mainwnd.h
// QdX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef MAINWND_H
#define MAINWND_H

#include "include/includes_common.h"
#include <qmainwindow.h>
#include <qdir.h>
#include <qpopupmenu.h>
#include "include/ogl_view.h"

//#include <qtimer.h>

class XSurface_main_wnd : public QMainWindow {

    Q_OBJECT
    
public:
    
    XSurface_main_wnd(QWidget* parent = 0, const char* name = 0, 
                         WFlags flags = WType_TopLevel);
    
    ~XSurface_main_wnd();
    
public slots:

    void load_base_surfaces(); //! displays dialog to load a set of base surfs
    void read_surface_set(); //! opens a set of surfaces intersecting base 
    void add_surface_set();  //! adds more surfaces   
    
    void compute_arr_on_surface();
    void compute_overlay_arr();
    void render_arr_mesh();

    void reset_view();
    void snapshot();

    void toggle_2d3d_view();
    void toggle_show_sweepline();
    
    void toggle_fill_mode();
    void toggle_bg_color();
    void toggle_lighting();
    void toggle_transparency();
    void toggle_bumpy();
    void toggle_flat_color();

    void toggle_show_grid();    
    void toggle_show_base_surf();
    void toggle_show_curves();
    void toggle_show_points();

    void toggle_show_outer_circle();
    void toggle_show_tube_circle();
    void toggle_show_pole();

    void screenshot();
    
    void about();   
    void aboutQt();
    void howto();   
    
    void set_resolution(int);

protected:
    
    int read_surf_dlg(int which, const char *msg, bool clear_flag);
   
    QDir data_dir;
    QPopupMenu *file;
    QPopupMenu *view;
    QPopupMenu *video;
    QPopupMenu *action;
    QPopupMenu *video4x3, *video16x9;

    int aspect_ratio; 
    int cur_resolution; //! index of the currently set resolution
};

#endif // MAINWND_H
