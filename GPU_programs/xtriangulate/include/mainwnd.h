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
// File          : demos/xtriangulate/include/mainwnd.h
// QdX_release   : $Name:  $
// Revision      : $Revision: 1.4 $
// Revision_date : $Date: 2009-07-24 13:21:30 $
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
//#include <qtimer.h>

class XTri_main_wnd : public QMainWindow {

    Q_OBJECT
    
public:
    
    XTri_main_wnd(QWidget* parent = 0, const char* name = 0, 
                         WFlags flags = WType_TopLevel);
    
    ~XTri_main_wnd();
    
public slots:

    void load_surf();
    void overlay_surf();
    void triangulate_surf();
    void write_STL();
    
    void reset_view();
    void snapshot();

    void toggle_2d3d_view();
    
    void toggle_fill_mode();
    void toggle_bg_color();
    void toggle_lighting();
    void toggle_transparency();
    void toggle_normals();
    void toggle_bumpy_surf();
    void toggle_show_grid();    
    void toggle_show_triangles();
    void toggle_show_normals();

    void about();    
    void aboutQt();
    void howto();   
    
    void set_resolution(int);

protected:
   
    QDir data_dir;
    QPopupMenu *file;
    QPopupMenu *view;
    QPopupMenu *video;
    QPopupMenu *action;
    QPopupMenu *video4x3, *video16x9;

    int aspect_ratio; 
    int cur_resolution; 
};

#endif // MAINWND_H
