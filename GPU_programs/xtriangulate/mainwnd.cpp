// ============================================================================
//
// Copyright (c) 2001-2007 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/);
// you may redistribute it under the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with EXACUS.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRcomANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : QdX
// File          : demos/xtriangulate/mainwnd.C
// QdX_release   : $Name:  $
// Revision      : $Revision: 1.7 $
// Revision_date : $Date: 2009-07-24 13:21:29 $
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#include <CGAL/basic.h>

#include "include/skeletonizer_interface.h"
#include "include/mainwnd.h"
#include "include/ogl_view.h"

#include <qapplication.h>
#include <qstatusbar.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qfiledialog.h>
#include <qlayout.h>
#include <qmessagebox.h>

#include <fstream>

//! moc: /KM/projects/ecg/software/linux/g++-4.1/qt-x11-free-3.3.6/bin/moc include/mainwnd.h -o mainwnd.moc

//! /KM/projects/ecg/software/linux/g++-4.1/qt-x11-free-3.3.6/bin/moc include/ogl_view.h -o ogl_view.moc

static OGL_view *ogl_view = NULL;

XTri_main_wnd *mainwnd = NULL;

XSkeletonizer *skeletonizer = NULL; 

enum {ASPECT_RATIO_4_3, ASPECT_RATIO_16_9};

static const struct { int w, h; } screen_res[] =  {
    {384, 288}, {768, 576}, {960, 720}, {1440, 1080},
    {512, 288}, {1024, 576}, {1280, 720}, {1920, 1080}
};

// ~/work/binaries/bin/xtriangulate --slices_x=200 --slices_y=200 --left_enlarge=2.5 --right_enlarge=2.5 --btm_enlarge=-7.0 --top_enlarge=7.0

static const char *ps[] = {"-sx=", "-sy=", 
      "-fpbits=", "-l=", "-r=",
        "-b=", "-t=", "-below=", "-above=",
        "-abs", "-help"};
static const int n_ps = sizeof(ps) / sizeof(const char *);

void usage(const char *name) {

    std::cout << "\nUSAGE: " << name;
    for(int i = 0; i < n_ps; i++)
         std::cout << " [" << ps[i] << "?]";  
    std::cout << "\n\n";
}

static unsigned slices_x = 1, slices_y = 1, fp_bits = 24;

static double left_enlarge = 0.1, right_enlarge = 0.1,
        btm_enlarge = 0.1, top_enlarge = 0.1,
        below_bnd = -2.0, above_bnd = 2.0;

static bool use_auto_bounds = true;

void parse_params(int argc, char **argv) {
    
    int *ls = new int[n_ps];
    for(int i = 0; i < n_ps; i++)
        ls[i] = strlen(ps[i]);
    
    for(int i = 1; i < argc; i++) {

        if(!strncmp(argv[i], ps[0], ls[0])) 
            slices_x = atoi(argv[i] + ls[0]);
        
        else if(!strncmp(argv[i], ps[1], ls[1])) 
            slices_y = atoi(argv[i] + ls[1]);

        else if(!strncmp(argv[i], ps[2], ls[2])) 
            fp_bits = atoi(argv[i] + ls[2]);

        else if(!strncmp(argv[i], ps[3], ls[3])) 
            left_enlarge = atof(argv[i] + ls[3]);

        else if(!strncmp(argv[i], ps[4], ls[4])) 
            right_enlarge = atof(argv[i] + ls[4]);

        else if(!strncmp(argv[i], ps[5], ls[5])) 
            btm_enlarge = atof(argv[i] + ls[5]);

        else if(!strncmp(argv[i], ps[6], ls[6])) 
            top_enlarge = atof(argv[i] + ls[6]);
 
        else if(!strncmp(argv[i], ps[7], ls[7])) 
            below_bnd = atof(argv[i] + ls[7]);

        else if(!strncmp(argv[i], ps[8], ls[8])) 
            above_bnd = atof(argv[i] + ls[8]);

        else if(!strncmp(argv[i], ps[9], ls[9]))
            use_auto_bounds = false;

        else {
            if(strcmp(argv[i], "-help"))
                std::cerr << "Unknown option: " << argv[i] << std::endl;
            usage(argv[0]); exit(1);
        }
    }
    delete []ls;
}

XTri_main_wnd::XTri_main_wnd(QWidget *parent, const char *name, 
        WFlags flags) : QMainWindow (parent, name, flags) {

    data_dir = "data";
    cur_resolution = -1;

    file = new QPopupMenu(this);
    menuBar()->insertItem("&File", file);

    file->insertItem("&Open surface..", this, 
        SLOT(load_surf()), Qt::CTRL+Qt::Key_O, 0);
    file->insertItem("O&verlay with..", this, 
        SLOT(overlay_surf()), Qt::CTRL+Qt::Key_V, 1);

    file->insertSeparator();
    file->insertItem("&Quit", qApp, SLOT(quit()), Qt::CTRL+Qt::Key_Q,4);

    file->setItemEnabled(1, false);
    //file->setItemEnabled(2, false);
   
    view = new QPopupMenu (this);	
    menuBar()->insertItem("&View", view);
    view->insertItem("Show &grid", this, 
                      SLOT(toggle_show_grid()), Qt::CTRL+Qt::Key_G,0);
    view->setItemChecked(0, false);
    view->insertItem("Show triangles", this, 
                      SLOT(toggle_show_triangles()), 0, 1);
    view->setItemChecked(1, true);
    view->insertItem("Show normals", this,
                      SLOT(toggle_show_normals()), 0, 13);
    view->setItemChecked(13, false);

    view->insertSeparator();
   
    view->insertItem("Draw &wireframe", this,
                       SLOT(toggle_fill_mode()), Qt::ALT+Qt::Key_W, 5);
    view->insertItem("White background", this,
                      SLOT(toggle_bg_color()), 0, 6);

    view->insertSeparator();
    view->insertItem("Enable &lighting", this,
            SLOT(toggle_lighting()), Qt::ALT+Qt::Key_L, 7);
    view->setItemChecked(7, true);
    view->insertItem("&Transparency", this,
            SLOT(toggle_transparency()), Qt::ALT+Qt::Key_T, 8);
    view->insertItem("&Inverse normals", this,
            SLOT(toggle_normals()), Qt::ALT+Qt::Key_I, 9);   
    view->insertItem("&Bumpy ;)", this,
            SLOT(toggle_bumpy_surf()), 0, 4);   

    view->insertSeparator();
    view->insertItem("Plain mode", this,
                      SLOT(toggle_2d3d_view()), 0, 11);
    view->insertItem("Reset view", this,
                      SLOT(reset_view()), 0, 10);

    action = new QPopupMenu(this);   
    menuBar()->insertItem("&Action", action);

    action->insertItem ("Skele&tonize body", this,
             SLOT(triangulate_surf()), Qt::CTRL+Qt::Key_T, 1);
    
    action->insertSeparator();
    action->insertItem ("&Save STL model", this,
             SLOT(write_STL()), Qt::CTRL+Qt::Key_S, 2);

    action->setItemEnabled(1, false);    
    action->setItemEnabled(2, false);    

    video = new QPopupMenu(this);
    menuBar()->insertItem("&Video", video);
    QPopupMenu *resolution = new QPopupMenu(this);
    video->insertItem("Resolution", resolution, 0, 2);
    
    video->insertItem ("&Take snapshot", this, SLOT(snapshot()),
             0, 3);

    video4x3 = new QPopupMenu(this);
    resolution->insertItem("Aspect ratio 4:3", video4x3, 0, 0);    
    video4x3->insertItem("384x288", 0);
    video4x3->insertItem("720x576", 1);
    video4x3->insertItem("960x720", 2);
    video4x3->insertItem("1440x1080", 3);
    QObject::connect(video4x3, SIGNAL(activated(int)), this,
                        SLOT(set_resolution(int)));

    video16x9 = new QPopupMenu(this);
    resolution->insertItem("Aspect ratio 16:9", video16x9, 1, 1);
    video16x9->insertItem("512x288", 4);
    video16x9->insertItem("1024x576",5);
    video16x9->insertItem("1280x720",6);
    video16x9->insertItem("1920x1080",7);
    QObject::connect(video16x9, SIGNAL(activated(int)), this,
                        SLOT(set_resolution(int)));

    QPopupMenu *help = new QPopupMenu(this);
    menuBar()->insertItem("&Help", help);
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), 0);
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );
 
    // Create a nice frame to put around the OpenGL widget
    QFrame* frame = new QFrame(this, "frame");
    frame->setFrameStyle(QFrame::Sunken | QFrame::Panel);
    frame->setLineWidth(2);
   
    mainwnd = this;
    
//     static unsigned slices_x = 1, slices_y = 1, fp_bits = 64; 
// static double left_enlarge = 1.25, right_enlarge = 1.25,
//         btm_enlarge = -2.0, top_enlarge = 2.0;

    parse_params(qApp->argc(), qApp->argv());

//     unsigned n_slices_x = 1, n_slices_y = 1, fp_bits = 128;
//     if(qApp->argc() >= 2)
//         n_slices_x = atoi(qApp->argv()[1]);
//     if(qApp->argc() >= 3)
//         n_slices_y = atoi(qApp->argv()[2]);
//     if(qApp->argc() >= 4)
//         fp_bits = atoi(qApp->argv()[3]);
 
    skeletonizer = new XSkeletonizer();
    skeletonizer->reset_parameters(slices_x, slices_y, left_enlarge,
         right_enlarge, btm_enlarge, top_enlarge, below_bnd, above_bnd);

    ogl_view = new OGL_view(frame, "glbox");

    QHBoxLayout* flayout = new QHBoxLayout(frame, 2, 2, "flayout");
    flayout->addWidget(ogl_view, 1);
    setCentralWidget(frame);

    set_resolution(2);
    
    statusBar()->message("Ready", 2000);
    qApp->processEvents();
}


void XTri_main_wnd::set_resolution(int i) {
    if(cur_resolution == i)
        return;
    
    QPopupMenu *mm;
    if(cur_resolution != -1) {
        mm = (aspect_ratio == ASPECT_RATIO_4_3 ? video4x3 : video16x9);
        mm->setItemChecked(cur_resolution, false);
    }
    if(i <= 3) {
        aspect_ratio = ASPECT_RATIO_4_3;
        mm = video4x3;
    } else {
        aspect_ratio = ASPECT_RATIO_16_9;
        mm = video16x9;
    }
    mm->setItemChecked(i, true);
    cur_resolution = i;

    QRect desk = QApplication::desktop()->screenGeometry();
    resize(screen_res[i].w, screen_res[i].h);
    move((desk.width() - screen_res[i].w) / 2,
         (desk.height() - screen_res[i].h) / 2);
}

void XTri_main_wnd::reset_view() {
    ogl_view->reset_view();    
    ogl_view->updateGL();
}

void XTri_main_wnd::toggle_2d3d_view() {

}

void XTri_main_wnd::toggle_fill_mode() {

    ogl_view->solid_fill ^= 1;
    if(ogl_view->solid_fill) 
        view->changeItem(5, "Draw wireframe");
    else
        view->changeItem(5, "Draw solid");
    ogl_view->updateGL();
}

void XTri_main_wnd::toggle_bg_color() {

    ogl_view->black_bg_color ^= 1;
    if(ogl_view->black_bg_color) 
        view->changeItem(6, "White background");
    else
        view->changeItem(6, "Black background");
    ogl_view->updateGL();
}


void XTri_main_wnd::toggle_show_grid() {
    ogl_view->show_grid ^= 1;
    view->setItemChecked(0, ogl_view->show_grid);
    ogl_view->updateGL();
}

void XTri_main_wnd::toggle_show_triangles() {
    ogl_view->show_triangles ^= 1;
    view->setItemChecked(1, ogl_view->show_triangles);
    ogl_view->updateGL();
}

void XTri_main_wnd::toggle_lighting() {
    ogl_view->lighting ^= 1;
    view->setItemChecked(7, ogl_view->lighting);
    ogl_view->updateGL();
}

void XTri_main_wnd::toggle_transparency() {
    ogl_view->transparency ^= 1;
    view->setItemChecked(8, ogl_view->transparency);
    ogl_view->updateGL();
}

void XTri_main_wnd::toggle_normals() {
    ogl_view->inverse_normals ^= 1;
    ogl_view->surf_mesh_computed = false;
    view->setItemChecked(9, ogl_view->inverse_normals);
    ogl_view->updateGL();
}

void XTri_main_wnd::toggle_show_normals() {
    ogl_view->show_normals ^= 1;
    ogl_view->surf_mesh_computed = false;
    view->setItemChecked(13, ogl_view->show_normals);
    ogl_view->updateGL();
}

void XTri_main_wnd::toggle_bumpy_surf() {
    ogl_view->bumpy_surf ^= 1;
    view->setItemChecked(4, ogl_view->bumpy_surf);
    ogl_view->updateGL();
}

void XTri_main_wnd::triangulate_surf() {

    skeletonizer->skeletonize(use_auto_bounds);
    action->setItemEnabled(2, true);
    ogl_view->surf_mesh_computed = false; 
    ogl_view->updateGL();
}

void XTri_main_wnd::load_surf() {
    
    QString filename = QFileDialog::getOpenFileName(
        data_dir.absPath(), QString::null, this, "dummy", "Open surface data");

    if(filename.isEmpty())
        return;

    bool res = skeletonizer->load_surface(filename);
    action->setItemEnabled(1, res);
    action->setItemEnabled(2, false);
    if(!res) {
        std::cerr << "Unable to load surface\n";
        return;
    } 
}

void XTri_main_wnd::overlay_surf() {

    QString filename = QFileDialog::getOpenFileName(
        data_dir.absPath(), QString::null, this, "dummy", 
                "Open overlay surface data");
    if(filename.isEmpty())
        return;

    bool res = skeletonizer->overlay_surface(filename);
    action->setItemEnabled(1, res);
    action->setItemEnabled(2, false);
    if(!res) {
        std::cerr << "Unable to load overlay surface\n";
        return;
    } 
}

void XTri_main_wnd::write_STL() {

    QString filename = QFileDialog::getSaveFileName(QDir::currentDirPath(),
             QString::null, this, "dummy", "Save STL model"); 
    if(filename.isEmpty())
        return;

    if(QFile::exists(filename) && QMessageBox::warning(this, "?!?",
        "Overwrite existing file ?", "&Yes", "&No", QString::null, 1, 1) != 0)
        return;
        
    std::ofstream os(filename);
    if(!os)
        return;
    ogl_view->write_STL(os);
}

void XTri_main_wnd::snapshot() {

}

void XTri_main_wnd::howto() {
    QMessageBox::about(this, "How to ?",
      "Got in trouble ? Don't know what to do ? Call 911 !");
}

void XTri_main_wnd::about() {
    QMessageBox::about(this, "About",
      "The greatest demo program ever seen by humanity..");
}

void XTri_main_wnd::aboutQt() {
    QMessageBox::aboutQt(this, "About Qt");
}

XTri_main_wnd::~XTri_main_wnd () {
    if(ogl_view != NULL)
        delete ogl_view;
    if(skeletonizer != NULL)
        delete skeletonizer;
}

#include "mainwnd.moc"
