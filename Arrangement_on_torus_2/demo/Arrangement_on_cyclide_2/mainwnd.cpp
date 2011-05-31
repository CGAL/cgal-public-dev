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
// File          : demos/xsurface/mainwnd.C
// QdX_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#include <CGAL/basic.h>

#include "include/arrangements.h"
#include "include/mainwnd.h"
#include "include/ogl_view.h"

#include <qapplication.h>
#include <qstatusbar.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qfiledialog.h>
#include <qlayout.h>
#include <qmessagebox.h>

static XSurface_view *ogl_view = NULL;

XSurface_main_wnd *mainwnd = NULL;

XSurface_arrangements *engine = NULL; 

enum {ASPECT_RATIO_4_3, ASPECT_RATIO_16_9};

static const struct { int w, h; } screen_res[] =  {
    {384, 288}, {768, 576}, {960, 720}, {1440, 1080},
    {512, 288}, {1024, 576}, {1280, 720}, {1920, 1080}
};

XSurface_main_wnd::XSurface_main_wnd(QWidget *parent, const char *name, 
        WFlags flags) : QMainWindow (parent, name, flags) {

    data_dir = "data";
    cur_resolution = -1;

    file = new QPopupMenu(this);
    menuBar()->insertItem("&File", file);

    file->insertItem("&Select base surface", this, 
        SLOT(load_base_surfaces()), Qt::CTRL+Qt::Key_O, 0);
    file->insertItem("&Open surface set", this, SLOT(read_surface_set()),
         Qt::CTRL+Qt::Key_A,1);
    file->insertItem("Add &more surfaces", this, SLOT(add_surface_set()),
         Qt::CTRL+Qt::Key_M,2);
   
    file->insertSeparator();
    file->insertItem("&Quit", qApp, SLOT(quit()), Qt::CTRL+Qt::Key_Q,4);

    file->setItemEnabled(1, false);
    file->setItemEnabled(2, false);
   
    view = new QPopupMenu (this);	
    menuBar()->insertItem("&View", view);
    view->insertItem("Show &grid", this, 
                      SLOT(toggle_show_grid()), Qt::CTRL+Qt::Key_G,0);
    view->setItemChecked(0, false);

    view->insertItem("Show &base surface", this,
                        SLOT(toggle_show_base_surf()), Qt::CTRL+Qt::Key_B, 2);
    view->insertItem("Show &intersection curves", this,
                        SLOT(toggle_show_curves()), Qt::CTRL+Qt::Key_I, 3);
    view->insertItem("Show &event points", this, 
                      SLOT(toggle_show_points()), Qt::CTRL+Qt::Key_E,4);
    view->setItemChecked(2, false);
    view->setItemChecked(3, true);
    view->setItemChecked(4, true);

    view->insertSeparator();

    view->insertItem("Outer circle", this,
                        SLOT(toggle_show_outer_circle()), 0, 12);
    view->insertItem("Tube circle", this, 
                      SLOT(toggle_show_tube_circle()), 0, 13);
    view->insertItem("Pole", this, 
                      SLOT(toggle_show_pole()), 0, 14);
    view->setItemChecked(12, false);
    view->setItemChecked(13, false);
    view->setItemChecked(14, false);   

    view->insertSeparator(); 
    
    view->insertItem("Draw &wireframe", this,
                       SLOT(toggle_fill_mode()), Qt::ALT+Qt::Key_W, 5);
    view->insertItem("White background", this,
                      SLOT(toggle_bg_color()), 0, 6);
    view->insertItem("Flat color", this,
                      SLOT(toggle_flat_color()), 0, 15);

    view->insertSeparator();
    view->insertItem("Enable &lighting", this,
            SLOT(toggle_lighting()), Qt::ALT+Qt::Key_L, 7);
    view->setItemChecked(7, true);
    view->insertItem("&Transparency", this,
            SLOT(toggle_transparency()), Qt::ALT+Qt::Key_T, 8);
    view->setItemChecked(8, false);
    view->insertItem("&Bumpy surface", this,
            SLOT(toggle_bumpy()), Qt::ALT+Qt::Key_B, 9);
    view->setItemChecked(8, false);
    
    view->insertSeparator();
    view->insertItem("Plain mode", this,
                      SLOT(toggle_2d3d_view()), 0, 11);
    view->insertItem("Reset view", this,
                      SLOT(reset_view()), 0, 10);
    view->insertItem("Screenshot", this, SLOT(screenshot()), 
         Qt::ALT+Qt::Key_S, 16);

    action = new QPopupMenu(this);   
    menuBar()->insertItem("&Action", action);
    action->insertItem ("CGAL s&weep on base surface", this,
             SLOT(compute_arr_on_surface()), Qt::CTRL+Qt::Key_W, 0);
    action->setItemEnabled(0, false);

    action->insertItem ("CGAL overlay arrangements", this,
             SLOT(compute_overlay_arr()), 0, 1);
    action->setItemEnabled(1, false);
    
    action->insertItem ("&Render arrangement mesh", this,
             SLOT(render_arr_mesh()), Qt::CTRL+Qt::Key_R, 2);
    action->setItemEnabled(2, false);
    
    video = new QPopupMenu(this);
    menuBar()->insertItem("&Video", video);
//     video->insertItem("Snap&hot", this,
//                       SLOT(snapshot()), Qt::Key_H,0);
//     video->insertSeparator();
//     video->insertItem ("Start &video grabbing", this,
//                       SLOT(toggle_video_grabbing()), Qt::Key_V,1);

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
 
//     QToolBar *layers_toolbar = 
//     new QToolBar("Tools", this, QMainWindow::DockTop, TRUE, "Tools");
//     QToolButton *axis_button = new QToolButton(QPixmap(axis_xpm),
//                "Show axis", 0, this, SLOT(axis_toggle()), layers_toolbar);
//     axis_button->setToggleButton(true);
//     axis_button->toggle();

    // Create a nice frame to put around the OpenGL widget
    QFrame* frame = new QFrame(this, "frame");
    frame->setFrameStyle(QFrame::Sunken | QFrame::Panel);
    frame->setLineWidth(2);
   
    mainwnd = this;
    ogl_view = new XSurface_view(frame, "glbox");
    engine = new XSurface_arrangements;

    QHBoxLayout* flayout = new QHBoxLayout(frame, 2, 2, "flayout");
    flayout->addWidget(ogl_view, 1);
    setCentralWidget(frame);

    set_resolution(2);
    
    statusBar()->message("Ready", 2000);
    qApp->processEvents();
}

void XSurface_main_wnd::compute_arr_on_surface() {
    
    engine->compute_arr_on_surface();
    action->setItemEnabled(1, true);
    action->setItemEnabled(2, true);
}

void XSurface_main_wnd::compute_overlay_arr() {

    int res = read_surf_dlg(1, "second surface set", true);
    if(res != 1) 
        return;
    
    engine->compute_overlay_arr();
    action->setItemEnabled(2, true);
}

void XSurface_main_wnd::render_arr_mesh() {
    double w = 10;
    CGAL::Bbox_2 bbox(-w, -w, w, w);
    engine->render(bbox, 600, 600);
    ogl_view->tube_mesh_computed = false;
    ogl_view->updateGL();
}

void XSurface_main_wnd::set_resolution(int i) {
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

void XSurface_main_wnd::reset_view() {
    ogl_view->reset_view();    
    ogl_view->updateGL();
}

void XSurface_main_wnd::toggle_2d3d_view() {
// 	if (ogl_view->in_plain_mode()) {
// 		ogl->set_3d_mode();
// 		view->changeItem(11, "Plain mode");
// 	} else {
// 		ogl->set_plain_mode();
// 		view->changeItem(11, "3D mode");
// 	}
//     view->setItemChecked(7, ogl->draw_base_quadric);
//     view->setItemChecked(8, ogl->draw_all_quadrics);
//     view->setItemChecked(9, ogl->draw_transparent!=0);
  //  ogl->updateGL();
}

void XSurface_main_wnd::toggle_fill_mode() {

    ogl_view->solid_fill ^= 1;
    if(ogl_view->solid_fill) 
        view->changeItem(5, "Draw wireframe");
    else
        view->changeItem(5, "Draw solid");
    ogl_view->updateGL();
}

void XSurface_main_wnd::toggle_bg_color() {

    ogl_view->black_bg_color ^= 1;
    if(ogl_view->black_bg_color) 
        view->changeItem(6, "White background");
    else
        view->changeItem(6, "Black background");
    ogl_view->updateGL();
}

void XSurface_main_wnd::toggle_flat_color() {

    ogl_view->flat_color ^= 1;
    if(ogl_view->flat_color) 
        view->changeItem(15, "Different colors");
    else
        view->changeItem(15, "Flat color");
    ogl_view->updateGL();
}

void XSurface_main_wnd::screenshot() {

    static unsigned ii = 0;
    char bbf[256];
    sprintf(bbf, "fluffy_bunny_%d.png", ii);

    ogl_view->updateGL();
    ogl_view->take_screenshot(bbf);
}

void XSurface_main_wnd::toggle_show_grid() {
    ogl_view->show_grid ^= 1;
    view->setItemChecked(0, ogl_view->show_grid);
    ogl_view->updateGL();
}

void XSurface_main_wnd::toggle_show_base_surf() {
    ogl_view->draw_base_surf ^= 1;
    view->setItemChecked(2, ogl_view->draw_base_surf);
    ogl_view->updateGL();
}

void XSurface_main_wnd::toggle_show_curves() {
    ogl_view->draw_curves ^= 1;
    view->setItemChecked(3, ogl_view->draw_curves);
    ogl_view->updateGL();
}

void XSurface_main_wnd::toggle_show_points() {
    ogl_view->draw_points ^= 1;
    view->setItemChecked(4, ogl_view->draw_points);
    ogl_view->updateGL();
}

void XSurface_main_wnd::toggle_show_outer_circle() {
    ogl_view->draw_outer_circle ^= 1;
    view->setItemChecked(12, ogl_view->draw_outer_circle);
    ogl_view->updateGL();
}

void XSurface_main_wnd::toggle_show_tube_circle() {
    ogl_view->draw_tube_circle ^= 1;
    view->setItemChecked(13, ogl_view->draw_tube_circle);
    ogl_view->updateGL();
}

void XSurface_main_wnd::toggle_show_pole() {
    ogl_view->draw_pole ^= 1;
    view->setItemChecked(14, ogl_view->draw_pole);
    ogl_view->updateGL();
}

void XSurface_main_wnd::toggle_lighting() {
    ogl_view->lighting ^= 1;
    view->setItemChecked(7, ogl_view->lighting);
    ogl_view->updateGL();
}

void XSurface_main_wnd::toggle_transparency() {
    ogl_view->transparency ^= 1;
    view->setItemChecked(8, ogl_view->transparency);
    ogl_view->updateGL();
}

void XSurface_main_wnd::toggle_bumpy() {
    ogl_view->bumpy_surf ^= 1;
    view->setItemChecked(9, ogl_view->bumpy_surf);
    ogl_view->updateGL();
}

void XSurface_main_wnd::toggle_show_sweepline() {
    //ogl_view->show_sweep_curve=!ogl_view->show_sweep_curve;
    //ogl_view->updateGL();
}

void XSurface_main_wnd::load_base_surfaces() {
    
    data_dir = "data";

    QFileDialog dlg(data_dir.absPath() + "/Cyclides/", NULL, this, 
                   "open base surface", "Choose a file" );

    dlg.setMode(QFileDialog::ExistingFile);
    if(dlg.exec() == QDialog::Accepted) {
        data_dir = *dlg.dir();
        bool res = engine->read_base_surfaces(dlg.selectedFile().latin1());
        file->setItemEnabled(1, res);
        file->setItemEnabled(2, res);     
        action->setItemEnabled(1, false);
        action->setItemEnabled(2, false);

        ogl_view->surf_mesh_computed = false;   
        ogl_view->aux_mesh_computed = false;
        ogl_view->tube_mesh_computed = false;
        engine->clear_approximations();
        ogl_view->updateGL();
    } 
}

void XSurface_main_wnd::read_surface_set() {
    
    int res = read_surf_dlg(0, "create a new surface set", true);
    if(res != -1) {
        action->setItemEnabled(0, (res == 1));
        action->setItemEnabled(1, false);
        action->setItemEnabled(2, false);
        ogl_view->updateGL();
    }
}

void XSurface_main_wnd::add_surface_set() {
    
    int res = read_surf_dlg(0, "add more surfaces", false);
    if(res != -1) {
        action->setItemEnabled(0, (res == 1));
        action->setItemEnabled(1, false);
        action->setItemEnabled(2, false);
        ogl_view->updateGL();
    }
}

//!\c which - which surface set to read in (0 or 1)
//!\c clear_flag - whether to clear existing surface
//! @return -1 dialog cancelled; 0 - accepted with errors; 1 - no errors
int XSurface_main_wnd::read_surf_dlg(int which, const char *msg,
                                     bool clear_flag) {
  
    data_dir = "data";

    std::cout << "dir: " << data_dir.absPath() << std::endl;

    QFileDialog dlg(data_dir.absPath() + "/Surfaces/", 
                    NULL, this, msg, "Choose file");

    dlg.setMode(QFileDialog::ExistingFile);
    if(dlg.exec() != QDialog::Accepted) 
        return -1;

    data_dir = *dlg.dir();    
    return engine->read_surface_set(dlg.selectedFile().latin1(), which,
             clear_flag);
}

void XSurface_main_wnd::snapshot() {

}

void XSurface_main_wnd::howto() {
    QMessageBox::about(this, "How to ?",
      "Got in trouble ? Don't know what to do ? Call 911 !");
}

void XSurface_main_wnd::about() {
    QMessageBox::about(this, "About",
      "The greatest demo program ever seen by humanity..");
}

void XSurface_main_wnd::aboutQt() {
    QMessageBox::aboutQt(this, "About Qt");
}

XSurface_main_wnd::~XSurface_main_wnd () {
    delete ogl_view;
    delete engine;
}

