// Copyright (c) 2019  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Kaimo Hu

#ifndef CGAL_VIEWER_H
#define CGAL_VIEWER_H

#include <QMap>
#include <CGAL/Qt/qglviewer.h>

// forward declarations
class QWidget;
class Scene;
class Viewer : public CGAL::QGLViewer {

  Q_OBJECT

 public:
  explicit Viewer(QWidget *parent);
  Viewer(const Viewer &) = delete;
  Viewer &operator = (const Viewer &) = delete;
  Viewer(const Viewer &&) = delete;
  Viewer &operator = (const Viewer &&) = delete;

  // overload several QGLViewer virtual functions
  void draw();
  void initializeGL();
  void setScene(Scene *pScene);

 protected:
  /*virtual void mousePressEvent(QMouseEvent *e);
  virtual void mouseReleaseEvent(QMouseEvent *e);*/

 private:
  Scene *m_pScene;
  bool m_custom_mouse;
}; // end class Viewer

#endif // CGAL_VIEWER_H
