// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
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
// Author(s) : Saar Katz <kats.saar@gmail.com>

#ifndef CGAL_MAINWINDOW_H
#define CGAL_MAINWINDOW_H

#include <QMainWindow>
#include <QGraphicsScene>
#include <QActionGroup>

#include <CGAL/Qt/GraphicsViewNavigation.h>

#include "EventFilterManagerGroup.h"
#include "PolygonTableModel.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
	explicit MainWindow(QWidget* parent = 0);
	~MainWindow();

public slots:
  void onAction_ImportPolygon();
  void onAction_DrawNewPolygon(bool checked);
  void onAction_NavigateScene(bool checked);
  void onAction_SelfMinkowskiSum();

private:
	Ui::MainWindow* ui;

	EventFilterManagerGroup* m_manager;
  PolygonTableModel* m_polygonList;
  QActionGroup* m_drawNavigateActionGroup;

  QList<PolygonWithHoles*>* readPolygonFile(QString aFilename);

  void setScene();
  void setEventFilterManagerGroup();
  void setViewControl();
  void setPolygonInputDevice();
  void setDrawNavigateActionGroup();
  void connectActions();
};

#endif // CGAL_MAINWINDOW_H
