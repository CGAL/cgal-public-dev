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
#include <QCheckBox>
#include <QLabel>
#include <QLineEdit>

#include <CGAL/Qt/GraphicsViewNavigation.h>

#include "EventFilterManagerGroup.h"
#include "PolygonTableModel.h"
#include "MinkowskiSumCalculator.h"

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
  void onAction_StartMinkowskiSumStar();
  void onAction_ZoomToFit();
  void onAction_IncrementViewLayer();
  void onAction_DecrementViewLayer();
  void onAction_SetViewLayer();

  void onAction_hidePolygon(PolygonGraphicsItem* polygon);
  void onAction_showPolygon(PolygonGraphicsItem* polygon);

  void onAction_CheckStayOnMax(int state);
  void setMaxLayer(int k);

  PolygonListItem* createNewPolygon(PolygonWithHoles* polygon,
    QString name, int layer);
  PolygonListItem* createNewPolygon(PolygonWithHoles* polygon,
    QString name, QColor color, bool isVisble, int layer);
  PolygonListItem* createNewPolygon(PolygonWithHoles* polygon,
    QString name, QColor fillColor, QColor edgeColor, QColor vertexColor, bool isVisble, int layer);

private:
	Ui::MainWindow* ui;

	EventFilterManagerGroup* m_manager;
  PolygonTableModel* m_polygonList;
  QActionGroup* m_drawNavigateActionGroup;
  int m_viewLayer;
  QLineEdit* m_layerLine;
  int m_maxLayer;
  QLabel* m_maxLayerText;
  QCheckBox* m_stayOnMaxCheckBox;
  bool m_stayOnMax;

  QString m_lastFolder;

  MinkowSkiSumStarController* m_minkSumStarController;

  QList<PolygonWithHoles*>* readPolygonFile(QString aFilename);
  QRectF* getRectOfPolygons(QList<PolygonGraphicsItem*>* polygons);
  QList<PolygonGraphicsItem*>* getAllVisiblePolygons();
  StdPolygonList* getAllVisiblePolygonObjects();
  void setPolygonTextBrowser(int numberVertices, qreal Area);
  void setTextToAllVisiblePolygons();
  void removePolygonFromScene(PolygonGraphicsItem* polygon);
  void addPolygonToScene(PolygonGraphicsItem* polygon);
  void setViewLayer(int layer);
  void setMaxLayerText(int maxLayer);
  void setStayOnMax(bool shouldStay);

  void setScene();
  void setEventFilterManagerGroup();
  void setViewControl();
  void setPolygonInputDevice();
  void setDrawNavigateActionGroup();
  void connectActions();

  int static getNumberVertices(PolygonWithHoles* polygon);
  double static getArea(PolygonWithHoles* polygon);
};

#endif // CGAL_MAINWINDOW_H
