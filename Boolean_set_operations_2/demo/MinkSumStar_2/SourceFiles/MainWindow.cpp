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

#include <fstream>

#include <QFileDialog>

#include "MainWindow.h"
#include "ui_MainWindow.h"

#include "Typedefs.h"
#include "ChoiseDialogPolygon.h"

#define CGAL_POLYGON_EDGE_DEFAULT_COLOR "blue"
#define CGAL_POLYGON_FILL_DEFAULT_COLOR "blue"
#define CGAL_VIEW_DEFAULT_MARGIN_PRECENT 0.1

MainWindow::MainWindow(QWidget* parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow)
{
  ui->setupUi(this);
  setAcceptDrops(true);

  qRegisterMetaType<StdPolygonList>("StdPolygonList");

  m_polygonList = new PolygonTableModel(this);
  m_layerLine = new QLineEdit();
  m_maxLayerText = new QLabel();
  m_stayOnMaxCheckBox = new QCheckBox();
  m_viewLayer = 0;
  m_lastFolder = "./";

  m_minkSumStarController = new MinkowSkiSumStarController();
  connect(m_minkSumStarController,
    SIGNAL(requestNewPolygonCreation(PolygonWithHoles*, QString, int)),
    this,
    SLOT(createNewPolygon(PolygonWithHoles*, QString, int)));

  //ui->PolygonsTableView->setModel(m_polygonList);
  setPolygonTextBrowser(0, 0);
  setMaxLayer(1);
  setStayOnMax(true);

  setScene();
  setEventFilterManagerGroup();
  setViewControl();
  //setPolygonInputDevice();
  setDrawNavigateActionGroup();
  connectActions();

  m_layerLine->setMaximumWidth(30);
  m_layerLine->setText(QString::number(m_viewLayer));
  m_layerLine->setValidator(new QRegExpValidator(QRegExp(tr("-?[0-9]*"))));
  ui->toolBar->insertWidget(ui->actionIncrementLayer, m_layerLine);

  ui->toolBar->addWidget(m_maxLayerText);

  m_stayOnMaxCheckBox->setText("Stay on max");
  ui->toolBar->addWidget(m_stayOnMaxCheckBox);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::onAction_ImportPolygon() {
  QFileInfo file = QFileInfo(QFileDialog::getOpenFileName(
    this,
    tr("Import polygon"),
    m_lastFolder,
    tr("All types (*.*)")));
  if (file.exists())
  {
    QList<PolygonWithHoles*>* pwhList = readPolygonFile(file.filePath());
    QList<PolygonGraphicsItem*>* polygonList = new QList<PolygonGraphicsItem*>();
    if (pwhList != NULL)
    {
      int numPolygons = pwhList->count();
      int i = 1;
      for (QList<PolygonWithHoles*>::iterator pwh = pwhList->begin(); pwh != pwhList->end(); ++pwh) {
        PolygonListItem* pli = createNewPolygon(
          *pwh,
          numPolygons > 1 ? file.baseName() + "_" + QString::number(i) : file.baseName(),
          m_viewLayer);
        polygonList->append(pli->getPolygonItem()->graphics());
        i++;
      }
      QRectF* zoomRect = getRectOfPolygons(polygonList);
      ui->SceneGraphicsView->setSceneRect(*zoomRect);
      ui->SceneGraphicsView->fitInView(*zoomRect, Qt::KeepAspectRatio);
    }
    setTextToAllVisiblePolygons();
    m_lastFolder = file.path();
  }
}
void MainWindow::onAction_DrawNewPolygon(bool checked)
{
  if (checked)
    m_manager->setCurrentFilterName("draw");
}
void MainWindow::onAction_NavigateScene(bool checked)
{
  if (checked)
    m_manager->setCurrentFilterName("navigate");
}
void MainWindow::onAction_StartMinkowskiSumStar() {
  StdPolygonList* allPolygons = getAllVisiblePolygonObjects();
  if (allPolygons->size() > 0) {
    m_minkSumStarController->setInput(*allPolygons);
    m_minkSumStarController->runCalculation(100);
    ui->actionSelfMinkowskiSum->setDisabled(true);
  }
}
void MainWindow::onAction_ZoomToFit() {
  QList<PolygonGraphicsItem*>* poligonsList = getAllVisiblePolygons();
  if (poligonsList->count() > 0) {
    QRectF* zoomRect = getRectOfPolygons(poligonsList);
    ui->SceneGraphicsView->setSceneRect(*zoomRect);
    ui->SceneGraphicsView->fitInView(*zoomRect, Qt::KeepAspectRatio);
  }
}
void MainWindow::onAction_IncrementViewLayer()
{
  setStayOnMax(false);
  setViewLayer(m_viewLayer + 1);
}
void MainWindow::onAction_DecrementViewLayer()
{
  setStayOnMax(false);
  setViewLayer(m_viewLayer - 1);
}
void MainWindow::onAction_SetViewLayer()
{
  setStayOnMax(false);
  int layer = m_layerLine->text().toInt();
  setViewLayer(layer);
}

PolygonListItem* MainWindow::createNewPolygon(PolygonWithHoles* polygon,
  QString name, int layer)
{
  return createNewPolygon(polygon, name, CGAL_POLYGON_EDGE_DEFAULT_COLOR, true, layer);
}
PolygonListItem* MainWindow::createNewPolygon(PolygonWithHoles* polygon,
  QString name, QColor color, bool isVisible, int layer)
{
  return createNewPolygon(polygon, name, color, color, Qt::transparent, isVisible, layer);
}
PolygonListItem* MainWindow::createNewPolygon(PolygonWithHoles* polygon,
  QString name, QColor fillColor, QColor edgeColor, QColor vertexColor, bool isVisble, int layer)
{

  PolygonListItem* pli = new PolygonListItem(
    name,
    edgeColor,
    false,
    new PolygonItem(polygon));
  m_polygonList->appendItem(pli);
  int height = pli->getPolygonItem()->graphics()->boundingRect().height();
  int width = pli->getPolygonItem()->graphics()->boundingRect().width();
  int maxSide = height > width ? height : width;
  pli->getPolygonItem()->graphics()->setBrush(QBrush(fillColor));
  pli->getPolygonItem()->graphics()->setEdgesPen(QPen(QBrush(edgeColor), 0.005 * maxSide));
  pli->getPolygonItem()->graphics()->setVerticesPen(QPen(vertexColor));

  connect(pli, SIGNAL(hidePolygon(PolygonGraphicsItem*)),
    this, SLOT(onAction_hidePolygon(PolygonGraphicsItem*)));
  connect(pli, SIGNAL(showPolygon(PolygonGraphicsItem*)),
    this, SLOT(onAction_showPolygon(PolygonGraphicsItem*)));

  pli->setLayer(layer);
  if (layer > m_maxLayer) {
    setMaxLayer(layer);
  }
  pli->setVisible(isVisble && pli->getLayer() == m_viewLayer);
  return pli;
}

void MainWindow::onAction_hidePolygon(PolygonGraphicsItem* polygon) {
  removePolygonFromScene(polygon);
}
void MainWindow::onAction_showPolygon(PolygonGraphicsItem* polygon) {
  addPolygonToScene(polygon);
}

void MainWindow::onAction_CheckStayOnMax(int state)
{
  switch (state) {
  case Qt::Unchecked:
    setStayOnMax(false);
    break;
  case Qt::Checked:
    setStayOnMax(true);
    break;
  }
}

void MainWindow::setMaxLayer(int k)
{
  if (k > 0) {
    m_maxLayer = k;
  }
  else {
    m_maxLayer = 1;
  }
  if (m_stayOnMax || m_viewLayer > m_maxLayer) {
    setViewLayer(m_maxLayer);
  }
  setMaxLayerText(m_maxLayer);
}

void MainWindow::setScene()
{
  QGraphicsScene* scene = new QGraphicsScene();
  scene->setItemIndexMethod(QGraphicsScene::NoIndex);
  scene->setSceneRect(-50, -50, 100, 100);
  ui->SceneGraphicsView->setScene(scene);
  ui->SceneGraphicsView->setMouseTracking(true);
  ui->SceneGraphicsView->scale(1, -1);
}
void MainWindow::setEventFilterManagerGroup()
{
  m_manager = new EventFilterManagerGroup();
  m_manager->addObjectToWatch("viewport", ui->SceneGraphicsView->viewport());
  m_manager->addObjectToWatch("scene", ui->SceneGraphicsView->scene());
}
void MainWindow::setViewControl()
{
  m_manager->addFilterWidget("viewport", "navigate", new CGAL::Qt::GraphicsViewNavigation());
}
void MainWindow::setPolygonInputDevice()
{
  //manager->addFilterWidget("scene", "draw",
  //	new PolylineInput(ui->SceneGraphicsView, ui->SceneGraphicsView->scene, 0, false));
}
void MainWindow::setDrawNavigateActionGroup()
{
  m_drawNavigateActionGroup = new QActionGroup(this);
  m_drawNavigateActionGroup->addAction(ui->actionNavigateScene);
  //m_drawNavigateActionGroup->addAction(ui->actionDrawNewPolygon);
}
void MainWindow::connectActions()
{
  connect(ui->actionNavigateScene, SIGNAL(toggled(bool)),
          this, SLOT(onAction_NavigateScene(bool)));
  //connect(ui->actionDrawNewPolygon, SIGNAL(toggled(bool)),
  //        this, SLOT(onAction_DrawNewPolygon(bool)));
  ui->actionNavigateScene->setChecked(true);

  connect(ui->actionImportPolygon, SIGNAL(triggered()),
          this, SLOT(onAction_ImportPolygon()));

  connect(ui->actionSelfMinkowskiSum, SIGNAL(triggered()),
          this, SLOT(onAction_StartMinkowskiSumStar()));

  connect(ui->actionZoomToFit, SIGNAL(triggered()), this, SLOT(onAction_ZoomToFit()));
  connect(ui->actionIncrementLayer, SIGNAL(triggered()), this, SLOT(onAction_IncrementViewLayer()));
  connect(ui->actionDecrementLayer, SIGNAL(triggered()), this, SLOT(onAction_DecrementViewLayer()));
  connect(m_layerLine, SIGNAL(returnPressed()), this, SLOT(onAction_SetViewLayer()));
  connect(m_stayOnMaxCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onAction_CheckStayOnMax(int)));
}




QList<PolygonWithHoles*>* MainWindow::readPolygonFile(QString aFilename)
{
  std::ifstream in_file(qPrintable(aFilename));

  if (in_file) {
    unsigned int n_regions;
    in_file >> n_regions;
    auto polygonList = new QList<PolygonWithHoles*>();
    for (unsigned int r = 0; r < n_regions; ++r) {
      unsigned int n_boundaries;
      in_file >> n_boundaries;
      Polygon outer;
      std::vector<Polygon> holes;
      for (unsigned int r = 0; r < n_boundaries; ++r) {
        Polygon p;
        in_file >> p;
        if (r == 0) {
          outer = p;
        }
        else {
          holes.push_back(p);
        }
      }
      polygonList->append(new PolygonWithHoles(outer, holes.begin(), holes.end()));
    }
    return polygonList;
  }
  return NULL;
}
QRectF* MainWindow::getRectOfPolygons(QList<PolygonGraphicsItem*>* polygons)
{
  boost::optional<QRectF> totalRect;

  for (QList<PolygonGraphicsItem*>::iterator polygon = polygons->begin(); polygon != polygons->end(); ++polygon)
  {
    QRectF rect = (*polygon)->boundingRect();
    if (totalRect)
      totalRect = *totalRect | rect;
    else totalRect = rect;
  }

  if (totalRect)
  {
    qreal left = CGAL_VIEW_DEFAULT_MARGIN_PRECENT * totalRect->width();
    qreal top = CGAL_VIEW_DEFAULT_MARGIN_PRECENT * totalRect->height();
    qreal right = CGAL_VIEW_DEFAULT_MARGIN_PRECENT * totalRect->width();
    qreal bottom = CGAL_VIEW_DEFAULT_MARGIN_PRECENT * totalRect->height();
    // QMarginsF margins = QMarginsF(left, top, right, bottom);
    // *totalRect = totalRect->marginsAdded(margins);
  }

  return new QRectF(*totalRect);
}
QList<PolygonGraphicsItem*>* MainWindow::getAllVisiblePolygons() {
  QList<PolygonGraphicsItem*>* allVisiblePolygons = new QList<PolygonGraphicsItem*>();
  for (int i = 0; i < m_polygonList->rowCount(); i++) {
    if (m_polygonList->getPolygonItem(i)->isVisible()) {
      allVisiblePolygons->append(m_polygonList->getPolygonItem(i)->getPolygonItem()->graphics());
    }
  }
  return allVisiblePolygons;
}
StdPolygonList* MainWindow::getAllVisiblePolygonObjects()
{
  StdPolygonList* allVisiblePolygonObjects = new StdPolygonList();
  for (int i = 0; i < m_polygonList->rowCount(); i++) {
    if (m_polygonList->getPolygonItem(i)->isVisible()) {
      allVisiblePolygonObjects->push_back(*m_polygonList->getPolygonItem(i)->getPolygonItem()->polygon());
    }
  }
  return allVisiblePolygonObjects;
}

void MainWindow::setPolygonTextBrowser(int numberVertices, qreal Area) {
  ui->PolygonInfoTextBrowser->setText((boost::format("Number of vertices: %1%\nArea: %2%\n")
    % numberVertices
    % Area).str().c_str());
}
void MainWindow::setTextToAllVisiblePolygons() {
  int totalVertices = 0;
  qreal totalArea = 0;
  for (int i = 0; i < m_polygonList->rowCount(); i++) {
    if (m_polygonList->getPolygonItem(i)->isVisible()) {
      totalVertices += getNumberVertices(m_polygonList->getPolygonItem(i)->getPolygonItem()->polygon());
      totalArea += getArea(m_polygonList->getPolygonItem(i)->getPolygonItem()->polygon());
    }
  }
  setPolygonTextBrowser(totalVertices, totalArea);
}

void MainWindow::removePolygonFromScene(PolygonGraphicsItem* polygon) {
  ui->SceneGraphicsView->scene()->removeItem(polygon);
}
void MainWindow::addPolygonToScene(PolygonGraphicsItem* polygon) {
  ui->SceneGraphicsView->scene()->addItem(polygon);
}

void MainWindow::setViewLayer(int layer) {
  if (0 < layer && layer <= m_maxLayer) {
    m_viewLayer = layer;
    for (int i = 0; i < m_polygonList->rowCount(); i++) {
      if (m_polygonList->getPolygonItem(i)->getLayer() == m_viewLayer) {
        m_polygonList->getPolygonItem(i)->setVisible(true);
      }
      else {
        m_polygonList->getPolygonItem(i)->setVisible(false);
      }
    }
    setTextToAllVisiblePolygons();
  }
  m_layerLine->setText(QString::number(m_viewLayer));
}
void MainWindow::setMaxLayerText(int maxLayer)
{
  m_maxLayerText->setText("Max: " + QString::number(maxLayer) + " ");
}

void MainWindow::setStayOnMax(bool shouldStay)
{
  m_stayOnMax = shouldStay;
  m_stayOnMaxCheckBox->setChecked(shouldStay);
  if (m_stayOnMax) {
    setViewLayer(m_maxLayer);
  }
}

int MainWindow::getNumberVertices(PolygonWithHoles* polygon) {
  int vertices = 0;
  vertices += polygon->outer_boundary().size();
  for (auto iHole = polygon->holes_begin(); iHole != polygon->holes_end(); ++iHole) {
    vertices += iHole->size();
  }
  return vertices;
}
double MainWindow::getArea(PolygonWithHoles* polygon) {
  Kernel::FT area = 0;
  area += polygon->outer_boundary().area();
  for (auto iHole = polygon->holes_begin(); iHole != polygon->holes_end(); ++iHole) {
    area += iHole->area();
  }
  return CGAL::to_double(area);
}
