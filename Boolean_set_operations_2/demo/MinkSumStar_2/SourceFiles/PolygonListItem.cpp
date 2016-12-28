#include "PolygonListItem.h"

PolygonListItem::PolygonListItem(QString name) {
  m_name = name;
  m_color = *new QColor();
  m_visibility = true;
  m_polygonItem = NULL;
  m_layer = 1;
}
PolygonListItem::PolygonListItem(QString name, QColor color, bool visibility, 
                                 PolygonItem* polygonItem) {
  m_name = name;
  m_color = color;
  m_visibility = visibility;
  m_polygonItem = polygonItem;
  m_layer = 1;
}

QString PolygonListItem::getName() {
  return m_name;
}
QColor PolygonListItem::getColor() { 
  return m_color;
}
bool PolygonListItem::isVisible() {
  return m_visibility;
}
PolygonItem* PolygonListItem::getPolygonItem() {
  return m_polygonItem;
}
int PolygonListItem::getLayer() {
  return m_layer;
}

void PolygonListItem::setName(QString name) {
  m_name = name;
}
void PolygonListItem::setColor(QColor color) {
  m_color = color; 
}
void PolygonListItem::setVisible(bool visibility) {
  if (m_visibility != visibility) {
    m_visibility = visibility;
    if (visibility == true) {
      emit(showPolygon(m_polygonItem->graphics()));
    }
    else {
      emit(hidePolygon(m_polygonItem->graphics()));
    }
  }
}
void PolygonListItem::setPolygonItem(PolygonItem* polygonItem) {
  m_polygonItem = polygonItem;
}
void PolygonListItem::setLayer(int layer) {
  m_layer = layer;
}
