#include <QObject>


#include "Scene_swept_volume_3_item.h"
#include "Scene_swept_volume_3_item_config.h"


Scene_swept_volume_3_item::Scene_swept_volume_3_item()
  : Scene_item_with_display_list(),
    m_swept_volume_3(new Swept_volume_3)
{
}

Scene_swept_volume_3_item::Scene_swept_volume_3_item(Swept_volume_3* const p)
  : Scene_item_with_display_list(),
    m_swept_volume_3(p)
{
}

Scene_swept_volume_3_item::Scene_swept_volume_3_item(const Swept_volume_3& p)
  : Scene_item_with_display_list(),
    m_swept_volume_3(new Swept_volume_3(p))
{
}

// Scene_swept_volume_3_item::Scene_swept_volume_3_item(const Scene_swept_volume_3_item& item)
//   : Scene_item_with_display_list(item),
//     m_swept_volume_3(new Swept_volume_3(*item.m_swept_volume_3))
// {
// }

Scene_swept_volume_3_item::~Scene_swept_volume_3_item()
{
  delete m_swept_volume_3;
}

Scene_swept_volume_3_item* 
Scene_swept_volume_3_item::clone() const {
  return new Scene_swept_volume_3_item(*m_swept_volume_3);
}

// Load swept_volume_3 from .SV3 file
bool
Scene_swept_volume_3_item::load(std::istream& in)
{
  in >> *m_swept_volume_3;
  return in && !isEmpty();
}

// Write swept_volume_3 to .SV3 file
bool 
Scene_swept_volume_3_item::save(std::ostream& out) const
{
  out << *m_swept_volume_3;
  return out;
}

QString 
Scene_swept_volume_3_item::toolTip() const
{
  if(!m_swept_volume_3)
    return QString();

  return QObject::tr("<p>Swept_volume_3 <b>%1</b> </p>")
    .arg(this->name());
}

// So far this just draws the BBox 
void Scene_swept_volume_3_item::direct_draw() const {
  const Bbox& b = bbox();
  
  ::glDisable(GL_LIGHTING);
  ::glColor3f(0.f,0.f,0.f);
  ::glBegin(GL_LINES);
  
  ::glVertex3d(b.xmin,b.ymin,b.zmin);
  ::glVertex3d(b.xmin,b.ymin,b.zmax);
  
  ::glVertex3d(b.xmin,b.ymin,b.zmin);
  ::glVertex3d(b.xmin,b.ymax,b.zmin);
  
  ::glVertex3d(b.xmin,b.ymin,b.zmin);
  ::glVertex3d(b.xmax,b.ymin,b.zmin);
  
  ::glVertex3d(b.xmax,b.ymin,b.zmin);
  ::glVertex3d(b.xmax,b.ymax,b.zmin);
  
  ::glVertex3d(b.xmax,b.ymin,b.zmin);
  ::glVertex3d(b.xmax,b.ymin,b.zmax);
  
  ::glVertex3d(b.xmin,b.ymax,b.zmin);
  ::glVertex3d(b.xmin,b.ymax,b.zmax);
  
  ::glVertex3d(b.xmin,b.ymax,b.zmin);
  ::glVertex3d(b.xmax,b.ymax,b.zmin);
  
  ::glVertex3d(b.xmax,b.ymax,b.zmin);
  ::glVertex3d(b.xmax,b.ymax,b.zmax);
  
  ::glVertex3d(b.xmin,b.ymin,b.zmax);
  ::glVertex3d(b.xmin,b.ymax,b.zmax);
  
  ::glVertex3d(b.xmin,b.ymin,b.zmax);
  ::glVertex3d(b.xmax,b.ymin,b.zmax);
  
  ::glVertex3d(b.xmax,b.ymax,b.zmax);
  ::glVertex3d(b.xmin,b.ymax,b.zmax);
  
  ::glVertex3d(b.xmax,b.ymax,b.zmax);
  ::glVertex3d(b.xmax,b.ymin,b.zmax);
  
  ::glEnd();  
}

Swept_volume_3* 
Scene_swept_volume_3_item::swept_volume_3()       { return m_swept_volume_3; }
const Swept_volume_3* 
Scene_swept_volume_3_item::swept_volume_3() const { return m_swept_volume_3; }

bool
Scene_swept_volume_3_item::isEmpty() const {
  return (m_swept_volume_3 == 0);
}

Scene_swept_volume_3_item::Bbox
Scene_swept_volume_3_item::bbox() const {
  return Bbox(-1.0,-1.0,-1.0,1.0,1.0,1.0);
}

#include "Scene_swept_volume_3_item.moc"
