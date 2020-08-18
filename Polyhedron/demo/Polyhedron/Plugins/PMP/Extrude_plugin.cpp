#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <QApplication>
#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QInputDialog>
#include <QMessageBox>

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Viewer_interface.h>

#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/constraint.h>

#include <CGAL/number_type_config.h>

#include "Messages_interface.h"
#ifdef USE_SURFACE_MESH
#include "Kernel_type.h"
#include "Scene_surface_mesh_item.h"
#else
#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"
#endif
#include "Scene_polyhedron_selection_item.h"
#include "Scene.h"
#ifdef USE_SURFACE_MESH
typedef Scene_surface_mesh_item Scene_face_graph_item;
#else
typedef Scene_polyhedron_item Scene_face_graph_item;
#endif

typedef Scene_face_graph_item::Face_graph Face_graph;
typedef CGAL::qglviewer::Vec Vec;
using namespace CGAL::Three;
//use frame to get dist and dir.
//fix frame in translation and use mousewheel to choose dist
//finding frame's position : first try at the center of the item's bbox
//maybe find intersection between surface and diag bbox.
//orientation : PCA : normal.
typedef Kernel::Plane_3 Plane;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
class Scene_arrow_item : public Scene_item
{
  Q_OBJECT
public :
  Scene_arrow_item(Vec center_, double r, double length_)
    :Scene_item(2, 1), center_(center_), length_(length_), R(r),
       frame(new Scene_item::ManipulatedFrame())
  {
    const CGAL::qglviewer::Vec offset = static_cast<Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
        frame->setPosition( center_+offset);
    nb_pos = 0;
    tick = length_/10.0f;
    ctrl_pressing = false;
  }
  // Indicates if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE {
    return (m == Gouraud);
  }
  //Displays the item
  void draw(Viewer_interface* viewer) const Q_DECL_OVERRIDE
  {
    if(!are_buffers_filled)
    {
      initializeBuffers(viewer);
    }
    vaos[0]->bind();
    program = getShaderProgram(PROGRAM_WITH_LIGHT);
    GLdouble d_mat[16];
    QMatrix4x4 mvp_mat;
    GLdouble matrix[16];
    QMatrix4x4 f_matrix;
    frame->getMatrix(matrix);
    for (int i=0; i<16; ++i)
      f_matrix.data()[i] = (float)matrix[i];
    viewer->camera()->getModelViewProjectionMatrix(d_mat);
    for (int i=0; i<16; ++i)
      mvp_mat.data()[i] = GLfloat(d_mat[i]);
    mvp_mat = mvp_mat*f_matrix;
    QMatrix4x4 mv_mat;
    viewer->camera()->getModelViewMatrix(d_mat);
    for (int i=0; i<16; ++i)
      mv_mat.data()[i] = GLfloat(d_mat[i]);
    mv_mat = mv_mat*f_matrix;
    attribBuffers(viewer, PROGRAM_WITH_LIGHT);
    program->bind();
    program->setUniformValue("mvp_matrix", mvp_mat);
    program->setUniformValue("mv_matrix", mv_mat);
    program->setUniformValue("is_clipbox_on", false);
    program->setAttributeValue("radius",0.01*diagonalBbox());
    program->setAttributeValue("colors", QColor(Qt::red));
    program->setAttributeValue("colors", this->color());
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(nb_pos/3));
    vaos[0]->release();
    program->release();
  }
  void invalidateOpenGLBuffers() Q_DECL_OVERRIDE
  {
    are_buffers_filled = false;
  }
  Scene_item* clone() const Q_DECL_OVERRIDE {return 0;}
  QString toolTip() const Q_DECL_OVERRIDE {return QString();}
  Vec center()const { return center_; }
  Scene_item::ManipulatedFrame* manipulatedFrame() Q_DECL_OVERRIDE { return frame; }
  bool manipulatable() const Q_DECL_OVERRIDE { return true; }
  bool eventFilter(QObject *, QEvent *event) Q_DECL_OVERRIDE
  {
    if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease)
    {
      ctrl_pressing = static_cast<QKeyEvent*>(event)->modifiers().testFlag(Qt::ControlModifier);
    }
    if(event->type() == QEvent::Wheel && ctrl_pressing)
    {
      QWheelEvent *mouseEvent = static_cast<QWheelEvent*>(event);
      int steps = mouseEvent->delta() / 120;
      if (steps > 0)
        length_+=tick;
      else
        length_-=tick;
      invalidateOpenGLBuffers();
      redraw();
      return true;
    }
    return false;
  }
  double length()const { return length_; }
private:
  //make an arrow showing the length and direction of the transformation for the extrusion. 
  void initializeBuffers(Viewer_interface *viewer)const
  {
    std::vector<float> vertices;
    std::vector<float> normals;
    int prec = 60;
    //Head
    const float Rf = static_cast<float>(R);
    for(int d = 0; d<360; d+= 360/prec)
    {
      float D = (float) (d * CGAL_PI / 180.);
      float a = (float) std::atan(Rf / 0.33);
      QVector4D pR(0., 1.*length_, 0, 1.);
      QVector4D nR(Rf*2.*sin(D), sin(a), Rf*2.*cos(D), 1.);
      
      //point A1
      vertices.push_back(pR.x());
      vertices.push_back(pR.y());
      vertices.push_back(pR.z());
      normals.push_back(nR.x());
      normals.push_back(nR.y());
      normals.push_back(nR.z());
      
      //point B1
      pR = QVector4D(Rf*2.*sin(D), 0.66f*length_, Rf*2.* cos(D), 1.f);
      nR = QVector4D(sin(D), sin(a), cos(D), 1.);
      vertices.push_back(pR.x());
      vertices.push_back(pR.y());
      vertices.push_back(pR.z());
      normals.push_back(nR.x());
      normals.push_back(nR.y());
      normals.push_back(nR.z());
      //point C1
      D = (d+360/prec)*CGAL_PI/180.0;
      pR = QVector4D(Rf*2.* sin(D), 0.66f*length_, Rf *2.* cos(D), 1.f);
      nR = QVector4D(sin(D), sin(a), cos(D), 1.0);
      
      vertices.push_back(pR.x());
      vertices.push_back(pR.y());
      vertices.push_back(pR.z());
      normals.push_back(nR.x());
      normals.push_back(nR.y());
      normals.push_back(nR.z());
    }
    
    //cylinder
    //body of the cylinder
    for(int d = 0; d<360; d+= 360/prec)
    {
      //point A1
      double D = d*CGAL_PI/180.0;
      QVector4D pR(Rf*sin(D), 0.66f*length_, Rf*cos(D), 1.f);
      QVector4D nR(sin(D), 0.f, cos(D), 1.f);
      
      vertices.push_back(pR.x());
      vertices.push_back(pR.y());
      vertices.push_back(pR.z());
      normals.push_back(nR.x());
      normals.push_back(nR.y());
      normals.push_back(nR.z());
      //point B1
      pR = QVector4D(Rf * sin(D),0,Rf*cos(D), 1.0);
      nR = QVector4D(sin(D), 0, cos(D), 1.0);
      
      vertices.push_back(pR.x());
      vertices.push_back(pR.y());
      vertices.push_back(pR.z());
      normals.push_back(nR.x());
      normals.push_back(nR.y());
      normals.push_back(nR.z());
      //point C1
      D = (d+360/prec)*CGAL_PI/180.0;
      pR = QVector4D(Rf * sin(D),0,Rf*cos(D), 1.0);
      nR = QVector4D(sin(D), 0, cos(D), 1.0);
      vertices.push_back(pR.x());
      vertices.push_back(pR.y());
      vertices.push_back(pR.z());
      normals.push_back(nR.x());
      normals.push_back(nR.y());
      normals.push_back(nR.z());
      //point A2
      D = (d+360/prec)*CGAL_PI/180.0;
      
      pR = QVector4D(Rf * sin(D),0,Rf*cos(D), 1.0);
      nR = QVector4D(sin(D), 0, cos(D), 1.0);

      vertices.push_back(pR.x());
      vertices.push_back(pR.y());
      vertices.push_back(pR.z());
      normals.push_back(nR.x());
      normals.push_back(nR.y());
      normals.push_back(nR.z());
      //point B2
      pR = QVector4D(Rf * sin(D), 0.66f*length_, Rf*cos(D), 1.f);
      nR = QVector4D(sin(D), 0, cos(D), 1.0);
      
      vertices.push_back(pR.x());
      vertices.push_back(pR.y());
      vertices.push_back(pR.z());
      normals.push_back(nR.x());
      normals.push_back(nR.y());
      normals.push_back(nR.z());
      //point C2
      D = d*CGAL_PI/180.0;
      pR = QVector4D(Rf * sin(D), 0.66f*length_, Rf*cos(D), 1.f);
      nR = QVector4D(sin(D), 0.f, cos(D), 1.f);
      
      vertices.push_back(pR.x());
      vertices.push_back(pR.y());
      vertices.push_back(pR.z());
      normals.push_back(nR.x());
      normals.push_back(nR.y());
      normals.push_back(nR.z());
    }
    
    //fill buffers
    //vao containing the data for the facets
    program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
    program->bind();
    vaos[0]->bind();
    buffers[0].bind();
    buffers[0].allocate(vertices.data(),
                        static_cast<GLsizei>(vertices.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[0].release();
    buffers[1].bind();
    buffers[1].allocate(normals.data(),
                        static_cast<GLsizei>(normals.size()*sizeof(float)));
    program->enableAttributeArray("normals");
    program->setAttributeBuffer("normals",GL_FLOAT,0,3);
    buffers[1].release();
    vaos[0]->release();
    program->release();
    //once the buffers are filled, we can empty the vectors to optimize memory consumption
    nb_pos = vertices.size();
    _bbox = Bbox(0,0,0,0,0,0);
    for(std::size_t i = 0; i< vertices.size(); i+=3)
    {
      _bbox += Point(vertices[i],
                     vertices[i+1],
                     vertices[i+2]).bbox();
    }
    
    are_buffers_filled = true;
  }
  
  mutable std::size_t nb_pos;
  Vec center_;
  double length_;
  double tick;
  double R;
  bool ctrl_pressing;
  mutable QOpenGLShaderProgram *program;
  Scene_item::ManipulatedFrame* frame;
}; //end of class Scene_arrow_item


template <typename TriangleMesh, typename OutputIterator>
CGAL::Bbox_3 triangles(const TriangleMesh& mesh,
                         OutputIterator out)
{
  CGAL::Bbox_3 bb;
  typename boost::property_map<TriangleMesh,CGAL::vertex_point_t>::const_type vpm =
      get(CGAL::vertex_point, mesh);
  BOOST_FOREACH(typename boost::graph_traits<TriangleMesh>::face_descriptor fd, faces(mesh)){
    typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd = halfedge(fd,mesh);
    Triangle t(get(vpm,source(hd,mesh)),
               get(vpm,target(hd,mesh)),
               get(vpm,target(next(hd,mesh),mesh)));
    *out++ = t;
    bb = bb + t.bbox();
  }
  return bb;
}

Vector estimate_normals(const std::vector<Triangle>& tris)
{
  Vector moy(0,0,0);
  
  BOOST_FOREACH(const Triangle& tri, tris)
  {
    Vector norm = CGAL::Polygon_mesh_processing::internal::triangle_normal(
          tri[0], tri[1], tri[2], Kernel());
    norm /= CGAL::sqrt(norm.squared_length());
    moy += norm;
  }
  return moy;
}


class ExtrudePlugin :
    public QObject,
    public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:
 
  bool applicable(QAction* action) const Q_DECL_OVERRIDE
  {
    if(action == actionCreateItem)
    {
      return !oliver_queen && 
          (qobject_cast<Scene_face_graph_item*>(scene->item(scene->mainSelectionIndex()))
           || qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex())));
    }
    else if(oliver_queen)
      return true;
    return false;
  }
 
  QList<QAction*> actions() const Q_DECL_OVERRIDE
  {
    return _actions;
  }

  void init(QMainWindow* mainWindow, Scene_interface* sc, Messages_interface* mi) Q_DECL_OVERRIDE
  {
    this->messageInterface = mi;
    this->scene = sc;
    this->mw = mainWindow;
    oliver_queen = NULL;
    target = NULL;
    actionCreateItem = new QAction(QString("Extrude Item"), mw);
    actionCreateItem->setProperty("submenuName", "Polygon Mesh Processing");
    connect(actionCreateItem, SIGNAL(triggered()),
            this, SLOT(createItem()));
    _actions << actionCreateItem;
    actionExtrude = new QAction(QString("Perform Extrusion"), mw);
    actionExtrude->setProperty("submenuName", "Polygon Mesh Processing");
    connect(actionExtrude, SIGNAL(triggered()),
            this, SLOT(do_extrude()));
    _actions << actionExtrude;
  }
private Q_SLOTS:
  void createItem()
  {
    Scene_item * item = scene->item(scene->mainSelectionIndex());
    Scene_polyhedron_selection_item* sel_item = qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
    Scene_face_graph_item* fg_item = qobject_cast<Scene_face_graph_item*>(item);
    Face_graph* pMesh = NULL;
    if(sel_item)
    {
      pMesh = new Face_graph();
      if(!sel_item->export_selected_facets_as_polyhedron(pMesh))
      {
        messageInterface->error("Face selection is not valid. Aborting.");
        
        return; 
      }
      fg_item = new Scene_facegraph_item(pMesh);
      fg_item->setName(QString("%1 selection").arg(sel_item->polyhedron_item()->name()));
      scene->addItem(fg_item);
      sel_item->polyhedron_item()->setWireframeMode();
      sel_item->polyhedron_item()->redraw();
    }
    if(fg_item)
      pMesh = fg_item->face_graph();
    else
      return;
    if(CGAL::is_closed(*pMesh))
    {
      messageInterface->error("The face graph must be open. Aborting.");
      return;
    }
    std::vector<Triangle> triangles;
    ::triangles(*pMesh,std::back_inserter(triangles));
    Plane plane;
    CGAL::linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<2>());
    
    // compute centroid
    Point c = CGAL::centroid(triangles.begin(),triangles.end());
    
    oliver_queen = new Scene_arrow_item(Vec(c.x(),c.y(),c.z()), fg_item->diagonalBbox() / 50.0f, 
                                        fg_item->diagonalBbox()/3.0f);
    Vec dir(plane.orthogonal_vector().x(),
            plane.orthogonal_vector().y(), 
            plane.orthogonal_vector().z());
    if(CGAL::scalar_product(Vector(dir.x, dir.y, dir.z), estimate_normals(triangles)) > 0)
      dir = -dir;
    
    CGAL::qglviewer::Quaternion orientation(CGAL::qglviewer::Vec(0,1,0), dir);
    oliver_queen->manipulatedFrame()->setOrientation(orientation);
    constraint.setRotationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FREE);
    constraint.setTranslationConstraintType(CGAL::qglviewer::AxisPlaneConstraint::FORBIDDEN);
    oliver_queen->manipulatedFrame()->setConstraint(&constraint);
    oliver_queen->setColor(QColor(Qt::green));
    oliver_queen->setName("Extrude item");
    CGAL::QGLViewer* viewer = *CGAL::QGLViewer::QGLViewerPool().begin();
    viewer->installEventFilter(oliver_queen);
    mw->installEventFilter(oliver_queen);
    scene->addItem(oliver_queen);
    target = fg_item;
    
    connect(oliver_queen, &Scene_arrow_item::aboutToBeDestroyed,
            [this](){
      oliver_queen = NULL;
    });
    
    //!@todo : add a way to track scene's bbox recomputation and reset frame's position when triggered.
  }
  void do_extrude()
  {
    if(!target)
      return;
    Face_graph pMesh = *target->face_graph();
    target->face_graph()->clear();
    double length = oliver_queen->length();
    double matrix[16];
    oliver_queen->manipulatedFrame()->getMatrix(matrix);
    QMatrix4x4 rotate_matrix;
    QMatrix4x4 transform_matrix;
    for(int i=0; i<16; ++i)
      transform_matrix.data()[i] = (float)matrix[i];
    rotate_matrix = transform_matrix;
    rotate_matrix.setColumn(3, QVector4D(0,0,0,1));
    QVector3D dir = rotate_matrix * QVector3D(0,1,0);
    dir.normalize();
    dir = length * dir;
    
    CGAL::Polygon_mesh_processing::extrude_mesh(pMesh, *target->face_graph(),
                                                Kernel::Vector_3(dir.x(), dir.y(), dir.z()));
    scene->erase(scene->item_id(oliver_queen));
    oliver_queen = NULL;
    target->invalidateOpenGLBuffers();
    target->itemChanged();
    target = NULL;
  }
  
private:
  
  QList<QAction*> _actions;
  Messages_interface* messageInterface;
  Scene_interface* scene;
  QMainWindow* mw;
  QAction *actionCreateItem;
  QAction *actionExtrude;
  Scene_arrow_item* oliver_queen;
  Scene_face_graph_item* target;
  CGAL::qglviewer::LocalConstraint constraint;
};
#include "Extrude_plugin.moc"
