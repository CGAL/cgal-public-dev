#include "config.h"

#include "Scene_implicit_function_item.h"
#include <QColor>
#include <map>
#include <CGAL/gl.h>
#include <CGAL/Simple_cartesian.h>

#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/qglviewer.h>

#include "Color_ramp.h"


Scene_implicit_function_item::
Scene_implicit_function_item(Implicit_function_interface* f)
  : function_(f)
  , frame_(new ManipulatedFrame())
  , grid_size_(SCENE_IMPLICIT_GRID_SIZE)
  , max_value_(0.)
  , min_value_(0.)
  , blue_color_ramp_()
  , red_color_ramp_()
{
  blue_color_ramp_.build_blue();
  red_color_ramp_.build_red();
  compute_min_max();
  compute_function_grid();
  
  connect(frame_, SIGNAL(modified()), this, SLOT(compute_function_grid()));
}


Scene_implicit_function_item::~Scene_implicit_function_item()
{
  delete frame_;
}


Scene_implicit_function_item::Bbox
Scene_implicit_function_item::bbox() const
{
  return function_->bbox();
}

void
Scene_implicit_function_item::direct_draw() const
{
  
  draw_function_grid(red_color_ramp_, blue_color_ramp_);
  draw_bbox();
}



QString
Scene_implicit_function_item::toolTip() const
{
  return tr("<p>Function <b>%1</b>")
    .arg(this->name());
}

bool
Scene_implicit_function_item::supportsRenderingMode(RenderingMode m) const
{ 
  switch ( m )
  {
    case Gouraud:
      return false;
      
    case Points:
    case Wireframe:
    case Flat:
    case FlatPlusEdges:
      return true;
      
    default:
      return false;
  }
  
  return false;
}

void
Scene_implicit_function_item::
draw_bbox() const
{
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

void 
Scene_implicit_function_item::
draw_function_grid(const Color_ramp& ramp_pos,
                   const Color_ramp& ramp_neg) const
{
  ::glDisable(GL_LIGHTING);
  ::glShadeModel(GL_SMOOTH);
  
  ::glBegin(GL_QUADS);
  const int nb_quads = grid_size_ - 1;
  for( int i=0 ; i < nb_quads ; i++ )
  {
    for( int j=0 ; j < nb_quads ; j++)
    {
      draw_grid_vertex(implicit_grid_[i][j], ramp_pos, ramp_neg);
      draw_grid_vertex(implicit_grid_[i][j+1], ramp_pos, ramp_neg);
      draw_grid_vertex(implicit_grid_[i+1][j+1], ramp_pos, ramp_neg);
      draw_grid_vertex(implicit_grid_[i+1][j], ramp_pos, ramp_neg);
    }
  }
  ::glEnd();
}


void
Scene_implicit_function_item::
draw_grid_vertex(const Point_value& pv,
                 const Color_ramp& ramp_positive,
                 const Color_ramp& ramp_negative) const
{
  const Point& p = pv.first;
  double v = pv.second;
  
  // determines grey level
  if ( v > 0 )
  {
    v = v/max_value_;
    ::glColor3d(ramp_positive.r(v),ramp_positive.g(v),ramp_positive.b(v));
  }
  else
  {
    v = v/min_value_;
    ::glColor3d(ramp_negative.r(v),ramp_negative.g(v),ramp_negative.b(v));
  }
  
  ::glVertex3d(p.x,p.y,p.z);
}


void
Scene_implicit_function_item::
compute_function_grid()
{
  typedef CGAL::Simple_cartesian<double>  K;
  typedef K::Aff_transformation_3         Aff_transformation;
  typedef K::Point_3                      Point_3;
  
  // Get transformation
  const ::GLdouble* m = frame_->matrix();
  
  // OpenGL matrices are row-major matrices
  Aff_transformation t (m[0], m[4], m[8], m[12],
                        m[1], m[5], m[9], m[13],
                        m[2], m[6], m[10], m[14]);
  
  double diag = bbox().diagonal_length() * .6;
  
  const double dx = diag;
  const double dy = diag;
  const double z (0);
  
  for(int i=0 ; i<grid_size_ ; ++i)
  {
    double x = -diag/2. + double(i)/double(grid_size_) * dx;
    
    for(int j=0 ; j<grid_size_ ; ++j)
    {
      double y = -diag/2. + double(j)/double(grid_size_) * dy;
      
      Point_3 query = t( Point_3(x,y,z) );
      double v = function_->operator()(query.x(), query.y(), query.z());
      
      implicit_grid_[i][j] = Point_value(Point(query.x(),query.y(),query.z()),v);
    }
  }
  
  // Update display list
  this->changed();
}

void
Scene_implicit_function_item::
compute_min_max()
{
  max_value_ = 0;
  min_value_ = 0;
  
  double probes_nb = double(grid_size_) / 2;
  
  // Probe bounding box
  const Bbox& b = bbox();
  
  for ( int i = 0 ; i <= probes_nb ; ++i )
  {
    double x = b.xmin + double(i) * (b.xmax - b.xmin) / probes_nb;
    
    for ( int j = 0 ; j <= probes_nb ; ++j )
    {
      double y = b.ymin + double(j) * (b.ymax - b.ymin) / probes_nb;
      
      for ( int k = 0 ; k <= probes_nb ; ++k )
      {
        double z = b.zmin + double(k) * (b.zmax - b.zmin) / probes_nb;
        
        double v = (*function_)(x,y,z);
        max_value_ = (std::max)(v,max_value_);
        min_value_ = (std::min)(v,min_value_);
      }
    }
  }
}


#include "Scene_implicit_function_item.moc"

