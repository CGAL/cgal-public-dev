// Copyright (c) 2005-2007 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id:  $
// 
//
// Author(s): Ophir Setter          <ophir.setter@post.tau.ac.il>
//

/*!
  \file   segment_voronoi_diagram.C
  \brief  A demo file for the linear Voronoi diagram.
  \todo   Move this into a more serious demo/example.
*/


#include <boost/program_options.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/bool.hpp>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h> // traits for Arr_2
#include <CGAL/Algebraic_curve_kernel_2_generator.h>

#include <CGAL/envelope_3.h>
#include <CGAL/Envelope_voronoi_traits_2/Linear_objects_Voronoi_traits_2.h>

#include <CGAL/Arr_linear_traits_2.h>

#include <CGAL/Timer.h>

#if !defined(MWA_NO_UI)
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_Curve_renderer_2.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Fig_stream_Curve_renderer_2.h>

#include <qapplication.h>
#include <qmainwindow.h>

#endif

namespace po = boost::program_options;

typedef CGAL::Arithmetic_kernel                       AK;
typedef AK::Integer                                   Integer;
typedef AK::Rational                                  Rational;
  


// Definition of Algebraic_kernel_2 (Algebraic_curve_kernel_2)
typedef Integer                                       Coefficient;
typedef CGAL::Algebraic_curve_kernel_2_generator<Coefficient>::
Algebraic_curve_kernel_with_qir_and_bitstream_2       Algebraic_curve_kernel_2;
typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_curve_kernel_2 > 
                                                      Curved_kernel_2;

typedef Curved_kernel_2::Point_2                      Point_2;
typedef Curved_kernel_2::Curve_2                      Curve_2;
typedef Curved_kernel_2::X_monotone_curve_2           X_monotone_curve_2;

typedef CGAL::Linear_objects_Voronoi_traits_2< 
  Curved_kernel_2 >                                   Segment_traits_2;


typedef Segment_traits_2::Site_2::Kernel              Linear_kernel_2;
typedef CGAL::Arr_linear_traits_2<Linear_kernel_2>    RTraits;
typedef CGAL::Arrangement_2<RTraits>                  Arrangement_2;
typedef CGAL::Envelope_diagram_2<Segment_traits_2>    Segment_diagram_2;

Segment_diagram_2                 segment_diagram;
Arrangement_2                     sites;

#if !defined(MWA_NO_UI)
typedef Linear_kernel_2           Fig_kernel;
CGAL::Fig_stream<Fig_kernel>      fig_stream;
std::string                       fig_file;

int W, H;


/*****************************************************************/
// Output operators section
/*****************************************************************/
// template <class Kernel>
// CGAL::Fig_stream<Kernel>& operator<< (CGAL::Fig_stream<Kernel>& fig_str, 
//                                       RTraits::Point_2& point)
// {
//   fig_str.write_point(point);
//   return fig_str;
// }

// template <class Kernel>
// CGAL::Fig_stream<Kernel>& operator<< (CGAL::Fig_stream<Kernel>& fig_str, 
//                                       RTraits::Segment_2& point)
// {
//   fig_str.write_segment(point);
//   return fig_str;
// }

// template <class Kernel>
// CGAL::Fig_stream<Kernel>& operator<< (CGAL::Fig_stream<Kernel>& fig_str, 
//                                       RTraits::Ray_2& point)
// {
//   fig_str.write_ray(point);
//   return fig_str;
// }

// template <class Kernel>
// CGAL::Fig_stream<Kernel>& operator<< (CGAL::Fig_stream<Kernel>& fig_str, 
//                                       RTraits::Line_2& point)
// {
//   fig_str.write_line(point);
//   return fig_str;
// }

// template <class Kernel>
// CGAL::Fig_stream<Kernel>& operator<< (CGAL::Fig_stream<Kernel>& fig_str, 
//                                       const char* c)
// {
//   return fig_str;
// }

CGAL::Qt_widget& operator<< (CGAL::Qt_widget& os,
                             const RTraits::X_monotone_curve_2& lobj)
{
  // Print a letter identifying the object type, then the object itself.
  if (lobj.is_segment())
    os << lobj.segment();
  else if (lobj.is_ray())
    os << lobj.ray();
  else
    os << lobj.line();

  return (os);
}

CGAL::Fig_stream<Fig_kernel>& operator<< (CGAL::Fig_stream<Fig_kernel>& os,
                             const RTraits::X_monotone_curve_2& lobj)
{
  // Print a letter identifying the object type, then the object itself.
  if (lobj.is_segment())
    os << lobj.segment();
  else if (lobj.is_ray())
    os << lobj.ray();
  else
    os << lobj.line();

  return (os);
}


// CGAL::Qt_widget& operator<< (CGAL::Qt_widget& os,
//                              const RTraits::Point_2& lobj)
// {
//   // Print a letter identifying the object type, then the object itself.
//   os << lobj;

//   return (os);
// }


/*****************************************************************/
// Qt widget & layers
/*****************************************************************/
template <class Arr>
void draw_arr(CGAL::Qt_widget &widget, const Arr& arr, 
              CGAL::Color c)
{
  widget << c;
  typename Arr::Halfedge_const_iterator eit;
  for (eit = arr.halfedges_begin(); eit != arr.halfedges_end(); ++eit)
  {
    widget << eit->curve();
  }
    
  typename Arr::Vertex_const_iterator vit;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    widget << vit->point();
  }
}


template <class Arr>
void draw_arr(CGAL::Qt_widget &widget, const Arr& arr, 
              CGAL::Fig_color fc, CGAL::Fig_line_style l_style)
{
  if (fig_stream.is_open())
    fig_stream.set_color(fc);
  
  fig_stream << l_style;
  
  typename Arr::Halfedge_const_iterator eit;
  for (eit = arr.halfedges_begin(); eit != arr.halfedges_end(); ++eit)
  {
    if (fig_stream.is_open())
      fig_stream << eit->curve();
  }
    
  typename Arr::Vertex_const_iterator vit;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    if (fig_stream.is_open())
      fig_stream << vit->point();
  }
}


template <class T, bool draw_fig = true>
class Arrangement_layer : public CGAL::Qt_widget_layer
{
public:
  typedef boost::mpl::bool_<draw_fig>               Draw_fig_tag;

  Arrangement_layer(T *t, 
                    CGAL::Color color = CGAL::ORANGE,
                    CGAL::Fig_color fig_color = CGAL::FIG_GREEN_2,
                    CGAL::Fig_line_style line_style = CGAL::FIG_SOLID)
    : m_t(t), m_color(color), m_fig_color(fig_color), m_line_style(line_style)
    {}

  void internal_draw(boost::mpl::bool_<true>)
  {
    draw_arr(*widget, *m_t, m_color);
    draw_arr(*widget, *m_t, m_fig_color, m_line_style);
  }

  void internal_draw(boost::mpl::bool_<false>)
  {
    draw_arr(*widget, *m_t, m_color);
  }


  void draw()
  {
    internal_draw(Draw_fig_tag());
  }
protected:
  T *m_t;
  CGAL::Color m_color;
  CGAL::Fig_color m_fig_color;
  CGAL::Fig_line_style m_line_style;
};

/*****************************************************************/
// Window 
/*****************************************************************/

class My_window : public QMainWindow
{
  Q_OBJECT
public:
  My_window(int x, int y)
    : sites_layer(&sites, CGAL::RED, CGAL::FIG_RED, CGAL::FIG_DASHED),
      vd_layer(&segment_diagram, CGAL::VIOLET, CGAL::FIG_BLUE)
  {
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);
    resize(W, H);
    widget->show();
    widget->set_window(-x, x, -y, y);

    //How to attach the standard toolbar
    std_toolbar = new CGAL::Qt_widget_standard_toolbar(widget, this,
						       "Standard Toolbar");
    widget->attach(&sites_layer);
    widget->attach(&vd_layer);

//     connect(widget, SIGNAL(redraw_on_back()),
// 	    this, SLOT(redraw_win()) );
  }

  void keyPressEvent(QKeyEvent *k) 
    {
      QMainWindow::keyPressEvent(k);
      
      if (k->key() == Qt::Key_S)
        sites_layer.toggle(!sites_layer.is_active());

      if (k->key() == Qt::Key_D)
        vd_layer.toggle(!vd_layer.is_active());

      if (k->key() == Qt::Key_P)
      {
        std::cout << "x_min " << widget->x_min() << std::endl;
        std::cout << "x_max " << widget->x_max() << std::endl;
        std::cout << "y_min " << widget->y_min() << std::endl;
        std::cout << "y_max " << widget->y_max() << std::endl;
        std::cout << "width " << widget->width() << std::endl;
        std::cout << "height " << widget->height() << std::endl;

        
        Fig_kernel::Iso_rectangle_2 irect (widget->x_min(), widget->y_min(), 
                                           widget->x_max(), widget->y_max());
        if (fig_file.size() > 0)
        {
          // the resolution of xfig is 1200ppi
          fig_stream.open (fig_file.c_str(), irect,
                           5*1200, 5*1200);
//          fig_stream.set_line_width(5);
          fig_stream.set_line_width(3);
          fig_stream.set_point_size(Rational(widget->x_max() - widget->x_min()) / 75);
        }
      }
      
      widget->redraw();

      if (k->key() == Qt::Key_P)
        fig_stream.close();
    }

private slots:	//functions
  void redraw_win()
  {
//       widget->lock();
//       widget->unlock();
  }

private:	//members
  CGAL::Qt_widget *widget;
  CGAL::Qt_widget_standard_toolbar *std_toolbar;
  Arrangement_layer<Arrangement_2> sites_layer;
  Arrangement_layer<Segment_diagram_2> vd_layer;
};

CGAL::Qt_widget&
operator <<(CGAL::Qt_widget & widget, const Point_2& p)
{
  widget.get_painter().drawPoint(widget.x_pixel(CGAL::to_double(p.xy().x())), 
                                 widget.x_pixel(CGAL::to_double(p.xy().y())));
  return widget;
}

// moc_source_file: min_width_annulus.C
#include "segment_voronoi_diagram.moc"

#endif

/*****************************************************************/
// Main 
/*****************************************************************/

int main( int argc, char **argv )
{

  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
#if !defined(MWA_NO_UI)
    ("show,S", "Opens QT window to watch the results.")
    ("width,W", po::value<int>()->default_value(500), 
     "Width of window")
    ("height,H", po::value<int>()->default_value(500), 
     "Height of window")
#endif
    ("fig_file_name,F", po::value<std::string>()->default_value(""), 
     "XFig filename")
    ;
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    
  
  if (vm.count("help")) 
  {
    std::cout << desc << "\n";
    return 1;
  }

#if !defined(MWA_NO_UI)
  W = vm["width"].as<int>();
  H = vm["height"].as<int>();
  fig_file = vm["fig_file_name"].as<std::string>();
#endif

  typedef Segment_traits_2::Surface_3           Surface_3;

  std::list<Surface_3>                          segments;

  int n;
  std::cin >> n;
  for (int i = 0; i < n; ++i)
  {  
    char type;
    std::cin >> type;

    Rational a, b, c, d;

    if (type == 'P' || type == 'p')
    {
      std::cin >> a >> b;
      segments.push_back(Surface_3(Linear_kernel_2::Point_2(a, b)));
    }
    else
    {
      std::cin >> a >> b >> c >> d;
      std::cout << "(" << a << ", " << b << "), (" << c << ", " << d << ")" << 
        std::endl;
      
      Linear_kernel_2::Point_2 s (a, b);
      Linear_kernel_2::Point_2 t (c, d);
      
      if (type == 'S' || type == 's')
        segments.push_back(Surface_3(Linear_kernel_2::Segment_2(s, t)));
      else if (type == 'R' || type == 'r')
        segments.push_back(Surface_3(Linear_kernel_2::Ray_2(s, t)));
      else if (type == 'L' || type == 'l')
        segments.push_back(Surface_3(Linear_kernel_2::Line_2(s, t)));
      else
        CGAL_error_msg("The curve type is not supported.");
    }
  }

  CGAL::Timer all_timer;
  all_timer.start();

  CGAL::lower_envelope_3 (segments.begin(), segments.end(),
                          segment_diagram);
  std::cout << "Segments VD:" << std::endl <<
    "T = " << all_timer.time() <<
    "\tV = " << segment_diagram.number_of_vertices() <<
    "\tE = " << segment_diagram.number_of_edges() <<
    "\tF = " << segment_diagram.number_of_faces() << std::endl;
  
  if (vm.count("show") == false)
    return 0;

  
  for (std::list<Surface_3>::iterator it = segments.begin(); 
       it != segments.end(); ++it)
  {
    if (it->is_point())
      insert_point (sites, RTraits::Point_2(it->point()));
    else if (it->is_segment())
      insert (sites, RTraits::Curve_2(it->segment()));
    else if (it->is_ray())
      insert (sites, RTraits::Curve_2(it->ray()));
    else
      insert (sites, RTraits::Curve_2(it->line()));
  }
          
#if !defined(MWA_NO_UI)
  // *************** UI SECTION ***************
  QApplication app( argc, argv );
  My_window W(5,5);
  app.setMainWidget( &W );
  W.show();
  W.setCaption("Using the Standard Toolbar");
#endif

  return app.exec();

}
