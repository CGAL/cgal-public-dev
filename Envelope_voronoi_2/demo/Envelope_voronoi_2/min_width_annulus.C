#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_impl.h> // traits for Arr_2
#include <CGAL/Algebraic_curve_kernel_2_generator.h>

#include <CGAL/envelope_3.h>
#include <CGAL/Envelope_voronoi_traits_2/Algebraic_apollonius_traits_2.h>
#include <CGAL/Envelope_voronoi_traits_2/Algebraic_farthest_point_farthest_site_traits_2.h>

#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Min_width_annulus_2/MWA_overlay_traits_2.h>
#include <CGAL/Min_width_annulus_2/MWA_arrangement_2.h>

#include <CGAL/Timer.h>

#if !defined(MWA_NO_UI)
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_Curve_renderer_2.h>
#include <CGAL/IO/Qt_widget_layer.h>

#include <qapplication.h>
#include <qmainwindow.h>

#include <CGAL/IO/Fig_stream_Curve_renderer_2.h>
#endif

#include <boost/program_options.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

namespace po = boost::program_options;

typedef CGAL::Arithmetic_kernel                       AK;
typedef AK::Integer                                   Integer;
typedef AK::Rational                                  Rational;
typedef AK::Field_with_sqrt                           Field_with_sqrt;
typedef boost::numeric::interval<Field_with_sqrt>     Interval;
  


// Definition of Algebraic_kernel_2 (Algebraic_curve_kernel_2)
typedef Integer                                       Coefficient;
typedef CGAL::Algebraic_curve_kernel_2_generator<Coefficient>::
Algebraic_curve_kernel_with_qir_and_bitstream_2       Algebraic_curve_kernel_2;
typedef CGAL::Curved_kernel_via_analysis_2< Algebraic_curve_kernel_2 > 
                                                      Curved_kernel_2;

typedef Curved_kernel_2::Point_2                      Point_2;
typedef Curved_kernel_2::Curve_2                      Curve_2;
typedef Curved_kernel_2::X_monotone_curve_2           X_monotone_curve_2;

typedef CGAL::Algebraic_apollonius_traits_2
<Curved_kernel_2>                                     Nearest_traits_2;
typedef CGAL::Algebraic_farthest_point_farthest_site_traits_2
<Curved_kernel_2>                                     Furthest_traits_2;


typedef CGAL::Arrangement_2<Curved_kernel_2>          Arrangement_2;
typedef CGAL::Envelope_diagram_2<Nearest_traits_2>    Nearest_diagram_2;
typedef CGAL::Envelope_diagram_2<Furthest_traits_2>   Furthest_diagram_2;
// Xy_monotone_surface_3 of both Nearest and Furthest should be the same.
typedef Nearest_traits_2::Xy_monotone_surface_3       Xy_monotone_surface_3;

// We can use Furthest_traits_2. We only use it to get the 
// Xy_monotone_surface_3 from the traits.
typedef CGAL::MWA_arrangement_2< Nearest_traits_2 >   MWA_arrangement_2;

typedef CGAL::MWA_overlay_traits_2< Nearest_diagram_2,
                                    Furthest_diagram_2,
                                    MWA_arrangement_2 > Overlay_traits_2;

Nearest_diagram_2                 nearest_diagram;
Furthest_diagram_2                furthest_diagram;
Arrangement_2                     sites;
Arrangement_2                     inner_outer_circles_arr;

#if !defined(MWA_NO_UI)
typedef CGAL::Cartesian<double>   Fig_kernel;
CGAL::Fig_stream<Fig_kernel>      fig_stream;
std::string                       fig_file;

int W, H;

template <class Arr>
void draw_arr(CGAL::Qt_widget &widget, const Arr& arr, 
              CGAL::Color c, CGAL::Fig_color fc)
{
  if (fig_stream.is_open())
    fig_stream.set_color(fc);

  widget << c;
  typename Arr::Halfedge_const_iterator eit;
  for (eit = arr.halfedges_begin(); eit != arr.halfedges_end(); ++eit)
  {
    widget << eit->curve();
    if (fig_stream.is_open())
      fig_stream << eit->curve();
  }
    
  typename Arr::Vertex_const_iterator vit;
  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    widget << vit->point();
    if (fig_stream.is_open())
      fig_stream << vit->point();
  }
}

template <class T>
class Min_annulus_layer : public CGAL::Qt_widget_layer
{
public:
  Min_annulus_layer(T *t, 
                    CGAL::Color color = CGAL::ORANGE,
                    CGAL::Fig_color fig_color = CGAL::FIG_GREEN_2)
    : m_t(t), m_color(color), m_fig_color(fig_color) {}

  void draw()
    {
      draw_arr(*widget, *m_t, m_color, m_fig_color);
    }
protected:
  T *m_t;
  CGAL::Color m_color;
  CGAL::Fig_color m_fig_color;
};


class My_window : public QMainWindow
{
  Q_OBJECT
public:
  My_window(int x, int y)
    : annulus_layer(&inner_outer_circles_arr, CGAL::ORANGE, CGAL::FIG_GREEN_2),
      sites_layer(&sites, CGAL::VIOLET, CGAL::FIG_MAGENTA_2), 
      nearest_layer(&nearest_diagram, CGAL::RED, CGAL::FIG_RED),
      furthest_layer(&furthest_diagram, CGAL::BLUE, CGAL::FIG_BLUE)
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
    widget->attach(&nearest_layer);
    widget->attach(&furthest_layer);
    widget->attach(&annulus_layer);

    connect(widget, SIGNAL(redraw_on_back()),
	    this, SLOT(redraw_win()) );
  }

  void keyPressEvent(QKeyEvent *k) 
    {
      QMainWindow::keyPressEvent(k);
      if (k->key() == Qt::Key_A)
        annulus_layer.toggle(!annulus_layer.is_active());
      
      if (k->key() == Qt::Key_S)
        sites_layer.toggle(!sites_layer.is_active());

      if (k->key() == Qt::Key_N)
        nearest_layer.toggle(!nearest_layer.is_active());

      if (k->key() == Qt::Key_F)
        furthest_layer.toggle(!furthest_layer.is_active());

      if (k->key() == Qt::Key_P)
      {
        Fig_kernel::Iso_rectangle_2 irect (widget->x_min(), widget->y_min(), 
                                           widget->x_max(), widget->y_max());
        if (fig_file.size() > 0)
        {
          fig_stream.open (fig_file.c_str(), irect,
                           W, H);
//          fig_stream.set_line_width(5);
          fig_stream.set_line_width(7);
          fig_stream.set_point_size(1);
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
  Min_annulus_layer<Arrangement_2> annulus_layer;
  Min_annulus_layer<Arrangement_2> sites_layer;
  Min_annulus_layer<Nearest_diagram_2> nearest_layer;
  Min_annulus_layer<Furthest_diagram_2> furthest_layer;
};

CGAL::Qt_widget&
operator<<(CGAL::Qt_widget & widget, const Point_2& p)
{
  widget.get_painter().drawPoint(widget.x_pixel(CGAL::to_double(p.xy().x())), 
                                 widget.x_pixel(CGAL::to_double(p.xy().y())));
  return widget;
}

// moc_source_file: min_width_annulus.C
#include "min_width_annulus.moc"

#endif

//! Functor to check if a vertex is a candidate for annulus center.
//  The annulus should intersect the disks in 4 points.
struct Is_annulus_candidate
{
  bool operator() (const MWA_arrangement_2::Vertex &v) const
    {
      return ((v.first.number_of_surfaces() + 
              v.second.number_of_surfaces()) >= 4);
    }
};

//! Functor to compare two vertices of the overlaid arrangements.
struct Comp_annulus_width_at_vertex
{
  
  template <typename NT>
  static NT annulus_width(const NT &px,
                          const NT &py, 
                          const Xy_monotone_surface_3& near_surf,
                          const Xy_monotone_surface_3& fur_surf)
    {
      NT f = Furthest_traits_2::distance(px, py, fur_surf);
      NT n = Nearest_traits_2::distance(px, py, near_surf);

      return f - n;
        
    }
  
  template <typename Iterator>
  static void print_sites(Iterator begin, Iterator end)
    {
      while (begin != end)
        std::cout << *begin++ << std::endl;
    }

  bool operator() (const MWA_arrangement_2::Vertex_const_iterator &v1,
                   const MWA_arrangement_2::Vertex_const_iterator &v2) const
    {
//       CGAL::set_pretty_mode(std::cout);
//       std::cout  << v1->point() << " " << v2->point() << std::endl;

      // if the nearest and the furthest contain the same site then
      // the distance is the radius of that site.
      typedef std::vector< Xy_monotone_surface_3 > Sites;
      Sites fs1(v1->first.surfaces_begin(), v1->first.surfaces_end());
      Sites ss1(v1->second.surfaces_begin(), v1->second.surfaces_end());
      std::sort (fs1.begin(), fs1.end());
      std::sort (ss1.begin(), ss1.end());
      Sites s1;
      std::set_intersection(fs1.begin(), fs1.end(), ss1.begin(), ss1.end(),
                            std::back_inserter(s1));
//       std::cout << "s1 : " << std::endl;
//       print_sites(s1.begin(), s1.end());
      if (s1.empty() == false)
      {
        Sites fs2(v2->first.surfaces_begin(), v2->first.surfaces_end());
        Sites ss2(v2->second.surfaces_begin(), v2->second.surfaces_end());
        std::sort (fs2.begin(), fs2.end());
        std::sort (ss2.begin(), ss2.end());
        Sites s2;
        std::set_intersection(fs2.begin(), fs2.end(), ss2.begin(), ss2.end(),
                              std::back_inserter(s2));
//         std::cout << "s2 : " << std::endl;
//         print_sites(s2.begin(), s2.end());

        for (Sites::iterator it = s2.begin(); it != s2.end(); ++it)
          if (std::find(s1.begin(), s1.end(), *it) != s1.end())
          {
            std::cout << "filtered!!1" << std::endl;
            return false;
          }
      }

      Curved_kernel_2 &ker = Curved_kernel_2::instance();

      const Point_2 &p1 = v1->point();
      const Point_2 &p2 = v2->point();

//       ker.kernel().refine_boundary_x_2_object() (p1.xy());
//       ker.kernel().refine_boundary_y_2_object() (p1.xy());
      Rational x1l = ker.kernel().lower_boundary_x_2_object() (p1.xy());
      Rational x1u = ker.kernel().upper_boundary_x_2_object() (p1.xy());

      Rational y1l = ker.kernel().lower_boundary_y_2_object() (p1.xy());
      Rational y1u = ker.kernel().upper_boundary_y_2_object() (p1.xy());
      
      Interval x1i(x1l, x1u);
      Interval y1i(y1l, y1u);

//       ker.kernel().refine_boundary_x_2_object() (p2.xy());
//       ker.kernel().refine_boundary_y_2_object() (p2.xy());
      Rational x2l = ker.kernel().lower_boundary_x_2_object() (p2.xy());
      Rational x2u = ker.kernel().upper_boundary_x_2_object() (p2.xy());
      
      Rational y2l = ker.kernel().lower_boundary_y_2_object() (p2.xy());
      Rational y2u = ker.kernel().upper_boundary_y_2_object() (p2.xy());
      
      Interval x2i(x2l, x2u);
      Interval y2i(y2l, y2u);

//       std::cout << "11111111111111" << std::endl;
//       std::cout << "nearest sites: " << std::endl;
//       print_sites(v1->first.surfaces_begin(), v1->first.surfaces_end());

//       std::cout << "furthest sites: " << std::endl;
//       print_sites(v1->second.surfaces_begin(), v1->second.surfaces_end());

      const Xy_monotone_surface_3& near1 = *v1->first.surfaces_begin();
      const Xy_monotone_surface_3& fur1 = *v1->second.surfaces_begin();

//       std::cout << "2222222222222" << std::endl;
//       std::cout << "nearest sites: " << std::endl;
//       print_sites(v2->first.surfaces_begin(), v2->first.surfaces_end());
      
//       std::cout << "furthest sites: " << std::endl;
//       print_sites(v2->second.surfaces_begin(), v2->second.surfaces_end());

      const Xy_monotone_surface_3& near2 = *v2->first.surfaces_begin();
      const Xy_monotone_surface_3& fur2 = *v2->second.surfaces_begin();

//      std::cout << "approx comp" << std::endl;

      Interval d1 = annulus_width(x1i, y1i, near1, fur1);
      Interval d2 = annulus_width(x2i, y2i, near2, fur2);
      
      {
        using namespace boost::numeric::interval_lib::compare::tribool;
        boost::tribool res = (d1 < d2);
        if (res)
          return true;
        if (!res)
          return false;
      }
      
//      std::cout << "conv to exact (CORE)" << std::endl;
      
      Field_with_sqrt p1x, p1y;
      Nearest_traits_2::convert_to_CORE(p1, p1x, p1y);
      Field_with_sqrt p2x, p2y;
      Nearest_traits_2::convert_to_CORE(p2, p2x, p2y);

//      std::cout << "exact comp" << std::endl;
      Field_with_sqrt de1 = annulus_width(p1x, p1y, near1, fur1);
      Field_with_sqrt de2 = annulus_width(p2x, p2y, near2, fur2);
//      std::cout << "the real comp" << std::endl;
      CGAL::Comparison_result res = CGAL::compare(de1, de2);
      if (res == CGAL::SMALLER)
        return true;
      return false;
    }
};

Curve_2 construct_circle(const Rational& x, const Rational& y,
                         const Rational& r)
{
  Rational A = -2*x;
  Rational B = -2*y;
  Rational C = x*x + y*y - r*r;

  typedef CGAL::Exponent_vector EV;
  typedef std::pair<EV, Coefficient> Monomial;
  std::vector<Monomial> monomials;
  
  typedef CGAL::Fraction_traits< Rational > FT;
  FT::Decompose decompose;
  Coefficient An, Ad;
  Coefficient Bn, Bd;
  Coefficient Cn, Cd;
  decompose(A, An, Ad);
  decompose(B, Bn, Bd);
  decompose(C, Cn, Cd);
  
  Coefficient d = Ad * Bd * Cd;
  // need to multiply each numenator with the multiplication of the
  // denominator (except the current denom).
  An *= Bd*Cd;
  Bn *= Cd*Ad;
  Cn *= Ad*Bd;

  
  monomials.push_back(Monomial(EV(2, 0), d));
  monomials.push_back(Monomial(EV(0, 2), d));
  monomials.push_back(Monomial(EV(1, 0), An));
  monomials.push_back(Monomial(EV(0, 1), Bn));
  monomials.push_back(Monomial(EV(0, 0), Cn));

  typedef CGAL::Polynomial< Coefficient > Polynomial_1;
  typedef CGAL::Polynomial< Polynomial_1 > Polynomial_2;
  CGAL::Polynomial_traits_d< Polynomial_2 >::Construct_polynomial 
    construct_polynomial;
  Polynomial_2 result = construct_polynomial(monomials.begin(),
                                             monomials.end());
  return Curved_kernel_2::instance().kernel().
    construct_curve_2_object()(result);
}

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

  BOOST_STATIC_ASSERT((boost::is_same<Nearest_traits_2::Surface_3,      \
                       Furthest_traits_2::Surface_3>::value));
  
  typedef Nearest_traits_2::Surface_3           Surface_3;

  std::list<Surface_3>                          circles;

  int n;
  std::cin >> n;
  for (int i = 0; i < n; ++i)
  {  
    Rational X, Y, R;
    std::cin >> X >> Y >> R;
    std::cout << X << " " << Y << " " << R << std::endl;
    
    CGAL_assertion_msg(R >= 0, "input invalid");

    circles.push_back(Surface_3(std::make_pair(X, Y), R));
  }

  CGAL::Timer all_timer;
  CGAL::Timer apollonius_timer;
  all_timer.start();
  apollonius_timer.start();
  CGAL::lower_envelope_3 (circles.begin(), circles.end(),
                          nearest_diagram);
  apollonius_timer.stop();
  std::cout << "Apollonius:" << std::endl <<
    "T = " << apollonius_timer.time() <<
    "\tV = " << nearest_diagram.number_of_vertices() <<
    "\tE = " << nearest_diagram.number_of_edges() <<
    "\tF = " << nearest_diagram.number_of_faces() << std::endl;
  
  CGAL::Timer fpfs_timer;
  fpfs_timer.start();
  CGAL::upper_envelope_3 (circles.begin(), circles.end(),
                          furthest_diagram);
  fpfs_timer.stop();
  std::cout << "FPFS:" << std::endl <<
    "T = " << fpfs_timer.time() <<
    "\tV = " << furthest_diagram.number_of_vertices() <<
    "\tE = " << furthest_diagram.number_of_edges() <<
    "\tF = " << furthest_diagram.number_of_faces() << std::endl;

  MWA_arrangement_2            overlay_arr;
  Overlay_traits_2             overlay_traits;

  CGAL::Timer overlay_timer;
  overlay_timer.start();
  CGAL::overlay(nearest_diagram, furthest_diagram, 
                overlay_arr, overlay_traits);
  overlay_timer.stop();
  std::cout << "Overlay:" << std::endl <<
    "T = " << overlay_timer.time() <<
    "\tV = " << overlay_arr.number_of_vertices() <<
    "\tE = " << overlay_arr.number_of_edges() <<
    "\tF = " << overlay_arr.number_of_faces() << std::endl;
  
  CGAL::Timer MWA_from_cand_timer;
  MWA_from_cand_timer.start();
  typedef std::list< MWA_arrangement_2::Vertex_const_iterator > Cand_list;
  Cand_list cand_verts;
  Is_annulus_candidate cand;
  MWA_arrangement_2::Vertex_const_iterator vit;
  for (vit = overlay_arr.vertices_begin();
       vit != overlay_arr.vertices_end(); ++vit)
  {
    if (cand(*vit))
      cand_verts.push_back(vit);
  }  
  
  Cand_list::iterator it = 
    std::min_element(cand_verts.begin(), cand_verts.end(), 
                     Comp_annulus_width_at_vertex());
  MWA_from_cand_timer.stop();
  all_timer.stop();
  std::cout << "Getting MWA from candidates: " << MWA_from_cand_timer.time() <<
    "s" << std::endl;
    
  std::cout << "All timer: " << all_timer.time() << std::endl;
  
  if (vm.count("show") == false)
    return 0;

  if (it != cand_verts.end())
  {
    const Point_2 &p = (*it)->point();

    Field_with_sqrt px, py;
    Nearest_traits_2::convert_to_CORE(p, px, py);

    Field_with_sqrt neari = Nearest_traits_2::distance(
      px, py, *(*it)->first.surfaces_begin());
    Field_with_sqrt furi = Furthest_traits_2::distance(
      px, py, *(*it)->second.surfaces_begin());
    
    long round_fact = 1000;

    long x = CGAL::to_double(px) * round_fact;
    long y = CGAL::to_double(py) * round_fact;
    long near_d = CGAL::to_double(neari) * round_fact;
    long fur_d = CGAL::to_double(furi) * round_fact;
    
    typedef CGAL::Fraction_traits< Rational > FT;
    FT::Compose compose;

    Rational r_x = compose(x, round_fact);
    Rational r_y = compose(y, round_fact);
    Rational r_near = compose(near_d, round_fact);
    Rational r_fur = compose(fur_d, round_fact);

    std::cout << "x: " << x << " y: " << y << " near: " << near_d << 
      " fur: " << fur_d << std::endl;

    Curve_2 circle1 = construct_circle(r_x, r_y, r_near);
    Curve_2 circle2 = construct_circle(r_x, r_y, r_fur);
    insert (inner_outer_circles_arr, circle2);
    insert (inner_outer_circles_arr, circle1);

    std::cout << "curves: " << circle1 << " " << circle2 << std::endl;
    std::cout << "number of faces: " << 
      inner_outer_circles_arr.number_of_faces() << std::endl;
  }
  
  for (std::list<Surface_3>::iterator it = circles.begin(); 
       it != circles.end(); ++it)
  {
    Curve_2 circle1 = construct_circle(it->center().first, 
                                       it->center().second,
                                       it->r());
    insert (sites, circle1);
  }

#if defined(CGAL_MEASURE_TIME)
  std::cout << "time of overlay: " << timer_overlay.time() << std::endl;
  std::cout << "time of resolve: " << timer_resolve.time() << std::endl;
  std::cout << "time of cleanup: " << timer_clean.time() << std::endl;
  std::cout << "time of bisector: " << timer_bisector.time() << std::endl;
  std::cout << "time of copy_face: " << timer_copy_face.time() << std::endl;
  std::cout << "time of insert: " << timer_insert_curve.time() << std::endl;
  std::cout << "time of decide over face: " << timer_decide_over_face.time() 
            << std::endl;
#endif

#if defined(CGAL_MEASURE_FACES)
  std::cout << "intersecting faces: " << intersect_faces << 
    " non-intersecting faces: " << non_intersect_faces << std::endl;
#endif
          
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
