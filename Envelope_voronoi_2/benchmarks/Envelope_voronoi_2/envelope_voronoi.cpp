#include "config.h"

#include <iostream>
#include <fstream>
#include <string>

#include <CGAL/Benchmark/config.hpp>
#include <CGAL/Benchmark/Benchmark.hpp>
#include <CGAL/Benchmark/Option_parser.hpp>

#include "Option_parser.hpp"
#include "bench_classes.h"

#include <CGAL/Envelope_voronoi_2.h>

enum MaxFilesNumber 
{
  MAX_FILES_NUMBER = 20
};

// PostScript support:
#if BENCH_TRAITS != LEDA_CONIC_TRAITS && BENCH_TRAITS != CORE_CONIC_TRAITS && \
  BENCH_TRAITS != EXACUS_CONIC_TRAITS && BENCH_TRAITS != CK_CIRCLE_TRAITS && \
  BENCH_TRAITS != CK_CONIC_TRAITS
#define POSTSCRIPT_SUPPORTED 1
#endif
#undef POSTSCRIPT_SUPPORTED





namespace po = boost::program_options;
namespace cb = CGAL::benchmark;


// Window stream:
#if defined(USE_CGAL_WINDOW)
typedef CGAL::Window_stream Window_stream;
#else
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Linear_object_2.h>
#include <qapplication.h>

typedef CGAL::Qt_widget Window_stream;
QApplication * App;
#endif

inline Window_stream & operator<<(Window_stream & ws, 
                                  Envelope_diagram_2 & arr)
{
  Envelope_diagram_2::Edge_iterator ei;
  ws << CGAL::BLUE;
  for (ei = arr.edges_begin(); ei != arr.edges_end(); ++ei)
    ws << (*ei).curve();
  Envelope_diagram_2::Vertex_iterator vi;
  ws << CGAL::RED;
  for (vi = arr.vertices_begin(); vi != arr.vertices_end(); ++vi)
    ws << (*vi).point();
  return ws;
}



class Display_voronoi : public Bench_envelope_voronoi
{
public:
  Display_voronoi() 
    : m_width(DEF_WIN_WIDTH), m_height(DEF_WIN_HEIGHT)
    {}

  void op()
    {
      Envelope_diagram_2 diagram;
      CGAL::voronoi_2(_sites.begin(), _sites.end(), diagram);

      if (m_verbose_level > 0) 
      {
        if (m_verbose_level > 1) 
        {
          if (!diagram.is_valid()) 
            std::cerr << "map invalid!" << std::endl;
        }
        
        std::cout << "# of vertices: " << diagram.number_of_vertices() << std::endl;
        std::cout << "# of halfedges: " << diagram.number_of_halfedges() << std::endl;
        std::cout << "# of faces: " << diagram.number_of_faces() << std::endl;
      }
      
      double x_max = -1000000, y_max = -1000000, 
        x_min = 1000000, y_min = 1000000;
      Envelope_diagram_2::Vertex_iterator vit;
      for (vit = diagram.vertices_begin(); vit != diagram.vertices_end();
           ++vit)
      {
        double x = CGAL::to_double(vit->point().x());
        double y = CGAL::to_double(vit->point().y());

        x_max = std::max(x_max, x);
        y_max = std::max(y_max, y);
        x_min = std::min(x_min, x);
        y_min = std::min(y_min, y);
      }

      float x_range = x_max - x_min;
      float y_range = y_max - y_min;
      float height = (y_range * m_width) / x_range;
      
      float min_range = (x_range < y_range) ? x_range : y_range;
      float x_margin = min_range / 4;
      float y_margin = (height * x_margin) / m_width;
      
      m_x0 = x_min - x_margin;
      m_x1 = x_max + x_margin;
      m_y0 = y_min - y_margin;
      
      m_height = static_cast<int>(height);

#if defined(USE_CGAL_WINDOW)
      m_window->init(m_x0, m_x1, m_y0);   // logical window size 

      m_window->set_redraw(&Display_arr::redraw);
      m_window->set_mode(leda_src_mode);
      m_window->set_node_width(3);
      m_window->set_point_style(leda_cross_point);
      m_window->set_line_width(1);
      m_window->display(leda_window::center, leda_window::center);
#else
      m_y1 = y_max + y_margin;

      App->setMainWidget(m_window);
      m_window->resize(m_width, m_height);
      m_window->set_window(m_x0, m_x1, m_y0, m_y1);   // logical window size 
      m_window->setLineWidth(1);
      m_window->setPointSize(3);
      m_window->show();
#endif

#if defined(USE_CGAL_WINDOW)
      m_window->set_flush(0);
      (*m_window) << arr;
      m_window->set_flush(1);
      m_window->flush();
#else
      m_window->lock();
      *m_window << CGAL::BackgroundColor(CGAL::WHITE) << CGAL::RED;
      (*m_window) << diagram;
      m_window->unlock();
      App->flush();
#endif

#if defined(POSTSCRIPT_SUPPORTED)
/*      if (m_postscript) 
      {
        // Print to Postscript file:
        std::cout << "Print to Postscript file" << std::endl;
        CGAL::Postscript_file_stream ps_stream(m_width, m_height ,"arr.ps");
        ps_stream.init(m_x0, m_x1, m_y0);
        // CGAL::cgalize(ps_stream);
        ps_stream.set_line_width(1);
        CGAL::Arr_drawer<Arr, CGAL::Postscript_file_stream> drawer(ps_stream);
        // drawer.draw_faces(arr.faces_begin(), arr.faces_end());
        ps_stream << CGAL::BLUE;
        drawer.draw_halfedges(arr.halfedges_begin(), arr.halfedges_end());
        ps_stream << CGAL::RED;
        drawer.draw_vertices(arr.vertices_begin(), arr.vertices_end());
        
        // draw_arr(arr, drawer, ps_stream);
        // ps_stream << arr;
      }
*/
#endif

    }
  
  int init()
    {
      int rc = Bench_envelope_voronoi::init();
      if (rc < 0) 
        return rc;

#if defined(USE_CGAL_WINDOW)
      m_window = new Window_stream(m_width, m_height);
      if (!m_window) return -1;
#else
      m_window = new Window_stream();
      if (!m_window) return -1;
#endif

      
      return 0;
    }

  void clean()
    {
      Bench_envelope_voronoi::clean();
      delete m_window;
    }
  
  void set_width(unsigned int width) { m_width = width; }
  void set_height(unsigned int height) { m_height = height; }
  
private:

#if defined(USE_CGAL_WINDOW)
  static void
  redraw(leda_window * wp, double x0, double y0, double x1, double y1) 
    { wp->flush_buffer(x0,y0,x1,y1); }
#endif
  
  Window_stream * m_window;

  int m_width, m_height;
  float m_x0, m_x1, m_y0, m_y1;
};

typedef cb::Benchmark<Display_voronoi>              Display_voronoi_bench;


template <class Bench_inst, class Benchable>
void run_bench(Bench_inst & bench_inst, Benchable & benchable,
               const char * fullname, 
               int samples, int iterations, unsigned int verbose_level,
               bool postscript = false, 
               const char * fullname2 = 0)
{
  //set some variable
  benchable.set_file_name(fullname);
  benchable.set_verbose_level(verbose_level);
  benchable.set_postscript(postscript);
  if (fullname2) benchable.set_file_name(fullname2, 1);
  
  if (samples > 0) 
    bench_inst.set_samples(samples);
  else if (iterations > 0) 
    bench_inst.set_iterations(iterations);
  
  //opertor () in the Bench - does all the work !
  bench_inst();
}

template <class BenchInstance>
void run_bench_for_type(Option_parser &option_parser,
                        Option_parser::Type_code type_code,
                        unsigned int type_mask,
                        const char *file_name,
                        const char * full_name,
                        int seconds,
                        unsigned int samples,
                        int iterations,
                        unsigned int verbose_level,
                        bool postscript)
{
  if (type_mask & (0x1 << type_code))
  {
    std::string name =
      std::string(option_parser.get_type_name(type_code)) + " " +
      + " (" + std::string(file_name) + ")";

    cb::Benchmark<BenchInstance>  bench_inst(name, seconds, 
                                             false);
    BenchInstance & benchable = bench_inst.get_benchable();
    run_bench(bench_inst, benchable,
              full_name,
              samples, iterations,
              verbose_level, postscript);
  }
}


int main(int argc, char * argv[])
{
  typedef Option_parser::Type_code                    Type_code;

  Option_parser option_parser;
  try 
  {
    option_parser(argc, argv);
  } 
  catch(Option_parser::Generic_option_exception & /* e */) 
  {
    return 0;
  } 
  catch(std::exception & e) 
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  
  unsigned int verbose_level = option_parser.get_verbose_level();
  unsigned int type_mask = option_parser.get_type_mask();
  unsigned int samples = option_parser.get_samples();
  int iterations = option_parser.get_iterations();
  int seconds = option_parser.get_seconds();
  bool print_header = option_parser.is_print_header();
  int nameLength = option_parser.get_name_length();
  unsigned int files_number = option_parser.get_number_files();
  for (unsigned int i = 0; i < files_number; ++i) 
  {
    if (verbose_level > 0) 
    {
      std::cout << "filename: " << option_parser.get_file_name(i).c_str()
                << std::endl
                << "fullname: " << option_parser.get_full_name(i).c_str()
                << std::endl;
    }
  }
  
  bool postscript = option_parser.get_postscript();
  if (postscript) 
    samples = 1;
  
  cb::Benchmark_base::set_name_length(nameLength);
  if (print_header) cb::Benchmark_base::print_header();
  //end parameters from command line
  
  //start bench
  if (verbose_level > 0) 
    std::cout << "start bench " << std::endl;
  
  //general definitions
  Type_code type_code;           // operation we want to apply
  
  if (verbose_level > 0) 
  {
    std::cout << "type_mask = "<< type_mask << std::endl;
  }
  
  for (unsigned int i = 0; i < files_number; ++i)
  {
    const char * file_name = option_parser.get_file_name(i).c_str();
    const char * full_name = option_parser.get_full_name(i).c_str();
    if (!file_name || !full_name) 
    {
      std::cerr << "error in getting file names" << std::endl;
      return -1;
    }

    // envelope voronoi
    run_bench_for_type<
    Bench_envelope_voronoi>(option_parser, 
                            Option_parser::TYPE_ENVELOPE_VORONOI,
                            type_mask, file_name, full_name, 
                            seconds, samples, iterations, 
                            verbose_level, postscript);
  
    // triangulation voronoi
    run_bench_for_type<
    Bench_triangulation_voronoi>(option_parser, 
                                 Option_parser::TYPE_TRIANGULATION_VORONOI,
                                 type_mask, file_name, full_name, 
                                 seconds, samples, iterations, 
                                 verbose_level, postscript);

    // envelope apollonius
    run_bench_for_type<
    Bench_envelope_apollonius>(option_parser, 
                               Option_parser::TYPE_ENVELOPE_APOLLONIUS,
                               type_mask, file_name, full_name, 
                               seconds, samples, iterations, 
                               verbose_level, postscript);

    // CGAL apollonius
    run_bench_for_type<
    Bench_apollonius>(option_parser, 
                      Option_parser::TYPE_CGAL_APOLLONIUS,
                      type_mask, file_name, full_name, 
                      seconds, samples, iterations, 
                      verbose_level, postscript);

#ifdef USE_EXACUS

    // EXACUS apollonius
    run_bench_for_type<
    Bench_EXACUS_apollonius>(option_parser, 
                      Option_parser::TYPE_EXACUS_APOLLONIUS,
                      type_mask, file_name, full_name, 
                      seconds, samples, iterations, 
                      verbose_level, postscript);

    // EXACUS apollonius using Qdx
    run_bench_for_type<
    Bench_EXACUS_Qdx_apollonius>(option_parser, 
                                 Option_parser::TYPE_EXACUS_QDX_APOLLONIUS,
                                 type_mask, file_name, full_name, 
                                 seconds, samples, iterations, 
                                 verbose_level, postscript);
#endif // #ifdef USE_EXACUS

    // Voronoi on Sphere
    run_bench_for_type<
    Bench_sphere>(option_parser, 
                  Option_parser::TYPE_SPHERE,
                  type_mask, file_name, full_name, 
                  seconds, samples, iterations, 
                  verbose_level, postscript);

    // Construct and Display
    type_code = Option_parser::TYPE_DISPLAY;
    if (type_mask & (0x1 << type_code)) 
    {
      unsigned int width = option_parser.get_width();
      unsigned int height = option_parser.get_height();
    
#if !defined(USE_CGAL_WINDOW)
      QApplication app(argc, argv);
      App = &app;
#endif
      if (verbose_level > 0) std::cout << "TYPE_DISPLAY " << std::endl;
      std::string name =
        std::string(option_parser.get_type_name(type_code)) + " " +
        " (" + std::string(file_name) + ")";
      Display_voronoi_bench bench_inst(name, seconds, false);
      Display_voronoi & benchable = bench_inst.get_benchable();
      benchable.set_width(width);
      benchable.set_height(height);
      run_bench<Display_voronoi_bench, Display_voronoi>
        (bench_inst, benchable, full_name, 
         samples, iterations,
         verbose_level, postscript);
    }
  }
  return 0;
}



#if 0

struct Help_exception {};

int main(int argc, char * argv[])
{
  po::options_description opts("Options");
  opts.add_options()("help,h", "print help message");
  cb::Option_parser bench_opts;
  opts.add(bench_opts.get_opts()); cops
  po::variables_map var_map;

  try 
  {
    po::store(po::command_line_parser(argc, argv).options(opts).run(), var_map);
    po::notify(var_map);
    if (var_map.count("help")) 
    {
      std::cout << opts << std::endl;
      return 0;
    }
  }
  catch(std::exception & e) 
  {
    std::cerr << e.what() << std::endl;
    return 1;
  }
  
  cb::Benchmark<Bench_envelope_voronoi> bench("Envelope voronoi", bench_opts.get_seconds());
  bench();
  return 0;
}

#endif

