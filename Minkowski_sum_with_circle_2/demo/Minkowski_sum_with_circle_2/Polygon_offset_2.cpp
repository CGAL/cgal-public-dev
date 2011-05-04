#include <fstream>
#include <iostream>

#include "include/typedefs.h"

#include "include/MainWindow.h"

typedef Minsum_with_circle_2 MSWC;

void construct_kgon(MSWC& mswc);
void construct_all(MSWC& mswc);    
void compare_offset_construction(MSWC& mswc, const std::string& fileName);
void construct_critical_epsilon(MSWC& mswc, const std::string& fileName);

int main(int argc, char **argv)
{
  // check arguments
  if(argc > 4)
  {
    std::cerr 
      << "Usage: Polygon_offset_2 (gui) or Polygon_offset_2 InputPolygonFile "
      << "[Radius [Epsilon [Delta]]]]" << std::endl;
    return 1;

  }

  if(argc == 1)
  {
    QApplication app(argc, argv);
  
    app.setOrganizationDomain("cgal.org");
    app.setOrganizationName("CGAL");
    app.setApplicationName("Polygon_offset_2 demo");
    
    // Import resources from libCGALQt4.
    // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
    Q_INIT_RESOURCE(File);
    Q_INIT_RESOURCE(CGAL);

    MainWindow mainWindow;
    mainWindow.show();
    return app.exec();
  }
  else
  {
    // the program receives up to 4 arguments:
    // [1] the name of file with input polygon points
    // [2] the radius of the offset (1 by default)
    // [3] the epsilon (1/10 by default)
    // [4] the delta (1/40 by default)
        
    std::string fileName(argv[1]);
    NT radius(1, 1);
    NT epsilon_ratio(1, 10);
    NT epsilon(1, 10);
    NT delta_ratio(1, 4);
    NT delta(1, 40);
    
    if (argc > 2)
    {
      std::istringstream str_radius(argv[2]);
      str_radius >> radius;
      std::clog << "radius: " << radius << std::endl;
    }

    // if eps is not given - test eps from 10^-1 to 10^-8
    bool no_eps = true;
    if (argc > 3)
    {
      no_eps = false;
      std::istringstream str_epsilon(argv[3]);
      str_epsilon >> epsilon;
      std::clog << "epsilon: " << epsilon << std::endl;
    }

    // if delta is not given - test delta from 2^-2 to 2^-10
    bool no_delta = true;
    if (argc > 4)
    {
      no_eps = false;
      std::istringstream str_delta(argv[4]);
      str_delta >> delta;
      std::clog << "delta: " << delta << std::endl;
    }
    


    MSWC mswc;

    Polygon_2 polygon;
    std::ifstream ifs(fileName.c_str());
    ifs >> polygon;
    std::clog << "polygon of size " << polygon.size() << std::endl;
    
    mswc.polygon(polygon);
    mswc.offset(radius);
    mswc.epsilon((no_eps? radius * epsilon_ratio: epsilon));
    mswc.delta((no_delta? epsilon * delta_ratio: delta));
    // test dependent type
    mswc.kgon_type(Circle_app_2::Kgon_regular);
//    mswc.kgon_type(Circle_app_2::Kgon_dependent);

    // compute all things with usual out/log
    if(!no_eps)
    {
      //construct_kgon(mswc);
      construct_all(mswc);
    }
    else
    // if eps is not given - test ofset construction with eps from 10^-1 to 10^-MAX_DEGREE
    {
      if(!no_delta)
      {
        compare_offset_construction(mswc, fileName);
      }
      else
      // if delta is given - find critical epsilon for a given radius
      // TODO: find critical epsilon for a bunch of radius values
      {
        construct_critical_epsilon(mswc, fileName);
      }
    }

    return 0;
  }
}


void construct_all(MSWC& mswc)
{
  mswc.update();
}

void construct_kgon(MSWC& mswc)
{
  // update kgon
  mswc.update(MSWC::Data_kgon);  
}

void compare_offset_construction(MSWC& mswc, const std::string& fileName)
{
    NT radius = mswc.offset();
    NT epsilon = mswc.epsilon();
    
    mswc.update(MSWC::Data_exact_offset_polygon);
    static const int MAX_DEGREE = 6; //6; // 8; // 10;
    double times[MAX_DEGREE][6];
    unsigned int sizes[MAX_DEGREE];
    NT epsilon_ratio(1, 1);
    for(int degree = 1; degree <= MAX_DEGREE; ++degree)
    {
      epsilon_ratio = epsilon_ratio / 10;
      epsilon = radius * epsilon_ratio;

      mswc.epsilon(epsilon);

      mswc.update(MSWC::Data_approximate_offset_polygon);

      // update kgon
      mswc.update(MSWC::Data_kgon);
      // update sum
      mswc.update(MSWC::Data_kgon_sum);
      // reconstruct arcs
      mswc.update(MSWC::Data_kgon_induced_circles);
      // check self-intersections and construct offset
      mswc.update(MSWC::Data_kgon_offset_polygon);

      // get statistics
      int index = degree - 1;
      mswc.get_times(times[index][0], times[index][1], times[index][2], sizes[index],
        times[index][3], times[index][4], times[index][5]);

    }

    // print statistcs to the predefined file (inputFileName_time)
    std::ostringstream fileNameStream;
    fileNameStream << fileName.c_str() << "_time";
    std::string timeFileName = fileNameStream.str();
    std::ofstream time_file(timeFileName.c_str());

    time_file << "file: " << fileName << " radius: " << radius << " polygon: " << mswc.polygon().size() <<std::endl;
    time_file << "eps \texact \tapp. \tkgon\t(size) \tk_sum \tarcs \tint. \ttotal_koff" << std::endl;

    //time_file << setprecision(5) << fixed;
    for(int degree = 1; degree <= MAX_DEGREE; ++degree)
    {
      int index = degree - 1;
      time_file << -degree
                << "\t" << times[index][0]
                << "\t" << times[index][1]
                << "\t" << times[index][2] << "\t(" << sizes[index] << ")"
                << "\t" << times[index][3]
                << "\t" << times[index][4]
                << "\t" << times[index][5]
                << "\t" << (times[index][2] + times[index][3] + times[index][4] + times[index][5])
                << std::endl;
    }
}

void construct_critical_epsilon(MSWC& mswc, const std::string& fileName)
{
    // print critical epsilon value(s) to the predefined file (inputFileName_criteps)
    std::ostringstream fileNameStream;
    fileNameStream << fileName.c_str() << "_criteps";
    std::string epsFileName = fileNameStream.str();
    std::ofstream eps_file(epsFileName.c_str());

    eps_file << "file: " << fileName << " delta: " << mswc.delta() << " polygon: " << mswc.polygon().size() <<std::endl;
  
}
