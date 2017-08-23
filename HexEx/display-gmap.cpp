#ifdef CGAL_LCC_USE_QT
#include "linear_cell_complex_3_viewer_qt.h"
#include "typedefs.h"

int main(int argc, char** argv)
{
  if (argc!=2)
  {
    std::cout<<"Usage: "<<argv[1]<<" filename"<<std::endl;
    std::cout<<"   will load the given filename, must be a GMap file, and display it."<<std::endl;
    return EXIT_FAILURE;
  }

  std::string str=argv[1];
  std::ifstream in(str);
  if (!in.is_open())
  {
    std::cout<<"Can't load file "<<str<<std::endl;
    return EXIT_FAILURE;
  }

  LCC_3 lcc;
  in>>lcc;
  in.close();
  display_lcc(lcc);

  return EXIT_SUCCESS;
}

#endif

