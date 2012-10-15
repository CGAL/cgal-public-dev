#include <CGAL/basic.h>
#include <fstream>

#ifndef CGAL_VDOL_IO_H
#define CGAL_VDOL_IO_H

namespace CGAL{
namespace VDOL_3{

void copy_file_into(const char *infilename, const char *outfilename){
  std::ifstream     in_file (infilename);
  std::ofstream     out_file (outfilename);
  std::cout << "save " << infilename << " to " << outfilename << "... " << std::flush;
  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << infilename << "!" << std::endl;
    return; 
  }
  if (! out_file.is_open()) {
    std::cerr << "Failed to open " << outfilename << "!" << std::endl;
    return; 
  }
  out_file << in_file.rdbuf();
  std::cout << "  done "<< std::endl;
}





template <class LinearKernel, typename OutputIterator>
OutputIterator read_lines_3(const LinearKernel& linear_kernel,const char *filename, OutputIterator oi){
  
  typedef typename LinearKernel::FT      FT; 
  typename LinearKernel::Construct_point_3 cpoint_3 = linear_kernel.construct_point_3_object(); 
  typename LinearKernel::Construct_line_3  cline_3  = linear_kernel.construct_line_3_object(); 

  // Open the input file.
  std::ifstream     in_file (filename);
  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << filename << "!" << std::endl;
    return oi;
  }
  // Read the lines from the file
  // The input file format should be (all coordinate values are integers):
  // <n>                                       // number of lines.
  // <a1_x> <a1_y> <a1_z> <b1_x> <b1_y> <b1_z> // line #1.
  // <a2_x> <a2_y> <a2_z> <b2_x> <b2_y> <b2_z> // line #2.
  //   :      :       :      :
  // <an_x> <an_y> <an_z> <bn_x> <bn_y> <bn_z> // line #n.
  
  // read number of lines   
  unsigned int n; in_file >> n;

  for (int k = 0; k < n; ++k) {
    int a,b,c,d,e,f;
    in_file >> a >> b >> c>> d >> e >> f;
    *oi++ = cline_3(cpoint_3(FT(a),FT(b),FT(c)),cpoint_3(FT(d),FT(e),FT(f)));
  }
  in_file.close();
}

} // namespace VDOL_3 
} // namespace CGAL 

#endif // CGAL_VDOL_IO_H 
