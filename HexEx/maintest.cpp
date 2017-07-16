#include<CGAL/Linear_cell_complex_for_generalized_map.h>
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<iostream>

typedef struct dartinfo
{
  double cell_no;
  CGAL::Point_3<CGAL::Exact_predicates_inexact_constructions_kernel> parameters;
} dart_info;

struct myItem
{
  template < class GMap >
  struct Dart_wrapper
  {
    typedef dart_info Dart_info;
    typedef CGAL::Cell_attribute_with_point< GMap, int, CGAL::Tag_true>   Vertex_attribute;
    
    typedef CGAL::cpp11::tuple<Vertex_attribute> Attributes;
  };
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel                    K; 
typedef CGAL::Linear_cell_complex_traits<3, K>	                               Traits;
typedef CGAL::Linear_cell_complex_for_generalized_map<3, 3, Traits, myItem>    LCC_3;
typedef LCC_3::Dart_handle                                                     Dart_handle;
typedef LCC_3::Dart_const_handle                                               Dart_const_handle;
typedef LCC_3::Point                                       Point;

int main(int argc, char** argv)
{
  std::string str;
  if (argc==1)
  {
    std::cout<<"Enter filename: (eg: data/elephant.off)"<<std::endl;
    std::cin>>str;
  }
  else
  {
    str=argv[1];
  }

  LCC_3 lcc;
  Dart_handle d1 = lcc.make_tetrahedron(Point(-1, 0, 0), Point(0, 2, 0), 
                                        Point(1, 0, 0), Point(1, 1, 2));
  Dart_handle d2 = lcc.make_tetrahedron(Point(0, 2, -1),
                                        Point(-1, 0, -1),
                                        Point(1, 0, -1),
                                        Point(1, 1, -3));
  lcc.sew<3>(d1, d2);

  std::ofstream ofile;
  ofile.open("triangulation");
  ofile<<lcc; //works
  ofile.close();
  std::ifstream in;
  lcc.clear();
  in.open("triangulation");
  if (in.is_open())
  {
    in>>lcc; //doesn't work: segfault - figure out why
    in.close();
  }
  lcc.display_characteristics(std::cout);
  
  return EXIT_SUCCESS;
}

