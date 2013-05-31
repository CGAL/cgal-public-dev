#include <CGAL/basic.h>
#include <iostream>

int main() {

#if CGAL_USE_NTL
  std::cout << "CGAL does use NTL" << std::endl;
#else
  std::cout << "CGAL does NOT use NTL" << std::endl;
#endif

#if CGAL_USE_GPU
  std::cout << "CGAL does use GPU" << std::endl;
#else
  std::cout << "CGAL does NOT use GPU" << std::endl;
#endif

  return 0;

}
