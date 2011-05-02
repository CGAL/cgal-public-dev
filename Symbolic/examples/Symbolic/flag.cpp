#include <iostream>

int main() {

#if CGAL_USE_NTL
  std::cout << "CGAL does use NTL" << std::endl;
#else
  std::cout << "CGAL does NOT use NTL" << std::endl;
#endif
  return 0;

}
