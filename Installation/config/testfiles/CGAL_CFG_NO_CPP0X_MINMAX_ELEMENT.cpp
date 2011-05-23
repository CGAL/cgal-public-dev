//| If a compiler does not support C++0x minmax_element
//| CGAL_CFG_NO_CPP0X_MINMAX_ELEMENT is set. 

#undef NDEBUG
#include <algorithm> 
#include <cassert>

int main()
{
  int f[5] = {1, 2, 3, 4, 5};
  std::pair<int*, int*> ret = std::minmax_element(f, f + 5);
  assert(*(ret.first) == 1);
  assert(*(ret.second) == 5);
}
