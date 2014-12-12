#include <vector>
#include <deque>
#include <CGAL/Memory_sizer.h>
#include <CGAL/array.h>
#include <CGAL/Block_list.h>
#include <CGAL/Timer.h>
#include <boost/container/deque.hpp>

struct X {
  double a,b,c,d,e,f,g,h,i,j,k,l;
};

int main()
{
  CGAL::Memory_sizer ms;
  CGAL::Timer t;
  char c;
  t.start();
#if 0

  boost::container::deque<X> de;
#elif 1
  std::deque<X> de;
#else
  std::vector<X> de;
  de.reserve(100000000);
#endif

  for(int i=0; i < 10000000; i++){
    de.push_back(X());
  }
  t.stop();
  std::cout << t.time() << "sec"<< std::endl;
  std::cout << ms.virtual_size() << " " << ms.resident_size() << std::endl;

  std::cout << "cont>" << std::endl;
  std::cin >> c;
  return 0;
}
